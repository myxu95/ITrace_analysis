"""Extract web-ready assets from raw GROMACS trajectories.

For each trajectory directory:
  - load xtc + pdb (pdb as topology)
  - strip water / ions (keep protein only)
  - downsample frames by STRIDE
  - write topology.pdb (first frame), traj.xtc (downsampled),
    rmsd.json (parsed xvg) and meta.json (derived metadata)

Usage:
    python -m pipeline.extract                # process all trajectories
    python -m pipeline.extract 1ao7_run2 ...  # process only the given ids
    python -m pipeline.extract --force        # re-process even if outputs exist
"""
from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path

import mdtraj as md
import numpy as np
from tqdm import tqdm

from . import config

# Three-letter -> one-letter amino acid codes.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # common protonation-state variants
    "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "CYM": "C",
    "ASH": "D", "GLH": "E", "LYN": "K",
}


def parse_xvg(path: Path) -> dict:
    """Parse a GROMACS .xvg (time, rmsd) file into {time: [...], rmsd: [...]}."""
    times, vals = [], []
    if not path.exists():
        return {"time": [], "rmsd": []}
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in "#@":
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    times.append(float(parts[0]))
                    vals.append(float(parts[1]))
                except ValueError:
                    continue
    return {"time": times, "rmsd": vals}


def chain_sequence(traj: md.Trajectory, chain_index: int) -> str:
    """One-letter sequence for a chain, in residue order."""
    chain = traj.topology.chain(chain_index)
    seq = []
    for res in chain.residues:
        seq.append(THREE_TO_ONE.get(res.name.upper(), "X"))
    return "".join(seq)


def summarize_chains(traj: md.Trajectory) -> list[dict]:
    """Per-chain composition: id letter, residue count, sequence."""
    top = traj.topology
    out = []
    for ci, chain in enumerate(top.chains):
        residues = list(chain.residues)
        # mdtraj loses original chain letters; reconstruct A, B, C... by order.
        letter = chr(ord("A") + ci) if ci < 26 else str(ci)
        out.append({
            "index": ci,
            "id": letter,
            "n_residues": len(residues),
            "n_atoms": chain.n_atoms,
            "sequence": "".join(
                THREE_TO_ONE.get(r.name.upper(), "X") for r in residues
            ),
        })
    return out


def guess_peptide_chain(chains: list[dict]) -> dict | None:
    """The peptide is the shortest protein chain (typically 8-15 residues)."""
    protein_chains = [c for c in chains if 5 <= c["n_residues"] <= 30]
    if not protein_chains:
        return None
    return min(protein_chains, key=lambda c: c["n_residues"])


def process_one(traj_dir: Path, force: bool = False) -> dict:
    traj_id = traj_dir.name
    out_dir = config.WEB_DATA / traj_id
    out_dir.mkdir(parents=True, exist_ok=True)

    topo_out = out_dir / config.OUT_TOPOLOGY
    traj_out = out_dir / config.OUT_TRAJ
    meta_out = out_dir / config.OUT_META

    if not force and topo_out.exists() and traj_out.exists() and meta_out.exists():
        with meta_out.open() as fh:
            return json.load(fh)

    xtc = traj_dir / config.XTC_NAME
    pdb = traj_dir / config.PDB_NAME

    # Select protein heavy atoms from the (cheap) topology first, then load ONLY
    # those atoms and only every STRIDE-th frame. This avoids reading the full
    # ~574MB solvated trajectory into memory (≈5x less IO than load-then-slice).
    # Hydrogens are dropped: they are ~half the atoms and invisible in the
    # cartoon view, so omitting them halves the served file size and the
    # per-frame work the browser does during playback.
    top = md.load_topology(str(pdb))
    protein_sel = top.select("protein and not element H")
    if protein_sel.size == 0:
        raise ValueError("no protein atoms selected")
    sub = md.load(
        str(xtc), top=str(pdb), stride=config.STRIDE, atom_indices=protein_sel
    )

    # Write derived assets.
    sub[0].save_pdb(str(topo_out))
    sub.save_xtc(str(traj_out))

    # RMSD time series + quality stats.
    rmsd = parse_xvg(traj_dir / config.RMSD_XVG)

    # Full frame count: from quality report / rmsd series (we no longer load all
    # frames). Fall back to inferring from the strided count.
    n_frames_full = None
    with (out_dir / config.OUT_RMSD).open("w") as fh:
        json.dump(rmsd, fh)

    quality = {}
    qpath = traj_dir / config.QUALITY_JSON
    if qpath.exists():
        with qpath.open() as fh:
            quality = json.load(fh)

    # Resolve full frame count without loading all frames.
    if quality.get("n_frames"):
        n_frames_full = int(quality["n_frames"])
    elif rmsd["rmsd"]:
        n_frames_full = len(rmsd["rmsd"])
    else:
        n_frames_full = (sub.n_frames - 1) * config.STRIDE + 1

    chains = summarize_chains(sub)
    peptide = guess_peptide_chain(chains)

    has_rmsd_png = (traj_dir / config.RMSD_PNG).exists()

    meta = {
        "traj_id": traj_id,
        "pdb_id": config.pdb_id_from_traj(traj_id),
        "n_atoms": sub.n_atoms,
        "n_residues": sub.topology.n_residues,
        "n_chains": sub.topology.n_chains,
        "n_frames_full": n_frames_full,
        "n_frames_web": sub.n_frames,
        "stride": config.STRIDE,
        "chains": chains,
        "peptide_chain": peptide["id"] if peptide else None,
        "peptide_seq": peptide["sequence"] if peptide else None,
        "peptide_length": peptide["n_residues"] if peptide else None,
        "rmsd_n_points": len(rmsd["rmsd"]),
        # Simulation length in ns (xvg time is ps). Frame counts vary with save
        # interval, so duration is the meaningful, comparable quantity.
        "duration_ns": round(rmsd["time"][-1] / 1000.0, 1) if rmsd["time"] else None,
        "quality": quality,
        "has_rmsd_png": has_rmsd_png,
        # RCSB annotation gets merged in later by enrich_rcsb.py
        "rcsb": None,
    }
    with meta_out.open("w") as fh:
        json.dump(meta, fh, indent=2)
    return meta


def main(argv: list[str]) -> int:
    force = "--force" in argv
    ids = [a for a in argv if not a.startswith("--")]

    config.WEB_DATA.mkdir(parents=True, exist_ok=True)

    if ids:
        dirs = [config.SOURCE_ROOT / i for i in ids]
    else:
        dirs = config.list_trajectory_dirs()

    print(f"Processing {len(dirs)} trajectories (stride={config.STRIDE}, force={force})")
    failures = []
    for d in tqdm(dirs, unit="traj"):
        try:
            process_one(d, force=force)
        except Exception as exc:  # noqa: BLE001 - keep batch going
            failures.append((d.name, repr(exc)))
            tqdm.write(f"FAILED {d.name}: {exc}")
            tqdm.write(traceback.format_exc())

    print(f"\nDone. {len(dirs) - len(failures)} ok, {len(failures)} failed.")
    if failures:
        log = config.WEB_DATA / "extract_failures.log"
        with log.open("w") as fh:
            for name, err in failures:
                fh.write(f"{name}\t{err}\n")
        print(f"Failures logged to {log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
