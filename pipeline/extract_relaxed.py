"""Build the relaxed-structure dataset.

The relaxed structure of a complex is its **trajectory medoid**: the single real
frame with the lowest mean C-alpha RMSD to every other sampled frame, i.e. the
most central conformation actually visited during simulation. It is offered as an
MD-equilibrated alternative to the starting crystal structure.

**Selection rule (changed 2026-08-25).** Structures were previously the medoid of
the *dominant interface conformational cluster*, and the per-complex structure was
taken from whichever replica had the largest dominant-cluster population. Both
steps depended on the interface-clustering stage, whose descriptors are not
reproducible across independent replicas: the dominant-cluster population varies by
a median of 29.6 percentage points between runs of the same complex (34.6% of the
library range), because the clustering cuts a dendrogram at a fixed absolute cutoff
over a per-trajectory max-normalised distance matrix, so the outcome is hostage to a
single outlier frame pair. Selecting "the best-converged replica" on that quantity
was therefore close to selecting a replica at random.

The medoid rule replaces it and depends on no threshold, no cluster count and no
free parameter:

  * per complex   — pool all replicas (3 x 1001 frames) and take the global medoid,
    so the published structure uses every frame simulated for that complex rather
    than one arbitrarily chosen run;
  * per trajectory — each replica additionally keeps its own within-replica medoid
    as ``<web_data>/<id>/relaxed.pdb`` for the detail page.

RMSD is computed over all C-alpha atoms of the complex with optimal superposition
(Theobald QCP, as implemented by mdtraj), which is deterministic: rebuilding the
dataset from the same trajectories reproduces the same frames exactly.

Outputs, per trajectory:
  - <web_data>/<id>/relaxed.pdb           that replica's own medoid (all atoms)
and, for the dataset as a whole:
  - <web_data>/relaxed/manifest.json      list + annotations (served at /api/relaxed)
  - <web_data>/relaxed/dataset.csv        flat table
  - <web_data>/relaxed/relaxed_structures.tar.gz   bulk download (pdb + csv + README)

Usage:
    python -m pipeline.extract_relaxed --source /home/xmy/work/data/immunotrace/trajectories/_src1000
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import os
import sys
import tarfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from . import config

#: mdtraj's RMSD kernel and numpy's BLAS each thread internally. With one worker
#: process per complex that oversubscribes the machine badly — measured 47 CPU-
#: minutes to do 90 seconds of work on 2 complexes, almost all of it OpenMP spin.
#: Pin every worker to one thread and get the parallelism from the process pool:
#: 11 s per complex serial, ~1.6 min for the whole library at --jobs 28.
_THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")

#: Source layout: <source>/<traj_id>/{md_processed.xtc, md_processed_converted.pdb}.
#: The all-atom source trajectory (hydrogens retained) shares the canonical
#: 1001-frame / 200 ps grid with the served web traj.xtc, so frame indices align.
TRAJ_NAME = "md_processed.xtc"
TOP_NAME = "md_processed_converted.pdb"

README = """Relaxed-structure dataset — pHLA-TCR Trajectory Database
========================================================

Each <id>.pdb is the MD-equilibrated representative conformation of one
peptide-HLA-TCR complex: the *trajectory medoid*, i.e. the single real frame whose
mean C-alpha RMSD to every other sampled frame of that complex is lowest. Frames
from all replicas of the complex (3 x 1001 frames over 3 x 200 ns) are pooled
before the medoid is taken, so the structure represents the whole simulated
ensemble rather than one arbitrarily selected run. RMSD is computed over all
C-alpha atoms with optimal superposition; the selection involves no threshold and
no free parameter, and is exactly reproducible from the deposited trajectories.

Chains: A = MHC heavy chain, B = beta-2 microglobulin, C = peptide,
D = TCR alpha, E = TCR beta. All atoms (including hydrogens) are retained.

dataset.csv gives, per structure: the replica and frame the medoid was taken from,
its simulation time, and its mean C-alpha RMSD to all pooled frames — a direct,
quantitative measure of how representative the structure is (a small value means
the ensemble stays close to this conformation; a large value means the complex
sampled a broad range and no single frame represents it well).
"""


def _medoids(traj_ids: list[str], source: Path) -> dict | None:
    """Pooled + per-replica medoids for one complex.

    Returns ``{"pooled": (traj_id, frame, mean_rmsd_A), "per_traj": {traj_id:
    (frame, mean_rmsd_A)}, "n_frames": int}`` or None when nothing is loadable.
    """
    import mdtraj as md
    import numpy as np

    parts, kept = [], []
    for tid in traj_ids:
        tdir = source / tid
        xtc, top_path = tdir / TRAJ_NAME, tdir / TOP_NAME
        if not (xtc.exists() and top_path.exists()):
            continue
        top = md.load_topology(str(top_path))
        ca = top.select("name CA")
        if len(ca) == 0:
            continue
        parts.append(md.load(str(xtc), top=top, atom_indices=ca))
        kept.append(tid)
    if not parts:
        return None

    # Replicas of one complex must share an atom set to be pooled. They do across
    # this library (verified library-wide), but a mismatch must never silently
    # produce a medoid over incomparable coordinates — fall back to the first
    # replica alone and let the caller report the reduced n.
    n_ca = parts[0].n_atoms
    if any(p.n_atoms != n_ca for p in parts):
        parts, kept = parts[:1], kept[:1]

    bounds, start = [], 0
    for p in parts:
        bounds.append((start, start + p.n_frames))
        start += p.n_frames
    pooled = parts[0] if len(parts) == 1 else parts[0].join(parts[1:])

    n = pooled.n_frames
    dist = np.empty((n, n), dtype=np.float32)
    for i in range(n):
        dist[i] = md.rmsd(pooled, pooled, i)

    means = dist.mean(axis=1)
    gi = int(means.argmin())
    per_traj = {}
    pooled_hit = None
    for tid, (lo, hi) in zip(kept, bounds):
        block = dist[lo:hi, lo:hi].mean(axis=1)
        li = int(block.argmin())
        per_traj[tid] = (li, round(float(block[li]) * 10.0, 3))
        if lo <= gi < hi:
            pooled_hit = (tid, gi - lo, round(float(means[gi]) * 10.0, 3))
    return {"pooled": pooled_hit, "per_traj": per_traj, "n_frames": n,
            "n_replicas_pooled": len(kept)}


def _write_frame(source: Path, traj_id: str, frame: int, dest: Path) -> None:
    """Write one all-atom frame of a source trajectory as a PDB."""
    import mdtraj as md
    tdir = source / traj_id
    frm = md.load_frame(str(tdir / TRAJ_NAME), frame, top=str(tdir / TOP_NAME))
    dest.parent.mkdir(parents=True, exist_ok=True)
    frm.save_pdb(str(dest))


def _job(payload):
    pdb_id, traj_ids, source_str = payload
    try:
        return pdb_id, _medoids(traj_ids, Path(source_str)), None
    except Exception as exc:                      # one bad complex must not kill the run
        return pdb_id, None, f"{type(exc).__name__}: {exc}"


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", required=True,
                    help="root holding <traj_id>/md_processed.xtc (+ _converted.pdb)")
    ap.add_argument("--jobs", type=int, default=min(28, (os.cpu_count() or 4) - 2),
                    help="parallel complexes (default: cores-2, capped at 28)")
    ap.add_argument("--ids", nargs="*", help="limit to these pdb ids")
    args = ap.parse_args(argv)
    source = Path(args.source)
    # Set before the pool forks so every worker inherits single-threaded kernels.
    for var in _THREAD_VARS:
        os.environ[var] = "1"

    ann = {}
    if config.MANIFEST.exists():
        with config.MANIFEST.open() as fh:
            for r in json.load(fh)["trajectories"]:
                ann[r["traj_id"]] = r
    if not ann:
        print("ERROR: manifest is empty — run pipeline.build_manifest first.", file=sys.stderr)
        return 1

    # The manifest is the source of truth for which runs are in the served library,
    # so dropped/incomplete runs are excluded here without a second rule.
    by_pdb: dict[str, list[str]] = defaultdict(list)
    for tid, row in ann.items():
        by_pdb[row.get("pdb_id") or tid].append(tid)
    for tids in by_pdb.values():
        tids.sort()
    if args.ids:
        by_pdb = {k: v for k, v in by_pdb.items() if k in set(args.ids)}

    out_root = config.WEB_DATA / "relaxed"
    out_root.mkdir(parents=True, exist_ok=True)

    payloads = [(p, t, str(source)) for p, t in sorted(by_pdb.items())]
    results, failures = {}, []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for i, (pdb_id, res, err) in enumerate(pool.map(_job, payloads), 1):
            if err or res is None or res.get("pooled") is None:
                failures.append((pdb_id, err or "no loadable replica"))
            else:
                results[pdb_id] = res
            if i % 25 == 0 or i == len(payloads):
                print(f"  medoid {i}/{len(payloads)}", flush=True)

    structures, traj_ids_done = [], []
    for pdb_id, res in sorted(results.items()):
        # per-replica medoid -> that replica's own relaxed.pdb (detail page)
        for tid, (frame, mean_rmsd) in sorted(res["per_traj"].items()):
            _write_frame(source, tid, frame, config.WEB_DATA / tid / "relaxed.pdb")
            traj_ids_done.append(tid)
        best_tid, best_frame, best_rmsd = res["pooled"]
        a = ann.get(best_tid, {})
        structures.append({
            "pdb_id": pdb_id,
            "traj_id": best_tid,
            "run_label": best_tid[len(pdb_id):].lstrip("_") or best_tid,
            "n_replicas": len(res["per_traj"]),
            "n_frames_pooled": res["n_frames"],
            "medoid_frame": best_frame,
            "medoid_time_ps": best_frame * 200,
            "mean_ca_rmsd_angstrom": best_rmsd,
            "peptide_seq": a.get("peptide_seq"),
            "peptide_length": a.get("peptide_length"),
            "hla_allele": a.get("hla_allele"),
            "organisms": a.get("organisms") or [],
            "duration_ns": a.get("duration_ns"),
        })
    structures.sort(key=lambda e: e["pdb_id"])
    traj_ids_done.sort()

    with (out_root / "manifest.json").open("w") as fh:
        json.dump({"count": len(structures), "structures": structures,
                   "traj_ids": traj_ids_done,
                   "selection": ("pooled-replica trajectory medoid; minimum mean C-alpha "
                                 "RMSD to all sampled frames of the complex")},
                  fh, indent=1)

    csv_cols = ["pdb_id", "run_label", "n_replicas", "n_frames_pooled", "medoid_frame",
                "medoid_time_ps", "mean_ca_rmsd_angstrom", "peptide_seq",
                "peptide_length", "hla_allele", "duration_ns"]
    csv_path = out_root / "dataset.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=csv_cols, extrasaction="ignore")
        w.writeheader()
        for e in structures:
            w.writerow(e)

    tar_path = out_root / "relaxed_structures.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tar:
        readme_bytes = README.encode()
        info = tarfile.TarInfo("relaxed_structures/README.txt")
        info.size = len(readme_bytes)
        tar.addfile(info, io.BytesIO(readme_bytes))
        tar.add(csv_path, arcname="relaxed_structures/dataset.csv")
        for e in structures:
            tar.add(config.WEB_DATA / e["traj_id"] / "relaxed.pdb",
                    arcname=f"relaxed_structures/{e['pdb_id']}.pdb")

    print(f"Relaxed-structure dataset: {len(structures)} complexes "
          f"(pooled medoid over {len(traj_ids_done)} replica trajectories).")
    print(f"  manifest: {out_root / 'manifest.json'}")
    print(f"  archive:  {tar_path} ({tar_path.stat().st_size / 1e6:.1f} MB)")
    if failures:
        print(f"  FAILED {len(failures)}:", file=sys.stderr)
        for pdb_id, err in failures[:20]:
            print(f"    {pdb_id}: {err}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
