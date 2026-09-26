"""Drop hydrogens from already-extracted web trajectories.

The first version of extract.py kept all protein atoms (including hydrogens).
Hydrogens are ~half the atoms and are invisible in the cartoon view, so they
only inflate the served file size and the per-frame work the browser does
during playback. This one-off rewrites each web_data/<id>/{traj.xtc,
topology.pdb} in place to heavy atoms only and updates meta.json's atom counts.

It reads the small, already-extracted protein trajectories (not the ~574 MB
source), so it is fast and does not touch the source disk.

Usage:
    python -m pipeline.lighten_web                 # all trajectories
    python -m pipeline.lighten_web 1ao7_run2 ...   # only the given ids
"""
from __future__ import annotations

import json
import sys
import traceback

import mdtraj as md
from tqdm import tqdm

from . import config


def lighten_one(traj_id: str) -> bool:
    out_dir = config.WEB_DATA / traj_id
    topo = out_dir / config.OUT_TOPOLOGY
    xtc = out_dir / config.OUT_TRAJ
    meta_path = out_dir / config.OUT_META
    if not (topo.exists() and xtc.exists()):
        return False

    t = md.load(str(xtc), top=str(topo))
    heavy = t.topology.select("not element H")
    if heavy.size == 0 or heavy.size == t.n_atoms:
        return False  # nothing to drop (already light, or selection failed)

    sub = t.atom_slice(heavy)
    sub[0].save_pdb(str(topo))
    sub.save_xtc(str(xtc))

    if meta_path.exists():
        with meta_path.open() as fh:
            meta = json.load(fh)
        meta["n_atoms"] = int(sub.n_atoms)
        for ci, chain in enumerate(sub.topology.chains):
            if ci < len(meta.get("chains", [])):
                meta["chains"][ci]["n_atoms"] = chain.n_atoms
        with meta_path.open("w") as fh:
            json.dump(meta, fh, indent=2)
    return True


def main(argv: list[str]) -> int:
    ids = [a for a in argv if not a.startswith("--")]
    if not ids:
        ids = [d.name for d in sorted(config.WEB_DATA.glob("*")) if d.is_dir()]

    done = skipped = 0
    failures = []
    for tid in tqdm(ids, unit="traj"):
        try:
            if lighten_one(tid):
                done += 1
            else:
                skipped += 1
        except Exception as exc:  # noqa: BLE001 - keep batch going
            failures.append((tid, repr(exc)))
            tqdm.write(f"FAILED {tid}: {exc}")
            tqdm.write(traceback.format_exc())

    print(f"\nLightened {done}, skipped {skipped}, failed {len(failures)}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
