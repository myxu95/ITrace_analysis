"""Build the uniform low-frame-count *viewer* trajectory (``traj_view.xtc``).

The web ``traj.xtc`` is the analysis substrate: protein-only, hydrogen-stripped,
on the canonical 0-200000 ps / 200 ps grid (1001 frames). It is what every
in-repo descriptor stage reads and what users download, so it must never be
resampled in place.

The 3-D player only needs a light copy, so this writes a separate
``traj_view.xtc`` = ``traj.xtc[::VIEW_STRIDE]`` (1001 -> 101 frames) and records
both counts in meta.json. ``traj.xtc`` is opened read-only.

HISTORY: before 2026-08-25 this module resampled ``traj.xtc`` itself down to 200
frames (``sub.save_xtc(str(xtc))``), which destroyed the analysis substrate and
left ``n_frames_web: 200`` in 700 meta.json files long after the trajectories had
been restored to 1001 frames. Do not reintroduce an in-place write here.

Usage:
    python -m pipeline.normalize_frames                 # all trajectories
    python -m pipeline.normalize_frames 1ao7_run2 ...   # only the given ids
"""
from __future__ import annotations

import json
import os
import sys

import mdtraj as md
from tqdm import tqdm

from . import config

VIEW_STRIDE = int(os.environ.get("IMMUNO_VIEW_STRIDE", "10"))
VIEW_TRAJ = "traj_view.xtc"


def normalize_one(traj_id: str, stride: int = VIEW_STRIDE) -> tuple[int, int] | None:
    """Write ``traj_view.xtc`` for one trajectory. Returns (n_web, n_view)."""
    out_dir = config.WEB_DATA / traj_id
    topo = out_dir / config.OUT_TOPOLOGY
    xtc = out_dir / config.OUT_TRAJ
    view = out_dir / VIEW_TRAJ
    meta_path = out_dir / config.OUT_META
    if not (topo.exists() and xtc.exists()):
        return None

    t = md.load(str(xtc), top=str(topo))          # read-only; never written back
    sub = t[::stride] if t.n_frames > stride else t
    tmp = view.with_suffix(".xtc.tmp")
    sub.save_xtc(str(tmp))
    os.replace(tmp, view)

    if meta_path.exists():
        with meta_path.open() as fh:
            meta = json.load(fh)
        meta["n_frames_web"] = int(t.n_frames)      # frames in traj.xtc
        meta["n_frames_view"] = int(sub.n_frames)   # frames in traj_view.xtc
        meta["view_stride"] = int(stride)
        meta.pop("web_resampled", None)             # stale: traj.xtc is not resampled
        with meta_path.open("w") as fh:
            json.dump(meta, fh, indent=2)
    return int(t.n_frames), int(sub.n_frames)


def main(argv: list[str]) -> int:
    ids = [a for a in argv if not a.startswith("--")]
    if not ids:
        ids = [d.name for d in sorted(config.WEB_DATA.glob("*")) if d.is_dir()]

    done = 0
    counts: dict[tuple[int, int], int] = {}
    for tid in tqdm(ids, unit="traj"):
        try:
            r = normalize_one(tid)
        except Exception as exc:  # noqa: BLE001
            tqdm.write(f"FAILED {tid}: {exc}")
            continue
        if r is None:
            continue
        done += 1
        counts[r] = counts.get(r, 0) + 1
    print(f"\nViewer trajectories written: {done} (stride {VIEW_STRIDE})")
    for (nw, nv), k in sorted(counts.items()):
        print(f"  traj.xtc {nw} -> traj_view.xtc {nv}: {k}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
