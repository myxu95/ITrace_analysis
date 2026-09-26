"""Per-trajectory equilibration QC from the whole-complex RMSD (bravo §1.2).

meta.quality was empty for every system, so the reserved QC slot (and the
manifest rmsd_* columns that read it) carried nothing, and genuine drift
outliers were served with no warning. This fills meta.quality from rmsd.json:

  * tail90_* — mean/std/max/variation (nm) over the last 90% of frames (the
    post-equilibration window the manifest already expects)
  * final_rmsd_nm, drift_slope_nm_per_ns — last-20% level and linear slope
  * equilibration — "ok" | "drifting" | "outlier", a single surfaced flag

All in nm (the unit of rmsd.json / GROMACS gmx rms). The frontend converts to Å
for display.

    IMMUNO_WEB_DATA=.../immuno-dyn python -m pipeline.compute_qc
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import config

OUTLIER_MEAN_NM = 1.0      # tail90 mean above this (10 Å) = drift outlier
OUTLIER_MAX_NM = 1.3       # tail90 max above this = outlier
DRIFT_SLOPE_NM_PER_NS = 0.0015  # positive last-20% slope above this = still rising


def assess(time_ps: list, rmsd_nm: list) -> dict | None:
    if not rmsd_nm or len(rmsd_nm) < 10:
        return None
    t = np.asarray(time_ps, float) / 1000.0  # ns
    r = np.asarray(rmsd_nm, float)
    n = len(r)
    tail = r[n // 10:]                         # last 90%
    last20 = r[int(n * 0.8):]
    t20 = t[int(n * 0.8):]
    slope = float(np.polyfit(t20, last20, 1)[0]) if len(t20) > 2 else 0.0

    q = {
        "tail90_mean_rmsd_nm": round(float(tail.mean()), 4),
        "tail90_std_rmsd_nm": round(float(tail.std()), 4),
        "tail90_max_rmsd_nm": round(float(tail.max()), 4),
        "tail90_variation_nm": round(float(tail.max() - tail.min()), 4),
        "final_rmsd_nm": round(float(last20.mean()), 4),
        "drift_slope_nm_per_ns": round(slope, 5),
    }
    if q["tail90_mean_rmsd_nm"] > OUTLIER_MEAN_NM or q["tail90_max_rmsd_nm"] > OUTLIER_MAX_NM:
        q["equilibration"] = "outlier"
    elif slope > DRIFT_SLOPE_NM_PER_NS:
        q["equilibration"] = "drifting"
    else:
        q["equilibration"] = "ok"
    return q


def main(argv=None) -> int:
    web = config.WEB_DATA
    counts = {"ok": 0, "drifting": 0, "outlier": 0}
    n_ok = n_skip = 0
    for meta_path in sorted(web.glob("*/" + config.OUT_META)):
        rmsd_path = meta_path.parent / "rmsd.json"
        if not rmsd_path.exists():
            n_skip += 1
            continue
        d = json.loads(rmsd_path.read_text())
        q = assess(d.get("time"), d.get("rmsd"))
        if q is None:
            n_skip += 1
            continue
        meta = json.loads(meta_path.read_text())
        # do-not-surface-QC: the RMSD-derived descriptors (tail90_*, final_rmsd,
        # drift_slope) and the categorical equilibration label are computed here for
        # INTERNAL QC/counting only. Nothing in the site, API or MCP reads meta.quality,
        # so the whole block (including any stale trim bookkeeping) is stripped from the
        # served meta rather than merged in.
        meta.pop("quality", None)
        meta_path.write_text(json.dumps(meta, indent=2))
        counts[q["equilibration"]] += 1
        n_ok += 1
    print(f"QC: {n_ok} written, {n_skip} skipped. {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
