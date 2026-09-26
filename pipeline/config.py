"""Shared configuration for the pHLA-TCR trajectory preprocessing pipeline."""
from __future__ import annotations

import os
from pathlib import Path

# --- Source data (raw GROMACS trajectories) ------------------------------
SOURCE_ROOT = Path(
    os.environ.get(
        "IMMUNO_SOURCE_ROOT",
        "/home/xmy/work/data/immunotrace/trajectories/run3_raw",
    )
)

# Per-trajectory source file names
XTC_NAME = "md_processed.xtc"
PDB_NAME = "md_processed_converted.pdb"
RMSD_XVG = Path("analysis") / "rmsd" / "rmsd.xvg"
RMSD_PNG = Path("plots") / "rmsd" / "rmsd.png"
QUALITY_JSON = Path("quality") / "preprocess_quality_report.json"

# --- Raw full-system trajectories (with water + ions) --------------------
# The "most original" coordinates: full-system, PBC-corrected, 1000-frame
# ``pbc_1000frames.xtc`` + its ``md.tpr`` topology. These live on the compute
# cluster and are fetched lazily on first download into a local cache (see
# ``backend.app._raw_cached_file``); the served site never bulk-copies them.
RAW_XTC_NAME = "pbc_1000frames.xtc"
RAW_TPR_NAME = "md.tpr"
RAW_CACHE = Path(os.environ.get("IMMUNO_RAW_CACHE", "/home/xmy/work/data/raw_cache"))
RAW_CLUSTER_HOST = os.environ.get("IMMUNO_RAW_HOST", "128")
RAW_CLUSTER_ROOT = os.environ.get(
    "IMMUNO_RAW_ROOT", "/public/home/xmy/dataset/parallel_run/pbc_processed")

# --- Derived web assets ---------------------------------------------------
# WEB_DATA holds the served assets. The project root lives on an exfat (USB)
# volume that cannot hold symlinks and is slow to stream, so the assets are
# kept on a fast internal disk under a project-specific subdirectory of the
# shared /home/xmy/work/data folder. Override with IMMUNO_WEB_DATA.
#
# Default is web_data_1000 (the 1001-frame set BOTH :8011 live and :8012 preview
# serve as of 2026-07-17). The old 200-frame ``web_data`` is retired/backup — do
# NOT default to it, or pipeline runs write to a directory nothing serves.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
WEB_DATA = Path(os.environ.get("IMMUNO_WEB_DATA", "/home/xmy/work/data/immunotrace/web_data_1000"))
MANIFEST = WEB_DATA / "manifest.json"
RCSB_CACHE = WEB_DATA / "rcsb_cache.json"

# Per-trajectory derived file names (inside web_data/<traj_id>/)
OUT_TOPOLOGY = "topology.pdb"   # protein-only single frame
OUT_TRAJ = "traj.xtc"           # protein-only, downsampled
OUT_RMSD = "rmsd.json"          # parsed time series
OUT_META = "meta.json"          # full per-trajectory metadata

# --- Processing parameters ------------------------------------------------
# Keep every Nth frame of the (typically 1001-frame) trajectory.
STRIDE = int(os.environ.get("IMMUNO_STRIDE", "5"))

# Dirs that are known to be empty / incomplete and must be skipped.
# 6vmc: the bad initial model was rebuilt with PDBFixer (parallel1_redo) and all
# three replicas now have valid 200 ns trajectories, so 6vmc is no longer skipped.
# 6vqo_run2: was a pending resubmission; a valid 200 ns trajectory is now supplied
# (pbc_processed), so it is no longer skipped.
SKIP_DIRS = {"7Q9B_sd_run2"}

# RCSB REST API
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_POLYMER_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
RCSB_TIMEOUT = 30


def list_trajectory_dirs() -> list[Path]:
    """Return sorted list of usable trajectory directories (non-empty, with data)."""
    dirs = []
    for d in sorted(SOURCE_ROOT.iterdir()):
        if not d.is_dir():
            continue
        if d.name in SKIP_DIRS:
            continue
        if (d / XTC_NAME).exists() and (d / PDB_NAME).exists():
            dirs.append(d)
    return dirs


def pdb_id_from_traj(traj_id: str) -> str:
    """`1ao7_run2` -> `1ao7`; `1G6R_sd_run2` -> `1G6R`."""
    return traj_id.split("_", 1)[0]
