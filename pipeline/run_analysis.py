"""Run the full ImmunoScope analysis suite over the pHLA-TCR trajectories.

The source trajectories under ``SOURCE_ROOT`` are already PBC-corrected and
converted (``md_processed.xtc`` + ``md_processed_converted.pdb``), so the
GROMACS-based preprocessing stage is skipped. We feed the converted PDB as both
the structure and the topology (no ``.tpr`` is available) and the processed XTC
as the trajectory, then run each analysis pipeline through ``BatchExecutor``
with a configurable number of parallel workers.

Output layout (mirrors the proven ``fulai_demo`` recipe)::

    <output_root>/
        <stage>/<traj_id>/...        # per-system, per-stage results
        manifests/<stage>_summary.json

Usage::

    python -m pipeline.run_analysis --workers 8                # all trajectories
    python -m pipeline.run_analysis --ids 1ao7_run2 --workers 1
    python -m pipeline.run_analysis --stages identity,rmsf,contact --limit 5
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from immunoscope.core.task_discovery import TaskDiscoverer
from immunoscope.pipeline import (
    AnnotatedRMSFPipeline,
    BatchExecutor,
    BiologicalIdentityPipeline,
    BSAPipeline,
    CationPiInteractionPipeline,
    ContactFrequencyPipeline,
    DockingAnglePipeline,
    HydrogenBondInteractionPipeline,
    HydrophobicInteractionPipeline,
    InterfaceClusteringPipeline,
    PiStackingInteractionPipeline,
    SaltBridgeInteractionPipeline,
)

from pipeline import config

# Analysis stages and their parameters, matching the validated fulai_demo run.
# Each entry is (stage_name, pipeline_factory). Pipelines are constructed lazily
# so an import/constructor error in one stage does not abort the others.
# ``stride`` defaults reproduce the original run exactly; --stride overrides the
# per-frame sampling for every stage that takes one (identity has no time axis).
STAGE_FACTORIES = {
    "identity": lambda stride=None: BiologicalIdentityPipeline(),
    "contact": lambda stride=10: ContactFrequencyPipeline(cutoff=4.5, stride=stride, min_frequency=0.0),
    "hbond": lambda stride=10: HydrogenBondInteractionPipeline(stride=stride, distance_cutoff=3.5, angle_cutoff=150.0),
    "saltbridge": lambda stride=10: SaltBridgeInteractionPipeline(stride=stride, distance_cutoff=4.0),
    "hydrophobic": lambda stride=10: HydrophobicInteractionPipeline(stride=stride, distance_cutoff=4.5),
    "pipi": lambda stride=10: PiStackingInteractionPipeline(stride=stride, distance_cutoff=6.5),
    "cationpi": lambda stride=10: CationPiInteractionPipeline(stride=stride, distance_cutoff=6.0, normal_angle_cutoff=60.0),
    "angle": lambda stride=10: DockingAnglePipeline(stride=stride, auto_identify_chains=True, print_each_frame=False),
    "bsa": lambda stride=10: BSAPipeline(stride=stride, probe_radius=1.4, time_unit="ps"),
    "rmsf": lambda stride=10: AnnotatedRMSFPipeline(stride=stride, time_unit="ps"),
    "inter_cluster": lambda stride=5: InterfaceClusteringPipeline(stride=stride, contact_cutoff_angstrom=4.5, distance_cutoff=0.4),
}

# Default execution order (identity first so chain roles are available; the
# interaction stages next; clustering last as it is the heaviest).
DEFAULT_STAGE_ORDER = [
    "identity",
    "contact",
    "hbond",
    "saltbridge",
    "hydrophobic",
    "pipi",
    "cationpi",
    "angle",
    "bsa",
    "rmsf",
    "inter_cluster",
]


def setup_logging(log_dir: Path) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        handlers=[
            logging.FileHandler(log_dir / "run_analysis.log", encoding="utf-8"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger("run_analysis")


def build_task_list(ids: list[str] | None, limit: int | None) -> list[dict]:
    """Build a processed-input task list from the source trajectory dirs."""
    dirs = config.list_trajectory_dirs()
    if ids:
        wanted = set(ids)
        dirs = [d for d in dirs if d.name in wanted]
    if limit:
        dirs = dirs[:limit]

    tasks = []
    for d in dirs:
        pdb = d / config.PDB_NAME
        xtc = d / config.XTC_NAME
        tasks.append(
            {
                "task_id": d.name,
                "task_root": str(d),
                # No .tpr exists; the converted PDB serves as both structure and
                # topology. MDAnalysis loads coordinates from the trajectory.
                "structure": str(pdb),
                "topology": str(pdb),
                "trajectory": str(xtc),
                "metadata": {"pdb_id": config.pdb_id_from_traj(d.name)},
            }
        )
    return tasks


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def run_trajectory(task: dict, stages: list[str], output_root_str: str, stride: int | None = None) -> dict:
    """Run all requested stages for a SINGLE trajectory, sequentially.

    Executed inside a worker process. Each stage opens its own MDAnalysis
    Universe, but because we keep one trajectory in one process the OS page
    cache keeps the ~600 MB XTC hot across stages, so the file is read from the
    (slow, USB) source disk only once instead of once per stage.
    """
    # Single-threaded BLAS so N trajectory workers do not oversubscribe cores.
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(var, "1")

    from immunoscope.core.task_discovery import TaskDiscoverer
    from immunoscope.pipeline import BatchExecutor

    output_root = Path(output_root_str)
    discoverer = TaskDiscoverer()
    report = discoverer.discover_tasks_from_list(
        [task],
        source_root=str(config.SOURCE_ROOT),
        required_files=["structure", "topology", "trajectory"],
    )

    result = {"task_id": task["task_id"], "stages": {}}
    for stage in stages:
        t0 = time.time()
        try:
            pipeline = (STAGE_FACTORIES[stage](stride) if stride is not None
                        else STAGE_FACTORIES[stage]())
            contexts = BatchExecutor(max_workers=1).execute_pipeline(
                report,
                pipeline,
                show_progress=False,
                output_base_dir=str(output_root / stage),
            )
            ctx = contexts[0] if contexts else None
            status = "success" if (ctx is not None and not ctx.has_errors()) else "failed"
            errors = list(ctx.errors) if ctx is not None else ["no context returned"]
        except Exception as exc:
            status, errors = "failed", [repr(exc)]
        result["stages"][stage] = {"status": status, "elapsed_seconds": round(time.time() - t0, 1), "errors": errors}
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Run the ImmunoScope analysis suite over the trajectory set.")
    parser.add_argument("--ids", default=None, help="Comma-separated trajectory ids (default: all).")
    parser.add_argument("--limit", type=int, default=None, help="Process at most N trajectories.")
    parser.add_argument("--workers", type=int, default=4, help="Parallel workers per stage.")
    parser.add_argument("--stages", default=None, help="Comma-separated subset of stages (default: all).")
    parser.add_argument("--stride", type=int, default=None,
                        help="Override the frame stride for every stage (default: each stage's own).")
    parser.add_argument(
        "--output",
        type=Path,
        default=config.PROJECT_ROOT / "analysis_output",
        help="Analysis output root directory.",
    )
    args = parser.parse_args(argv)

    output_root: Path = args.output
    manifest_dir = output_root / "manifests"
    logger = setup_logging(output_root / "logs")

    ids = [s.strip() for s in args.ids.split(",")] if args.ids else None
    stages = [s.strip() for s in args.stages.split(",")] if args.stages else DEFAULT_STAGE_ORDER
    unknown = [s for s in stages if s not in STAGE_FACTORIES]
    if unknown:
        parser.error(f"Unknown stage(s): {', '.join(unknown)}. Choose from: {', '.join(STAGE_FACTORIES)}")

    task_list = build_task_list(ids, args.limit)
    if not task_list:
        logger.error("No trajectories matched the selection.")
        return 1

    logger.info(
        "Selected %d trajectories; stages=%s; workers=%d stride=%s source=%s (trajectory-parallel)",
        len(task_list), ",".join(stages), args.workers,
        args.stride if args.stride is not None else "per-stage default", config.SOURCE_ROOT,
    )
    write_json(
        manifest_dir / "selected_tasks.json",
        {"count": len(task_list), "ids": [t["task_id"] for t in task_list], "tasks": task_list},
    )

    # Trajectory-outer parallelism: each worker process runs ALL stages for one
    # trajectory (reading its XTC once), N trajectories in flight at a time.
    results: dict[str, dict] = {}
    t_start = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(run_trajectory, task, stages, str(output_root), args.stride): task["task_id"]
            for task in task_list
        }
        for fut in as_completed(futures):
            tid = futures[fut]
            try:
                res = fut.result()
            except Exception as exc:
                res = {"task_id": tid, "stages": {}, "fatal": repr(exc)}
            results[tid] = res
            done += 1
            n_ok = sum(1 for s in res.get("stages", {}).values() if s.get("status") == "success")
            n_tot = len(res.get("stages", {}))
            logger.info("[%d/%d] %s: %d/%d stages ok", done, len(task_list), tid, n_ok, n_tot)
            # Persist incrementally so a crash/stop keeps prior progress.
            write_json(manifest_dir / "progress.json", {"done": done, "total": len(task_list), "results": results})

    # Per-stage roll-up across trajectories.
    overall = {}
    for stage in stages:
        ok = sum(1 for r in results.values() if r.get("stages", {}).get(stage, {}).get("status") == "success")
        overall[stage] = {"successful": ok, "total": len(results)}
    write_json(manifest_dir / "overall_summary.json", {
        "stages": overall,
        "elapsed_seconds": round(time.time() - t_start, 1),
        "workers": args.workers,
        "trajectories": len(results),
    })
    logger.info("All trajectories complete in %.1f min.", (time.time() - t_start) / 60.0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
