"""Trajectories router - manages processed trajectory assets for reuse."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException

from immunoscope.web.services.job_store import job_store

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/trajectories", tags=["trajectories"])


def _read_json(path: Path) -> dict[str, Any]:
    """Safely read a JSON file, return empty dict if not found."""
    if not path.exists():
        return {}
    try:
        with path.open("r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        logger.warning(f"Failed to read {path}: {e}")
        return {}


def _extract_chain_info(prep_summary: dict[str, Any]) -> dict[str, Any]:
    """Extract chain standardization info."""
    chain_std = prep_summary.get("chain_standardization", {})
    return {
        "num_chains": chain_std.get("num_chains", 0),
        "chain_mapping": chain_std.get("chain_mapping", {}),
        "status": chain_std.get("status", "unknown"),
    }


def _extract_quality_info(quality_report_path: Path) -> dict[str, Any]:
    """Extract trajectory quality metrics from preprocess report."""
    quality = _read_json(quality_report_path)
    if not quality:
        return {}

    return {
        "n_frames": quality.get("n_frames", 0),
        "duration_ps": quality.get("tail90_end_time_ps", 0),
        "duration_ns": round(quality.get("tail90_end_time_ps", 0) / 1000, 1),
        "tail90_mean_rmsd_nm": round(quality.get("tail90_mean_rmsd_nm", 0), 3),
        "tail90_std_rmsd_nm": round(quality.get("tail90_std_rmsd_nm", 0), 3),
        "convergence_quality": _assess_convergence(quality),
    }


def _assess_convergence(quality: dict[str, Any]) -> str:
    """Assess trajectory convergence quality based on RMSD."""
    std = quality.get("tail90_std_rmsd_nm", 0)
    variation = quality.get("tail90_variation_nm", 0)
    if std == 0:
        return "unknown"
    if std < 0.05 and variation < 0.2:
        return "excellent"
    if std < 0.1 and variation < 0.3:
        return "good"
    if std < 0.15:
        return "moderate"
    return "unstable"


def _get_file_size(path_str: str | None) -> int:
    """Get file size in bytes, 0 if not found."""
    if not path_str:
        return 0
    path = Path(path_str)
    if not path.exists():
        return 0
    return path.stat().st_size


@router.get("")
async def list_trajectories() -> list[dict[str, Any]]:
    """List all processed trajectories from jobs."""
    all_jobs = await job_store.list_all()
    trajectories = []

    for job in all_jobs:
        # Need full job details to access run_dir
        detail = await job_store.get(job.id)
        if not detail or not detail.run_dir:
            continue

        run_dir = Path(detail.run_dir)
        prep_dir = run_dir / "preparation"
        prep_summary_path = prep_dir / "preparation_summary.json"
        quality_report_path = prep_dir / "quality" / "preprocess_quality_report.json"

        # Skip if preparation hasn't completed
        if not prep_summary_path.exists():
            continue

        prep_summary = _read_json(prep_summary_path)
        chain_info = _extract_chain_info(prep_summary)
        quality_info = _extract_quality_info(quality_report_path)

        # Trajectory paths
        analysis_trajectory = prep_summary.get("analysis_trajectory")
        analysis_structure = prep_summary.get("analysis_structure")
        processed_trajectory = prep_summary.get("full_system_trajectory")

        trajectories.append({
            "job_id": detail.id,
            "job_name": detail.name,
            "job_status": detail.status.value if hasattr(detail.status, "value") else str(detail.status),
            "created_at": detail.created_at.isoformat() if detail.created_at else None,
            "finished_at": detail.finished_at.isoformat() if detail.finished_at else None,
            "input_files": detail.input_files or [],
            "analysis_structure": analysis_structure,
            "analysis_trajectory": analysis_trajectory,
            "processed_trajectory": processed_trajectory,
            "trajectory_size_bytes": _get_file_size(analysis_trajectory),
            "structure_size_bytes": _get_file_size(analysis_structure),
            "chain_info": chain_info,
            "quality": quality_info,
            "preprocess_method": detail.config.preprocess_method if detail.config else "unknown",
            "stride": detail.config.stride if detail.config else 5,
            "can_reuse": bool(analysis_trajectory and analysis_structure and
                              Path(analysis_trajectory).exists() and
                              Path(analysis_structure).exists()),
        })

    return trajectories


@router.get("/{job_id}")
async def get_trajectory(job_id: str) -> dict[str, Any]:
    """Get detailed info for a specific trajectory."""
    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")

    if not job.run_dir:
        raise HTTPException(status_code=404, detail=f"No trajectory for job: {job_id}")

    run_dir = Path(job.run_dir)
    prep_dir = run_dir / "preparation"
    prep_summary_path = prep_dir / "preparation_summary.json"
    quality_report_path = prep_dir / "quality" / "preprocess_quality_report.json"

    if not prep_summary_path.exists():
        raise HTTPException(status_code=404, detail=f"Trajectory not prepared yet for job: {job_id}")

    prep_summary = _read_json(prep_summary_path)
    chain_info = _extract_chain_info(prep_summary)
    quality_info = _extract_quality_info(quality_report_path)

    return {
        "job_id": job.id,
        "job_name": job.name,
        "preparation_summary": prep_summary,
        "chain_info": chain_info,
        "quality": quality_info,
    }
