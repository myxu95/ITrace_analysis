"""Compare job result browser API.

Exposes the structured artifacts produced by `ims compare systems` so the
frontend can render a multi-tab detail view (Quality / Interface / Flexibility
/ Hotspots / FEL / Cluster) without forcing the user to open the HTML report.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

from immunoscope.web.services.job_store import job_store

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/compare/jobs", tags=["compare"])


def _comparison_dir(job) -> Path | None:
    """Resolve the comparison artifacts directory inside a compare job."""
    if not job or not job.run_dir:
        return None
    candidate = Path(job.run_dir) / "analysis" / "comparison"
    return candidate if candidate.exists() else None


def _read_json_safe(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as e:
        logger.warning(f"Failed to read {path}: {e}")
        return {}


def _read_csv_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    try:
        with path.open(encoding="utf-8") as f:
            return list(csv.DictReader(f))
    except OSError as e:
        logger.warning(f"Failed to read {path}: {e}")
        return []


@router.get("/{job_id}/summary")
async def get_compare_summary(job_id: str) -> dict[str, Any]:
    """Return the top-level comparison summary (case A/B, takeaways, comparability)."""
    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")

    comp_dir = _comparison_dir(job)
    if not comp_dir:
        return {
            "status": "no_artifacts",
            "message": "Comparison artifacts not yet generated (job may still be running or has failed).",
            "job_status": job.status.value if hasattr(job.status, "value") else str(job.status),
        }

    summary = _read_json_safe(comp_dir / "comparison_summary.json")
    return {
        "status": "ready",
        "job_id": job_id,
        "job_status": job.status.value if hasattr(job.status, "value") else str(job.status),
        "summary": summary,
    }


@router.get("/{job_id}/table")
async def get_compare_table(job_id: str) -> dict[str, Any]:
    """Return the flat metric comparison table grouped by category."""
    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")

    comp_dir = _comparison_dir(job)
    if not comp_dir:
        return {"status": "no_artifacts", "rows": []}

    rows = _read_csv_rows(comp_dir / "comparison_table.csv")

    # Group by category
    grouped: dict[str, list[dict]] = {}
    for r in rows:
        cat = r.get("category", "other")
        grouped.setdefault(cat, []).append(r)

    return {
        "status": "ready",
        "categories": grouped,
        "n_metrics": len(rows),
    }


@router.get("/{job_id}/tables/{name}")
async def get_compare_sub_table(job_id: str, name: str) -> dict[str, Any]:
    """Return one of the detailed comparison CSVs.

    Available names:
      - identity_comparison
      - rmsf_region_comparison
      - rrcs_region_comparison
      - interaction_family_comparison
    """
    allowed = {
        "identity_comparison",
        "rmsf_region_comparison",
        "rrcs_region_comparison",
        "interaction_family_comparison",
    }
    if name not in allowed:
        raise HTTPException(status_code=400, detail=f"Unknown table: {name}. Allowed: {sorted(allowed)}")

    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")

    comp_dir = _comparison_dir(job)
    if not comp_dir:
        return {"status": "no_artifacts", "rows": []}

    rows = _read_csv_rows(comp_dir / f"{name}.csv")
    return {"status": "ready", "name": name, "rows": rows}


@router.get("/{job_id}/plots")
async def list_compare_plots(job_id: str) -> dict[str, Any]:
    """List all PNG plots produced by the comparison."""
    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")

    comp_dir = _comparison_dir(job)
    if not comp_dir:
        return {"status": "no_artifacts", "plots": []}

    plots = []
    for path in sorted(comp_dir.glob("*.png")):
        # Derive a friendly label from the filename
        stem = path.stem
        label = stem.replace("_", " ").replace("comparison", "").strip().title()
        plots.append({
            "name": stem,
            "label": label or stem,
            "url": f"/api/compare/jobs/{job_id}/plots/{path.name}",
        })

    # Also surface enhanced HTML report if it exists
    html_paths = []
    for h in ("comparison_report_enhanced.html", "comparison_report.html"):
        candidate = Path(job.run_dir) / h
        if candidate.exists():
            html_paths.append({
                "name": h,
                "label": "Full HTML report" if "enhanced" not in h else "Enhanced HTML report",
                "url": f"/api/compare/jobs/{job_id}/report?name={h}",
            })

    return {
        "status": "ready",
        "plots": plots,
        "html_reports": html_paths,
    }


@router.get("/{job_id}/plots/{filename}")
async def get_compare_plot(job_id: str, filename: str):
    """Serve a comparison PNG plot."""
    if "/" in filename or ".." in filename or not filename.endswith(".png"):
        raise HTTPException(status_code=400, detail="Invalid filename")

    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")

    comp_dir = _comparison_dir(job)
    if not comp_dir:
        raise HTTPException(status_code=404, detail="Comparison artifacts not found")

    path = comp_dir / filename
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"Plot not found: {filename}")
    return FileResponse(str(path), media_type="image/png")


@router.get("/{job_id}/report")
async def get_compare_report(job_id: str, name: str = "comparison_report_enhanced.html"):
    """Serve the comparison HTML report."""
    if "/" in name or ".." in name or not name.endswith(".html"):
        raise HTTPException(status_code=400, detail="Invalid filename")

    job = await job_store.get(job_id)
    if not job or not job.run_dir:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")

    path = Path(job.run_dir) / name
    if not path.exists():
        # Try the alternative name
        alt = "comparison_report.html" if "enhanced" in name else "comparison_report_enhanced.html"
        alt_path = Path(job.run_dir) / alt
        if alt_path.exists():
            path = alt_path
        else:
            raise HTTPException(status_code=404, detail=f"Report not found: {name}")
    return FileResponse(str(path), media_type="text/html")
