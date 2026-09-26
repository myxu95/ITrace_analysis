"""Reports router for viewing generated analysis reports."""

import os
from pathlib import Path
from typing import List, Dict, Any
from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse, HTMLResponse
import logging

from immunoscope.web.config import settings

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/reports", tags=["reports"])


def scan_reports(base_dir: Path) -> List[Dict[str, Any]]:
    """Scan for HTML reports in the output directory."""
    reports = []

    if not base_dir.exists():
        logger.warning(f"Output directory does not exist: {base_dir}")
        return reports

    # Find all HTML files
    for html_file in base_dir.rglob("*.html"):
        try:
            stat = html_file.stat()
            relative_path = html_file.relative_to(base_dir)

            # Determine report type
            report_type = "comparison" if "comparison" in html_file.name.lower() else "analysis"

            # Extract a clean name
            name = html_file.stem.replace("_", " ").title()
            if "report" in name.lower():
                name = name.replace("Report", "").strip()

            reports.append({
                "name": name or html_file.stem,
                "path": str(html_file),
                "relative_path": str(relative_path),
                "type": report_type,
                "size": stat.st_size,
                "modified": stat.st_mtime,
            })
        except Exception as e:
            logger.error(f"Error processing {html_file}: {e}")
            continue

    # Sort by modification time (newest first)
    reports.sort(key=lambda x: x["modified"], reverse=True)

    return reports


@router.get("/list")
async def list_reports():
    """List all available reports."""
    try:
        output_dir = Path(settings.DATA_DIR) / "output"

        # Also check if there's a global output directory
        if not output_dir.exists():
            # Try the project root output directory
            project_root = Path(__file__).parent.parent.parent.parent
            output_dir = project_root / "output"

        reports = scan_reports(output_dir)

        return {
            "success": True,
            "reports": reports,
            "count": len(reports),
            "output_dir": str(output_dir),
        }
    except Exception as e:
        logger.error(f"Failed to list reports: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/view")
async def view_report(path: str = Query(..., description="Full path to the report file")):
    """Serve a report file for viewing."""
    try:
        report_path = Path(path)

        # Security check: ensure the path is within allowed directories
        output_dir = Path(settings.DATA_DIR) / "output"
        project_root = Path(__file__).parent.parent.parent.parent
        project_output = project_root / "output"

        allowed_dirs = [output_dir, project_output]

        is_allowed = False
        for allowed_dir in allowed_dirs:
            try:
                report_path.resolve().relative_to(allowed_dir.resolve())
                is_allowed = True
                break
            except ValueError:
                continue

        if not is_allowed:
            raise HTTPException(status_code=403, detail="Access denied")

        if not report_path.exists():
            raise HTTPException(status_code=404, detail="Report not found")

        if not report_path.is_file():
            raise HTTPException(status_code=400, detail="Path is not a file")

        # Serve the HTML file
        return FileResponse(
            path=report_path,
            media_type="text/html",
            filename=report_path.name,
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to serve report: {e}")
        raise HTTPException(status_code=500, detail=str(e))
