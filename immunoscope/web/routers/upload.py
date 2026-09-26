"""File upload router for trajectory files."""

from __future__ import annotations
import shutil
import uuid
from pathlib import Path
from typing import List

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import JSONResponse

from immunoscope.web.config import settings


router = APIRouter(prefix="/agent", tags=["agent"])


@router.post("/upload")
async def upload_trajectory_files(
    trajectory_id: str = Form(...),
    trajectory_name: str = Form(...),
    files: List[UploadFile] = File(...),
):
    """Upload trajectory files for analysis.

    Accepts multiple files (topology, trajectory, structure) and stores them
    in a session-specific directory for agent access.
    """
    # Create upload directory
    upload_base = Path(settings.DATA_DIR) / "agent_uploads"
    upload_base.mkdir(parents=True, exist_ok=True)

    # Create trajectory-specific directory
    trajectory_dir = upload_base / trajectory_id
    trajectory_dir.mkdir(exist_ok=True)

    uploaded_files = []

    try:
        for file in files:
            if not file.filename:
                continue

            # Validate file extension
            allowed_extensions = {'.xtc', '.tpr', '.pdb', '.gro', '.trr', '.dcd'}
            file_ext = Path(file.filename).suffix.lower()

            if file_ext not in allowed_extensions:
                raise HTTPException(
                    status_code=400,
                    detail=f"File type {file_ext} not allowed. Allowed: {allowed_extensions}"
                )

            # Save file
            file_path = trajectory_dir / file.filename

            with file_path.open("wb") as buffer:
                shutil.copyfileobj(file.file, buffer)

            uploaded_files.append({
                "name": file.filename,
                "path": str(file_path),
                "size": file_path.stat().st_size,
                "type": _classify_file_type(file.filename)
            })

    except Exception as e:
        # Cleanup on error
        if trajectory_dir.exists():
            shutil.rmtree(trajectory_dir)
        raise HTTPException(status_code=500, detail=f"Upload failed: {str(e)}")

    return JSONResponse({
        "status": "success",
        "trajectory_id": trajectory_id,
        "trajectory_name": trajectory_name,
        "path": str(trajectory_dir),
        "files": uploaded_files,
        "message": f"Uploaded {len(uploaded_files)} files successfully"
    })


@router.delete("/upload/{trajectory_id}")
async def delete_trajectory_files(trajectory_id: str):
    """Delete uploaded trajectory files."""
    trajectory_dir = Path(settings.DATA_DIR) / "agent_uploads" / trajectory_id

    if not trajectory_dir.exists():
        raise HTTPException(status_code=404, detail="Trajectory not found")

    try:
        shutil.rmtree(trajectory_dir)
        return JSONResponse({
            "status": "success",
            "message": f"Trajectory {trajectory_id} deleted"
        })
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Delete failed: {str(e)}")


@router.get("/upload/{trajectory_id}")
async def get_trajectory_info(trajectory_id: str):
    """Get information about uploaded trajectory files."""
    trajectory_dir = Path(settings.DATA_DIR) / "agent_uploads" / trajectory_id

    if not trajectory_dir.exists():
        raise HTTPException(status_code=404, detail="Trajectory not found")

    files = []
    for file_path in trajectory_dir.iterdir():
        if file_path.is_file():
            files.append({
                "name": file_path.name,
                "path": str(file_path),
                "size": file_path.stat().st_size,
                "type": _classify_file_type(file_path.name)
            })

    return JSONResponse({
        "trajectory_id": trajectory_id,
        "path": str(trajectory_dir),
        "files": files
    })


@router.get("/uploads")
async def list_all_uploads():
    """List all uploaded trajectories."""
    upload_base = Path(settings.DATA_DIR) / "agent_uploads"

    if not upload_base.exists():
        return JSONResponse({"trajectories": []})

    trajectories = []
    for trajectory_dir in upload_base.iterdir():
        if trajectory_dir.is_dir():
            files = []
            for file_path in trajectory_dir.iterdir():
                if file_path.is_file():
                    files.append({
                        "name": file_path.name,
                        "size": file_path.stat().st_size,
                        "type": _classify_file_type(file_path.name)
                    })

            trajectories.append({
                "trajectory_id": trajectory_dir.name,
                "path": str(trajectory_dir),
                "files": files,
                "file_count": len(files)
            })

    return JSONResponse({"trajectories": trajectories})


def _classify_file_type(filename: str) -> str:
    """Classify file type based on extension."""
    ext = Path(filename).suffix.lower()

    type_map = {
        '.xtc': 'trajectory',
        '.trr': 'trajectory',
        '.dcd': 'trajectory',
        '.tpr': 'topology',
        '.pdb': 'structure',
        '.gro': 'structure'
    }

    return type_map.get(ext, 'unknown')
