"""Job management API for the ImmunoScope web app."""

from __future__ import annotations

import shutil
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, RedirectResponse

from immunoscope.analysis.reporter.result_schema import answer_from_agent_schema, build_agent_result_schema
from immunoscope.web.config import settings
from immunoscope.web.models import (
    AnalysisConfig,
    AnalysisProfile,
    AssistantSkillInvoke,
    AssistantQuery,
    CompareJobRequest,
    InputMode,
    JobDetail,
    JobInputPaths,
    JobStatus,
    JobSummary,
)
from immunoscope.web.services.job_runner import start_job
from immunoscope.web.services.job_store import job_store
from immunoscope.web.services.result_index import build_job_result_index

router = APIRouter(tags=["jobs"])

ALLOWED_INPUT_EXTENSIONS = {".pdb", ".gro", ".tpr", ".xtc", ".trr", ".ndx", ".csv", ".json", ".zip"}


def _validate_upload(filename: str) -> None:
    suffix = Path(filename).suffix.lower()
    if suffix not in ALLOWED_INPUT_EXTENSIONS:
        raise HTTPException(status_code=400, detail=f"Unsupported file type: {suffix}")


async def _save_upload(upload: UploadFile, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as handle:
        shutil.copyfileobj(upload.file, handle)


def _split_modules(value: str) -> list[str]:
    if value.strip().lower() in {"all", ""}:
        return ["quality", "identity", "contact", "rmsf", "bsa", "rrcs", "cluster", "landscape", "report"]
    if value.strip().lower() == "none":
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def _safe_filename(filename: str) -> str:
    name = Path(filename).name
    if not name or name in {".", ".."}:
        raise HTTPException(status_code=400, detail="Invalid filename")
    return name


def _classify_input(paths: list[Path]) -> dict[str, Path]:
    classified: dict[str, Path] = {}
    for path in paths:
        suffix = path.suffix.lower()
        if suffix in {".pdb", ".gro"} and "structure" not in classified:
            classified["structure"] = path
        elif suffix == ".tpr" and "topology" not in classified:
            classified["topology"] = path
        elif suffix in {".xtc", ".trr"} and "trajectory" not in classified:
            classified["trajectory"] = path
    return classified


async def _save_named_upload(upload: UploadFile | None, input_dir: Path, label: str, input_files: list[str]) -> Path | None:
    if not upload or not upload.filename:
        return None
    _validate_upload(upload.filename)
    filename = _safe_filename(upload.filename)
    destination = input_dir / f"{label}{Path(filename).suffix.lower()}"
    await _save_upload(upload, destination)
    input_files.append(filename)
    return destination


@router.post("/jobs", response_model=JobSummary)
async def create_job(
    structure: UploadFile | None = File(default=None),
    topology: UploadFile | None = File(default=None),
    trajectory: UploadFile | None = File(default=None),
    prepared_structure: UploadFile | None = File(default=None),
    processed_trajectory: UploadFile | None = File(default=None),
    files: list[UploadFile] | None = File(default=None),
    job_name: str = Form(""),
    profile: AnalysisProfile = Form(AnalysisProfile.STANDARD),
    input_mode: InputMode = Form(InputMode.RAW),
    modules: str = Form("quality,identity,contact,cluster,report"),
    notes: str = Form(""),
    auto_start: bool = Form(True),
    preprocess_method: str = Form("2step"),
    stride: int = Form(5),
    landscape_reducer: str = Form("pca"),
    tica_lag: int = Form(10),
    continue_on_error: bool = Form(True),
    reuse_from_job_id: str = Form(""),
) -> JobSummary:
    job_id = uuid.uuid4().hex[:12]
    job_dir = settings.JOBS_DIR / job_id
    input_dir = job_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)

    input_files: list[str] = []
    saved_paths: list[Path] = []

    # Handle trajectory reuse: copy prepared trajectory from another job
    reuse_from_job_id = reuse_from_job_id.strip()
    if reuse_from_job_id:
        source_run_dir = settings.JOBS_DIR / reuse_from_job_id / "run"
        source_prep_dir = source_run_dir / "preparation" / "analysis_input"
        source_structure = source_prep_dir / "analysis_structure.pdb"
        source_trajectory = source_prep_dir / "analysis_trajectory.xtc"
        source_topology = settings.JOBS_DIR / reuse_from_job_id / "run" / "input" / "raw_topology.tpr"

        if not source_structure.exists() or not source_trajectory.exists():
            raise HTTPException(
                status_code=400,
                detail=f"Source job {reuse_from_job_id} has no reusable trajectory",
            )

        target_structure = input_dir / "prepared_structure.pdb"
        target_trajectory = input_dir / "processed_trajectory.xtc"
        shutil.copy2(source_structure, target_structure)
        shutil.copy2(source_trajectory, target_trajectory)
        input_files.extend([source_structure.name, source_trajectory.name])
        saved_paths.extend([target_structure, target_trajectory])

        # Also copy topology if available (needed for some downstream analyses)
        if source_topology.exists():
            target_topology = input_dir / "topology.tpr"
            shutil.copy2(source_topology, target_topology)
            input_files.append(source_topology.name)
            saved_paths.append(target_topology)

        # Force prepared mode
        input_mode = InputMode.PREPARED
    explicit = {
        "structure": await _save_named_upload(structure, input_dir, "structure", input_files),
        "topology": await _save_named_upload(topology, input_dir, "topology", input_files),
        "trajectory": await _save_named_upload(trajectory, input_dir, "trajectory", input_files),
        "prepared_structure": await _save_named_upload(prepared_structure, input_dir, "prepared_structure", input_files),
        "processed_trajectory": await _save_named_upload(processed_trajectory, input_dir, "processed_trajectory", input_files),
    }
    saved_paths.extend(path for path in explicit.values() if path is not None)

    for upload in files or []:
        if not upload.filename:
            continue
        _validate_upload(upload.filename)
        filename = _safe_filename(upload.filename)
        destination = input_dir / filename
        await _save_upload(upload, destination)
        input_files.append(filename)
        saved_paths.append(destination)

    inferred = _classify_input(saved_paths)
    topology_path = explicit["topology"] or inferred.get("topology")
    if input_mode == InputMode.RAW:
        input_paths = JobInputPaths(
            structure=str(explicit["structure"] or inferred.get("structure")) if (explicit["structure"] or inferred.get("structure")) else None,
            topology=str(topology_path) if topology_path else None,
            trajectory=str(explicit["trajectory"] or inferred.get("trajectory")) if (explicit["trajectory"] or inferred.get("trajectory")) else None,
        )
    else:
        input_paths = JobInputPaths(
            topology=str(topology_path) if topology_path else None,
            prepared_structure=str(explicit["prepared_structure"] or inferred.get("structure")) if (explicit["prepared_structure"] or inferred.get("structure")) else None,
            processed_trajectory=str(explicit["processed_trajectory"] or inferred.get("trajectory")) if (explicit["processed_trajectory"] or inferred.get("trajectory")) else None,
        )

    if input_mode == InputMode.RAW:
        # In RAW mode, structure is optional - will be extracted from TPR if not provided
        missing = [name for name in ["topology", "trajectory"] if not getattr(input_paths, name)]
    else:
        missing = [name for name in ["topology", "prepared_structure", "processed_trajectory"] if not getattr(input_paths, name)]
    if missing:
        raise HTTPException(status_code=400, detail=f"Missing required input(s): {', '.join(missing)}")

    module_list = _split_modules(modules)
    config = AnalysisConfig(
        profile=profile,
        modules=module_list,
        notes=notes,
        input_mode=input_mode,
        auto_start=auto_start,
        preprocess_method=preprocess_method,
        stride=stride,
        landscape_reducer=landscape_reducer,
        tica_lag=tica_lag,
        continue_on_error=continue_on_error,
    )
    name = job_name.strip() or f"immunoscope_{profile.value}_{job_id}"
    job = JobSummary(
        id=job_id,
        name=name,
        status=JobStatus.QUEUED,
        config=config,
        input_files=input_files,
        created_at=datetime.now(timezone.utc),
        progress=0,
        job_dir=str(job_dir),
        run_dir=str(job_dir / "run"),
    )
    detail = JobDetail(**job.model_dump(), inputs=input_paths)
    await job_store.save(detail)
    (job_dir / "job.log").write_text("Job created.\n", encoding="utf-8")
    if auto_start:
        start_job(job_id)
    return job


@router.get("/jobs", response_model=list[JobSummary])
async def list_jobs() -> list[JobSummary]:
    return await job_store.list_all()


@router.post("/compare/jobs", response_model=JobSummary)
async def create_compare_job(payload: CompareJobRequest) -> JobSummary:
    case_a = Path(payload.case_a).expanduser().resolve()
    case_b = Path(payload.case_b).expanduser().resolve()
    missing = [label for label, path in [("case_a", case_a), ("case_b", case_b)] if not path.exists()]
    if missing:
        raise HTTPException(status_code=400, detail=f"Missing required compare input(s): {', '.join(missing)}")

    job_id = uuid.uuid4().hex[:12]
    job_dir = settings.JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    config = AnalysisConfig(
        profile=AnalysisProfile.SAMPLING_COMPARE,
        modules=["compare"],
        notes=payload.comparison_context,
        input_mode=InputMode.PREPARED,
        auto_start=payload.auto_start,
    )
    name = payload.job_name.strip() or f"compare_{job_id}"
    job = JobSummary(
        id=job_id,
        name=name,
        status=JobStatus.QUEUED,
        config=config,
        input_files=[],
        created_at=datetime.now(timezone.utc),
        progress=0,
        job_dir=str(job_dir),
        run_dir=str(job_dir / "run"),
    )
    inputs = JobInputPaths(
        case_a=str(case_a),
        case_b=str(case_b),
        label_a=payload.label_a,
        label_b=payload.label_b,
        comparison_mode=payload.comparison_mode,
        comparison_scope=payload.comparison_scope,
        alignment_selection=payload.alignment_selection,
        residue_mapping=payload.residue_mapping,
    )
    detail = JobDetail(**job.model_dump(), inputs=inputs)
    await job_store.save(detail)
    (job_dir / "job.log").write_text("Compare job created.\n", encoding="utf-8")
    if payload.auto_start:
        start_job(job_id)
    return job


@router.get("/jobs/{job_id}", response_model=JobDetail)
async def get_job(job_id: str) -> JobDetail:
    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return job


@router.post("/jobs/{job_id}/start")
async def start_existing_job(job_id: str) -> dict[str, str]:
    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.status == JobStatus.RUNNING:
        return {"status": "already_running"}
    result = start_job(job_id)
    return {"status": result["status"]}


@router.delete("/jobs/{job_id}")
async def delete_job(job_id: str) -> dict[str, str]:
    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    job_dir = settings.JOBS_DIR / job_id
    if job_dir.exists():
        shutil.rmtree(job_dir)
    await job_store.delete(job_id)
    return {"status": "deleted"}


@router.get("/jobs/{job_id}/files")
async def list_job_files(job_id: str) -> dict[str, list[str]]:
    output_dir = settings.JOBS_DIR / job_id / "output"
    if not output_dir.exists():
        return {"files": []}
    return {"files": sorted(path.name for path in output_dir.iterdir() if path.is_file())}


@router.get("/jobs/{job_id}/files/{filename}")
async def download_job_file(job_id: str, filename: str):
    safe_name = _safe_filename(filename)
    file_path = settings.JOBS_DIR / job_id / "output" / safe_name
    if not file_path.exists() or not file_path.is_file():
        raise HTTPException(status_code=404, detail="File not found")
    return FileResponse(file_path, filename=safe_name)


@router.get("/jobs/{job_id}/result-index")
async def get_result_index(job_id: str) -> dict[str, Any]:
    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    index_path = Path(job.job_dir) / "job_result_index.json"
    if not index_path.exists() and Path(job.run_dir or "").exists():
        return build_job_result_index(Path(job.job_dir))
    if not index_path.exists():
        return {}
    return build_job_result_index(Path(job.job_dir))


@router.get("/jobs/{job_id}/report")
async def open_job_report(job_id: str):
    return RedirectResponse(url=f"/api/jobs/{job_id}/report/index.html")


@router.get("/jobs/{job_id}/report/{asset_path:path}")
async def open_job_report_asset(job_id: str, asset_path: str):
    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    report_html = job.report_html or job.run_summary.get("report_html")
    if not report_html:
        raise HTTPException(status_code=404, detail="Report not available")
    report_path = Path(report_html)
    if not report_path.exists() or not report_path.is_file():
        raise HTTPException(status_code=404, detail="Report file not found")
    report_root = report_path.parent.resolve()
    requested = report_path if asset_path in {"", "index.html"} else (report_root / asset_path).resolve()
    if report_root not in requested.parents and requested != report_root:
        raise HTTPException(status_code=400, detail="Invalid report asset path")
    if not requested.exists() or not requested.is_file():
        raise HTTPException(status_code=404, detail="Report asset not found")
    media_type = "text/html" if requested.suffix.lower() in {".html", ".htm"} else None
    return FileResponse(requested, media_type=media_type)


@router.post("/jobs/{job_id}/assistant/query")
async def query_job_assistant(job_id: str, payload: AssistantQuery) -> dict[str, Any]:
    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if not payload.question.strip():
        raise HTTPException(status_code=400, detail="Question is required")
    report_html = job.report_html or job.run_summary.get("report_html")
    if report_html and Path(report_html).exists():
        case_dir = Path(report_html).parent.parent
    else:
        case_dir = Path(job.run_dir or job.job_dir)

    # Use the new Agent run_query entry point instead of ReporterSkillOrchestrator.
    from immunoscope.agent.query import run_query

    query_result = await run_query(
        payload.question,
        tool_whitelist=["query_analysis_results"],
        max_turns=2,
        case_hint={
            "case_dir": str(case_dir),
            "job_id": job_id,
            "system_id": job.name,
        },
    )

    # Map QueryResult to the legacy response shape expected by the UI.
    result = {
        "answer": query_result.answer,
        "sources": [e.content[:200] + "..." if len(e.content) > 200 else e.content
                    for e in query_result.evidence],
        "tool_calls": query_result.tool_calls,
        "cost_usd": query_result.cost_usd,
        "truncated": query_result.truncated,
    }
    if query_result.error:
        result["error"] = query_result.error

    # Preserve the agent_result_schema for UI compatibility.
    index = build_job_result_index(Path(job.job_dir))
    schema = build_agent_result_schema(index)
    result["agent_result_schema"] = schema

    return _json_safe(result)


@router.get("/jobs/{job_id}/assistant/skills")
async def list_job_assistant_skills(job_id: str) -> dict[str, Any]:
    job = await job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    index = build_job_result_index(Path(job.job_dir))
    schema = build_agent_result_schema(index)
    return _json_safe(
        {
            "job_id": job_id,
            "scope": "single-job",
            "available_modules": schema.get("available_modules", []),
            "failed_modules": schema.get("failed_modules", []),
            "skills": schema.get("suggested_skills", []),
            "evidence_cards": schema.get("evidence_cards", []),
        }
    )


@router.post("/jobs/{job_id}/assistant/invoke")
async def invoke_job_assistant_skill(job_id: str, payload: AssistantSkillInvoke) -> dict[str, Any]:
    question = payload.question.strip() or _default_skill_question(payload.skill)
    query = AssistantQuery(
        question=question,
        skill=payload.skill,
        intent=payload.intent,
        case_b=payload.parameters.get("case_b"),
        label_a=payload.parameters.get("label_a"),
        label_b=payload.parameters.get("label_b"),
    )
    return await query_job_assistant(job_id, query)


def _default_skill_question(skill: str) -> str:
    text = str(skill or "").strip().lower().replace("_", "-")
    defaults = {
        "summarize-existing-job": "Summarize this trajectory",
        "summarize_existing_job": "Summarize this trajectory",
        "identify-hotspots": "Identify hotspots",
        "identify_hotspots": "Identify hotspots",
        "diagnose-stability": "Diagnose stability",
        "diagnose_stability": "Diagnose stability",
        "explain-fel": "Explain FEL result",
        "explain_fel": "Explain FEL result",
        "suggest-mutation-sites": "Suggest mutation sites",
        "suggest_mutation_sites": "Suggest mutation sites",
    }
    return defaults.get(text, "Summarize this trajectory")


def _should_use_schema_answer(result: dict[str, Any], schema_answer: dict[str, Any]) -> bool:
    skill = schema_answer.get("skill")
    if skill in {"explain_fel", "summarize_existing_job"}:
        return True
    confidence = str(result.get("confidence") or "").lower()
    answer = str(result.get("answer") or "")
    if confidence in {"low", "provisional"}:
        return True
    if "insufficient" in answer.lower() or "not available" in answer.lower():
        return True
    return False


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float):
        if value != value or value in {float("inf"), float("-inf")}:
            return None
    return value
