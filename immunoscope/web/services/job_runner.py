"""Background execution bridge from ImmunoScope Web jobs to the CLI run layer."""

from __future__ import annotations

import argparse
import asyncio
import contextlib
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from immunoscope.cli.commands.compare import main as compare_main
from immunoscope.cli.commands.run import SingleRunOrchestrator
from immunoscope.web.models import AnalysisProfile, InputMode, JobStatus
from immunoscope.web.services.job_store import job_store
from immunoscope.web.services.result_index import build_job_result_index


_RUNNING_TASKS: dict[str, asyncio.Task] = {}


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _job_logger(job_id: str, log_path: Path) -> logging.Logger:
    logger = logging.getLogger(f"immunoscope.web.job.{job_id}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(log_path, encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
    logger.addHandler(handler)
    return logger


def _build_run_args(job) -> argparse.Namespace:
    inputs = job.inputs
    input_mode = job.config.input_mode
    run_dir = Path(job.run_dir or Path(job.job_dir) / "run")
    raw_mode = input_mode == InputMode.RAW
    return argparse.Namespace(
        action="single",
        structure=Path(inputs.structure) if raw_mode and inputs.structure else None,
        topology=Path(inputs.topology) if inputs.topology else None,
        trajectory=Path(inputs.trajectory) if raw_mode and inputs.trajectory else None,
        prepared_structure=Path(inputs.prepared_structure) if not raw_mode and inputs.prepared_structure else None,
        processed_trajectory=Path(inputs.processed_trajectory) if not raw_mode and inputs.processed_trajectory else None,
        skip_preparation=not raw_mode,
        output=run_dir,
        system_id=job.id,
        profile="standard",
        modules=",".join(job.config.modules) if job.config.modules else "none",
        preprocess_method=job.config.preprocess_method,
        analysis_group="Protein",
        dt=None,
        gmx="gmx",
        stride=job.config.stride,
        contact_cutoff=4.5,
        min_contact_frequency=0.0,
        landscape_reducer=job.config.landscape_reducer,
        tica_lag=job.config.tica_lag,
        continue_on_error=job.config.continue_on_error,
        dry_run=False,
        verbose=False,
    )


def _run_compare_sync(job, log_path: Path, logger: logging.Logger) -> int:
    output_dir = Path(job.run_dir)
    inputs = job.inputs
    argv = [
        "systems",
        "--case-a", str(inputs.case_a),
        "--case-b", str(inputs.case_b),
        "--label-a", inputs.label_a or "Condition A",
        "--label-b", inputs.label_b or "Condition B",
        "--comparison-mode", inputs.comparison_mode or "sampling",
        "--comparison-scope", inputs.comparison_scope or "same-system",
        "--alignment-selection", inputs.alignment_selection or "phla_core_ca",
        "--residue-mapping", inputs.residue_mapping or "auto",
        "-o", str(output_dir),
    ]
    if job.config.notes:
        argv.extend(["--comparison-context", job.config.notes])
    logger.info("Starting ImmunoScope compare for web job %s", job.id)
    with log_path.open("a", encoding="utf-8") as log_handle:
        with contextlib.redirect_stdout(log_handle), contextlib.redirect_stderr(log_handle):
            exit_code = compare_main(argv)
    report_html = output_dir / "comparison_report.html"
    summary = {
        "job_id": job.id,
        "status": "completed" if exit_code == 0 else "failed",
        "profile": job.config.profile,
        "modules": ["compare"],
        "module_results": {
            "compare": {
                "status": "completed" if exit_code == 0 else "failed",
                "root": str(output_dir),
                "html": str(report_html) if report_html.exists() else None,
            }
        },
        "report_html": str(report_html) if report_html.exists() else None,
        "error": None if exit_code == 0 else f"ims compare systems failed with exit code {exit_code}",
        "updated_at": _now().isoformat(),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "run_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    logger.info("Finished ImmunoScope compare for web job %s with exit code %s", job.id, exit_code)
    return exit_code


def _run_sync(job_id: str) -> int:
    job = asyncio.run(job_store.get(job_id))
    if not job:
        return 1
    job_dir = Path(job.job_dir)
    log_path = job_dir / "job.log"
    logger = _job_logger(job_id, log_path)
    if job.config.profile == AnalysisProfile.SAMPLING_COMPARE:
        return _run_compare_sync(job, log_path, logger)
    args = _build_run_args(job)
    logger.info("Starting ImmunoScope run for web job %s", job_id)
    with log_path.open("a", encoding="utf-8") as log_handle:
        with contextlib.redirect_stdout(log_handle), contextlib.redirect_stderr(log_handle):
            exit_code = SingleRunOrchestrator(args, logger).run()
    logger.info("Finished ImmunoScope run for web job %s with exit code %s", job_id, exit_code)
    return exit_code


async def run_job(job_id: str) -> None:
    job = await job_store.get(job_id)
    if not job:
        return
    await job_store.update(job_id, status=JobStatus.RUNNING, started_at=_now(), progress=5, error=None)
    try:
        exit_code = await asyncio.to_thread(_run_sync, job_id)
        index = build_job_result_index(Path(job.job_dir))
        status = JobStatus.COMPLETED if exit_code == 0 else JobStatus.FAILED
        error = None if exit_code == 0 else f"ims run single failed with exit code {exit_code}"
        await job_store.update(
            job_id,
            status=status,
            finished_at=_now(),
            progress=100,
            error=error,
            report_html=index.get("report_html"),
            result_index=str(Path(job.job_dir) / "job_result_index.json"),
        )
    except Exception as exc:  # pragma: no cover - runtime guard
        await job_store.update(job_id, status=JobStatus.FAILED, finished_at=_now(), progress=100, error=str(exc))
    finally:
        _RUNNING_TASKS.pop(job_id, None)


def start_job(job_id: str) -> dict[str, Any]:
    task = _RUNNING_TASKS.get(job_id)
    if task and not task.done():
        return {"status": "already_running"}
    task = asyncio.create_task(run_job(job_id))
    _RUNNING_TASKS[job_id] = task
    return {"status": "started"}


def get_running_job_ids() -> list[str]:
    return [job_id for job_id, task in _RUNNING_TASKS.items() if not task.done()]
