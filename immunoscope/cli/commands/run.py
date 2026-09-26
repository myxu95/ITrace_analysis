#!/usr/bin/env python3
"""IMS Run - product-level guided analysis entry points."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
import importlib.util
import json
import logging
from pathlib import Path
import shutil
import subprocess
from typing import Any, Callable

from immunoscope.analysis.structure import PDBChainStandardizer
from immunoscope.core.task_discovery import discover_tasks
from immunoscope.core.context import PipelineContext
from immunoscope.pipeline.analysis_pipelines import ContactFrequencyPipeline


STANDARD_MODULES = ["quality", "identity", "contact", "rmsf", "bsa", "rrcs", "cluster", "landscape", "angles", "dihedrals", "report"]


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ims run",
        description="Run guided ImmunoScope workflows from raw or prepared MD inputs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ims run single --structure raw.pdb --topology md.tpr --trajectory md.xtc -o ./output/job_001

  ims run single --structure raw.pdb --topology md.tpr --trajectory md.xtc \
    -o ./output/job_001 --modules identity,rrcs,bsa,report

  ims run single --prepared-structure processed.pdb --topology md.tpr \
    --processed-trajectory md_processed.xtc -o ./output/job_001 --skip-preparation --modules none

  ims run single --prepared-structure processed.pdb --topology md.tpr \
    --processed-trajectory md_processed.xtc -o ./output/job_001 --skip-preparation

  ims run validate --structure raw.pdb --topology md.tpr --trajectory md.xtc

  ims run validate --job-dir ./output/job_001

  ims run batch --input-root ./input/md_tasks -o ./output/batch_jobs --modules contact,report --workers 4
        """,
    )
    subparsers = parser.add_subparsers(dest="action", help="Run action")

    single = subparsers.add_parser("single", help="Run a single-trajectory guided workflow")
    input_group = single.add_argument_group("raw input mode")
    input_group.add_argument("--structure", type=Path, default=None, help="Raw or representative structure PDB (optional, will extract from TPR if not provided)")
    input_group.add_argument("--topology", type=Path, required=True, help="Topology file, usually md.tpr")
    input_group.add_argument("--trajectory", type=Path, default=None, help="Raw trajectory file, usually md.xtc")

    prepared_group = single.add_argument_group("prepared input mode")
    prepared_group.add_argument("--prepared-structure", type=Path, default=None, help="Prepared/standardized structure PDB")
    prepared_group.add_argument("--processed-trajectory", type=Path, default=None, help="Already preprocessed trajectory")
    prepared_group.add_argument("--skip-preparation", action="store_true", help="Use prepared inputs directly")

    single.add_argument("-o", "--output", type=Path, required=True, help="Job output directory")
    single.add_argument("--system-id", default=None, help="Optional system/job identifier")
    single.add_argument("--profile", choices=["standard"], default="standard", help="Workflow profile")
    single.add_argument(
        "--modules",
        default="all",
        help="Comma-separated modules to run: quality,identity,contact,rmsf,bsa,rrcs,cluster,landscape,report,all,none",
    )
    single.add_argument("--preprocess-method", choices=["2step", "3step"], default="2step", help="PBC preprocessing method")
    single.add_argument("--analysis-group", default="Protein", help="GROMACS group exported for downstream analysis")
    single.add_argument("--dt", type=float, default=None, help="Optional preprocessing frame interval in ps")
    single.add_argument("--gmx", default="gmx", help="GROMACS executable")
    single.add_argument("--stride", type=int, default=5, help="Default stride for heavier downstream analyses")
    single.add_argument("--contact-cutoff", type=float, default=4.5, help="Contact cutoff in Angstrom")
    single.add_argument("--min-contact-frequency", type=float, default=0.0, help="Minimum contact frequency to keep")
    single.add_argument("--landscape-reducer", choices=["pca", "umap", "tica"], default="pca", help="Landscape reducer")
    single.add_argument("--tica-lag", type=int, default=10, help="TICA lag when --landscape-reducer tica")
    single.add_argument("--continue-on-error", action="store_true", help="Continue later stages after a module failure")
    single.add_argument("--dry-run", action="store_true", help="Write manifest plan without executing stages")
    single.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")

    validate = subparsers.add_parser("validate", help="Validate raw inputs or an existing run job")
    validate_input_group = validate.add_argument_group("raw input mode")
    validate_input_group.add_argument("--structure", type=Path, default=None, help="Raw or representative structure PDB")
    validate_input_group.add_argument("--topology", type=Path, default=None, help="Topology file, usually md.tpr")
    validate_input_group.add_argument("--trajectory", type=Path, default=None, help="Raw trajectory file, usually md.xtc")

    validate_prepared_group = validate.add_argument_group("prepared input mode")
    validate_prepared_group.add_argument("--prepared-structure", type=Path, default=None, help="Prepared/standardized structure PDB")
    validate_prepared_group.add_argument("--processed-trajectory", type=Path, default=None, help="Already preprocessed trajectory")
    validate_prepared_group.add_argument("--skip-preparation", action="store_true", help="Validate prepared inputs directly")

    validate.add_argument("--job-dir", type=Path, default=None, help="Existing ims run single job directory to validate")
    validate.add_argument("--modules", default="all", help="Comma-separated modules to validate for a planned run")
    validate.add_argument("--output", type=Path, default=None, help="Optional validation JSON output path")
    validate.add_argument("--gmx", default="gmx", help="GROMACS executable")
    validate.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")

    batch = subparsers.add_parser("batch", help="Run guided single-job workflow over multiple tasks")
    batch.add_argument("--input-root", type=Path, required=True, help="Directory containing task subdirectories")
    batch.add_argument("-o", "--output", type=Path, required=True, help="Batch output directory")
    batch.add_argument("--workers", type=int, default=1, help="Number of parallel single-job workers")
    batch.add_argument("--task-depth", type=int, default=1, help="Task directory depth below input root")
    batch.add_argument("--system-prefix", default="", help="Optional prefix for generated job IDs")
    batch.add_argument("--profile", choices=["standard"], default="standard", help="Workflow profile")
    batch.add_argument("--modules", default="all", help="Comma-separated modules to run for each task")
    batch.add_argument("--skip-preparation", action="store_true", help="Treat discovered structure/trajectory as already prepared")
    batch.add_argument("--preprocess-method", choices=["2step", "3step"], default="2step", help="PBC preprocessing method")
    batch.add_argument("--analysis-group", default="Protein", help="GROMACS group exported for downstream analysis")
    batch.add_argument("--dt", type=float, default=None, help="Optional preprocessing frame interval in ps")
    batch.add_argument("--gmx", default="gmx", help="GROMACS executable")
    batch.add_argument("--stride", type=int, default=5, help="Default stride for heavier downstream analyses")
    batch.add_argument("--contact-cutoff", type=float, default=4.5, help="Contact cutoff in Angstrom")
    batch.add_argument("--min-contact-frequency", type=float, default=0.0, help="Minimum contact frequency to keep")
    batch.add_argument("--landscape-reducer", choices=["pca", "umap", "tica"], default="pca", help="Landscape reducer")
    batch.add_argument("--tica-lag", type=int, default=10, help="TICA lag when --landscape-reducer tica")
    batch.add_argument("--continue-on-error", action="store_true", help="Continue later stages after a module failure")
    batch.add_argument("--dry-run", action="store_true", help="Write per-task manifest plans without executing stages")
    batch.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return parser


def setup_logging(verbose: bool = False) -> None:
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format="%(levelname)s: %(message)s")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _json_safe(value: Any) -> Any:
    if is_dataclass(value):
        return _json_safe(asdict(value))
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_json_safe(payload), ensure_ascii=False, indent=2), encoding="utf-8")


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


class RunValidator:
    def __init__(self, args: argparse.Namespace, logger: logging.Logger):
        self.args = args
        self.logger = logger
        self.checks: list[dict[str, Any]] = []

    def run(self) -> int:
        if self.args.job_dir:
            report = self._validate_job_dir(self.args.job_dir.resolve())
        else:
            report = self._validate_inputs()

        output_path = self.args.output
        if output_path is None and self.args.job_dir:
            output_path = self.args.job_dir / "run_validation.json"
        if output_path is not None:
            _write_json(output_path, report)

        status = report["status"]
        passed = report["summary"]["passed"]
        failed = report["summary"]["failed"]
        warnings = report["summary"]["warnings"]
        self.logger.info(f"Validation status: {status} ({passed} passed, {failed} failed, {warnings} warnings)")
        for check in report["checks"]:
            level = check["level"].upper()
            self.logger.info(f"[{level}] {check['name']}: {check['message']}")
        return 0 if status == "passed" else 1

    def _add_check(self, name: str, passed: bool, message: str, level: str = "error", details: dict[str, Any] | None = None) -> None:
        self.checks.append(
            {
                "name": name,
                "passed": passed,
                "level": "info" if passed else level,
                "message": message,
                "details": details or {},
            }
        )

    def _finalize(self, mode: str, target: str | None = None) -> dict[str, Any]:
        failed = [check for check in self.checks if not check["passed"] and check["level"] == "error"]
        warnings = [check for check in self.checks if not check["passed"] and check["level"] == "warning"]
        return {
            "schema_version": "immunoscope.run_validation.v1",
            "mode": mode,
            "target": target,
            "status": "failed" if failed else "passed",
            "summary": {
                "total": len(self.checks),
                "passed": sum(1 for check in self.checks if check["passed"]),
                "failed": len(failed),
                "warnings": len(warnings),
            },
            "checks": self.checks,
            "updated_at": _now(),
        }

    def _parse_modules(self, value: str) -> list[str]:
        return SingleRunOrchestrator._parse_modules(self, value)

    def _validate_inputs(self) -> dict[str, Any]:
        mode = "prepared" if self.args.skip_preparation else "raw"
        required = (
            [("topology", self.args.topology), ("prepared_structure", self.args.prepared_structure), ("processed_trajectory", self.args.processed_trajectory)]
            if self.args.skip_preparation
            else [("topology", self.args.topology), ("structure", self.args.structure), ("trajectory", self.args.trajectory)]
        )
        for label, path in required:
            exists = bool(path and path.exists())
            self._add_check(
                f"input:{label}",
                exists,
                f"{label} exists: {path}" if exists else f"Missing required input: {label}",
                details={"path": str(path) if path else None},
            )

        try:
            modules = self._parse_modules(self.args.modules)
            self._add_check("modules", True, f"Selected modules: {', '.join(modules) if modules else 'none'}", details={"modules": modules})
            if "report" in modules and len([module for module in modules if module != "report"]) == 0:
                self._add_check("module_dependencies", False, "report is selected without upstream analysis modules", level="warning")
        except Exception as exc:
            self._add_check("modules", False, str(exc))

        if not self.args.skip_preparation:
            gmx_path = shutil.which(self.args.gmx)
            self._add_check(
                "gromacs",
                bool(gmx_path),
                f"GROMACS executable found: {gmx_path}" if gmx_path else f"GROMACS executable not found: {self.args.gmx}",
                details={"executable": self.args.gmx, "path": gmx_path},
            )

        return self._finalize(mode=f"{mode}_input")

    def _validate_job_dir(self, job_dir: Path) -> dict[str, Any]:
        self._add_check("job_dir", job_dir.exists() and job_dir.is_dir(), f"Job directory exists: {job_dir}", details={"path": str(job_dir)})
        manifest_path = job_dir / "run_manifest.json"
        summary_path = job_dir / "run_summary.json"
        manifest = self._load_required_json("manifest", manifest_path)
        summary = self._load_required_json("summary", summary_path)
        if manifest:
            self._validate_manifest(manifest, job_dir)
        if summary:
            self._validate_summary(summary, job_dir, manifest)
        return self._finalize(mode="job", target=str(job_dir))

    def _load_required_json(self, label: str, path: Path) -> dict[str, Any] | None:
        if not path.exists():
            self._add_check(label, False, f"Missing {path.name}", details={"path": str(path)})
            return None
        try:
            payload = _read_json(path)
        except Exception as exc:
            self._add_check(label, False, f"Invalid JSON in {path.name}: {exc}", details={"path": str(path)})
            return None
        self._add_check(label, True, f"Loaded {path.name}", details={"path": str(path)})
        return payload

    def _validate_manifest(self, manifest: dict[str, Any], job_dir: Path) -> None:
        for key in ["job_id", "profile", "status", "modules", "stages", "prepared_input"]:
            self._add_check(f"manifest:{key}", key in manifest, f"manifest contains {key}")
        modules = manifest.get("modules", [])
        stages = manifest.get("stages", [])
        stage_names = [stage.get("name") for stage in stages if isinstance(stage, dict)]
        for expected in ["preflight", "input", "preparation"]:
            self._add_check(f"stage:{expected}", expected in stage_names, f"stage present: {expected}")
        for module in modules:
            self._add_check(f"stage:module:{module}", f"module:{module}" in stage_names, f"module stage present: {module}")

        for stage in stages:
            if not isinstance(stage, dict):
                self._add_check("stage_schema", False, "Stage entry is not an object")
                continue
            status = stage.get("status")
            name = stage.get("name", "<unnamed>")
            self._add_check(f"stage_status:{name}", status in {"planned", "running", "completed", "failed"}, f"{name} status is {status}")
            if status == "failed":
                has_error = bool(stage.get("error") or (isinstance(stage.get("details"), dict) and stage["details"].get("error")))
                self._add_check(f"stage_error:{name}", has_error, f"{name} has failure detail")

        for label, raw_path in manifest.get("prepared_input", {}).items():
            if label == "mode" or raw_path is None:
                continue
            path = Path(raw_path)
            self._add_check(f"prepared_input:{label}", path.exists(), f"prepared input exists: {label}", details={"path": str(path)})

    def _validate_summary(self, summary: dict[str, Any], job_dir: Path, manifest: dict[str, Any] | None) -> None:
        for key in ["job_id", "status", "manifest", "modules", "module_results", "report_html"]:
            self._add_check(f"summary:{key}", key in summary, f"summary contains {key}")
        manifest_ref = summary.get("manifest")
        if manifest_ref:
            self._add_check("summary:manifest_ref", Path(manifest_ref).exists(), "summary manifest reference exists", details={"path": str(manifest_ref)})

        module_results = summary.get("module_results", {})
        for module, result in module_results.items():
            if not isinstance(result, dict):
                self._add_check(f"module_result:{module}", False, f"module result is not an object: {module}")
                continue
            status = result.get("status")
            self._add_check(f"module_status:{module}", status in {"completed", "failed"}, f"{module} status is {status}")
            if status == "completed":
                self._validate_module_artifacts(module, result)
            if status == "failed":
                self._add_check(f"module_error:{module}", bool(result.get("error")), f"{module} has failure detail")

        report_html = summary.get("report_html")
        if report_html:
            self._add_check("report_html", Path(report_html).exists(), "report HTML exists", details={"path": str(report_html)})
        elif "report" in summary.get("modules", []):
            report_result = module_results.get("report", {})
            self._add_check(
                "report_html",
                report_result.get("status") != "completed",
                "report_html is missing although report completed" if report_result.get("status") == "completed" else "report_html unavailable because report did not complete",
            )

        if manifest:
            self._add_check("job_id_match", summary.get("job_id") == manifest.get("job_id"), "summary job_id matches manifest job_id")
            self._add_check("module_list_match", summary.get("modules") == manifest.get("modules"), "summary modules match manifest modules")

    def _validate_module_artifacts(self, module: str, result: dict[str, Any]) -> None:
        root = result.get("root")
        if root:
            self._add_check(f"module_root:{module}", Path(root).exists(), f"{module} root exists", details={"path": str(root)})
        artifact_keys = {
            "contact": ["contact_report"],
            "report": ["html"],
        }
        for key in artifact_keys.get(module, []):
            value = result.get(key)
            self._add_check(
                f"module_artifact:{module}:{key}",
                bool(value and Path(value).exists()),
                f"{module} artifact exists: {key}",
                details={"path": str(value) if value else None},
            )


class BatchRunOrchestrator:
    def __init__(self, args: argparse.Namespace, logger: logging.Logger):
        self.args = args
        self.logger = logger
        self.input_root = args.input_root.resolve()
        self.batch_root = args.output.resolve()
        self.jobs_dir = self.batch_root / "jobs"
        self.manifest_path = self.batch_root / "batch_manifest.json"
        self.summary_json_path = self.batch_root / "batch_summary.json"
        self.summary_csv_path = self.batch_root / "batch_summary.csv"

    def run(self) -> int:
        self.batch_root.mkdir(parents=True, exist_ok=True)
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        try:
            discovery = discover_tasks(
                self.input_root,
                task_depth=self.args.task_depth,
                required_files=["structure", "topology", "trajectory"],
            )
        except Exception as exc:
            self.logger.error(f"Task discovery failed: {exc}")
            return 1

        manifest = self._build_manifest(discovery)
        _write_json(self.manifest_path, manifest)
        if discovery.num_valid == 0:
            self.logger.error("No valid tasks found for guided batch run")
            self._write_batch_summary(discovery, [])
            return 1

        self.logger.info(f"Guided batch run: {discovery.num_valid} valid task(s), {self.args.workers} worker(s)")
        results: list[dict[str, Any]] = []
        with ThreadPoolExecutor(max_workers=max(1, self.args.workers)) as executor:
            future_map = {executor.submit(self._run_task, task): task for task in discovery.valid_tasks}
            for future in as_completed(future_map):
                task = future_map[future]
                try:
                    result = future.result()
                except Exception as exc:
                    result = {
                        "task_id": task.task_id,
                        "job_id": self._job_id(task.task_id),
                        "status": "failed",
                        "exit_code": 1,
                        "error": str(exc),
                    }
                results.append(result)
                self.logger.info("[%s] %s", result["status"], result["task_id"])

        summary = self._write_batch_summary(discovery, results)
        self.logger.info("Batch summary JSON: %s", self.summary_json_path)
        self.logger.info("Batch summary CSV: %s", self.summary_csv_path)
        return 0 if summary["failed"] == 0 else 1

    def _run_task(self, task) -> dict[str, Any]:
        job_id = self._job_id(task.task_id)
        job_dir = self.jobs_dir / job_id
        input_files = task.input_files
        single_args = argparse.Namespace(
            action="single",
            structure=None if self.args.skip_preparation else Path(input_files.structure_path),
            topology=Path(input_files.topology_path),
            trajectory=None if self.args.skip_preparation else Path(input_files.trajectory_path),
            prepared_structure=Path(input_files.structure_path) if self.args.skip_preparation else None,
            processed_trajectory=Path(input_files.trajectory_path) if self.args.skip_preparation else None,
            skip_preparation=self.args.skip_preparation,
            output=job_dir,
            system_id=job_id,
            profile=self.args.profile,
            modules=self.args.modules,
            preprocess_method=self.args.preprocess_method,
            analysis_group=self.args.analysis_group,
            dt=self.args.dt,
            gmx=self.args.gmx,
            stride=self.args.stride,
            contact_cutoff=self.args.contact_cutoff,
            min_contact_frequency=self.args.min_contact_frequency,
            landscape_reducer=self.args.landscape_reducer,
            tica_lag=self.args.tica_lag,
            continue_on_error=self.args.continue_on_error,
            dry_run=self.args.dry_run,
            verbose=self.args.verbose,
        )
        exit_code = SingleRunOrchestrator(single_args, self.logger).run()
        validation_code = RunValidator(
            argparse.Namespace(
                job_dir=job_dir,
                output=job_dir / "run_validation.json",
                structure=None,
                topology=None,
                trajectory=None,
                prepared_structure=None,
                processed_trajectory=None,
                skip_preparation=False,
                modules=self.args.modules,
                gmx=self.args.gmx,
            ),
            self.logger,
        ).run()
        summary_path = job_dir / "run_summary.json"
        run_summary = _read_json(summary_path) if summary_path.exists() else {}
        return {
            "task_id": task.task_id,
            "job_id": job_id,
            "status": "completed" if exit_code == 0 and validation_code == 0 else "failed",
            "exit_code": exit_code,
            "validation_exit_code": validation_code,
            "job_dir": str(job_dir),
            "manifest": str(job_dir / "run_manifest.json"),
            "summary": str(summary_path),
            "validation": str(job_dir / "run_validation.json"),
            "report_html": run_summary.get("report_html"),
            "error": run_summary.get("error"),
        }

    def _job_id(self, task_id: str) -> str:
        return f"{self.args.system_prefix}{task_id}" if self.args.system_prefix else task_id

    def _build_manifest(self, discovery) -> dict[str, Any]:
        return {
            "schema_version": "immunoscope.run_batch_manifest.v1",
            "input_root": str(self.input_root),
            "output_root": str(self.batch_root),
            "created_at": _now(),
            "parameters": {
                "profile": self.args.profile,
                "modules": self.args.modules,
                "workers": self.args.workers,
                "task_depth": self.args.task_depth,
                "skip_preparation": self.args.skip_preparation,
                "dry_run": self.args.dry_run,
            },
            "discovery": {
                "total_tasks": discovery.total_tasks,
                "valid_tasks": discovery.num_valid,
                "invalid_tasks": discovery.num_invalid,
                "ambiguous_tasks": discovery.num_ambiguous,
            },
            "tasks": [
                {
                    "task_id": task.task_id,
                    "task_root": task.task_root,
                    "validation_status": task.validation_status,
                    "validation_messages": list(task.validation_messages),
                    "job_id": self._job_id(task.task_id),
                    "job_dir": str(self.jobs_dir / self._job_id(task.task_id)),
                }
                for task in discovery.all_tasks
            ],
        }

    def _write_batch_summary(self, discovery, results: list[dict[str, Any]]) -> dict[str, Any]:
        successful = [item for item in results if item.get("status") == "completed"]
        failed = [item for item in results if item.get("status") != "completed"]
        summary = {
            "schema_version": "immunoscope.run_batch_summary.v1",
            "input_root": str(self.input_root),
            "output_root": str(self.batch_root),
            "manifest": str(self.manifest_path),
            "total_discovered": discovery.total_tasks,
            "valid_tasks": discovery.num_valid,
            "invalid_tasks": discovery.num_invalid,
            "ambiguous_tasks": discovery.num_ambiguous,
            "successful": len(successful),
            "failed": len(failed),
            "results": sorted(results, key=lambda item: item.get("task_id", "")),
            "updated_at": _now(),
        }
        _write_json(self.summary_json_path, summary)
        self._write_summary_csv(summary["results"])
        return summary

    def _write_summary_csv(self, results: list[dict[str, Any]]) -> None:
        fieldnames = ["task_id", "job_id", "status", "exit_code", "validation_exit_code", "job_dir", "report_html", "error"]
        with self.summary_csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for result in results:
                writer.writerow({key: result.get(key, "") for key in fieldnames})


class SingleRunOrchestrator:
    def __init__(self, args: argparse.Namespace, logger: logging.Logger):
        self.args = args
        self.logger = logger
        self.job_root = args.output.resolve()
        self.system_id = args.system_id or self.job_root.name
        self.input_dir = self.job_root / "input"
        self.preparation_dir = self.job_root / "preparation"
        self.analysis_dir = self.job_root / "analysis"
        self.report_dir = self.job_root / "report"
        self.manifest_path = self.job_root / "run_manifest.json"
        self.summary_path = self.job_root / "run_summary.json"
        self.module_results: dict[str, dict[str, Any]] = {}
        self.manifest: dict[str, Any] = {
            "job_id": self.system_id,
            "profile": args.profile,
            "created_at": _now(),
            "updated_at": _now(),
            "status": "created",
            "input": {},
            "prepared_input": {},
            "modules": self._parse_modules(args.modules),
            "stages": [],
        }

    def run(self) -> int:
        self.job_root.mkdir(parents=True, exist_ok=True)
        self._write_manifest(status="planning")
        if self.args.dry_run:
            self._populate_dry_run_plan()
            self._write_manifest(status="planned")
            self._write_summary(status="planned")
            self.logger.info(f"Dry-run manifest written: {self.manifest_path}")
            return 0

        try:
            self._stage("preflight", self._run_preflight)
            self._stage("input", self._prepare_input_files)
            if self.args.skip_preparation:
                self._stage("preparation", self._use_prepared_inputs)
            else:
                self._stage("preparation", self._run_preparation)
            self._run_analysis_modules()
            self._write_manifest(status="completed")
            self._write_summary(status="completed")
            self.logger.info("Guided single run completed")
            self.logger.info(f"Run manifest: {self.manifest_path}")
            self.logger.info(f"Run summary: {self.summary_path}")
            return 0
        except Exception as exc:
            self.logger.exception("Guided single run failed")
            self._write_manifest(status="failed")
            self._write_summary(status="failed", error=str(exc))
            return 1

    def _parse_modules(self, value: str) -> list[str]:
        selected = [item.strip() for item in value.split(",") if item.strip()]
        if not selected or "all" in selected:
            return list(STANDARD_MODULES)
        if "none" in selected:
            return []
        invalid = [item for item in selected if item not in STANDARD_MODULES]
        if invalid:
            raise ValueError(f"Unknown module(s): {', '.join(invalid)}")
        return selected

    def _populate_dry_run_plan(self) -> None:
        self.manifest["input"] = {
            "structure": str(self.args.structure) if self.args.structure else None,
            "topology": str(self.args.topology) if self.args.topology else None,
            "trajectory": str(self.args.trajectory) if self.args.trajectory else None,
            "prepared_structure": str(self.args.prepared_structure) if self.args.prepared_structure else None,
            "processed_trajectory": str(self.args.processed_trajectory) if self.args.processed_trajectory else None,
        }
        self.manifest["stages"] = [
            {"name": "preflight", "status": "planned", "details": {}},
            {"name": "input", "status": "planned", "details": {"mode": "prepared" if self.args.skip_preparation else "raw"}},
            {"name": "preparation", "status": "planned", "details": {"skip": bool(self.args.skip_preparation)}},
            *[
                {"name": f"module:{module_name}", "status": "planned", "details": {"module": module_name}}
                for module_name in self.manifest["modules"]
            ],
        ]

    def _stage(self, name: str, func: Callable[[], dict[str, Any] | None]) -> None:
        record = {"name": name, "status": "running", "started_at": _now(), "finished_at": None, "details": {}}
        self.manifest["stages"].append(record)
        self._write_manifest(status="running")
        self.logger.info("=" * 60)
        self.logger.info(f"Stage: {name}")
        self.logger.info("=" * 60)
        try:
            details = func() or {}
            record["details"] = details
            record["status"] = "failed" if isinstance(details, dict) and details.get("status") == "failed" else "completed"
            record["finished_at"] = _now()
            self._write_manifest(status="running")
        except Exception as exc:
            record["status"] = "failed"
            record["error"] = str(exc)
            record["finished_at"] = _now()
            self._write_manifest(status="failed")
            raise

    def _run_preflight(self) -> dict[str, Any]:
        checks: dict[str, Any] = {
            "mode": "prepared" if self.args.skip_preparation else "raw",
            "modules": list(self.manifest["modules"]),
            "gmx": None,
            "output_dir": str(self.job_root),
            "inputs": {},
        }

        self.job_root.mkdir(parents=True, exist_ok=True)
        probe = self.job_root / ".write_test"
        probe.write_text("ok", encoding="utf-8")
        probe.unlink()

        input_paths = {
            "topology": self.args.topology,
            "structure": self.args.structure,
            "trajectory": self.args.trajectory,
            "prepared_structure": self.args.prepared_structure,
            "processed_trajectory": self.args.processed_trajectory,
        }
        for label, path in input_paths.items():
            if path is not None:
                checks["inputs"][label] = {"path": str(path), "exists": path.exists()}

        if self.args.skip_preparation:
            required = [("topology", self.args.topology), ("prepared_structure", self.args.prepared_structure), ("processed_trajectory", self.args.processed_trajectory)]
        else:
            required = [("topology", self.args.topology), ("structure", self.args.structure), ("trajectory", self.args.trajectory)]

        missing = [label for label, path in required if path is None or not path.exists()]
        if missing:
            raise FileNotFoundError(f"Missing required input(s): {', '.join(missing)}")

        if not self.args.skip_preparation:
            gmx_path = shutil.which(self.args.gmx)
            if not gmx_path:
                raise FileNotFoundError(f"GROMACS executable not found: {self.args.gmx}")
            checks["gmx"] = {"executable": self.args.gmx, "path": gmx_path}
            try:
                completed = subprocess.run(
                    [self.args.gmx, "--version"],
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=15,
                )
                checks["gmx"]["returncode"] = completed.returncode
                checks["gmx"]["version_line"] = next(
                    (line.strip() for line in completed.stdout.splitlines() if "GROMACS version" in line),
                    "",
                )
            except Exception as exc:
                checks["gmx"]["warning"] = str(exc)

        return checks

    def _prepare_input_files(self) -> dict[str, Any]:
        self.input_dir.mkdir(parents=True, exist_ok=True)
        topology = self._copy_input(self.args.topology, "raw_topology")
        self.manifest["input"]["topology"] = str(topology)

        if self.args.skip_preparation:
            if not self.args.prepared_structure or not self.args.processed_trajectory:
                raise ValueError("--skip-preparation requires --prepared-structure and --processed-trajectory")
            prepared_structure = self._copy_input(self.args.prepared_structure, "prepared_structure")
            processed_trajectory = self._copy_input(self.args.processed_trajectory, "processed_trajectory")
            self.manifest["input"]["prepared_structure"] = str(prepared_structure)
            self.manifest["input"]["processed_trajectory"] = str(processed_trajectory)
        else:
            if not self.args.trajectory:
                raise ValueError("raw input mode requires --trajectory")

            # Structure handling:
            # - No structure: extract PDB from TPR
            # - .gro structure: extract PDB from TPR (GRO lacks chain info)
            # - .pdb structure: use directly
            if not self.args.structure:
                logging.info("No structure provided, extracting PDB from TPR...")
                structure = self._extract_structure_from_tpr(topology)
            elif self.args.structure.suffix.lower() == ".gro":
                logging.info("GRO format detected (no chain info), extracting PDB from TPR instead...")
                structure = self._extract_structure_from_tpr(topology)
            else:
                structure = self._copy_input(self.args.structure, "raw_structure")

            trajectory = self._copy_input(self.args.trajectory, "raw_trajectory")
            self.manifest["input"]["structure"] = str(structure)
            self.manifest["input"]["trajectory"] = str(trajectory)

        return dict(self.manifest["input"])

    def _copy_input(self, source: Path, stem: str) -> Path:
        source = source.resolve()
        if not source.exists():
            raise FileNotFoundError(f"Input file not found: {source}")
        target = self.input_dir / f"{stem}{source.suffix}"
        if source != target.resolve():
            shutil.copy2(source, target)
        return target

    def _extract_structure_from_tpr(self, tpr_path: Path) -> Path:
        """Extract PDB structure from TPR file using gmx editconf."""
        output_pdb = self.input_dir / "raw_structure.pdb"
        cmd = [self.args.gmx, "editconf", "-f", str(tpr_path), "-o", str(output_pdb)]

        logging.info(f"Extracting structure from TPR: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(f"Failed to extract structure from TPR: {result.stderr}")

        if not output_pdb.exists():
            raise RuntimeError(f"Structure extraction succeeded but output file not found: {output_pdb}")

        logging.info(f"Successfully extracted structure to {output_pdb}")
        return output_pdb

    def _use_prepared_inputs(self) -> dict[str, Any]:
        structure = Path(self.manifest["input"]["prepared_structure"]).resolve()
        trajectory = Path(self.manifest["input"]["processed_trajectory"]).resolve()
        topology = Path(self.manifest["input"]["topology"]).resolve()
        self.preparation_dir.mkdir(parents=True, exist_ok=True)
        self.manifest["prepared_input"] = {
            "structure": str(structure),
            "topology": str(topology),
            "trajectory": str(trajectory),
            "mode": "prepared",
        }
        summary = {
            "status": "completed",
            "mode": "prepared",
            "structure": str(structure),
            "topology": str(topology),
            "trajectory": str(trajectory),
        }
        _write_json(self.preparation_dir / "preparation_summary.json", summary)
        return summary

    def _run_preparation(self) -> dict[str, Any]:
        self.preparation_dir.mkdir(parents=True, exist_ok=True)
        raw_structure = Path(self.manifest["input"]["structure"]).resolve()
        raw_trajectory = Path(self.manifest["input"]["trajectory"]).resolve()
        topology = Path(self.manifest["input"]["topology"]).resolve()

        standardized_structure = self.preparation_dir / "standardized_structure.pdb"
        standardizer = PDBChainStandardizer()
        std_result = standardizer.process_single(
            raw_structure,
            standardized_structure,
            task_name=self.system_id,
            skip_if_standard=False,
        )
        _write_json(self.preparation_dir / "chain_standardization_summary.json", asdict(std_result))
        if std_result.status not in {"OK", "ALREADY_STANDARD"} or not standardized_structure.exists():
            raise RuntimeError(f"Chain standardization failed: {std_result.status} {std_result.error_message or ''}".strip())

        from immunoscope.cli.commands import preprocess as preprocess_command

        processed_trajectory = self.preparation_dir / "processed_trajectory.xtc"
        preprocess_args = [
            "-f", str(raw_trajectory),
            "-s", str(topology),
            "-o", str(processed_trajectory),
            "-m", self.args.preprocess_method,
            "--gmx", self.args.gmx,
        ]
        if self.args.dt is not None:
            preprocess_args.extend(["--dt", str(self.args.dt)])
        if self.args.verbose:
            preprocess_args.append("-v")
        ret = preprocess_command.main(preprocess_args)
        if ret != 0:
            raise RuntimeError(f"Trajectory preprocessing failed with exit code {ret}")

        processed_structure = processed_trajectory.with_name(f"{processed_trajectory.stem}_converted.pdb")
        analysis_input = self._extract_analysis_inputs(
            processed_trajectory=processed_trajectory,
            topology=topology,
            group_name=self.args.analysis_group,
        )
        semantic_structure = analysis_input["structure"]
        viewer_structure = analysis_input["structure"]
        self.manifest["prepared_input"] = {
            "structure": str(semantic_structure.resolve()),
            "viewer_structure": str(viewer_structure.resolve()),
            "topology": str(analysis_input["topology"].resolve()),
            "trajectory": str(analysis_input["trajectory"].resolve()),
            "full_system_topology": str(topology),
            "full_system_trajectory": str(processed_trajectory.resolve()),
            "processed_structure": str(processed_structure.resolve()) if processed_structure.exists() else None,
            "analysis_group": self.args.analysis_group,
            "mode": "raw_prepared",
        }
        summary = {
            "status": "completed",
            "mode": "raw_prepared",
            "chain_standardization": asdict(std_result),
            "standardized_structure": str(standardized_structure.resolve()),
            "analysis_structure": str(analysis_input["structure"].resolve()),
            "analysis_topology": str(analysis_input["topology"].resolve()),
            "analysis_trajectory": str(analysis_input["trajectory"].resolve()),
            "analysis_group": self.args.analysis_group,
            "full_system_trajectory": str(processed_trajectory.resolve()),
            "processed_structure": str(processed_structure.resolve()) if processed_structure.exists() else None,
        }
        _write_json(self.preparation_dir / "preparation_summary.json", summary)
        return summary

    def _extract_analysis_inputs(self, processed_trajectory: Path, topology: Path, group_name: str) -> dict[str, Path]:
        analysis_dir = self.preparation_dir / "analysis_input"
        analysis_dir.mkdir(parents=True, exist_ok=True)
        analysis_trajectory = analysis_dir / "analysis_trajectory.xtc"
        analysis_structure = analysis_dir / "analysis_structure.pdb"

        self._run_trjconv_group(
            topology=topology,
            input_trajectory=processed_trajectory,
            output_path=analysis_trajectory,
            group_name=group_name,
        )
        self._run_trjconv_group(
            topology=topology,
            input_trajectory=processed_trajectory,
            output_path=analysis_structure,
            group_name=group_name,
            dump_time=0,
        )
        return {
            "structure": analysis_structure,
            "topology": analysis_structure,
            "trajectory": analysis_trajectory,
        }

    def _run_trjconv_group(
        self,
        *,
        topology: Path,
        input_trajectory: Path,
        output_path: Path,
        group_name: str,
        dump_time: int | None = None,
    ) -> None:
        command = [
            self.args.gmx,
            "trjconv",
            "-s",
            str(topology),
            "-f",
            str(input_trajectory),
            "-o",
            str(output_path),
        ]
        if dump_time is not None:
            command.extend(["-dump", str(dump_time)])

        completed = subprocess.run(
            command,
            input=f"{group_name}\n",
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0 or not output_path.exists():
            raise RuntimeError(
                "Failed to export analysis group "
                f"{group_name!r} with GROMACS trjconv: {completed.stderr or completed.stdout}".strip()
            )

    def _run_analysis_modules(self) -> dict[str, Any]:
        self.module_results = {}
        for module_name in self.manifest["modules"]:
            self._stage(f"module:{module_name}", lambda module_name=module_name: self._execute_module_stage(module_name))
        _write_json(self.analysis_dir / "module_status.json", self.module_results)
        return self.module_results

    def _execute_module_stage(self, module_name: str) -> dict[str, Any]:
        handler = getattr(self, f"_run_module_{module_name}")
        try:
            result = handler()
        except Exception as exc:
            result = {"status": "failed", "error": str(exc)}
            self.module_results[module_name] = result
            _write_json(self.analysis_dir / "module_status.json", self.module_results)
            if not self.args.continue_on_error:
                raise
            return result
        self.module_results[module_name] = result
        _write_json(self.analysis_dir / "module_status.json", self.module_results)
        return result

    @property
    def prepared(self) -> dict[str, Any]:
        return self.manifest["prepared_input"]

    def _module_output(self, name: str) -> Path:
        path = self.analysis_dir / name
        path.mkdir(parents=True, exist_ok=True)
        return path

    def _call_module(self, name: str, argv: list[str]) -> dict[str, Any]:
        self.logger.info(f"Running module: {name}")
        module = __import__(f"immunoscope.cli.commands.{name}", fromlist=["main"])
        ret = module.main(argv)
        status = "completed" if ret == 0 else "failed"
        result = {"status": status, "exit_code": ret, "argv": argv}
        if ret != 0:
            raise RuntimeError(f"Module {name} failed with exit code {ret}")
        return result

    def _run_module_quality(self) -> dict[str, Any]:
        output = self._module_output("quality")
        result = self._call_module("quality", [
            "-f", self.prepared["trajectory"],
            "-s", self.prepared["topology"],
            "-o", str(output),
            "--stride", str(self.args.stride),
        ])
        result["root"] = str(output)
        return result

    def _run_module_identity(self) -> dict[str, Any]:
        output = self._module_output("identity")
        result = self._call_module("identity", [
            "--structure", self.prepared["structure"],
            "-o", str(output),
        ])
        result["root"] = str(output)
        return result

    def _run_module_rmsf(self) -> dict[str, Any]:
        output = self._module_output("rmsf")
        result = self._call_module("rmsf", [
            "-f", self.prepared["trajectory"],
            "-s", self.prepared["topology"],
            "--structure", self.prepared["structure"],
            "--annotated",
            "-d", str(output),
            "--stride", str(self.args.stride),
        ])
        result["root"] = str(output)
        return result

    def _run_module_bsa(self) -> dict[str, Any]:
        output = self._module_output("bsa")
        result = self._call_module("bsa", [
            "-f", self.prepared["trajectory"],
            "-s", self.prepared["topology"],
            "--structure", self.prepared["structure"],
            "-o", str(output),
            "--stride", str(self.args.stride),
        ])
        result["root"] = str(output)
        return result

    def _run_module_contact(self) -> dict[str, Any]:
        output = self._module_output("contact")
        context = PipelineContext(
            system_id=self.system_id,
            topology=self.prepared["topology"],
            trajectory_raw=self.prepared["trajectory"],
            trajectory_processed=self.prepared["trajectory"],
            structure_pdb=self.prepared["structure"],
            output_dir=str(output),
        )
        pipeline = ContactFrequencyPipeline(
            cutoff=self.args.contact_cutoff,
            stride=self.args.stride,
            min_frequency=self.args.min_contact_frequency,
        )
        result_context = pipeline.execute(context)
        result_context.save(str(output / "pipeline_context.json"))
        if result_context.has_errors():
            raise RuntimeError("; ".join(result_context.errors))
        contact_report = output / "analysis" / "contacts" / "contact_report.csv"
        return {
            "status": "completed",
            "root": str(output),
            "contact_report": str(contact_report) if contact_report.exists() else None,
            "n_contact_pairs": result_context.results.get("contact_frequency", {}).get("n_contact_pairs"),
            "cutoff_angstrom": self.args.contact_cutoff,
            "min_frequency": self.args.min_contact_frequency,
            "stride": self.args.stride,
        }

    def _run_module_rrcs(self) -> dict[str, Any]:
        output = self._module_output("rrcs")
        result = self._call_module("rrcs", [
            "--structure", self.prepared["structure"],
            "--topology", self.prepared["topology"],
            "--trajectory", self.prepared["trajectory"],
            "-o", str(output),
            "--stride", str(self.args.stride),
            "--pair-scope", "interface",
        ])
        result["root"] = str(output)
        return result

    def _run_module_cluster(self) -> dict[str, Any]:
        output = self._module_output("inter_cluster")
        result = self._call_module("inter_cluster", [
            "-f", self.prepared["trajectory"],
            "-s", self.prepared["topology"],
            "--structure", self.prepared["structure"],
            "-o", str(output),
            "--stride", str(self.args.stride),
        ])
        result["root"] = str(output)
        return result

    def _run_module_landscape(self) -> dict[str, Any]:
        output = self._module_output("landscape")
        args = [
            "--structure", self.prepared["structure"],
            "--topology", self.prepared["topology"],
            "--trajectory", self.prepared["trajectory"],
            "-o", str(output),
            "--stride", str(self.args.stride),
            "--reducer", self.args.landscape_reducer,
        ]
        if self.args.landscape_reducer == "tica":
            args.extend(["--tica-lag", str(self.args.tica_lag)])
        if importlib.util.find_spec("freesasa") is None:
            args.append("--no-bsa")
        result = self._call_module("landscape", args)
        result["root"] = str(output)
        return result

    def _run_module_angles(self) -> dict[str, Any]:
        """Compute TCR-pMHC docking angles (crossing / incident / tilt)."""
        import numpy as np
        from immunoscope.analysis.angles import DockingAngleAnalyzer, DockingAngleInput

        output = self._module_output("angles")
        output.mkdir(parents=True, exist_ok=True)

        angle_input = DockingAngleInput(
            topology=self.prepared["topology"],
            trajectory=self.prepared["trajectory"],
            auto_identify_chains=True,
            use_anarci=True,
            stride=self.args.stride,
            output_dir=str(output),
        )

        try:
            analyzer = DockingAngleAnalyzer()
            result = analyzer.analyze(angle_input)
        except Exception as e:
            return {"status": "failed", "error": str(e), "root": str(output)}

        # Persist time series in the format the report builder expects
        # (Time(ps), Crossing(deg), Incident(deg), [Tilt(deg)])
        csv_path = output / "docking_angles_timeseries.csv"
        try:
            import pandas as pd
            crossing = np.asarray(result.crossing_angles)
            incident = np.asarray(result.incident_angles) if result.incident_angles is not None else None
            n = len(crossing)
            # Build time axis. Stride awareness: we don't have the actual frame
            # interval here, but most MD pipelines step in ps. Use stride * 10ps
            # as a sensible default (matches typical 10ps frame rate).
            dt = 10.0 * float(self.args.stride)
            time_ps = np.arange(n) * dt
            data = {"Time(ps)": time_ps, "Crossing(deg)": crossing}
            if incident is not None:
                data["Incident(deg)"] = incident
            if hasattr(result, "tilt_angles") and result.tilt_angles is not None:
                data["Tilt(deg)"] = np.asarray(result.tilt_angles)
            df = pd.DataFrame(data)
            df.to_csv(csv_path, index=False)
        except Exception as e:
            logging.warning(f"Failed to write docking_angles_timeseries.csv: {e}")
            csv_path = None

        # Persist summary JSON
        def stats(arr):
            if arr is None:
                return None
            arr = np.asarray(arr)
            if arr.size == 0:
                return None
            return {
                "mean": float(np.mean(arr)),
                "std": float(np.std(arr)),
                "min": float(np.min(arr)),
                "max": float(np.max(arr)),
                "n": int(arr.size),
            }

        summary = {
            "n_frames": int(np.asarray(result.crossing_angles).size) if result.crossing_angles is not None else 0,
            "stride": self.args.stride,
            "crossing_angle": stats(result.crossing_angles),
            "incident_angle": stats(result.incident_angles),
        }
        # Tilt may not always be present
        if hasattr(result, "tilt_angles") and result.tilt_angles is not None:
            summary["tilt_angle"] = stats(result.tilt_angles)

        summary_path = output / "docking_angles_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

        return {
            "status": "completed",
            "root": str(output),
            "summary_json": str(summary_path),
            "timeseries_csv": str(csv_path) if csv_path else "",
        }

    def _run_module_dihedrals(self) -> dict[str, Any]:
        """Compute per-residue backbone phi/psi (Ramachandran).

        Delegates to `analysis.geometry.backbone_dihedrals` so the computation
        (correct residue↔angle alignment via `ag1`, circular φ/ψ statistics in
        degrees) lives in one reusable place rather than being duplicated here.
        """
        from immunoscope.analysis.geometry.backbone_dihedrals import (
            compute_backbone_dihedrals, write_backbone_dihedrals,
        )

        output = self._module_output("dihedrals")
        output.mkdir(parents=True, exist_ok=True)
        try:
            rows, summary = compute_backbone_dihedrals(
                self.prepared["topology"], self.prepared["trajectory"],
                stride=self.args.stride,
            )
        except ImportError as e:
            return {"status": "failed", "root": str(output),
                    "error": f"MDAnalysis dihedrals not available: {e}"}
        except Exception as e:  # noqa: BLE001
            return {"status": "failed", "root": str(output),
                    "error": f"Dihedral analysis failed: {e}"}

        write_backbone_dihedrals(output, rows, summary)
        csv_path = output / "residue_dihedrals.csv"

        return {
            "status": "completed",
            "root": str(output),
            "summary_json": str(output / "dihedrals_summary.json"),
            "residue_csv": str(csv_path) if rows else "",
        }

    def _run_module_report(self) -> dict[str, Any]:
        base_dir = self.report_dir / f"interaction_case_{self.system_id}"
        base_dir.mkdir(parents=True, exist_ok=True)
        argv = [
            "interaction",
            "--base-dir", str(base_dir),
            "--system-id", self.system_id,
            "--source-pdb", self.prepared.get("viewer_structure") or self.prepared["structure"],
        ]
        roots = {
            "contact": ("--contact-root", self.analysis_dir / "contact"),
            "bsa": ("--bsa-root", self.analysis_dir / "bsa"),
            "rmsf": ("--rmsf-root", self.analysis_dir / "rmsf"),
            "identity": ("--identity-root", self.analysis_dir / "identity"),
            "rrcs": ("--rrcs-root", self.analysis_dir / "rrcs"),
            "cluster": ("--cluster-root", self.analysis_dir / "inter_cluster"),
            "landscape": ("--landscape-root", self.analysis_dir / "landscape"),
            "angles": ("--angles-root", self.analysis_dir / "angles"),
            "dihedrals": ("--dihedrals-root", self.analysis_dir / "dihedrals"),
        }
        for module_name, (flag, root) in roots.items():
            if module_name in self.manifest["modules"] and root.exists():
                argv.extend([flag, str(root)])
        result = self._call_module("report", argv)
        result["root"] = str(base_dir)
        html_path = base_dir / "overview" / "interaction_report_demo.html"
        result["html"] = str(html_path) if html_path.exists() else None
        return result

    def _write_manifest(self, status: str) -> None:
        self.manifest["status"] = status
        self.manifest["updated_at"] = _now()
        _write_json(self.manifest_path, self.manifest)

    def _write_summary(self, status: str, error: str | None = None) -> None:
        report_result = self.module_results.get("report", {})
        report_html = report_result.get("html")
        if report_html and not Path(report_html).exists():
            report_html = None
        summary = {
            "job_id": self.system_id,
            "status": status,
            "profile": self.args.profile,
            "manifest": str(self.manifest_path),
            "prepared_input": self.manifest.get("prepared_input", {}),
            "modules": self.manifest.get("modules", []),
            "module_results": self.module_results,
            "report_html": report_html,
            "error": error,
            "updated_at": _now(),
        }
        _write_json(self.summary_path, summary)


def handle_single(args: argparse.Namespace, logger: logging.Logger) -> int:
    return SingleRunOrchestrator(args, logger).run()


def handle_validate(args: argparse.Namespace, logger: logging.Logger) -> int:
    return RunValidator(args, logger).run()


def handle_batch(args: argparse.Namespace, logger: logging.Logger) -> int:
    return BatchRunOrchestrator(args, logger).run()


def main(argv=None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    if not args.action:
        parser.print_help()
        return 1
    if args.action == "single":
        return handle_single(args, logger)
    if args.action == "validate":
        return handle_validate(args, logger)
    if args.action == "batch":
        return handle_batch(args, logger)
    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
