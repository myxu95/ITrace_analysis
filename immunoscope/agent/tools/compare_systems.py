"""Tool for comparing two MD simulation systems by job ID.

Wraps SystemComparisonPipeline. Accepts ImmunoScope job IDs (preferred) or
raw case directory paths. On completion emits a `comparison_ready` event
so the frontend can show a compact "comparison ready" card in the chat.
"""

from __future__ import annotations

import asyncio
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class CompareSystemsInput(BaseModel):
    """Input schema for system comparison.

    Either provide job IDs (job_a_id, job_b_id) OR raw paths (case_a, case_b).
    Job IDs are preferred when the user is referring to ImmunoScope jobs.
    """

    job_a_id: str = Field(
        default="",
        description=(
            "ImmunoScope job ID for Case A (the reference system). "
            "If provided, case_a is auto-resolved to that job's analysis directory."
        ),
    )
    job_b_id: str = Field(
        default="",
        description="ImmunoScope job ID for Case B (the comparison system).",
    )
    case_a: str = Field(
        default="",
        description="Fallback: raw path to Case A directory. Only used if job_a_id is empty.",
    )
    case_b: str = Field(
        default="",
        description="Fallback: raw path to Case B directory. Only used if job_b_id is empty.",
    )
    label_a: str = Field(
        default="Condition A",
        description="Display label for Case A",
    )
    label_b: str = Field(
        default="Condition B",
        description="Display label for Case B",
    )
    comparison_mode: str = Field(
        default="generic",
        description="One of: 'generic' | 'mutation' | 'sampling' | 'replicate'",
    )
    comparison_scope: str = Field(
        default="same-system",
        description="One of: 'same-system' | 'cross-system' | 'auto'",
    )
    alignment_selection: str = Field(
        default="phla_core_ca",
        description="Selection used for structure alignment",
    )
    comparison_context: str = Field(
        default="",
        description="Free-text note describing what is being compared and why",
    )


@register_tool
class CompareSystemsTool(Tool):
    """Run a structured comparison between two analyzed MD jobs.

    The tool creates a new ImmunoScope compare job (so it appears in the
    Jobs list and gets a proper detail page), runs the comparison pipeline,
    and emits a `comparison_ready` event with a link to the detail view.
    """

    name: ClassVar[str] = "compare_systems"
    description: ClassVar[str] = (
        "Compare two completed ImmunoScope analysis jobs. Generates "
        "structured comparison artifacts (quality, flexibility, RRCS, "
        "interactions, FEL). Prefer job_a_id / job_b_id over raw paths."
    )
    Input: ClassVar[type[BaseModel]] = CompareSystemsInput

    is_read_only: ClassVar[bool] = False
    is_concurrency_safe: ClassVar[bool] = False
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 1200.0  # 20 minutes

    async def call(self, args: CompareSystemsInput, ctx: ToolContext) -> ToolResult:
        # Resolve case directories
        case_a_path, case_b_path = await self._resolve_paths(args)
        if isinstance(case_a_path, str):
            return ToolResult(is_error=True, content=case_a_path)  # error string
        if isinstance(case_b_path, str):
            return ToolResult(is_error=True, content=case_b_path)

        # Validate mode / scope
        if args.comparison_mode not in ("generic", "mutation", "sampling", "replicate"):
            return ToolResult(is_error=True, content=f"Invalid comparison_mode: {args.comparison_mode}")
        if args.comparison_scope not in ("auto", "same-system", "cross-system"):
            return ToolResult(is_error=True, content=f"Invalid comparison_scope: {args.comparison_scope}")

        # Create a new compare job so the result has a stable home
        from immunoscope.web.config import settings as web_settings
        from immunoscope.web.models import (
            AnalysisConfig, AnalysisProfile, InputMode, JobDetail,
            JobInputPaths, JobStatus, JobSummary,
        )
        from immunoscope.web.services.job_store import job_store

        job_id = uuid.uuid4().hex[:12]
        job_dir = web_settings.JOBS_DIR / job_id
        run_dir = job_dir / "run"
        run_dir.mkdir(parents=True, exist_ok=True)

        name = f"{args.label_a}_vs_{args.label_b}"
        job_summary = JobSummary(
            id=job_id,
            name=name,
            status=JobStatus.RUNNING,
            config=AnalysisConfig(
                profile=AnalysisProfile.SAMPLING_COMPARE,
                modules=["compare"],
                notes=args.comparison_context,
                input_mode=InputMode.PREPARED,
                auto_start=True,
            ),
            input_files=[],
            created_at=datetime.now(timezone.utc),
            progress=0,
            job_dir=str(job_dir),
            run_dir=str(run_dir),
        )
        detail = JobDetail(
            **job_summary.model_dump(),
            inputs=JobInputPaths(
                case_a=str(case_a_path),
                case_b=str(case_b_path),
                label_a=args.label_a,
                label_b=args.label_b,
                comparison_mode=args.comparison_mode,
                comparison_scope=args.comparison_scope,
                alignment_selection=args.alignment_selection,
                residue_mapping="auto",
            ),
        )
        await job_store.save(detail)

        if ctx.on_event:
            try:
                await ctx.on_event({
                    "type": "progress",
                    "message": f"Comparison job {job_id} started: {args.label_a} vs {args.label_b}",
                })
            except Exception:
                pass

        # Run the pipeline in a thread
        try:
            await asyncio.to_thread(
                self._run_pipeline,
                case_a_path, case_b_path, run_dir, args,
            )
        except Exception as e:
            await job_store.update(job_id, status=JobStatus.FAILED, error=str(e))
            return ToolResult(is_error=True, content=f"Comparison pipeline failed: {e}")

        # Mark complete + emit ready event
        await job_store.update(job_id, status=JobStatus.COMPLETED, progress=100)
        if ctx.on_event:
            try:
                await ctx.on_event({
                    "type": "comparison_ready",
                    "compare_job_id": job_id,
                    "name": name,
                    "label_a": args.label_a,
                    "label_b": args.label_b,
                })
            except Exception:
                pass

        return ToolResult(
            is_error=False,
            content=(
                f"Comparison complete. Created compare job '{name}' (ID: {job_id}). "
                f"Artifacts written to {run_dir}/analysis/comparison/. "
                f"The user can open the detail view at #/jobs/{job_id}."
            ),
        )

    async def _resolve_paths(self, args: CompareSystemsInput):
        """Resolve case directories from job IDs or raw paths.

        Returns (case_a_path, case_b_path) on success, or
        (error_str, None) on error.
        """
        from immunoscope.web.services.job_store import job_store

        async def resolve(job_id: str, raw_path: str, label: str):
            if job_id:
                job = await job_store.get(job_id)
                if not job:
                    return f"Job not found: {job_id} (for {label})"
                if not job.run_dir:
                    return f"Job {job_id} has no run_dir (for {label})"
                p = Path(job.run_dir)
                if not p.exists():
                    return f"Job {job_id}'s run_dir does not exist: {p}"
                return p
            if raw_path:
                p = Path(raw_path).expanduser().resolve()
                if not p.exists():
                    return f"Path does not exist: {raw_path} (for {label})"
                return p
            return f"Neither job_id nor path provided for {label}"

        a = await resolve(args.job_a_id, args.case_a, "Case A")
        b = await resolve(args.job_b_id, args.case_b, "Case B")
        if isinstance(a, str):
            return a, None
        if isinstance(b, str):
            return b, None
        return a, b

    def _run_pipeline(self, case_a_path, case_b_path, run_dir, args):
        """Blocking comparison execution. Runs in a thread."""
        from immunoscope.core import PipelineContext
        from immunoscope.pipeline.analysis_pipelines import SystemComparisonPipeline

        context = PipelineContext(
            system_id=f"{args.label_a}_vs_{args.label_b}",
            topology="",
            trajectory_raw="",
            output_dir=str(run_dir),
        )
        pipeline = SystemComparisonPipeline(
            case_a_root=str(case_a_path),
            case_b_root=str(case_b_path),
            label_a=args.label_a,
            label_b=args.label_b,
            comparison_mode=args.comparison_mode,
            comparison_context=args.comparison_context,
            comparison_scope=args.comparison_scope,
            alignment_selection=args.alignment_selection,
            residue_mapping="auto",
        )
        result_context = pipeline.execute(context)
        if result_context.has_errors():
            raise RuntimeError("; ".join(result_context.errors))

    @classmethod
    def system_prompt_section(cls) -> str:
        return (
            "## compare_systems\n"
            "Compare two ImmunoScope analysis jobs. Prefer passing job IDs "
            "(job_a_id, job_b_id) when the user references existing jobs in "
            "the system. The tool creates a new compare job, runs the full "
            "pipeline (5-15 min), and emits a `comparison_ready` event with "
            "a link to the detail page. Use 'mutation' mode for WT vs mutant, "
            "'sampling' for standard vs enhanced sampling."
        )
