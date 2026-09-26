"""Show a comparison summary card inline in the chat.

Lightweight tool for when the user wants to *see* an existing comparison
job without re-running it. Emits a `show_comparison_card` event with the
key plots + takeaways so the frontend can render a compact diff card.
"""

from __future__ import annotations

from typing import ClassVar
from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class ShowComparisonCardInput(BaseModel):
    """Input schema for displaying a comparison card."""

    compare_job_id: str = Field(
        ...,
        description="ImmunoScope compare job ID (the result of a previous compare_systems call).",
    )
    focus: str = Field(
        default="overview",
        description=(
            "Which aspect to highlight: 'overview' | 'quality' | 'flexibility' "
            "| 'interface' | 'hotspots'. Default 'overview' renders the takeaways "
            "+ 4 key plots."
        ),
    )
    note: str = Field(
        default="",
        description="Optional commentary to display alongside the card.",
    )


@register_tool
class ShowComparisonCardTool(Tool):
    """Render an inline comparison card in the chat for an existing compare job."""

    name: ClassVar[str] = "show_comparison_card"
    description: ClassVar[str] = (
        "Render a compact comparison card inline in the chat for an existing "
        "compare job. Shows key takeaways and plots without re-running the "
        "comparison. Use after compare_systems completes, or when the user "
        "refers to a previously-created compare job by ID."
    )
    Input: ClassVar[type[BaseModel]] = ShowComparisonCardInput

    is_read_only: ClassVar[bool] = True
    timeout_seconds: ClassVar[float] = 10.0

    # Design-only: the comparison card UI lives inside Design Copilot.
    allowed_modes: ClassVar[frozenset[str]] = frozenset({"design"})

    async def call(self, args: ShowComparisonCardInput, ctx: ToolContext) -> ToolResult:
        from immunoscope.web.services.job_store import job_store

        job = await job_store.get(args.compare_job_id)
        if not job:
            return ToolResult(
                is_error=True,
                content=f"Compare job not found: {args.compare_job_id}",
            )

        # Check it's actually a compare job
        modules = (job.config.modules if job.config else []) or []
        if "compare" not in modules:
            return ToolResult(
                is_error=True,
                content=f"Job {args.compare_job_id} is not a compare job (modules: {modules}).",
            )

        # Emit event for frontend to render
        if ctx.on_event:
            try:
                await ctx.on_event({
                    "type": "show_comparison_card",
                    "compare_job_id": args.compare_job_id,
                    "name": job.name,
                    "focus": args.focus,
                    "note": args.note,
                })
            except Exception:
                pass

        return ToolResult(
            is_error=False,
            content=(
                f"Displayed inline comparison card for '{job.name}' "
                f"(focus: {args.focus}). The user can now see the diff summary in the chat. "
                f"Full detail at #/jobs/{args.compare_job_id}."
            ),
        )

    @classmethod
    def system_prompt_section(cls) -> str:
        return (
            "## show_comparison_card\n"
            "After compare_systems completes (or when the user refers to a "
            "previous compare job), call this to render an inline diff card "
            "in the chat: takeaways + key plots + a link to the full detail "
            "page. Lightweight — does not re-run the comparison."
        )
