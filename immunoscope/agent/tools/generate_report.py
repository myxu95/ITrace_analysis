"""Tool for querying analysis results using the Agent's run_query entry point."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class GenerateReportInput(BaseModel):
    """Input schema for querying analysis results."""

    case_dir: str = Field(
        description="Path to the case directory containing analysis results"
    )
    question: str = Field(
        description=(
            "Natural language question about the analysis results. Examples: "
            "'What are the key interaction hotspots?', "
            "'Which residues have high flexibility?', "
            "'Show me the buried surface area breakdown'"
        )
    )


@register_tool
class GenerateReportTool(Tool):
    """Query analysis results and generate answers using the Agent's query system.

    This tool uses the same LLM client, prompt, and tool registry as the full
    Agent chat, but runs in constrained single-turn mode. It's designed for
    quick factual queries about a specific analysis case.

    For multi-system comparisons or deeper reasoning, recommend the user open
    the full Agent chat instead.
    """

    name: ClassVar[str] = "generate_report"
    description: ClassVar[str] = (
        "Query analysis results for a specific case directory. "
        "Answers factual questions about hotspots, flexibility, interactions, "
        "buried surface area, and other analysis outputs. "
        "Use this after running analyses to interpret results."
    )
    Input: ClassVar[type[BaseModel]] = GenerateReportInput

    is_read_only: ClassVar[bool] = True
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 120.0

    async def call(self, args: GenerateReportInput, ctx: ToolContext) -> ToolResult:
        """Execute query."""
        case_path = Path(args.case_dir)
        if not case_path.exists():
            return ToolResult(
                content=f"Error: Case directory not found: {args.case_dir}",
                is_error=True,
            )

        try:
            from immunoscope.agent.query import run_query

            query_result = await run_query(
                args.question,
                tool_whitelist=["query_analysis_results"],
                max_turns=2,
                case_hint={
                    "case_dir": str(case_path),
                    "system_id": case_path.name,
                },
            )

            # Format output for the Agent.
            output_lines = []
            output_lines.append(f"Question: {args.question}")
            output_lines.append("")
            output_lines.append("Answer:")
            output_lines.append(query_result.answer)

            if query_result.evidence:
                output_lines.append("")
                output_lines.append("Evidence:")
                for i, ev in enumerate(query_result.evidence, 1):
                    # Truncate long evidence for readability.
                    content = ev.content[:300] + "..." if len(ev.content) > 300 else ev.content
                    output_lines.append(f"  {i}. [{ev.tool}] {content}")

            if query_result.error:
                output_lines.append("")
                output_lines.append(f"Note: {query_result.error}")

            if query_result.truncated:
                output_lines.append("")
                output_lines.append("(Query was truncated due to turn limit)")

            output_text = "\n".join(output_lines)

            return ToolResult(
                content=output_text,
                is_error=bool(query_result.error),
            )

        except Exception as e:
            return ToolResult(
                content=f"Error: Query failed: {e}",
                is_error=True,
            )
