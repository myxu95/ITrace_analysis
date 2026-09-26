"""Visualize a residue / mutation site as an inline 3D viewer in the chat.

This tool is for ad-hoc structural inspection inside Design Copilot — when
the user (or the Agent itself) wants to bring up a 3D view of a residue
without committing it as a saved recommendation. The frontend listens for
the `visualize_residue` event and renders a 3D viewer card inline in the
conversation stream.
"""

from __future__ import annotations

from typing import ClassVar, Optional
from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class VisualizeResidueInput(BaseModel):
    """Input schema for visualizing a residue in the chat."""

    residue: str = Field(
        ...,
        description="Residue label, e.g. 'TYR99'",
    )
    chain: str = Field(
        default="",
        description="Chain identifier, e.g. 'TCR_alpha' or 'alpha'",
    )
    region: str = Field(
        default="",
        description="CDR region, e.g. 'CDR3_alpha', 'CDR1_beta', 'Framework'",
    )
    suggested_mutations: list[str] = Field(
        default_factory=list,
        description=(
            "Optional: amino acid codes the residue might be mutated to. "
            "If provided, the viewer card shows the proposed substitutions."
        ),
    )
    contact_partners: list[dict] = Field(
        default_factory=list,
        description=(
            "Optional list of contact partners to highlight, each like "
            "{'residue': 'MET4', 'chain': 'peptide'}"
        ),
    )
    note: str = Field(
        default="",
        description="Optional one-sentence commentary shown beneath the viewer",
    )


@register_tool
class VisualizeResidueTool(Tool):
    """Pop up a 3D structural view of a residue in the conversation.

    Use this when the user asks to *see* or *visualize* a residue, or when
    you want to give them a structural intuition before discussing a
    candidate. The view is rendered inline in the chat using NGL.js with:

    - Protein backbone as gray cartoon
    - Target residue highlighted in red (stick + transparent sphere)
    - Contact partners highlighted in orange
    - Hydrogen bonds shown as yellow dashes
    - Camera auto-centered on the residue

    Unlike `save_recommendation`, this tool does NOT commit anything to the
    user's draft panel. It's purely for structural inspection / discussion.

    Good times to use it:
    - "show me TYR99 in 3D"
    - "let me look at the contact environment around this residue"
    - After you identify a strong hotspot but before the user has decided
      to save it as a recommendation
    """

    name: ClassVar[str] = "visualize_residue"
    description: ClassVar[str] = (
        "Show a residue's 3D structural context inline in the chat. "
        "Use when the user wants to *see* a residue or its contact "
        "environment, not when committing a recommendation."
    )
    Input: ClassVar[type[BaseModel]] = VisualizeResidueInput

    is_read_only: ClassVar[bool] = True
    timeout_seconds: ClassVar[float] = 10.0

    # Design-only: the residue viewer is part of the Design Copilot UI.
    allowed_modes: ClassVar[frozenset[str]] = frozenset({"design"})

    async def call(
        self,
        args: VisualizeResidueInput,
        ctx: ToolContext,
    ) -> ToolResult:
        """Emit a `visualize_residue` event so the frontend renders a viewer."""
        from immunoscope.web.services.agent_store import get_agent_store

        store = get_agent_store()
        state = await store.get_session(ctx.session_id)
        if not state:
            return ToolResult(
                is_error=True,
                content="No active session.",
            )

        # Pull design task ID + source job ID from session's design_context
        design_context = getattr(state.session, "design_context", None) or {}
        task_id = design_context.get("task_id")
        if not task_id:
            return ToolResult(
                is_error=True,
                content=(
                    "This tool is available in Design Copilot sessions only. "
                    "No design task is associated with the current session."
                ),
            )

        # Build payload for the frontend
        payload = args.model_dump(exclude_none=True)
        payload["task_id"] = task_id

        # Emit event to frontend WebSocket
        if ctx.on_event:
            try:
                await ctx.on_event({
                    "type": "visualize_residue",
                    **payload,
                })
            except Exception:
                pass  # best-effort

        return ToolResult(
            is_error=False,
            content=(
                f"Displayed 3D view of {args.residue}"
                + (f" ({args.region})" if args.region else "")
                + (f" with {len(args.contact_partners)} contact partner(s)" if args.contact_partners else "")
                + " inline in the chat. The user can now see the structural context."
            ),
        )

    @classmethod
    def system_prompt_section(cls) -> str:
        return (
            "## visualize_residue\n"
            "Use to pop up an inline 3D viewer of a residue in the chat. "
            "Good after you discover a strong hotspot and want the user to "
            "*see* it before discussing further. Required: residue label "
            "(e.g. 'TYR99'). Optional: chain, region, suggested_mutations "
            "list, contact_partners list, and a short note explaining what "
            "to look at. The viewer auto-loads the source job's PDB, "
            "highlights the residue in red, contact partners in orange, "
            "and shows H-bonds in yellow."
        )
