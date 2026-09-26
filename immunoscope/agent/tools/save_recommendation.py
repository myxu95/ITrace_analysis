"""Save mutation recommendation tool for Design Copilot mode.

Called by the Agent when it commits to a mutation recommendation during a
Design Copilot conversation. Appends the recommendation to the design task's
draft list so the user sees it in real time on the right-side panel.
"""

from __future__ import annotations

from typing import ClassVar, Optional
from pydantic import BaseModel, Field

from immunoscope.agent.exceptions import ToolValidationError
from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool

# Minimum bar before save_recommendation will accept a submission. The numbers
# are calibrated to be hard to satisfy by an unevidenced LLM hallucination,
# while not blocking a legitimate recommendation backed by even one MD query.
_MIN_RATIONALE_CHARS = 80
_MD_EVIDENCE_TYPES = frozenset({"md_data", "md", "analysis", "rrcs", "occupancy", "rmsf"})


class SaveRecommendationInput(BaseModel):
    """Input schema for committing a mutation recommendation."""

    residue: str = Field(..., description="Residue label, e.g. 'TYR99'")
    chain: str = Field(..., description="Chain identifier, e.g. 'TCR_alpha' or 'alpha'")
    region: str = Field(
        ...,
        description="CDR region or domain, e.g. 'CDR3_alpha', 'CDR1_beta', 'Framework'",
    )
    current_aa: str = Field(..., description="Current amino acid (single letter), e.g. 'Y'")
    suggested_mutations: list[str] = Field(
        ...,
        description="List of suggested mutations as single-letter codes, e.g. ['W', 'F']",
    )
    priority: str = Field(
        ...,
        description="Priority level: 'high' | 'medium' | 'low'",
    )
    confidence: str = Field(
        ...,
        description="Confidence level based on evidence strength: 'high' | 'medium' | 'low'",
    )
    rationale: str = Field(
        ...,
        description=(
            "Why this residue + these mutations. Must cite specific MD metrics "
            "(RRCS, occupancy, RMSF) and ideally PMIDs. 2-4 sentences."
        ),
    )
    expected_effects: dict[str, str] = Field(
        default_factory=dict,
        description=(
            "Expected effects keyed by dimension: "
            "{'affinity': 'increased ~2x', 'specificity': 'maintained', 'stability': 'slight gain'}"
        ),
    )
    risks: list[str] = Field(
        default_factory=list,
        description="List of potential risks or caveats",
    )
    validation_experiments: list[str] = Field(
        default_factory=list,
        description="Suggested validation experiments (SPR, DSF, cell binding, etc.)",
    )
    supporting_evidence: list[dict] = Field(
        default_factory=list,
        description=(
            "List of evidence items, each like "
            "{'type': 'md_data', 'data': 'RRCS=15.19, occupancy=0.85'} or "
            "{'type': 'literature', 'pmid': '12345', 'title': '...'}"
        ),
    )
    notes: Optional[str] = Field(
        default=None,
        description="Optional free-text notes for the user",
    )


@register_tool
class SaveRecommendationTool(Tool):
    """Save a mutation recommendation to the user's draft panel.

    Use this in Design Copilot mode AFTER you have investigated a residue
    with tools (e.g. `query_analysis_results` view=residue or hotspots),
    synthesized the evidence, and decided this is a strong recommendation
    worth committing.

    The saved recommendation appears immediately in the user's draft panel
    on the right side of the workspace. The user can review, discard, or
    finalize the draft list later.

    DO NOT call this tool:
    - Before you have collected actual MD evidence with other tools
    - To pad a list of recommendations toward a target number
    - For residues that conflict with user-stated constraints
    - With fabricated PMIDs (cite only PMIDs from your literature retrieval)
    """

    name: ClassVar[str] = "save_recommendation"
    description: ClassVar[str] = (
        "Save a mutation recommendation to the user's draft panel. "
        "Only use AFTER investigating with query_analysis_results and "
        "synthesizing evidence. The recommendation appears live on the "
        "right side of the user's workspace."
    )
    Input: ClassVar[type[BaseModel]] = SaveRecommendationInput

    # This tool writes to disk (design draft file), not read-only
    is_read_only: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 10.0

    # Design-only: the draft panel only exists in Design Copilot sessions.
    # Exposing this in generic agent mode would let the LLM call into a
    # task store that has no active design_context.
    allowed_modes: ClassVar[frozenset[str]] = frozenset({"design"})

    async def validate(
        self,
        args: SaveRecommendationInput,
        ctx: ToolContext,
    ) -> None:
        """Reject submissions that look LLM-fabricated rather than evidence-backed.

        The system prompt already tells the LLM not to call save_recommendation
        until it has actual MD evidence, but a prompt is a request, not a
        constraint. This validate() turns the request into a hard gate:

        - rationale must be substantive (>= 80 chars; 2-4 sentences as documented)
        - supporting_evidence must contain at least one item
        - at least one evidence item must be MD-derived (type in _MD_EVIDENCE_TYPES)
          OR the rationale must cite a concrete metric (RRCS=, occupancy=, RMSF=)

        Literature-only recommendations are still allowed if the rationale cites
        a numeric MD metric inline — covers the legitimate case where the LLM
        synthesizes prior tool output into the rationale text without re-listing
        it in supporting_evidence.
        """
        rationale = (args.rationale or "").strip()
        if len(rationale) < _MIN_RATIONALE_CHARS:
            raise ToolValidationError(
                f"rationale is too short ({len(rationale)} chars; need at least "
                f"{_MIN_RATIONALE_CHARS}). Cite specific MD metrics (RRCS, occupancy, "
                "RMSF) and explain why this residue + these mutations."
            )

        evidence = args.supporting_evidence or []
        if not evidence:
            raise ToolValidationError(
                "supporting_evidence is empty. Investigate the residue with "
                "query_analysis_results / calculate_rrcs / analyze_hbonds etc. "
                "before calling save_recommendation."
            )

        has_md_evidence = any(
            isinstance(item, dict) and str(item.get("type", "")).lower() in _MD_EVIDENCE_TYPES
            for item in evidence
        )
        rationale_lower = rationale.lower()
        has_inline_metric = any(
            token in rationale_lower for token in ("rrcs", "occupancy", "rmsf", "bsa")
        )
        if not (has_md_evidence or has_inline_metric):
            raise ToolValidationError(
                "No MD-derived evidence found. supporting_evidence must include "
                "at least one item with type in {'md_data', 'analysis', 'rrcs', "
                "'occupancy', 'rmsf'}, OR the rationale must cite a concrete "
                "metric (RRCS, occupancy, RMSF, BSA) by name."
            )

    async def call(
        self,
        args: SaveRecommendationInput,
        ctx: ToolContext,
    ) -> ToolResult:
        """Append the recommendation to the active design task's drafts."""

        # Resolve the design task from session context.
        # We retrieve the design task_id from the session's design_context.
        from immunoscope.web.services.agent_store import get_agent_store

        store = get_agent_store()
        state = await store.get_session(ctx.session_id)
        if not state:
            return ToolResult(
                is_error=True,
                content="Cannot save recommendation: no active session.",
            )

        design_context = getattr(state.session, "design_context", None)
        if not design_context or not design_context.get("task_id"):
            return ToolResult(
                is_error=True,
                content=(
                    "This tool can only be used in Design Copilot sessions. "
                    "There is no design task associated with the current session."
                ),
            )

        task_id = design_context["task_id"]

        # Build the recommendation dict
        recommendation = args.model_dump(exclude_none=True)

        # Append to drafts file via design router helper
        from immunoscope.web.routers.design import append_draft
        try:
            n_drafts = append_draft(task_id, recommendation)
        except Exception as e:
            return ToolResult(
                is_error=True,
                content=f"Failed to save recommendation: {e}",
            )

        # Emit a custom event so the frontend can refresh its draft panel
        if ctx.on_event:
            try:
                await ctx.on_event({
                    "type": "draft_saved",
                    "task_id": task_id,
                    "n_drafts": n_drafts,
                    "residue": args.residue,
                    "chain": args.chain,
                    "priority": args.priority,
                })
            except Exception:
                pass  # Event emission is best-effort

        return ToolResult(
            is_error=False,
            content=(
                f"Saved recommendation #{n_drafts}: {args.residue} ({args.chain} "
                f"{args.region}) → {'/'.join(args.suggested_mutations)} "
                f"[priority={args.priority}, confidence={args.confidence}]. "
                f"The user can now see it in their draft panel."
            ),
        )

    @classmethod
    def system_prompt_section(cls) -> str:
        """Section appended to the system prompt explaining this tool."""
        return (
            "## save_recommendation\n"
            "When in Design Copilot mode and you decide a residue + mutation "
            "is a STRONG, evidence-backed recommendation, call this tool to "
            "commit it to the user's draft panel. Required: residue, chain, "
            "region, current_aa, suggested_mutations, priority, confidence, "
            "rationale (with MD metrics + PMIDs), expected_effects, risks, "
            "validation_experiments. Do NOT use this tool to dump candidates; "
            "use it deliberately for picks you're confident about."
        )
