"""
Analysis results query tool for ImmunoScope Agent.

Provides structured access to MD analysis results with intelligent filtering
and formatting to minimize token usage while maximizing information value.
"""

from pathlib import Path
from typing import ClassVar, Optional, Dict, Any
from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool
from immunoscope.agent.presenters import (
    CaseLocator,
    PresenterError,
    SPATIAL_LAYERS,
    SUB_FLAVORS,
    VIEW_REGISTRY,
    hierarchy_summary,
    render,
    views_for_layer,
    views_for_sub_flavor,
)


class QueryAnalysisResultsInput(BaseModel):
    """Input for querying analysis results."""

    case_dir: str = Field(
        ...,
        description="Path to analysis case directory (e.g., output/5c0a_run2_full_analysis)"
    )

    view: Optional[str] = Field(
        default=None,
        description=(
            "View to render: 'overview', 'hotspots', 'interface', 'flexibility', "
            "'quality', 'clustering', 'residue', 'pair', 'fingerprint', "
            "'angles', 'dihedrals', 'interface_comparison'. "
            "Use 'overview' first to get system summary. The legacy "
            "'candidates' view (design_priority_score) has been removed per "
            "D-B5 (2026-05-26); for mutation-target ranking call Pure-RRCS-"
            "rank as a baseline or query 'hotspots' (RRCS-only)."
        )
    )

    spatial_layer: Optional[str] = Field(
        default=None,
        description=(
            "D-B1 (2026-05-26) spatial-hierarchy filter. One of "
            "'complex', 'interface', 'region', 'residue', 'pair'. When "
            "set, the response is restricted to views whose primary "
            "spatial_layer matches. Use this when you want to make an "
            "explicit ablation cut (e.g. 'give me only Pair-layer "
            "evidence'). Optional; default behavior is unchanged."
        )
    )

    sub_flavor: Optional[str] = Field(
        default=None,
        description=(
            "D-B1 (2026-05-26) sub-flavor filter: 'static' or 'dynamic'. "
            "When set together with `spatial_layer`, returns only views "
            "that serve that sub-flavor at that layer. Optional."
        )
    )

    query_type: Optional[str] = Field(
        default=None,
        description=(
            "DEPRECATED: Use 'view' instead. Legacy aliases: "
            "'summary'→'overview', 'interactions'→'hotspots'"
        )
    )

    filters: Optional[Dict[str, Any]] = Field(
        default=None,
        description=(
            "Optional filters: "
            "{'region': 'CDR3', 'min_rrcs': 3.0, 'min_occupancy': 0.7, "
            "'top_n': 10, 'residue': 'ASP92', 'residue1': 'ASP92', "
            "'residue2': 'LYS66', 'case_b_dir': 'output/case_b'}"
        )
    )


@register_tool
class QueryAnalysisResultsTool(Tool):
    """
    Query MD analysis results with standardized Markdown views.

    This tool provides structured access to analysis results via standardized
    views optimized for LLM consumption. Each view returns ≤5000 chars of
    formatted Markdown with consistent field naming and provenance tracking.

    Available views:
    - overview: System TL;DR with metadata, module status, and key findings
    - hotspots: Top RRCS contact pairs, ranked and grouped by partner
    - interface: BSA statistics and interface composition
    - flexibility: RMSF by region, highlighting flexible residues
    - quality: Trajectory quality metrics and convergence
    - clustering: Conformational states and transitions
    - residue: Single-residue drill-down (requires filter: residue)
    - pair: Single-pair dynamics fingerprint (requires: residue1, residue2)
    - fingerprint: Interface occupancy distribution and binding mode
    - angles: TCR-pMHC docking geometry (crossing / incident / tilt)
    - dihedrals: Ramachandran (φ/ψ) per residue + backbone flexibility ranking
    - interface_comparison: Residue-level RRCS diff between two cases
      (requires filter: case_b_dir)

    Workflow:
    1. Start with 'overview' to understand system and available modules
    2. Query specific views based on analysis goals
    3. Use filters to narrow results (region, min_rrcs, top_n, etc.)
    4. Use drill-down views (residue, pair) for detailed investigation

    NOTE (D-B5, 2026-05-26): the legacy `candidates` view, which ranked
    residues by a composite `design_priority_score`, has been deleted
    from the LLM-facing surface. The Pure-RRCS-rank baseline (Methods
    §6.2) replaces the mono-signal-heuristic role; the LLM does its own
    integration over the spatial-hierarchy evidence.

    NOTE (D-B1, 2026-05-26): each registered view declares a
    `spatial_layer` and `sub_flavors_served` on its presenter class
    (Complex / Interface / Region / Residue / Pair × static / dynamic).
    `query_analysis_results` accepts an optional `spatial_layer` and
    `sub_flavor` filter so callers can request, e.g., only Pair-layer
    dynamic evidence; this also drives §6.5.2 ablation cuts.
    """

    name: ClassVar[str] = "query_analysis_results"
    description: ClassVar[str] = (
        "Query MD analysis results via standardized Markdown views. "
        "Returns formatted, LLM-optimized summaries (≤5000 chars). "
        "Start with view='overview' to get system summary, then query specific views."
    )
    Input: ClassVar[type[BaseModel]] = QueryAnalysisResultsInput

    is_read_only: ClassVar[bool] = True
    timeout_seconds: ClassVar[float] = 30.0

    # Legacy query_type → view mapping
    QUERY_TYPE_ALIASES = {
        "summary": "overview",
        "interactions": "hotspots",
    }

    async def call(
        self,
        args: QueryAnalysisResultsInput,
        ctx: ToolContext
    ) -> ToolResult:
        """Query analysis results."""

        case_dir = Path(args.case_dir)

        # Validate case directory
        if not case_dir.exists():
            return ToolResult(
                is_error=True,
                content=f"Case directory not found: {case_dir}"
            )

        # D-B1 (2026-05-26): validate spatial_layer / sub_flavor up front.
        layer = (args.spatial_layer or "").lower() or None
        flavor = (args.sub_flavor or "").lower() or None
        if layer and layer not in SPATIAL_LAYERS:
            return ToolResult(
                is_error=True,
                content=(
                    f"Invalid spatial_layer={layer!r}. "
                    f"Must be one of {SPATIAL_LAYERS}."
                ),
            )
        if flavor and flavor not in SUB_FLAVORS:
            return ToolResult(
                is_error=True,
                content=(
                    f"Invalid sub_flavor={flavor!r}. "
                    f"Must be one of {SUB_FLAVORS}."
                ),
            )

        # Determine view (handle legacy query_type)
        view = args.view
        if not view and args.query_type:
            # Map legacy query_type to view
            view = self.QUERY_TYPE_ALIASES.get(args.query_type, args.query_type)

        # D-B5 (2026-05-26): the `candidates` view was removed from the
        # LLM-facing surface. Fail fast with an actionable message so the
        # caller can pick a still-supported view instead of getting a
        # generic "view not found" error.
        if view == "candidates":
            return ToolResult(
                is_error=True,
                content=(
                    "The 'candidates' view was removed per D-B5 (2026-05-26): "
                    "the composite design_priority_score is no longer "
                    "surfaced to the LLM, since pre-ranking by a hand-tuned "
                    "heuristic would short-circuit the agent's reasoning. "
                    "Use 'hotspots' for pure-RRCS-ranked contact pairs, "
                    "'residue' for single-residue drill-down, or filter via "
                    "spatial_layer='residue' to discover other residue-layer "
                    "views. For benchmark comparison, Pure-RRCS-rank is a "
                    "separate baseline (Methods §6.2)."
                ),
            )

        if not view:
            # D-B1 discovery path: if the caller passed only spatial_layer /
            # sub_flavor, return the matching view directory instead of an
            # error — this is the documented way to ask "what's available at
            # the Pair layer?".
            if layer or flavor:
                if layer and flavor:
                    matches = views_for_sub_flavor(flavor, layer=layer)
                elif layer:
                    matches = views_for_layer(layer)
                else:
                    matches = views_for_sub_flavor(flavor)
                summary = hierarchy_summary()
                lines = ["# Views matching your D-B1 filter\n"]
                if not matches:
                    lines.append("_No registered view matches this filter._")
                else:
                    lines.append(
                        "| View | spatial_layer | sub_flavors |\n"
                        "|------|---------------|-------------|"
                    )
                    for name in matches:
                        info = summary[name]
                        lines.append(
                            f"| {name} | {info['spatial_layer']} | "
                            f"{', '.join(info['sub_flavors'])} |"
                        )
                lines.append(
                    "\nPass one of these as `view=` (plus any filters) to "
                    "render it."
                )
                return ToolResult(is_error=False, content="\n".join(lines))
            return ToolResult(
                is_error=True,
                content=(
                    "Missing required parameter: 'view'. "
                    f"Available views: {', '.join(sorted(VIEW_REGISTRY.keys()))}. "
                    "Alternatively pass `spatial_layer=` (and optionally "
                    "`sub_flavor=`) to discover matching views per D-B1."
                )
            )

        # D-B1: if a layer / flavor filter was supplied alongside a view,
        # verify the requested view actually serves that layer/flavor.
        if view in VIEW_REGISTRY and (layer or flavor):
            cls = VIEW_REGISTRY[view]
            if layer and cls.spatial_layer != layer:
                return ToolResult(
                    is_error=True,
                    content=(
                        f"View '{view}' belongs to spatial_layer="
                        f"'{cls.spatial_layer}', not '{layer}'. "
                        f"Views at '{layer}': "
                        f"{', '.join(views_for_layer(layer)) or '(none)'}."
                    ),
                )
            if flavor and flavor not in cls.sub_flavors_served:
                return ToolResult(
                    is_error=True,
                    content=(
                        f"View '{view}' does not serve sub_flavor "
                        f"'{flavor}' (serves: "
                        f"{', '.join(cls.sub_flavors_served)}). "
                        f"Views serving '{flavor}'"
                        + (f" at '{layer}'" if layer else "")
                        + f": {', '.join(views_for_sub_flavor(flavor, layer=layer)) or '(none)'}."
                    ),
                )

        # Create locator
        try:
            locator = CaseLocator(case_dir)
        except Exception as e:
            return ToolResult(
                is_error=True,
                content=f"Failed to initialize case locator: {e}"
            )

        # Render view
        try:
            markdown = render(view, locator, filters=args.filters or {})

            return ToolResult(
                is_error=False,
                content=markdown
            )

        except PresenterError as e:
            return ToolResult(
                is_error=True,
                content=str(e)
            )

        except Exception as e:
            return ToolResult(
                is_error=True,
                content=f"Failed to render view '{view}': {e}"
            )
