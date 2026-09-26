"""Design Critic — round-2.5 grounded adversarial reviewer (v1, 2026-06-02).

After the RecommendationAgent (round 2) produces its site picks, this agent
reviews them. It is NOT a second opinion that re-does the design — it is an
adversary whose job is to find material flaws: recommendations not supported
by the briefs, hard-constraint violations (a germline-conserved framework
position or an HLA anchor recommended for mutation without a flag), strong
signals the recommender ignored, conflicts left unresolved, or
misclassifications.

It emits a critique `AgentMessage` whose `structured_payload.verdict` is
`approve` or `revise`. On `revise`, the orchestrator runs the
RecommendationAgent once more in revision mode with this critique injected.

Pure reviewer: `enabled_tools = ()`. Everything it needs (the four briefs,
the candidate table, the recommendation) is in the injected payload; it
grounds every finding with a `MessageRef` back to a round-1 brief or the
recommendation.
"""

from __future__ import annotations

from ..base import BaseAgent
from ..context import MutationDesignContext
from ..knowledge import load_knowledge
from ._payload import (
    base_intent_header,
    format_recommendation,
    format_round1_briefs,
)


class DesignCriticAgent(BaseAgent):
    """Round-2.5 adversarial reviewer of the recommendation.

    May spot-check a doubted quantitative claim against the MD views via
    `query_analysis_results`, but is budgeted (max_tool_calls) so it verifies
    a handful of specific numbers rather than re-deriving the analysis.
    """

    name = "design_critic"
    enabled_tools = ("query_analysis_results",)
    # The critic reviews ALL four round-1 briefs + the recommendation, so it
    # needs at least the readers' budget (14). With max_turns=8 it sometimes
    # exhausted the loop right after the "emit now" nudge fires (which lands
    # only once max_tool_calls is hit), leaving too few turns to finalize and
    # raising AgentTimeoutError -> the whole run is marked not_succeeded even
    # though the recommendation is fine. max_tool_calls still caps how much it
    # spot-checks; this only widens the room to emit after the nudge.
    max_turns = 14
    # Spot-check budget: a few targeted verifications, not a full re-read.
    max_tool_calls = 6
    max_tokens_per_turn = 4096

    def build_system_prompt(self, ctx: MutationDesignContext) -> str:
        return load_knowledge("design_critic")

    def build_input_payload(self, ctx: MutationDesignContext) -> str:
        candidate_table = (
            ctx.recommendation_preamble
            or "(orchestrator produced no candidate baseline table)"
        )
        sections = [
            "# Multi-agent run — Design Critic stage (round 2.5)",
            "",
            base_intent_header(ctx),
            "",
            "## Round-1 briefs (the evidence the recommendation must rest on)",
            "",
            "Cite a brief in a finding with a `MessageRef` whose `msg_id` "
            "matches the heading.",
            "",
            format_round1_briefs(ctx),
            "",
            "## Deterministic candidate table (what the recommender started from)",
            "",
            candidate_table,
            "",
            "## The recommendation under review",
            "",
            format_recommendation(ctx),
            "",
            "## Deterministic numeric-fidelity flags (machine check)",
            "",
            "A model-free check flagged these numbers in the round-1 briefs as "
            "NOT tracing to any tool output that reader received — i.e. likely "
            "fabricated. If the recommendation relies on a flagged number, that "
            "is a `blocking` `unsupported_claim` finding.",
            "",
            (ctx.numeric_audit_summary
             or "(none — all round-1 numbers traced to tool output)"),
            "",
            "## Deterministic constraint / coverage findings (machine check)",
            "",
            "Model-free checks already confirmed these against the structured "
            "payloads. A `blocking` entry WILL force a revision regardless of "
            "your verdict — fold it into your findings so the reviser sees one "
            "coherent critique.",
            "",
            (ctx.deterministic_findings_summary
             or "(none — no hard-constraint or coverage violations detected)"),
            "",
            "## Your task",
            "",
            "Adversarially review the recommendation above against the briefs "
            "and candidate table. Find material flaws (do NOT re-design, do "
            "NOT flag mere preferences). You MAY call `query_analysis_results` "
            "a few times to verify a specific suspicious number (pass the "
            "`case_dir` from the header). Cite each flaw's brief via a "
            "`MessageRef`. Then call `emit_message` exactly once with the "
            "critique `structured_payload`. Set `verdict='revise'` only if "
            "there is ≥1 blocking finding or ≥2 major findings; otherwise "
            "`verdict='approve'` (advisory findings may still be listed).",
        ]
        return "\n".join(sections)


__all__ = ["DesignCriticAgent"]
