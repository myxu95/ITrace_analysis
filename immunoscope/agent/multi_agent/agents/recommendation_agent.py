"""Recommendation Agent — round-2 integration / site picker (v1, 2026-05-28).

Round 1 fans out four parallel readers (Bio / Interaction / Conformation /
Interface & Exposure) that each emit one `AgentMessage`. This agent is
round 2: a single integrator that consumes all four briefs plus a
deterministically pre-aggregated candidate table (built by the
orchestrator from the four `structured_payload`s) and produces ONE
recommendation `AgentMessage`.

Scope boundary (deliberate):
  - It recommends *which site* to mutate and *which direction* to push
    (a `design_hint_class` + constraints), NOT a specific amino acid.
    Choosing `V→F` is deferred to a future round-3 PhysicoChemicalAgent,
    for which `handoff_to_physchem` is the interface slot.
  - It is a pure integrator: `enabled_tools = ()`. All evidence is in
    the injected payload; it grounds claims with `MessageRef`s pointing
    back at the round-1 msg_ids.

It does not run any MD/analysis tool. Reading a residue's structural
neighbourhood is a future skill (not wired here).
"""

from __future__ import annotations

from ..base import BaseAgent
from ..context import MutationDesignContext
from ..knowledge import load_knowledge
from ._payload import (
    base_intent_header,
    format_critique,
    format_recommendation,
    format_round1_briefs,
    has_critique,
)


class RecommendationAgent(BaseAgent):
    """Round-2 integrator: folds the four round-1 briefs into site picks."""

    name = "recommendation_agent"
    enabled_tools = ()  # pure integrator — evidence comes from the payload
    # No tools, so every turn is a structured-emit attempt (not tool looping).
    # The revision pass (when the critic asks for changes) is harder to get
    # schema-valid on the first try than the initial pass, and max_turns=6
    # occasionally exhausted retries -> AgentTimeoutError -> the run is marked
    # not_succeeded despite sound candidates. 10 gives emit headroom; healthy
    # cases still emit on turn 1-2, so this only helps the hard ones.
    max_turns = 10
    max_tokens_per_turn = 4096

    def build_system_prompt(self, ctx: MutationDesignContext) -> str:
        return load_knowledge("recommendation_agent")

    def build_input_payload(self, ctx: MutationDesignContext) -> str:
        if has_critique(ctx):
            return self._build_revision_payload(ctx)
        candidate_table = (
            ctx.recommendation_preamble
            or "(orchestrator produced no candidate baseline table)"
        )
        sections = [
            "# Multi-agent run — Recommendation Agent stage (round 2)",
            "",
            base_intent_header(ctx),
            "",
            "## Round-1 briefs (the four parallel readers)",
            "",
            "Each section below is one reader's `AgentMessage`. Cite a "
            "brief in your recommendations with a `MessageRef` whose "
            "`msg_id` matches the heading.",
            "",
            format_round1_briefs(ctx),
            "",
            "## Pre-aggregated candidate baseline (deterministic)",
            "",
            "The orchestrator extracted these candidate sites from the "
            "briefs' structured_payloads by fixed rules. Treat it as a "
            "starting point: you MAY add sites the briefs support but the "
            "table missed, and you MAY drop sites whose evidence is thin. "
            "Each row carries the source msg_id(s) you should cite.",
            "",
            candidate_table,
            "",
            "## Your task",
            "",
            "Work in three phases, then finalize:",
            "",
            "  A. Review the candidate baseline against the four briefs — "
            "add or drop sites with justification.",
            "  B. For each surviving candidate, weigh the supporting vs. "
            "opposing evidence across the briefs and record any conflicts "
            "(e.g. one reader says preserve, another says mobile).",
            "  C. Call `emit_message` exactly once with the recommendation "
            "`structured_payload` described in your system prompt.",
            "",
            "Stay at site + direction granularity — do NOT name a specific "
            "amino acid substitution. No analysis tools are available; "
            "ground every claim in the briefs via `MessageRef`.",
        ]
        return "\n".join(sections)

    def _build_revision_payload(self, ctx: MutationDesignContext) -> str:
        """Revision pass: the critic flagged issues with the prior pick set."""
        candidate_table = (
            ctx.recommendation_preamble
            or "(orchestrator produced no candidate baseline table)"
        )
        sections = [
            "# Multi-agent run — Recommendation Agent REVISION (round 2)",
            "",
            base_intent_header(ctx),
            "",
            "A design critic reviewed your previous recommendation and asked "
            "for a revision. Produce a NEW, corrected recommendation that "
            "addresses every blocking/major finding while keeping the picks "
            "that were sound.",
            "",
            "## Round-1 briefs (the evidence)",
            "",
            format_round1_briefs(ctx),
            "",
            "## Pre-aggregated candidate baseline",
            "",
            candidate_table,
            "",
            "## Your previous recommendation",
            "",
            format_recommendation(ctx),
            "",
            "## Critic findings to address",
            "",
            format_critique(ctx),
            "",
            "## Machine-confirmed findings (deterministic — MUST be fixed)",
            "",
            "These were verified by model-free checks, not opinion. Any "
            "`blocking` entry is non-negotiable: a germline-conserved framework "
            "position must not be recommended for mutation; a flagged "
            "missed-signal residue must be recommended or explicitly skipped "
            "with a reason.",
            "",
            (ctx.deterministic_findings_summary or "(none)"),
            "",
            "## Your task",
            "",
            "Emit ONE revised recommendation via `emit_message`, using the same "
            "`structured_payload` schema. For each critic finding: fix it "
            "(drop/repriortise/reclassify the site, add the missing flag, "
            "surface the conflict, or move it to `skipped_candidates` with a "
            "reason). Do not name specific amino acids; cite briefs via "
            "`MessageRef`.",
        ]
        return "\n".join(sections)


__all__ = ["RecommendationAgent"]
