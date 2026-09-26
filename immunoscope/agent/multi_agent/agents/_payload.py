"""Shared helpers for building per-agent input payloads.

Each round-1 agent builds a payload from the design intent + the
pre-rendered system preamble (overview + quality) the orchestrator
attaches to `MutationDesignContext`. Round-1 agents do not see each
other's `AgentMessage`s, so they only use `base_intent_header`.

The round-2 RecommendationAgent *does* consume the four round-1 briefs;
`format_round1_briefs` renders them (narrative + structured_payload)
for injection into that agent's payload.
"""

from __future__ import annotations

import json

from ..context import MutationDesignContext


def base_intent_header(ctx: MutationDesignContext) -> str:
    """The case_id / case_dir / design_intent block used by every agent.

    When `ctx.system_preamble` is set, the pre-rendered system overview
    + trajectory-quality views are appended so every parallel agent
    shares a baseline picture without each calling those views.
    """
    sections = [
        f"**case_id**: `{ctx.case_id}`",
        f"**case_dir** (pass this to `query_analysis_results`): "
        f"`{ctx.case_dir}`",
        "",
        "## Design intent (from the user)",
        "",
        ctx.design_intent,
    ]
    if ctx.system_preamble:
        sections += [
            "",
            "## System preamble (pre-rendered overview + quality)",
            "",
            ctx.system_preamble.strip(),
        ]
    return "\n".join(sections) + "\n"


# The four round-1 reader senders, so brief formatting never picks up the
# round-2 recommendation or the round-2.5 critique that also live in history.
_ROUND1_SENDERS = (
    "bio_agent",
    "interaction_reader",
    "conformation_reader",
    "interface_exposure_reader",
)


def _render_message(msg) -> str:
    block = [f"### {msg.sender}  (msg_id: `{msg.msg_id}`)", "", msg.narrative.strip()]
    if msg.structured_payload:
        block += [
            "",
            "structured_payload:",
            "```json",
            json.dumps(msg.structured_payload, ensure_ascii=False, indent=2),
            "```",
        ]
    return "\n".join(block)


def format_round1_briefs(ctx: MutationDesignContext) -> str:
    """Render the four round-1 reader briefs for a downstream payload.

    Filters to the round-1 reader senders so it is stable even after the
    round-2 recommendation and round-2.5 critique have been appended to
    `message_history`. Each brief carries its msg_id for `MessageRef` citation.
    """
    briefs = [m for m in ctx.message_history if m.sender in _ROUND1_SENDERS]
    if not briefs:
        return "(no round-1 briefs available — round-1 produced no messages)"
    return "\n\n".join(_render_message(m) for m in briefs)


def _latest_from(ctx: MutationDesignContext, sender: str):
    msgs = [m for m in ctx.message_history if m.sender == sender]
    return msgs[-1] if msgs else None


def format_recommendation(ctx: MutationDesignContext) -> str:
    """Render the most recent RecommendationAgent message (for the critic /
    a revision pass). Returns a marker if none exists yet."""
    msg = _latest_from(ctx, "recommendation_agent")
    return _render_message(msg) if msg else "(no recommendation produced yet)"


def format_critique(ctx: MutationDesignContext) -> str:
    """Render the most recent DesignCritic message (for a revision pass)."""
    msg = _latest_from(ctx, "design_critic")
    return _render_message(msg) if msg else "(no critique available)"


def has_critique(ctx: MutationDesignContext) -> bool:
    return _latest_from(ctx, "design_critic") is not None


__all__ = [
    "base_intent_header",
    "format_round1_briefs",
    "format_recommendation",
    "format_critique",
    "has_critique",
]
