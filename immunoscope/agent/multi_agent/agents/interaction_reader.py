"""Interaction Reader — contact / hotspot analysis (v1, 2026-05-28).

Reads the residue / pair / region views via `query_analysis_results` and
surfaces which residues drive the TCR-pHLA interface. Runs in parallel
with the other three round-1 agents; does not see peer messages.
"""

from __future__ import annotations

from ..base import BaseAgent
from ..context import MutationDesignContext
from ..knowledge import load_knowledge
from ._payload import base_intent_header


class InteractionReader(BaseAgent):
    """Hotspot / pair / fingerprint reader."""

    name = "interaction_reader"
    enabled_tools = ("query_analysis_results",)
    max_turns = 14
    max_tokens_per_turn = 3072

    def build_system_prompt(self, ctx: MutationDesignContext) -> str:
        return load_knowledge("interaction_reader")

    def build_input_payload(self, ctx: MutationDesignContext) -> str:
        sections = [
            "# Multi-agent run — Interaction Reader stage",
            "",
            base_intent_header(ctx),
            "",
            "## Your task",
            "",
            "Use `query_analysis_results` against `case_dir` to read the "
            "`hotspots`, `pair`, and `fingerprint` views. Then call "
            "`emit_message` exactly once with the ranked hotspots and "
            "evidence described in your system prompt. You run in parallel "
            "with the other round-1 agents; no peer messages are available.",
        ]
        return "\n".join(sections)


__all__ = ["InteractionReader"]
