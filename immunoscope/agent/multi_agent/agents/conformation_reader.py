"""Conformation Reader — backbone / sidechain conformation reader (v1, 2026-05-28).

Reads the `flexibility` and `dihedrals` views via `query_analysis_results`
and surfaces which residues are rigid, which are mobile, and how the
backbone partitions across Ramachandran basins. Runs in parallel with
the other three round-1 agents, so it does not depend on any prior
agent's message — every claim is grounded in the views it reads itself.

Clustering (FES / TICA) is intentionally not on this reader's beat in
v1 — it belongs to the deferred 'iteration mode' that re-models and
re-runs MD on prior recommendations.
"""

from __future__ import annotations

from ..base import BaseAgent
from ..context import MutationDesignContext
from ..knowledge import load_knowledge
from ._payload import base_intent_header


class ConformationReader(BaseAgent):
    """Backbone (Ramachandran) and per-residue flexibility (RMSF) reader."""

    name = "conformation_reader"
    enabled_tools = ("query_analysis_results",)
    max_turns = 14
    max_tokens_per_turn = 3072

    def build_system_prompt(self, ctx: MutationDesignContext) -> str:
        return load_knowledge("conformation_reader")

    def build_input_payload(self, ctx: MutationDesignContext) -> str:
        sections = [
            "# Multi-agent run — Conformation Reader stage",
            "",
            base_intent_header(ctx),
            "",
            "## Your task",
            "",
            "Use `query_analysis_results` against `case_dir` to read the "
            "`flexibility` and `dihedrals` views. Identify rigid vs. "
            "mobile residues, flag positions sampling multiple Ramachandran "
            "basins, and call out conformational-risk residues that a "
            "mutation designer should avoid (or deliberately target). "
            "Then call `emit_message` exactly once — no peer messages "
            "are available because you run in parallel with the other "
            "round-1 agents.",
        ]
        return "\n".join(sections)


__all__ = ["ConformationReader"]
