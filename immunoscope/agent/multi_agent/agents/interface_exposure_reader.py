"""Interface & Exposure Reader — interface size + per-residue SASA (v1, 2026-05-28).

Reads the `interface` and `exposure` views via `query_analysis_results`
and surfaces:

  * The chain-level interface footprint (BSA mean / std / composition).
  * Per-residue burial-state distribution.
  * **Cavity candidates** — `interface_core` residues whose
    `relative_exposure` is anomalously high, suggesting a packing
    defect at the buried face that a bulkier / more complementary
    residue could fill.
  * Load-bearing burial events (top |delta_sasa|) for orientation.

Runs in parallel with the other three round-1 readers; does not see
their messages.
"""

from __future__ import annotations

from ..base import BaseAgent
from ..context import MutationDesignContext
from ..knowledge import load_knowledge
from ._payload import base_intent_header


class InterfaceExposureReader(BaseAgent):
    """Interface BSA + per-residue SASA / cavity-candidate reader."""

    name = "interface_exposure_reader"
    enabled_tools = ("query_analysis_results",)
    max_turns = 14
    max_tokens_per_turn = 3072

    def build_system_prompt(self, ctx: MutationDesignContext) -> str:
        return load_knowledge("interface_exposure_reader")

    def build_input_payload(self, ctx: MutationDesignContext) -> str:
        sections = [
            "# Multi-agent run — Interface & Exposure Reader stage",
            "",
            base_intent_header(ctx),
            "",
            "## Your task",
            "",
            "Use `query_analysis_results` against `case_dir` to read the "
            "`interface` view (BSA + composition) and the `exposure` view "
            "(per-residue SASA + cavity candidates). Focus on residues "
            "classified as `interface_core` whose `relative_exposure` is "
            "anomalously high — these are packing defects that a bulkier "
            "or chemically complementary residue could fill. Then call "
            "`emit_message` exactly once. You run in parallel with the "
            "other round-1 agents, so no peer messages are available.",
        ]
        return "\n".join(sections)


__all__ = ["InterfaceExposureReader"]
