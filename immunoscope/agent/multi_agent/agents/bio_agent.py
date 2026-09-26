"""Bio Agent — biological context provider (v1, 2026-05-28).

Runs first in the pipeline. Reads only the design intent (no MD data)
and consults the literature RAG to produce a biological briefing for
the three readers and the Writer.

Outputs one broadcast `AgentMessage`. The downstream agents cite this
message via `[ref:<bio_msg_id>]` / `MessageRef` when their reasoning
leans on the biology surfaced here.
"""

from __future__ import annotations

from typing import Optional

from immunoscope.agent.llm import LLMClient

from ..base import BaseAgent
from ..context import MutationDesignContext
from ..knowledge import load_knowledge
from ..tools import (
    SearchFactBlocksTool,
    SearchLiteratureTool,
    TcrConservationTool,
)


class BioAgent(BaseAgent):
    """Biological context provider — literature RAG + facts + germline conservation.

    Has no MD trajectory tools; its three extra_tools are all sequence /
    knowledge grounded: `search_literature` (abstract-level curated RAG),
    `search_fact_blocks` (precise, citable biochemical facts mined from
    full text), and `query_tcr_conservation` (per-residue IMGT germline
    conservation).
    """

    name = "bio_agent"
    enabled_tools = ()  # uses extra_tools instead
    max_turns = 14
    # Bio has three (relatively heavy) retrieval tools and a prompt that
    # invites several calls each; cap tool calls lower so it is nudged to
    # finalize before exhausting its turn budget.
    max_tool_calls = 8
    # Bio emits the largest payload (prior_engineering entries + quotes +
    # conservation block + many LiteratureRefs). A 2048 cap truncated the
    # emit_message tool call mid-generation, so it never finalized and timed
    # out. Give it room to emit in one turn.
    max_tokens_per_turn = 4096

    def __init__(
        self,
        llm: LLMClient,
        *,
        literature_db_path: Optional[str] = None,
    ) -> None:
        super().__init__(
            llm,
            extra_tools={
                "search_literature": SearchLiteratureTool(literature_db_path),
                "search_fact_blocks": SearchFactBlocksTool(),
                "query_tcr_conservation": TcrConservationTool(),
            },
        )

    def build_system_prompt(self, ctx: MutationDesignContext) -> str:
        return load_knowledge("bio_agent")

    def build_input_payload(self, ctx: MutationDesignContext) -> str:
        return (
            "# Multi-agent run — Bio Agent stage\n"
            "\n"
            f"**case_id**: `{ctx.case_id}`\n"
            f"**case_dir**: `{ctx.case_dir}`\n"
            "\n"
            "## Design intent (from the user)\n"
            "\n"
            f"{ctx.design_intent}\n"
            "\n"
            "## Your task\n"
            "\n"
            "1. Call `search_literature` 1–3 times for background on the\n"
            "   peptide / HLA / TCR combination implied by the design intent.\n"
            "2. Call `search_fact_blocks` 1–3 times for precise, citable\n"
            "   precedents (mutation→effect outcomes, binding measurements)\n"
            "   to populate `prior_engineering` — each entry should trace to\n"
            "   a returned fact block's PMID + verbatim quote.\n"
            "3. Call `query_tcr_conservation` once with this exact\n"
            f"   `case_dir`: `{ctx.case_dir}` to get the per-residue germline\n"
            "   conservation profile (which framework positions are\n"
            "   load-bearing, which CDR1/CDR2 positions are variable).\n"
            "Then call `emit_message` exactly once with the biological brief\n"
            "described in your system prompt. Stop when emit_message is\n"
            "accepted — do not call any tool after that.\n"
        )


__all__ = ["BioAgent"]
