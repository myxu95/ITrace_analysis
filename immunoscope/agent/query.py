"""
Lightweight single-query entry point for the Agent engine.

`run_query()` is a simplified, constrained flavor of the full `run_turn`
loop, designed for the Assistant drawer inside the report UI:

  - Single-turn (default max_turns=2: one LLM + optional tool + one LLM summary)
  - Tool whitelist (caller opts in to exactly which tools are allowed)
  - No hooks, no permissions (not interactive — this runs non-interactively
    against read-only query tools by default)
  - Returns a structured QueryResult with the answer text, evidence
    extracted from tool calls, and cost/token info

This shares the full Agent's LLM client, prompt layer, and tool registry.
The intent is that there is exactly ONE LLM client, ONE prompt, ONE tool
system for the project — Assistant and Agent chat both route through the
same components, just with different constraints.

Design rationale: rather than threading a "lightweight mode" flag through
the full run_turn, we have a separate coroutine that composes the same
building blocks (llm.stream, execute_tools) with a tighter loop. Keeps
both paths easy to read and prevents Assistant behavior from accidentally
drifting inside the bigger Agent loop.
"""

from __future__ import annotations

import asyncio
import json
import logging
import time
import uuid
from contextlib import suppress
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

import anyio

from immunoscope.agent.config import Settings, get_settings
from immunoscope.agent.context import build_system_prompt
from immunoscope.agent.cost_tracker import CostTracker
from immunoscope.agent.engine import _run_stream, execute_tools
from immunoscope.agent.exceptions import LLMStreamError
from immunoscope.agent.llm import make_llm_client
from immunoscope.agent.session import Session
from immunoscope.agent.stream_events import MessageStop, ToolUseEnd
from immunoscope.agent.tool import ToolContext
from immunoscope.agent.tools import TOOL_REGISTRY

log = logging.getLogger("immunoscope.agent.query")


# ---------------------------------------------------------------------- #
# Result types
# ---------------------------------------------------------------------- #

@dataclass
class EvidenceItem:
    """One piece of structured evidence extracted from a tool call.

    An EvidenceItem is meant to be renderable as a UI card next to the
    Assistant's answer. It references WHICH tool produced the data and
    the raw content returned (as a string — presenter Markdown is common).
    """

    tool: str                         # e.g. "query_analysis_results"
    view: Optional[str] = None        # e.g. "hotspots"
    content: str = ""                 # raw tool output (usually Markdown)
    is_error: bool = False


@dataclass
class QueryResult:
    """Structured response returned by run_query().

    Consumed by the /assistant/query Web endpoint and the Assistant drawer
    in the report UI.
    """

    answer: str                            # final LLM text (user-facing)
    evidence: list[EvidenceItem] = field(default_factory=list)
    tool_calls: list[dict] = field(default_factory=list)   # {name, input} records
    turns: int = 0
    cost_usd: float = 0.0
    input_tokens: int = 0
    output_tokens: int = 0
    truncated: bool = False                # True if hit max_turns
    error: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "answer": self.answer,
            "evidence": [
                {"tool": e.tool, "view": e.view, "content": e.content, "is_error": e.is_error}
                for e in self.evidence
            ],
            "tool_calls": list(self.tool_calls),
            "turns": self.turns,
            "cost_usd": round(self.cost_usd, 6),
            "tokens": {
                "input": self.input_tokens,
                "output": self.output_tokens,
            },
            "truncated": self.truncated,
            "error": self.error,
        }


# ---------------------------------------------------------------------- #
# Helpers
# ---------------------------------------------------------------------- #

def _filter_tools(
    whitelist: Optional[list[str]],
) -> list[dict]:
    """Build the tool schema list for the LLM, honoring an optional whitelist.

    An empty/None whitelist means "no tools" (pure text Q&A).
    A non-empty whitelist means "only these tools, in this order".
    """
    if whitelist is None:
        return []
    schemas: list[dict] = []
    for name in whitelist:
        tool = TOOL_REGISTRY.get(name)
        if tool is None:
            log.warning("run_query: whitelisted tool %r not registered — skipping", name)
            continue
        schemas.append(tool.to_anthropic_schema())
    return schemas


def _build_query_system_prompt(
    session: Session,
    case_hint: Optional[dict] = None,
) -> str:
    """Build a system prompt for single-turn query mode.

    Starts from the full Agent system prompt (so domain knowledge is
    identical), then appends a short directive that keeps the assistant
    terse and evidence-linked. Also injects the case hint if provided.
    """
    base = build_system_prompt(session)

    extras = [
        "",
        "# Query Mode",
        "",
        "You are answering a single factual question about a specific analysis case.",
        "Use the available tools to fetch the needed data, then produce a concise,",
        "evidence-linked answer:",
        "",
        "- Prefer calling a tool once to get the facts, then summarize.",
        "- Cite concrete numbers from tool output whenever possible.",
        "- Keep the answer short (2-5 sentences or a small table).",
        "- Do not speculate or suggest multi-step reasoning — if the question needs",
        "  deeper analysis, say so and recommend the user open the full Agent chat.",
    ]

    if case_hint:
        extras.append("")
        extras.append("# Current Case Context")
        for key, value in case_hint.items():
            if value is None or value == "":
                continue
            extras.append(f"- {key}: {value}")

    return base + "\n".join(extras)


def _extract_evidence(
    tool_calls: list[ToolUseEnd],
    tool_results: list[dict],
) -> list[EvidenceItem]:
    """Pair up tool_uses with their tool_results into evidence items."""
    # Build id → result content map.
    result_by_id: dict[str, dict] = {r["tool_use_id"]: r for r in tool_results}
    out: list[EvidenceItem] = []

    for tu in tool_calls:
        r = result_by_id.get(tu.id)
        if r is None:
            continue

        content = r.get("content", "")
        if isinstance(content, list):
            # Tool returned blocks — concatenate text blocks.
            content = "\n".join(
                blk.get("text", "") if isinstance(blk, dict) else str(blk)
                for blk in content
            )
        elif not isinstance(content, str):
            content = str(content)

        view = None
        if isinstance(tu.input, dict):
            view = tu.input.get("view")

        out.append(EvidenceItem(
            tool=tu.name,
            view=view,
            content=content,
            is_error=bool(r.get("is_error")),
        ))

    return out


# ---------------------------------------------------------------------- #
# Main entry point
# ---------------------------------------------------------------------- #

async def run_query(
    question: str,
    *,
    tool_whitelist: Optional[list[str]] = None,
    max_turns: int = 2,
    case_hint: Optional[dict] = None,
    settings: Optional[Settings] = None,
    turn_timeout: float = 60.0,
) -> QueryResult:
    """
    Run a single-turn constrained Agent query.

    This is the lightweight entry point used by the Assistant drawer and
    the /assistant/query Web endpoint. It reuses the full Agent's LLM
    client, prompt layer, and tool registry, but applies strict
    constraints:

      - At most `max_turns` LLM calls (default 2 — one to call a tool,
        one to summarize).
      - Only tools named in `tool_whitelist` are advertised to the LLM.
      - No hooks, no permissions, no abort flow — this is meant to be
        invoked synchronously from an HTTP handler.

    Args:
        question: The user's question.
        tool_whitelist: Tool names the LLM may call. Default (None) uses
            ["query_analysis_results"], which is the canonical read-only
            query tool. Pass [] to disable tool use entirely.
        max_turns: Hard cap on LLM calls. Default 2.
        case_hint: Optional dict injected into the system prompt to tell
            the LLM which case is in scope (e.g. {"case_dir": "...",
            "job_id": "...", "system_id": "..."}).
        settings: Agent settings; defaults to get_settings().
        turn_timeout: Per-LLM-call timeout in seconds.

    Returns:
        QueryResult with answer, evidence, tool_calls, cost, and tokens.
        On error, .error is set and .answer is a user-safe fallback
        message.
    """
    settings = settings or get_settings()

    # Default whitelist: the canonical read-only query tool.
    if tool_whitelist is None:
        tool_whitelist = ["query_analysis_results"]

    tools_schema = _filter_tools(tool_whitelist)

    # Ephemeral session: run_query is stateless; we don't persist message
    # history. This session exists only to satisfy build_system_prompt()
    # and to thread through _run_stream.
    session = Session(id=f"query-{uuid.uuid4().hex[:8]}")
    session.cost_tracker = CostTracker(model=settings.LLM_MODEL)
    session.messages.append({"role": "user", "content": question})

    system_prompt = _build_query_system_prompt(session, case_hint=case_hint)

    # Minimal ToolContext: no audit, no permissions, no on_event streaming
    # (this entry point is non-streaming by design).
    abort_event = anyio.Event()

    async def _noop_event(_ev: dict) -> None:  # noqa: ARG001
        return None

    ctx = ToolContext(
        session_id=session.id,
        abort_event=abort_event,
        on_event=_noop_event,
        settings=settings,
        db_path=getattr(settings, "AGENT_AUDIT_DB_PATH", ""),
    )

    llm = make_llm_client(settings)

    total_input_tokens = 0
    total_output_tokens = 0
    total_cost = 0.0
    all_tool_calls_record: list[dict] = []
    all_evidence: list[EvidenceItem] = []
    answer_text = ""
    error_msg: Optional[str] = None
    truncated = False

    try:
        for turn_idx in range(max_turns):
            session.turn_idx = turn_idx
            tool_use_blocks: list[ToolUseEnd] = []

            # Run one LLM stream, collecting tool_uses and the final content.
            try:
                stop_ev: MessageStop | None = await asyncio.wait_for(
                    _run_stream(
                        llm,
                        system=system_prompt,
                        messages=session.messages,
                        tools=tools_schema,
                        tool_use_blocks=tool_use_blocks,
                        on_event=_noop_event,
                    ),
                    timeout=turn_timeout,
                )
            except asyncio.TimeoutError:
                error_msg = f"LLM stream timed out after {turn_timeout}s"
                break
            except LLMStreamError as e:
                error_msg = f"LLM stream error: {e}"
                break

            # Accumulate token / cost usage.
            usage = stop_ev.usage if stop_ev else None
            turn_cost, _total = session.cost_tracker.add_turn(usage)
            total_cost += turn_cost
            if usage:
                total_input_tokens += int(usage.get("input_tokens") or 0)
                total_output_tokens += int(usage.get("output_tokens") or 0)

            # Capture the assistant content (so subsequent turns can reference
            # the tool_use blocks if we go another round).
            assistant_content = (stop_ev.content if stop_ev else None) or []
            if not assistant_content:
                break
            session.messages.append({"role": "assistant", "content": assistant_content})

            # If no tools were requested, extract text and finish.
            if not tool_use_blocks:
                answer_text = _extract_text_from_content(assistant_content)
                break

            # Record tool calls for the result.
            for tu in tool_use_blocks:
                all_tool_calls_record.append({"name": tu.name, "input": tu.input})

            # Execute tools, build user message with tool_result blocks.
            tool_results = await execute_tools(tool_use_blocks, ctx, turn_idx)
            all_evidence.extend(_extract_evidence(tool_use_blocks, tool_results))

            session.messages.append({"role": "user", "content": list(tool_results)})
        else:
            # max_turns exhausted without a plain-text answer
            truncated = True
            # Best-effort: pull last text from the final assistant message.
            if session.messages and session.messages[-1].get("role") == "assistant":
                answer_text = _extract_text_from_content(
                    session.messages[-1].get("content") or []
                )

    except Exception as e:
        log.exception("run_query failed: %s", e)
        error_msg = str(e)

    finally:
        with suppress(Exception):
            if hasattr(llm, "aclose"):
                await llm.aclose()

    # User-facing fallback if something went wrong and we have no answer.
    if not answer_text and error_msg:
        answer_text = (
            "I couldn't retrieve an answer for that question. "
            "Try again, or open the Agent for deeper analysis."
        )

    return QueryResult(
        answer=answer_text.strip(),
        evidence=all_evidence,
        tool_calls=all_tool_calls_record,
        turns=session.turn_idx + 1,
        cost_usd=total_cost,
        input_tokens=total_input_tokens,
        output_tokens=total_output_tokens,
        truncated=truncated,
        error=error_msg,
    )


def _extract_text_from_content(content: list[Any]) -> str:
    """Concatenate text blocks from an assistant content list."""
    parts: list[str] = []
    for blk in content or []:
        if isinstance(blk, dict):
            if blk.get("type") == "text":
                parts.append(blk.get("text", ""))
        elif hasattr(blk, "type") and getattr(blk, "type", "") == "text":
            parts.append(getattr(blk, "text", "") or "")
    return "\n".join(p for p in parts if p)


__all__ = [
    "run_query",
    "QueryResult",
    "EvidenceItem",
]