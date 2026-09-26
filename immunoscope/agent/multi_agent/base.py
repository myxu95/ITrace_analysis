"""BaseAgent skeleton for the multi-agent mutation-design pipeline.

A concrete agent (Bio / Interaction Reader / Dynamics Reader / Geometry
Reader / Writer) subclasses `BaseAgent` and supplies four things:

  - `name`           : sender id used in `AgentMessage.sender` and audit
                       trails ("bio_agent", "interaction_reader", ...).
  - `system_prompt`  : prompt body. May be a class attribute or computed
                       by overriding `build_system_prompt(ctx)`.
  - `enabled_tools`  : list of tool names from `agent.tools.TOOL_REGISTRY`
                       this agent is allowed to call. The Writer needs
                       none; the readers need `query_analysis_results`;
                       Bio needs a literature-RAG retriever (added in a
                       later patch).
  - `build_input_payload(ctx)`: returns the user-role content that the
                       LLM sees on turn 1 — typically the design intent,
                       the relevant prior agent messages (bio brief +
                       any peer reader narratives), and the per-reader
                       view excerpts.

The base class drives a bounded tool loop:

    1. Build the system prompt + input payload.
    2. Stream from the LLM, dispatching tool calls back into the registry
       and feeding their results back in.
    3. Loop ends when the LLM calls the special `emit_message` tool
       (`finalize=True`) — its arguments populate an `AgentMessage` which
       is validated and appended to `ctx.message_history`.
    4. If `max_turns` runs out before `emit_message` fires, raise
       `AgentTimeoutError`. The orchestrator decides whether to retry.

`emit_message` is exposed as an in-process pseudo-tool — it never
runs analysis, it just captures the agent's final structured output.
We model it as a tool (rather than a parsing step on free-form text)
because the Anthropic / OpenAI tool-use surface enforces schema
compliance for us, which is exactly what we want for a contract
that downstream validators (`validate_message`) depend on.

v1 keeps this synchronous + sequential. v1.1 will add concurrent
reader execution by making `BaseAgent.run` a coroutine over a shared
read-only context snapshot.
"""

from __future__ import annotations

import json
import logging
import time
import uuid
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

from immunoscope.agent.llm import LLMClient
from immunoscope.agent.stream_events import (
    MessageStop,
    TextDelta,
    ToolUseEnd,
    ToolUseStart,
)
from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import TOOL_REGISTRY

from .context import MutationDesignContext, SkillCallRecord
from .messaging import (
    AgentMessage,
    EvidenceRef,
    MessageValidationError,
    evidence_ref_from_dict,
    validate_message,
)

log = logging.getLogger("immunoscope.agent.multi_agent")


class AgentTimeoutError(RuntimeError):
    """Raised when an agent does not emit a final message within max_turns."""


class AgentToolError(RuntimeError):
    """Raised when an enabled tool is missing from TOOL_REGISTRY."""


# ---------------------------------------------------------------------------
# `emit_message` pseudo-tool schema
# ---------------------------------------------------------------------------

_EMIT_MESSAGE_SCHEMA: Dict[str, Any] = {
    "name": "emit_message",
    "description": (
        "Finalize this agent's contribution. Call EXACTLY ONCE at the end "
        "of your reasoning. Every non-common-knowledge claim in `narrative` "
        "MUST be backed by an inline `[ref:<id>]` token whose id appears in "
        "`evidence_refs`. Do not call any analysis tool after `emit_message`."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "narrative": {
                "type": "string",
                "description": (
                    "The natural-language briefing or analysis you are "
                    "handing to downstream agents. Cite evidence inline as "
                    "[ref:<id>]."
                ),
            },
            "evidence_refs": {
                "type": "array",
                "minItems": 1,
                "description": (
                    "Tagged-union list of evidence references. Every "
                    "[ref:<id>] in `narrative` must match a ref_id here."
                ),
                "items": {
                    "type": "object",
                    "oneOf": [
                        {
                            "properties": {
                                "ref_id": {"type": "string"},
                                "kind": {"const": "context_field"},
                                "path": {"type": "string"},
                                "quote": {"type": "string"},
                            },
                            "required": ["ref_id", "kind", "path", "quote"],
                        },
                        {
                            "properties": {
                                "ref_id": {"type": "string"},
                                "kind": {"const": "literature"},
                                "pmid": {"type": "string"},
                                "title": {"type": "string"},
                                "relevant_quote": {"type": "string"},
                            },
                            "required": [
                                "ref_id", "kind", "pmid", "title", "relevant_quote",
                            ],
                        },
                        {
                            "properties": {
                                "ref_id": {"type": "string"},
                                "kind": {"const": "skill_call"},
                                "skill_name": {"type": "string"},
                                "call_id": {"type": "string"},
                                "result_summary": {"type": "string"},
                            },
                            "required": [
                                "ref_id", "kind", "skill_name", "call_id",
                                "result_summary",
                            ],
                        },
                        {
                            "properties": {
                                "ref_id": {"type": "string"},
                                "kind": {"const": "message"},
                                "msg_id": {"type": "string"},
                                "excerpt": {"type": "string"},
                            },
                            "required": ["ref_id", "kind", "msg_id", "excerpt"],
                        },
                    ],
                },
            },
            "structured_payload": {
                "type": ["object", "null"],
                "description": (
                    "Optional dict for downstream numerical consumers. The "
                    "Writer fills this with the 5-field recommendation JSON; "
                    "readers usually leave it null."
                ),
            },
            "recipients": {
                "type": ["array", "string"],
                "description": (
                    "Either the literal string 'broadcast' or a list of agent "
                    "names that should consume this message."
                ),
            },
        },
        "required": ["narrative", "evidence_refs"],
    },
}


# ---------------------------------------------------------------------------
# Per-turn LLM call state
# ---------------------------------------------------------------------------


@dataclass
class _PendingMessage:
    """Populated when the LLM invokes `emit_message`."""

    narrative: str
    evidence_refs: List[EvidenceRef]
    structured_payload: Optional[Dict[str, Any]] = None
    recipients: Any = "broadcast"


@dataclass
class AgentRunResult:
    """What `BaseAgent.run` returns to the orchestrator."""

    msg_id: str
    message: AgentMessage
    turns_used: int
    tool_calls: List[SkillCallRecord] = field(default_factory=list)
    # Post-hoc numeric-fidelity report (None when the agent has no tools to
    # check against). `.ok is False` means a decimal in the message did not
    # trace to any tool output the agent received.
    numeric_audit: Optional[Any] = None


# ---------------------------------------------------------------------------
# Base agent
# ---------------------------------------------------------------------------


class BaseAgent(ABC):
    """Skeleton for one agent in the multi-agent pipeline.

    Subclass attributes (REQUIRED):
        name          : sender identifier
        enabled_tools : tool names from TOOL_REGISTRY this agent may call
                        (excluding `emit_message`, which is always exposed)

    Subclass overrides (REQUIRED):
        build_system_prompt(ctx) -> str
        build_input_payload(ctx) -> str

    Tunables (override per subclass if needed):
        max_turns          : hard ceiling on tool-use loops
        max_tokens_per_turn: per-turn LLM cap
    """

    name: str = ""
    enabled_tools: Sequence[str] = ()
    max_turns: int = 8
    max_tokens_per_turn: int = 4096
    # When True, a failed numeric-fidelity check (a decimal in the message that
    # does not trace to any tool output the agent received) is fed back to the
    # model for one self-correction instead of being merely logged. Off by
    # default so a single false positive cannot wedge the pipeline; readers
    # that emit measurements (Interaction / Conformation / Interface) opt in.
    numeric_audit_strict: bool = False
    # Soft cap on analysis/search tool calls before the agent is told to
    # stop gathering and emit. Guards against the model looping on tools
    # (e.g. re-querying thin views) and never finalizing within max_turns.
    max_tool_calls: int = 12

    def __init__(
        self,
        llm: LLMClient,
        *,
        extra_tools: Optional[Dict[str, Tool]] = None,
    ) -> None:
        if not self.name:
            raise ValueError(
                f"{type(self).__name__} must declare a non-empty `name` "
                "class attribute."
            )
        self._llm = llm
        self._tools: Dict[str, Tool] = self._resolve_tools(self.enabled_tools)
        # `extra_tools` lets a subclass inject tools that are not in the
        # global TOOL_REGISTRY — e.g. Bio Agent's `search_literature`,
        # which only makes sense inside the multi-agent pipeline and is
        # not exposed to the main conversational agent.
        if extra_tools:
            for name, tool in extra_tools.items():
                if name in self._tools:
                    raise AgentToolError(
                        f"extra_tools conflicts with enabled_tools on {name!r}"
                    )
                self._tools[name] = tool

    # ------------------------------------------------------------------
    # Subclass hooks
    # ------------------------------------------------------------------

    @abstractmethod
    def build_system_prompt(self, ctx: MutationDesignContext) -> str:
        """Return the LLM system prompt for this agent given the run context."""

    @abstractmethod
    def build_input_payload(self, ctx: MutationDesignContext) -> str:
        """Return the user-role input that frames the agent's task.

        Typically a markdown blob containing: the design intent, relevant
        prior agent messages (with [ref:<msg_id>] hooks), and any
        per-agent view excerpts the agent should reason over.
        """

    # ------------------------------------------------------------------
    # Public entry
    # ------------------------------------------------------------------

    async def run(self, ctx: MutationDesignContext) -> AgentRunResult:
        """Execute the agent and append its emitted message to `ctx`."""
        system_prompt = self.build_system_prompt(ctx)
        user_payload = self.build_input_payload(ctx)
        full_tools_schema = self._build_tools_schema()
        emit_only_schema = [_EMIT_MESSAGE_SCHEMA]
        messages: List[Dict[str, Any]] = [
            {"role": "user", "content": user_payload}
        ]

        pending: Optional[_PendingMessage] = None
        pending_msg: Optional[AgentMessage] = None
        pending_audit: Any = None
        tool_calls_made: List[SkillCallRecord] = []

        for turn_idx in range(self.max_turns):
            assistant_content: List[Dict[str, Any]] = []
            tool_uses: List[Dict[str, Any]] = []
            text_buffer: List[str] = []
            stop_reason: Optional[str] = None

            # Once the agent has used its tool-call budget, hide the analysis
            # / search tools entirely so the only remaining action is
            # `emit_message` — a hard guarantee that it finalizes rather than
            # looping on tools until it times out.
            over_budget = len(tool_calls_made) >= self.max_tool_calls
            tools_schema = emit_only_schema if over_budget else full_tools_schema

            async for event in self._llm.stream(
                system=system_prompt,
                messages=messages,
                tools=tools_schema,
                max_tokens=self.max_tokens_per_turn,
            ):
                if isinstance(event, TextDelta):
                    text_buffer.append(event.text)
                elif isinstance(event, ToolUseStart):
                    tool_uses.append({"id": event.id, "name": event.name, "input": {}})
                elif isinstance(event, ToolUseEnd):
                    for tu in tool_uses:
                        if tu["id"] == event.id:
                            tu["input"] = event.input or {}
                            break
                    else:
                        tool_uses.append(
                            {"id": event.id, "name": event.name, "input": event.input or {}}
                        )
                elif isinstance(event, MessageStop):
                    stop_reason = event.stop_reason
                    if event.content:
                        assistant_content = event.content
                # Other event types (input deltas, etc.) ignored — we only
                # need final structured assistant content for replay.

            if not assistant_content:
                # Fall back to a synthesized assistant turn for providers
                # that don't echo content in MessageStop.
                assistant_content = self._synthesize_assistant_content(
                    text_buffer, tool_uses
                )
            messages.append({"role": "assistant", "content": assistant_content})

            if not tool_uses:
                # No tool calls + no emit_message ⇒ malformed turn. Push a
                # nudge and let the agent try again next turn.
                if stop_reason != "tool_use":
                    messages.append(
                        {
                            "role": "user",
                            "content": (
                                "You did not call any tool this turn. When you "
                                "are ready to finalize, call `emit_message`."
                            ),
                        }
                    )
                continue

            tool_results: List[Dict[str, Any]] = []
            finalize_now = False
            for tu in tool_uses:
                tool_name = tu["name"]
                args = tu["input"] or {}
                if tool_name == "emit_message":
                    try:
                        pending = self._parse_emit(args)
                        # Build + validate the message HERE so that a
                        # narrative<->evidence_refs inconsistency (e.g. an
                        # inline [ref:X] that was never declared) is fed back
                        # to the model for self-correction instead of killing
                        # the agent after the loop.
                        candidate_msg = self._build_agent_message(pending)
                        validate_message(candidate_msg)
                    except (MessageValidationError, ValueError, KeyError) as e:
                        tool_results.append(
                            self._tool_result_block(
                                tu["id"],
                                f"emit_message validation failed: {e}. "
                                "Every inline [ref:<id>] token must match a "
                                "ref_id you declared in evidence_refs. Fix the "
                                "refs and call emit_message again.",
                                is_error=True,
                            )
                        )
                        continue

                    # Numeric-fidelity check: every decimal in the narrative /
                    # payload must trace to a tool result this agent received.
                    audit_report = self._audit_numbers(candidate_msg, tool_calls_made)
                    if (
                        audit_report is not None
                        and self.numeric_audit_strict
                        and not audit_report.ok
                    ):
                        bad = ", ".join(
                            f"{c.raw} (in {c.where})"
                            for c in audit_report.ungrounded_decimals[:8]
                        )
                        tool_results.append(
                            self._tool_result_block(
                                tu["id"],
                                "emit_message rejected: these numbers do not "
                                f"appear in any tool output you called: {bad}. "
                                "Every RRCS / SASA / occupancy / BSA / RMSF value "
                                "must be copied from a tool result. Remove or "
                                "correct the unsupported numbers and call "
                                "emit_message again.",
                                is_error=True,
                            )
                        )
                        continue

                    pending_msg = candidate_msg
                    pending_audit = audit_report
                    tool_results.append(
                        self._tool_result_block(
                            tu["id"],
                            "Message accepted. Stop now — do not call any "
                            "further tools.",
                            is_error=False,
                        )
                    )
                    finalize_now = True
                else:
                    result, record = await self._dispatch_tool(
                        tool_name, args, tu["id"], ctx
                    )
                    if record is not None:
                        tool_calls_made.append(record)
                        ctx.add_skill_call(record)
                    tool_results.append(
                        self._tool_result_block(
                            tu["id"], result.content, is_error=result.is_error
                        )
                    )

            messages.append({"role": "user", "content": tool_results})

            if finalize_now and pending_msg is not None:
                ctx.add_message(pending_msg)
                if pending_audit is not None:
                    log.info("%s", pending_audit.summary())
                return AgentRunResult(
                    msg_id=pending_msg.msg_id,
                    message=pending_msg,
                    turns_used=turn_idx + 1,
                    tool_calls=tool_calls_made,
                    numeric_audit=pending_audit,
                )

            # Stop the agent from spending its whole budget querying tools
            # and never finalizing. Nudge to emit once it is near the turn
            # ceiling OR has already made plenty of tool calls (looping on
            # thin views). Either condition triggers the same instruction.
            turns_left = self.max_turns - turn_idx - 1
            over_tool_budget = len(tool_calls_made) >= self.max_tool_calls
            if turns_left > 0 and (turns_left <= 3 or over_tool_budget):
                reason = (
                    "You have gathered enough data "
                    f"({len(tool_calls_made)} tool calls)."
                    if over_tool_budget
                    else f"You have {turns_left} turn(s) left."
                )
                messages.append(
                    {
                        "role": "user",
                        "content": (
                            f"{reason} Stop calling analysis/search tools now "
                            "and call `emit_message` with your brief based on "
                            "the evidence you already have. Do not request "
                            "more data."
                        ),
                    }
                )

        raise AgentTimeoutError(
            f"Agent {self.name!r} did not call emit_message within "
            f"{self.max_turns} turns."
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _resolve_tools(self, names: Sequence[str]) -> Dict[str, Tool]:
        resolved: Dict[str, Tool] = {}
        for name in names:
            tool = TOOL_REGISTRY.get(name)
            if tool is None:
                raise AgentToolError(
                    f"Agent {self.name!r} requested unknown tool {name!r}; "
                    f"available: {sorted(TOOL_REGISTRY.keys())}"
                )
            resolved[name] = tool
        return resolved

    def _build_tools_schema(self) -> List[Dict[str, Any]]:
        schema = [_EMIT_MESSAGE_SCHEMA]
        for name, tool in self._tools.items():
            schema.append(
                {
                    "name": name,
                    "description": getattr(tool, "description", ""),
                    "input_schema": self._tool_input_schema(tool),
                }
            )
        return schema

    @staticmethod
    def _tool_input_schema(tool: Any) -> Dict[str, Any]:
        """Resolve a tool's JSON-Schema parameters as a proper object schema.

        Tools expose their input as a pydantic `Input` model rather than a
        pre-built `input_schema`. An empty `{}` is rejected by strict
        OpenAI-compatible providers (DeepSeek) as `type: null`, so derive a
        real `{"type": "object", ...}` schema from the model when possible.
        """
        sch = getattr(tool, "input_schema", None)
        if isinstance(sch, dict) and sch.get("type"):
            return sch
        model = getattr(tool, "Input", None)
        if model is not None and hasattr(model, "model_json_schema"):
            js = model.model_json_schema()
            js.setdefault("type", "object")
            js.setdefault("properties", {})
            return js
        return {"type": "object", "properties": {}}

    async def _dispatch_tool(
        self,
        tool_name: str,
        args: Dict[str, Any],
        tool_use_id: str,
        ctx: MutationDesignContext,
    ) -> tuple[ToolResult, Optional[SkillCallRecord]]:
        tool = self._tools.get(tool_name)
        if tool is None:
            return (
                ToolResult(
                    content=(
                        f"Tool {tool_name!r} is not enabled for agent "
                        f"{self.name!r}. Allowed: {sorted(self._tools)}."
                    ),
                    is_error=True,
                ),
                None,
            )
        tool_ctx = self._make_tool_context(ctx)
        t0 = time.monotonic()
        try:
            validated = tool.Input.model_validate(args)
            result = await tool.call(validated, tool_ctx)
        except Exception as e:
            log.exception("tool %s raised", tool_name)
            return (
                ToolResult(content=f"Tool {tool_name} raised: {e}", is_error=True),
                None,
            )
        duration_ms = int((time.monotonic() - t0) * 1000)
        record = SkillCallRecord(
            call_id=tool_use_id,
            agent=self.name,
            skill_name=tool_name,
            arguments=args,
            result_summary=self._summarize_tool_result(result, duration_ms),
            raw_result_text=self._stringify_tool_content(result.content),
            is_error=bool(result.is_error),
        )
        return result, record

    def _make_tool_context(self, ctx: MutationDesignContext) -> ToolContext:
        """Tools expect a `ToolContext`; we synthesize a minimal one.

        Multi-agent v1 routes through tools that read `case_dir` from
        their own arguments (e.g. `query_analysis_results`), so the
        session-level fields here only need to be valid enough to satisfy
        the dataclass. Tools that touch the SQLite audit DB or emit UI
        events are not enabled for any v1 reader/writer.
        """
        import anyio

        from immunoscope.agent.config import get_settings

        async def _noop_event(_payload: dict) -> None:
            return None

        return ToolContext(
            session_id=f"multi_agent::{ctx.case_id}",
            abort_event=anyio.Event(),
            on_event=_noop_event,
            settings=get_settings(),
            db_path="",
        )

    def _tool_result_block(
        self, tool_use_id: str, content: Any, is_error: bool
    ) -> Dict[str, Any]:
        if isinstance(content, (dict, list)):
            content_str = json.dumps(content, ensure_ascii=False, default=str)
        else:
            content_str = "" if content is None else str(content)
        return {
            "type": "tool_result",
            "tool_use_id": tool_use_id,
            "content": content_str,
            "is_error": is_error,
        }

    @staticmethod
    def _stringify_tool_content(content: Any) -> str:
        """Full stringified tool output for the numeric-fidelity audit.

        Unlike `_summarize_tool_result` (a 200-char digest for the audit log),
        this keeps the whole content so the post-hoc check can confirm a
        number an agent cited actually appeared in a tool result."""
        if isinstance(content, str):
            return content
        if content is None:
            return ""
        try:
            return json.dumps(content, ensure_ascii=False, default=str)
        except Exception:  # noqa: BLE001
            return str(content)

    def _summarize_tool_result(self, result: ToolResult, duration_ms: int) -> str:
        if isinstance(result.content, str):
            head = result.content[:200]
        else:
            try:
                head = json.dumps(result.content, ensure_ascii=False, default=str)[:200]
            except Exception:
                head = repr(result.content)[:200]
        tag = "ERR" if result.is_error else "OK"
        return f"[{tag} {duration_ms}ms] {head}"

    def _parse_emit(self, args: Dict[str, Any]) -> _PendingMessage:
        narrative = args.get("narrative")
        if not isinstance(narrative, str) or not narrative.strip():
            raise ValueError("`narrative` must be a non-empty string.")
        raw_refs = args.get("evidence_refs")
        if not isinstance(raw_refs, list) or not raw_refs:
            raise ValueError("`evidence_refs` must be a non-empty list.")
        refs: List[EvidenceRef] = []
        for i, raw in enumerate(raw_refs):
            if not isinstance(raw, dict):
                raise ValueError(f"evidence_refs[{i}] must be an object.")
            refs.append(evidence_ref_from_dict(raw))
        recipients = args.get("recipients", "broadcast")
        if recipients != "broadcast" and not isinstance(recipients, list):
            raise ValueError(
                "`recipients` must be 'broadcast' or a list of agent names."
            )
        payload = args.get("structured_payload")
        if payload is not None and not isinstance(payload, dict):
            raise ValueError("`structured_payload` must be a dict or null.")
        return _PendingMessage(
            narrative=narrative,
            evidence_refs=refs,
            structured_payload=payload,
            recipients=recipients,
        )

    def _audit_numbers(
        self, msg: AgentMessage, tool_calls: Sequence[SkillCallRecord]
    ) -> Any:
        """Run the numeric-fidelity check for agents that produce measurements.

        Returns a `NumericAuditReport`, or None for tool-less integrators
        (Recommendation / Critic) whose numbers come from upstream briefs, not
        their own tool calls — grounding those is a different check.
        """
        if not tool_calls:
            return None
        from .numeric_audit import audit_agent_message

        try:
            return audit_agent_message(
                self.name, msg.narrative, msg.structured_payload, tool_calls
            )
        except Exception:  # noqa: BLE001 — audit must never break a run
            log.exception("numeric audit raised for agent %s", self.name)
            return None

    def _build_agent_message(self, pending: _PendingMessage) -> AgentMessage:
        return AgentMessage(
            msg_id=self._mint_msg_id(),
            sender=self.name,
            recipients=pending.recipients,
            narrative=pending.narrative,
            evidence_refs=pending.evidence_refs,
            structured_payload=pending.structured_payload,
        )

    def _mint_msg_id(self) -> str:
        return f"{self.name}:{uuid.uuid4().hex[:8]}"

    def _synthesize_assistant_content(
        self, text_chunks: Sequence[str], tool_uses: Sequence[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """Build a content list when the LLM client did not return one.

        Anthropic's client supplies `MessageStop.content` (full structured
        blocks) but the OpenAI-compatible adapter only fills it with the
        tool_use blocks. We reconstruct text blocks from streamed deltas
        so the conversation replay stays consistent.
        """
        content: List[Dict[str, Any]] = []
        joined = "".join(text_chunks).strip()
        if joined:
            content.append({"type": "text", "text": joined})
        for tu in tool_uses:
            content.append(
                {
                    "type": "tool_use",
                    "id": tu["id"],
                    "name": tu["name"],
                    "input": tu["input"],
                }
            )
        return content


__all__ = [
    "AgentRunResult",
    "AgentTimeoutError",
    "AgentToolError",
    "BaseAgent",
]
