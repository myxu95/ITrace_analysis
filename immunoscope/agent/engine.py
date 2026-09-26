from __future__ import annotations
import asyncio
import json
import logging
import time
from contextlib import suppress
from pathlib import Path
from typing import Any

from pydantic import ValidationError

from immunoscope.agent import audit, hooks
from immunoscope.agent import hooks_builtin  # noqa: F401  -- side-effect: registers BudgetGuardHook
from immunoscope.agent.budget import BudgetTracker
from immunoscope.agent.cost_tracker import CostTracker
from immunoscope.agent.context import build_system_prompt
from immunoscope.agent.exceptions import (
    AgentError,
    LLMStreamError,
    MaxTurnsExceeded,
    PermissionDenied,
    ToolError,
    ToolFatal,
    ToolTimeout,
    ToolValidationError,
    WallClockExceeded,
)
from immunoscope.agent.llm import LLMClient, make_llm_client
from immunoscope.agent.permissions import check_permission
from immunoscope.agent.session import Session
from immunoscope.agent.stream_events import (
    MessageStop,
    StreamError,
    TextDelta,
    ToolUseEnd,
    ToolUseInputDelta,
    ToolUseStart,
)
from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import TOOL_REGISTRY, get_enabled_tools

log = logging.getLogger("immunoscope.agent.engine")


def _redact(d: Any) -> Any:
    if isinstance(d, dict):
        return {k: _redact(v) for k, v in d.items()}
    if isinstance(d, list):
        return [_redact(v) for v in d]
    if isinstance(d, str) and len(d) > 500:
        return d[:500] + f"...[+{len(d) - 500} chars]"
    return d


def _summary_for_event(result: ToolResult, max_chars: int = 240) -> str:
    if isinstance(result.content, str):
        s = result.content
    elif result.content is None:
        s = ""
    else:
        try:
            s = json.dumps(result.content, ensure_ascii=False, default=str)
        except Exception:
            s = repr(result.content)
    if len(s) > max_chars:
        return s[:max_chars] + "…"
    return s


def _maybe_truncate(result: ToolResult, max_chars: int, session_id: str, tool_use_id: str,
                    upload_dir: Path) -> ToolResult:
    """If the serialized content exceeds max_chars, write the full payload to disk
    and replace the model-facing content with a head + a hint pointing at the
    spill file. Mirrors CLI's `maxResultSizeChars`."""
    if isinstance(result.content, str):
        s = result.content
    elif result.content is None:
        return result
    else:
        try:
            s = json.dumps(result.content, ensure_ascii=False, default=str)
        except Exception:
            s = repr(result.content)
    if len(s) <= max_chars:
        return result

    upload_dir.mkdir(parents=True, exist_ok=True)
    spill = upload_dir / f"{session_id}_{tool_use_id}.txt"
    try:
        spill.write_text(s, encoding="utf-8")
    except Exception:
        log.warning("spill write failed for %s", spill, exc_info=True)
        spill = None  # type: ignore[assignment]

    head = s[:max_chars]
    suffix = f"\n\n[truncated: {len(s) - max_chars} more chars"
    if spill:
        suffix += f"; full output saved to {spill}"
    suffix += "]"
    return ToolResult(
        content=head + suffix,
        is_error=result.is_error,
        truncated=True,
        saved_to=str(spill) if spill else None,
    )


async def _call_with_interrupt(tool: Tool, args, ctx: ToolContext, timeout: float) -> ToolResult:
    """Run tool.call honoring `interrupt_behavior` against ctx.abort_event.

    - "cancel": abort interrupts the running task → CancelledError.
    - "block":  abort is acknowledged but the task is allowed to complete
                so transactional state stays clean.
    Per-call timeout still applies in both modes."""
    call_task = asyncio.create_task(tool.call(args, ctx))
    timeout_task = asyncio.create_task(asyncio.sleep(timeout))
    abort_task = (
        asyncio.create_task(ctx.abort_event.wait())
        if tool.interrupt_behavior == "cancel"
        else None
    )

    waitset = {call_task, timeout_task}
    if abort_task is not None:
        waitset.add(abort_task)

    try:
        done, _pending = await asyncio.wait(waitset, return_when=asyncio.FIRST_COMPLETED)
        if call_task in done:
            return call_task.result()
        if timeout_task in done:
            call_task.cancel()
            with suppress(asyncio.CancelledError, Exception):
                await call_task
            raise asyncio.TimeoutError()
        # abort fired (only reachable when interrupt_behavior == "cancel")
        call_task.cancel()
        with suppress(asyncio.CancelledError, Exception):
            await call_task
        raise asyncio.CancelledError("aborted by user")
    finally:
        for t in (timeout_task, abort_task):
            if t is not None and not t.done():
                t.cancel()
                with suppress(asyncio.CancelledError, Exception):
                    await t


async def execute_tool(
    *,
    tool: Tool,
    tool_use_id: str,
    raw_input: dict,
    ctx: ToolContext,
    turn_idx: int,
) -> dict:
    """Validate → permission → call → audit → serialize. All exceptions are
    converted to tool_result(is_error=True) so the model can recover."""
    started = time.time()
    err: BaseException | None = None
    result: ToolResult

    try:
        try:
            args = tool.Input.model_validate(raw_input)
        except ValidationError as ve:
            raise ToolValidationError(str(ve))

        await tool.validate(args, ctx)

        decision = await check_permission(tool, args, ctx, tool_use_id)
        from immunoscope.agent.tool import Permission as _P
        if decision.behavior == _P.DENY:
            raise PermissionDenied(decision.reason or "denied")

        timeout = tool.timeout_seconds or ctx.settings.AGENT_TOOL_TIMEOUT_DEFAULT
        result = await _call_with_interrupt(tool, args, ctx, timeout)

    except ToolError as e:
        err = e
        result = ToolResult(content=f"<tool_error>{type(e).__name__}: {e}</tool_error>", is_error=True)
    except asyncio.TimeoutError:
        err = ToolTimeout(f"tool {tool.name} timed out")
        result = ToolResult(content=f"<tool_error>ToolTimeout: tool {tool.name} timed out</tool_error>", is_error=True)
    except asyncio.CancelledError as e:
        # Two cases:
        #   1. _call_with_interrupt raised CancelledError("aborted by user")
        #      because the user fired abort and tool.interrupt_behavior=="cancel".
        #      We turn this into an error tool_result so siblings in the same
        #      turn can still complete and the transcript stays consistent.
        #   2. The whole task was cancelled from outside (e.g. WS teardown).
        #      We can detect by inspecting ctx.aborted: if abort_event is set,
        #      this is case 1; otherwise propagate.
        if ctx.aborted:
            err = e
            result = ToolResult(
                content=f"<tool_error>Aborted: tool {tool.name} cancelled by user</tool_error>",
                is_error=True,
            )
        else:
            raise
    except Exception as e:
        log.exception("unexpected tool failure: %s", tool.name)
        err = ToolFatal(repr(e))
        result = ToolResult(content=f"<tool_error>InternalError: {e}</tool_error>", is_error=True)

    # Attach this turn's token/cost deltas only to the FIRST tool of the
    # turn so SUM aggregations don't double count when a turn fires several
    # tool_use blocks. ctx.last_turn_usage is set by run_turn after each
    # LLM stream end; we consume-and-clear it here.
    audit_usage = ctx.last_turn_usage
    audit_cost = ctx.last_turn_cost
    ctx.last_turn_usage = None
    ctx.last_turn_cost = None

    # Audit BEFORE truncation so we keep the full original size info, but with
    # a summary string. Audit failure must not block the turn.
    if ctx.settings.AGENT_AUDIT_ENABLED:
        try:
            await asyncio.wait_for(
                audit.write_tool_call_audit(
                    ctx.db_path,
                    session_id=ctx.session_id,
                    turn_idx=turn_idx,
                    tool_use_id=tool_use_id,
                    tool_name=tool.name,
                    args=_redact(raw_input),
                    result=result,
                    duration_seconds=time.time() - started,
                    error=err,
                    input_tokens=int(audit_usage.get("input_tokens", 0)) if audit_usage else None,
                    output_tokens=int(audit_usage.get("output_tokens", 0)) if audit_usage else None,
                    cache_read_tokens=int(audit_usage.get("cache_read_input_tokens", 0)) if audit_usage else None,
                    cache_creation_tokens=int(audit_usage.get("cache_creation_input_tokens", 0)) if audit_usage else None,
                    cost_usd=audit_cost,
                ),
                timeout=0.5,
            )
        except Exception:
            log.warning("audit write failed (non-fatal)", exc_info=True)

    # Spill oversize results to disk and replace with a head + hint.
    upload_dir = Path(ctx.settings.DATA_DIR) / "agent_tool_outputs"
    result = _maybe_truncate(
        result,
        ctx.settings.AGENT_TOOL_RESULT_MAX_CHARS,
        ctx.session_id,
        tool_use_id,
        upload_dir,
    )

    await ctx.on_event({
        "type": "tool_use_end",
        "id": tool_use_id,
        "name": tool.name,
        "ok": not result.is_error,
        "summary": _summary_for_event(result),
        "truncated": result.truncated,
    })
    return tool.serialize_result(result, tool_use_id)


async def execute_tools(
    tool_uses: list[ToolUseEnd],
    ctx: ToolContext,
    turn_idx: int,
) -> list[dict]:
    """Concurrency-safe tools (read-only by default) run in parallel; unsafe
    tools run serially. The order in the returned list matches `tool_uses`
    (Anthropic API requires tool_result order to match tool_use order in the
    paired assistant message)."""
    safe: list[tuple[Tool, ToolUseEnd]] = []
    unsafe: list[tuple[Tool, ToolUseEnd]] = []
    for tu in tool_uses:
        tool = TOOL_REGISTRY.get(tu.name)
        if tool is None:
            # unknown tool: synthesize an error tool_result so the model can recover
            log.warning("model called unknown tool: %s", tu.name)
            results_by_id: dict[str, dict] = {}
            results_by_id[tu.id] = {
                "type": "tool_result",
                "tool_use_id": tu.id,
                "content": f"<tool_error>UnknownTool: {tu.name!r} is not available</tool_error>",
                "is_error": True,
            }
            continue
        (safe if tool.is_concurrency_safe else unsafe).append((tool, tu))

    results_by_id: dict[str, dict] = {}
    if safe:
        coros = [
            execute_tool(tool=t, tool_use_id=tu.id, raw_input=tu.input, ctx=ctx, turn_idx=turn_idx)
            for t, tu in safe
        ]
        outcomes = await asyncio.gather(*coros)
        for (_t, tu), r in zip(safe, outcomes):
            results_by_id[tu.id] = r
    for t, tu in unsafe:
        results_by_id[tu.id] = await execute_tool(
            tool=t, tool_use_id=tu.id, raw_input=tu.input, ctx=ctx, turn_idx=turn_idx,
        )

    # Re-add the unknown-tool errors in the right place (rare path)
    ordered: list[dict] = []
    for tu in tool_uses:
        if tu.id in results_by_id:
            ordered.append(results_by_id[tu.id])
        else:
            ordered.append({
                "type": "tool_result",
                "tool_use_id": tu.id,
                "content": f"<tool_error>UnknownTool: {tu.name!r} is not available</tool_error>",
                "is_error": True,
            })
    return ordered


async def _run_stream(
    llm: LLMClient,
    *,
    system: str,
    messages: list[dict],
    tools: list[dict],
    tool_use_blocks: list[ToolUseEnd],
    on_event,
) -> MessageStop | None:
    """Inner stream consumer. Designed to live inside an asyncio.Task so it can
    be cancelled cleanly when the abort_event fires (cancellation propagates
    into the Anthropic SDK's underlying httpx read).

    Returns the MessageStop event so the caller has the full final.content
    (including provider-specific blocks like `thinking` that must be echoed
    verbatim on the next API call)."""
    async for ev in llm.stream(system=system, messages=messages, tools=tools):
        if isinstance(ev, TextDelta):
            await on_event({"type": "text_delta", "text": ev.text})
        elif isinstance(ev, ToolUseStart):
            await on_event({"type": "tool_use_start", "id": ev.id, "name": ev.name})
        elif isinstance(ev, ToolUseInputDelta):
            await on_event({
                "type": "tool_use_input_delta",
                "id": ev.id,
                "partial": ev.partial_json,
            })
        elif isinstance(ev, ToolUseEnd):
            tool_use_blocks.append(ev)
        elif isinstance(ev, MessageStop):
            return ev
        elif isinstance(ev, StreamError):
            raise LLMStreamError() from ev.error
    return None


async def run_turn(session: Session, user_msg: str, ctx: ToolContext) -> None:
    """Execute one user turn: append message → loop {LLM stream → tool exec}
    until no tool_use blocks come back. Caller handles errors emitted via
    on_event."""
    session.messages.append({"role": "user", "content": user_msg})
    deadline = time.time() + ctx.settings.AGENT_SESSION_WALL_CLOCK
    llm = make_llm_client(ctx.settings)
    tools_schema = [
        t.to_anthropic_schema()
        for t in get_enabled_tools(mode=ctx.settings.AGENT_MODE)
    ]
    tracker = BudgetTracker(
        budget=ctx.settings.AGENT_TOKEN_BUDGET,
        diminishing_threshold=ctx.settings.AGENT_DIMINISHING_THRESHOLD,
        diminishing_turns=ctx.settings.AGENT_DIMINISHING_TURNS,
    )
    if session.cost_tracker is None:
        session.cost_tracker = CostTracker(model=ctx.settings.LLM_MODEL)

    try:
        for turn_idx in range(ctx.settings.AGENT_MAX_TURNS):
            session.turn_idx = turn_idx
            if time.time() > deadline:
                raise WallClockExceeded()
            if ctx.aborted:
                return

            tool_use_blocks: list[ToolUseEnd] = []
            system_prompt = build_system_prompt(session)

            stream_task = asyncio.create_task(
                asyncio.wait_for(
                    _run_stream(
                        llm,
                        system=system_prompt,
                        messages=session.messages,
                        tools=tools_schema,
                        tool_use_blocks=tool_use_blocks,
                        on_event=ctx.on_event,
                    ),
                    timeout=ctx.settings.AGENT_TURN_TIMEOUT,
                )
            )
            abort_task = asyncio.create_task(ctx.abort_event.wait())

            done, _pending = await asyncio.wait(
                {stream_task, abort_task},
                return_when=asyncio.FIRST_COMPLETED,
            )

            if abort_task in done and stream_task not in done:
                stream_task.cancel()
                with suppress(asyncio.CancelledError, Exception):
                    await stream_task
                return

            abort_task.cancel()
            with suppress(asyncio.CancelledError):
                await abort_task

            try:
                stop_ev: MessageStop | None = await stream_task
            except asyncio.TimeoutError:
                raise LLMStreamError("turn timeout exceeded") from None

            # Update budget + cost trackers before any decision logic.
            # `usage` may be None on provider-side errors; trackers tolerate that.
            usage = stop_ev.usage if stop_ev else None
            tracker.update_after_turn(usage)
            turn_cost, total_cost = session.cost_tracker.add_turn(usage)
            await ctx.on_event({
                "type": "cost_update",
                "turn_idx": turn_idx,
                "turn_cost_usd": round(turn_cost, 6),
                **session.cost_tracker.snapshot(),
            })
            # Stash this turn's deltas so execute_tools can attach them to
            # the FIRST tool's audit row (avoids double-counting on SUM).
            ctx.last_turn_usage = usage  # type: ignore[attr-defined]
            ctx.last_turn_cost = turn_cost  # type: ignore[attr-defined]

            # Echo the SDK's final.content verbatim — this preserves
            # provider-specific blocks (e.g. DeepSeek's `thinking`) which the
            # API requires on subsequent requests. We do NOT reconstruct from
            # text_delta + ToolUseEnd; that drops blocks the model expects to see.
            assistant_content = (stop_ev.content if stop_ev else None) or []
            if not assistant_content:
                # Defensive: API may rarely produce no content at all (e.g.
                # provider-side error). Treat as turn end.
                await ctx.on_event({"type": "turn_done"})
                return
            session.messages.append({"role": "assistant", "content": assistant_content})

            # Hook dispatch: gives BudgetGuardHook (and future hooks) a say
            # before we either tool-execute or finish the turn.
            decision = await hooks.dispatch(
                "turn_end", ctx,
                tracker=tracker,
                session=session,
                turn_idx=turn_idx,
                stop_event=stop_ev,
                pending_tool_uses=tool_use_blocks,
            )

            if decision.action == "stop":
                log.info("hook stopped turn at idx=%d reason=%s", turn_idx, decision.reason)
                await ctx.on_event({"type": "turn_done"})
                return

            nudge_text = decision.content if decision.action == "inject_user_message" else None

            if not tool_use_blocks:
                if nudge_text:
                    # Give the model one more turn to honor the nudge.
                    session.messages.append({"role": "user", "content": nudge_text})
                    await ctx.on_event({
                        "type": "budget_notice",
                        "reason": decision.reason,
                        "message": nudge_text,
                    })
                    continue
                await ctx.on_event({"type": "turn_done"})
                return

            # Tool execution path — Anthropic requires the next user message
            # to start with tool_result blocks paired with the tool_uses we
            # just got back. We append the nudge text inside the SAME user
            # message (as a trailing text block) so it lands in chronological
            # order without producing two consecutive user messages.
            tool_results = await execute_tools(tool_use_blocks, ctx, turn_idx)
            user_content: list[dict] = list(tool_results)
            if nudge_text:
                user_content.append({"type": "text", "text": nudge_text})
                await ctx.on_event({
                    "type": "budget_notice",
                    "reason": decision.reason,
                    "message": nudge_text,
                })
            session.messages.append({"role": "user", "content": user_content})

        raise MaxTurnsExceeded(f"exceeded {ctx.settings.AGENT_MAX_TURNS} turns")
    finally:
        await llm.aclose() if hasattr(llm, "aclose") else None  # type: ignore[func-returns-value]
