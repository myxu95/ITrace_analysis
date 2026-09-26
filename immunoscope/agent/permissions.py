from __future__ import annotations
import asyncio
import logging
import uuid

from pydantic import BaseModel

from immunoscope.agent.tool import Permission, PermissionDecision, Tool, ToolContext

log = logging.getLogger("prism.agent.permissions")


def _summarize_args(args: BaseModel) -> dict:
    """Trim args for display in a permission_request event. Strings are clipped
    to 200 chars so we never ship a 1MB blob to the frontend."""
    raw = args.model_dump(mode="json")
    return {k: (v[:200] + "…") if isinstance(v, str) and len(v) > 200 else v
            for k, v in raw.items()}


def check_mode(tool: Tool, ctx: ToolContext) -> PermissionDecision | None:
    """Defense-in-depth: deny if the active AGENT_MODE isn't in tool.allowed_modes.

    The primary gate is `tools.get_enabled_tools(mode=...)`, which filters the
    schema before it ever reaches the LLM. This check covers the case where a
    stale tool_use_id from an earlier session, or a hand-crafted client, tries
    to invoke a tool not registered for the current mode.

    Returns None when the call is allowed (so the caller falls through to the
    normal needs_permission path)."""
    mode = getattr(ctx.settings, "AGENT_MODE", "agent")
    if mode not in tool.allowed_modes:
        return PermissionDecision(
            Permission.DENY,
            reason=f"tool {tool.name!r} not available in agent mode {mode!r}",
        )
    return None


async def check_permission(
    tool: Tool,
    args: BaseModel,
    ctx: ToolContext,
    tool_use_id: str,
) -> PermissionDecision:
    """Three-state gate run after validate() and before call().

    For first delivery all read-only tools have needs_permission=False, so this
    short-circuits to ALLOW. The ASK path is fully implemented so action tools
    in the next phase only need to set the flag."""
    mode_decision = check_mode(tool, ctx)
    if mode_decision is not None:
        return mode_decision
    if not tool.needs_permission:
        return PermissionDecision(Permission.ALLOW)
    return await _ask_user(tool, args, ctx, tool_use_id)


async def _ask_user(
    tool: Tool,
    args: BaseModel,
    ctx: ToolContext,
    tool_use_id: str,
) -> PermissionDecision:
    request_id = uuid.uuid4().hex
    loop = asyncio.get_running_loop()
    fut: asyncio.Future = loop.create_future()
    ctx.pending_permissions[request_id] = fut

    try:
        await ctx.on_event({
            "type": "permission_request",
            "request_id": request_id,
            "tool_use_id": tool_use_id,
            "tool": tool.name,
            "args_summary": _summarize_args(args),
            "is_destructive": tool.is_destructive,
        })
        decision_msg = await asyncio.wait_for(
            fut, timeout=ctx.settings.AGENT_PERMISSION_TIMEOUT
        )
    except asyncio.TimeoutError:
        log.info("permission ask timed out for tool=%s id=%s", tool.name, tool_use_id)
        return PermissionDecision(Permission.DENY, reason="user did not respond in time")
    except asyncio.CancelledError:
        # router cancels pending futures on disconnect/abort
        return PermissionDecision(Permission.DENY, reason="session aborted")
    finally:
        ctx.pending_permissions.pop(request_id, None)

    decision = decision_msg.get("decision") if isinstance(decision_msg, dict) else None
    if decision == "allow":
        return PermissionDecision(Permission.ALLOW)
    if decision == "deny":
        return PermissionDecision(Permission.DENY, reason="user denied")
    return PermissionDecision(Permission.DENY, reason=f"unrecognized response: {decision_msg!r}")
