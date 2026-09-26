"""Hook registry + dispatch.

This is the foundational scaffolding Sprint A delivers; full lifecycle
(PreToolUse / PostToolUse / SessionStart / SessionEnd) is Sprint B.

For now we expose ONE event — `turn_end` — fired immediately after the
LLM has produced its assistant message and before tool execution (so a
hook can short-circuit a runaway loop) and again at the very end of a
turn that produced no tool_use (so a hook can append a synthetic final
message before sending `turn_done`).

A hook returns a HookDecision:
  - {"action": "continue"}                                     — proceed
  - {"action": "stop", "reason": str}                          — emit
        turn_done immediately and exit run_turn
  - {"action": "inject_user_message", "content": str}          — append a
        synthetic user message and continue the loop. Used by
        BudgetGuardHook to nudge the model to summarize.

Decisions are processed in registration order; first non-continue wins.
"""
from __future__ import annotations
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Literal

log = logging.getLogger("prism.agent.hooks")


HookEvent = Literal["turn_end"]
HookAction = Literal["continue", "stop", "inject_user_message"]


@dataclass
class HookDecision:
    action: HookAction
    reason: str | None = None
    content: str | None = None  # only used when action == inject_user_message

    @classmethod
    def carry_on(cls) -> "HookDecision":
        return cls(action="continue")


class Hook(ABC):
    name: str

    @abstractmethod
    async def run(self, event: HookEvent, ctx: Any, /, **kw) -> HookDecision: ...


_HOOKS: list[Hook] = []


def register_hook(hook: Hook) -> Hook:
    """Append a hook to the registry. Order matters: earlier entries get
    first say at each event."""
    _HOOKS.append(hook)
    return hook


def reset_hooks() -> None:
    """Test-only: clear the registry."""
    _HOOKS.clear()


def list_hooks() -> list[Hook]:
    return list(_HOOKS)


async def dispatch(event: HookEvent, ctx: Any, /, **kw) -> HookDecision:
    for h in _HOOKS:
        try:
            d = await h.run(event, ctx, **kw)
        except Exception:
            log.exception("hook %s raised on %s; treating as continue", h.name, event)
            continue
        if d.action != "continue":
            log.info("hook %s decided action=%s reason=%s", h.name, d.action, d.reason)
            return d
    return HookDecision.carry_on()
