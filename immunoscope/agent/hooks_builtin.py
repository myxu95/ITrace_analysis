"""Bundled hooks. Registered at import time."""
from __future__ import annotations
import logging

from immunoscope.agent.budget import BudgetTracker
from immunoscope.agent.hooks import Hook, HookDecision, register_hook

log = logging.getLogger("prism.agent.hooks_builtin")


class BudgetGuardHook(Hook):
    """Watches the per-turn BudgetTracker and asks the model to summarize
    when budget is nearly exhausted or when the conversation has stalled
    (diminishing returns). This hook is the single source of truth for
    Sprint A's #1 (token budget) feature.

    Two-stage enforcement:
    1. First over-budget detection → inject nudge (gives the model one
       turn to wrap up).
    2. Second consecutive over-budget detection → hard stop.

    Some models (notably DeepSeek-v4-pro in our testing) ignore the
    nudge and keep calling tools, so the second-strike hard stop is
    what actually protects against runaway costs."""

    name = "budget_guard"

    def __init__(self) -> None:
        # Per-session strike counter, keyed by Session.id. Reset when
        # the tracker reports continue.
        self._strikes: dict[str, int] = {}

    async def run(self, event, ctx, /, **kw) -> HookDecision:
        if event != "turn_end":
            return HookDecision.carry_on()
        tracker: BudgetTracker | None = kw.get("tracker")
        session = kw.get("session")
        if tracker is None or session is None:
            return HookDecision.carry_on()

        action, msg = tracker.check()
        sid = session.id
        if action == "continue":
            self._strikes.pop(sid, None)
            return HookDecision.carry_on()

        n = self._strikes.get(sid, 0) + 1
        self._strikes[sid] = n
        if n == 1:
            return HookDecision(
                action="inject_user_message",
                reason=action,
                content=(msg or "") + "\n\nThis is your final turn — produce the answer now without any tool_use blocks.",
            )
        # Second strike: stop unconditionally. Engine emits turn_done.
        log.warning(
            "BudgetGuard hard-stopping session=%s after %d strikes (action=%s, tokens=%d, budget=%d)",
            sid, n, action, tracker.cumulative_tokens, tracker.budget,
        )
        return HookDecision(
            action="stop",
            reason=f"budget_guard_hard_stop_after_{n}_strikes",
        )


register_hook(BudgetGuardHook())
