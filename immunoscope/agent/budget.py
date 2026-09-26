"""Per-session token-budget tracker.

Mirrors `query/tokenBudget.ts:L3-92` from the CLI: the agent stops itself
when it has spent most of its budget, or when the last few turns have
each contributed very little new output (diminishing returns).

The check returns a structured signal that engine.run_turn forwards to
the hook registry. The hook decides whether to continue, inject a final
"please summarize and stop" user message, or hard-stop.
"""
from __future__ import annotations
import logging
from dataclasses import dataclass, field
from typing import Literal

log = logging.getLogger("prism.agent.budget")


# What the hook layer sees. Mirrors CLI tokenBudget action enum.
BudgetAction = Literal["continue", "near_limit", "diminishing"]


@dataclass
class BudgetSnapshot:
    cumulative_tokens: int
    cumulative_input: int
    cumulative_output: int
    cache_read: int
    cache_creation: int
    last_turn_delta: int
    consecutive_diminishing: int

    def as_dict(self) -> dict:
        return {
            "tokens": self.cumulative_tokens,
            "input": self.cumulative_input,
            "output": self.cumulative_output,
            "cache_read": self.cache_read,
            "cache_creation": self.cache_creation,
            "last_turn_delta": self.last_turn_delta,
            "consecutive_diminishing": self.consecutive_diminishing,
        }


@dataclass
class BudgetTracker:
    """Lifetime: one per `run_turn` call. The session itself is not the
    right scope — budget is per-task, not per-conversation."""

    budget: int  # AGENT_TOKEN_BUDGET; 0 disables limit
    diminishing_threshold: int = 500
    diminishing_turns: int = 3

    cumulative_input: int = 0
    cumulative_output: int = 0
    cache_read: int = 0
    cache_creation: int = 0
    _consecutive_diminishing: int = 0
    _last_total: int = 0

    @property
    def cumulative_tokens(self) -> int:
        return self.cumulative_input + self.cumulative_output

    def update_after_turn(self, usage: dict | None) -> int:
        """Record a turn's usage delta. Returns the per-turn delta in tokens."""
        if not usage:
            return 0
        self.cumulative_input += int(usage.get("input_tokens", 0) or 0)
        self.cumulative_output += int(usage.get("output_tokens", 0) or 0)
        self.cache_read += int(usage.get("cache_read_input_tokens", 0) or 0)
        self.cache_creation += int(usage.get("cache_creation_input_tokens", 0) or 0)
        delta = self.cumulative_tokens - self._last_total
        self._last_total = self.cumulative_tokens
        if delta < self.diminishing_threshold:
            self._consecutive_diminishing += 1
        else:
            self._consecutive_diminishing = 0
        return delta

    def snapshot(self) -> BudgetSnapshot:
        return BudgetSnapshot(
            cumulative_tokens=self.cumulative_tokens,
            cumulative_input=self.cumulative_input,
            cumulative_output=self.cumulative_output,
            cache_read=self.cache_read,
            cache_creation=self.cache_creation,
            last_turn_delta=self.cumulative_tokens - 0,  # not exposed; for snapshot only
            consecutive_diminishing=self._consecutive_diminishing,
        )

    def check(self) -> tuple[BudgetAction, str | None]:
        """Return (action, optional continuation message)."""
        if self.budget <= 0:
            return ("continue", None)
        if self.cumulative_tokens >= int(self.budget * 0.9):
            return (
                "near_limit",
                (
                    f"You have used {self.cumulative_tokens} of your "
                    f"{self.budget}-token budget. Wrap up: produce a final "
                    f"answer based on what you already know — do NOT call "
                    f"more tools."
                ),
            )
        if self._consecutive_diminishing >= self.diminishing_turns:
            return (
                "diminishing",
                (
                    f"The last {self._consecutive_diminishing} turns added "
                    f"very little new content. If you have enough to answer, "
                    f"give a concise final response now without calling tools."
                ),
            )
        return ("continue", None)
