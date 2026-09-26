"""USD-cost estimator.

Pricing comes from a static table — values are placeholder estimates
that match each provider's published rates closely enough for visibility
in the UI. The numbers are NOT used for billing; they're a budgeting aid
for users so they know roughly what a long session is costing.

Update PRICING when a provider adjusts rates. Missing models log once
and return 0.
"""
from __future__ import annotations
import logging
from dataclasses import dataclass, field
from typing import TypedDict

log = logging.getLogger("prism.agent.cost")


class _ModelPricing(TypedDict, total=False):
    input: float            # USD per 1M input tokens
    output: float           # USD per 1M output tokens
    cache_read: float       # USD per 1M cache-read tokens (lower than input)
    cache_creation: float   # USD per 1M cache-creation tokens (often higher)


# USD per million tokens. Conservative placeholders for DeepSeek; verified
# rates for Claude. Replace freely.
PRICING: dict[str, _ModelPricing] = {
    "deepseek-v4-pro":   {"input": 0.30, "output": 1.00, "cache_read": 0.07},
    "deepseek-v4-flash": {"input": 0.10, "output": 0.30, "cache_read": 0.03},
    "claude-sonnet-4-5": {"input": 3.00, "output": 15.00,
                          "cache_read": 0.30, "cache_creation": 3.75},
    "claude-opus-4-1":   {"input": 15.00, "output": 75.00,
                          "cache_read": 1.50, "cache_creation": 18.75},
}

_WARNED: set[str] = set()


def _warn_once(model: str) -> None:
    if model in _WARNED:
        return
    _WARNED.add(model)
    log.warning("no pricing entry for model=%r; cost will be reported as 0", model)


def estimate_turn_cost(usage: dict | None, model: str) -> float:
    """Return USD cost for one turn's `usage` dict. Tolerates missing
    fields and unknown models (returns 0 in either case)."""
    if not usage:
        return 0.0
    rates = PRICING.get(model)
    if rates is None:
        _warn_once(model)
        return 0.0
    inp = int(usage.get("input_tokens", 0) or 0)
    out = int(usage.get("output_tokens", 0) or 0)
    cr = int(usage.get("cache_read_input_tokens", 0) or 0)
    cc = int(usage.get("cache_creation_input_tokens", 0) or 0)
    cost = 0.0
    cost += inp * rates.get("input", 0.0) / 1_000_000.0
    cost += out * rates.get("output", 0.0) / 1_000_000.0
    cost += cr * rates.get("cache_read", 0.0) / 1_000_000.0
    cost += cc * rates.get("cache_creation", 0.0) / 1_000_000.0
    return cost


@dataclass
class CostTracker:
    """Cumulative USD across a session. Lifetime: per session, not per turn,
    so the sidebar meter shows total spent for as long as the WS lives."""
    model: str
    total_usd: float = 0.0
    total_input: int = 0
    total_output: int = 0
    total_cache_read: int = 0
    total_cache_creation: int = 0
    turns: int = 0

    def add_turn(self, usage: dict | None) -> tuple[float, float]:
        """Add a turn's usage. Returns (turn_cost_usd, total_usd)."""
        if not usage:
            return (0.0, self.total_usd)
        delta = estimate_turn_cost(usage, self.model)
        self.total_usd += delta
        self.total_input += int(usage.get("input_tokens", 0) or 0)
        self.total_output += int(usage.get("output_tokens", 0) or 0)
        self.total_cache_read += int(usage.get("cache_read_input_tokens", 0) or 0)
        self.total_cache_creation += int(usage.get("cache_creation_input_tokens", 0) or 0)
        self.turns += 1
        return (delta, self.total_usd)

    def snapshot(self) -> dict:
        return {
            "model": self.model,
            "turns": self.turns,
            "total_usd": round(self.total_usd, 6),
            "total_input": self.total_input,
            "total_output": self.total_output,
            "total_cache_read": self.total_cache_read,
            "total_cache_creation": self.total_cache_creation,
            "total_tokens": self.total_input + self.total_output,
        }
