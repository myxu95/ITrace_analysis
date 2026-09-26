"""Numeric-fidelity post-hoc check (paper §4.5, 2026-06-02).

The protocol's `validate_message` only checks that every inline `[ref:<id>]`
token resolves in `evidence_refs` — it does NOT check whether the *numbers* a
reader puts in its narrative / `structured_payload` actually came back from a
tool it called. That is the largest credibility hole in the pipeline: an LLM
can cite a real ref id next to a hallucinated RRCS/SASA/occupancy value and
pass validation.

This module closes the hole with a deterministic, model-free check:

  1. Collect every numeric token an agent emitted (narrative + payload).
  2. Collect the full text of every tool result THAT agent actually received
     this run (from `SkillCallRecord.raw_result_text`).
  3. A number is *grounded* if it appears — within a small relative tolerance,
     accounting for rounding / unit suffixes — somewhere in that tool text.
  4. Report ungrounded numbers (the ones a reviewer must treat as suspect).

Design choices that keep false positives low:
  - "Structural" integers that are not measurements (residue numbers, list
    indices, years, PMIDs, ref ids) are common and usually DO appear in the
    tool output anyway; when they don't, they are almost always identifiers,
    not fabricated data. We therefore SUPPRESS bare small/large integers and
    only hard-check *decimals* (RRCS, SASA, occupancy, BSA, RMSF — the values
    the readers are forbidden to invent). Integers are still reported, but in
    a separate, advisory bucket.
  - Matching is tolerant: a payload `4.2` matches a tool `4.21` (rounding),
    and `0.87` matches `87%` (occupancy rendered as a percentage).

The check is advisory by default: `BaseAgent` logs the report and attaches it
to the run result. A strict mode (feed back to the model for self-correction)
is available but off by default so one over-eager false positive cannot wedge
the pipeline.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

# A decimal number, optionally signed, with optional exponent. We capture the
# raw string so we can re-derive both the float and its printed precision.
_NUMBER = re.compile(r"(?<![\w.])([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)(?![\w])")

# Relative tolerance for float matching (covers 2-3 sig-fig rounding).
_REL_TOL = 0.02
_ABS_TOL = 0.01

# Payload paths whose values are legitimately DERIVED by the agent (sums,
# aggregates) and so will not appear verbatim in any tool output. A
# substring match against the json-ish path. These are exempted from the
# hard decimal check (routed to advisory) so a correct computed value like
# `cross_partner_profile.cdr3a_vs_peptide_rrcs_sum` cannot wedge the run.
_DERIVED_PATH_SUBSTRINGS = ("_sum", "_total", "_count", "n_ref", "n_frames")


def _is_derived_path(where: str) -> bool:
    if not where.startswith("payload:"):
        return False
    low = where.lower()
    return any(s in low for s in _DERIVED_PATH_SUBSTRINGS)


@dataclass
class NumberCheck:
    """One emitted number and whether it was found in tool output."""

    value: float
    raw: str
    where: str           # "narrative" | "payload:<path>"
    grounded: bool
    is_integer: bool


@dataclass
class NumericAuditReport:
    """Outcome of auditing one agent message against its tool calls."""

    agent: str
    checked: int = 0
    ungrounded_decimals: List[NumberCheck] = field(default_factory=list)
    ungrounded_integers: List[NumberCheck] = field(default_factory=list)
    had_tool_text: bool = True

    @property
    def ok(self) -> bool:
        """A run passes the *hard* check when no decimal is ungrounded.

        Integers are advisory only (identifiers vs. measurements are hard to
        tell apart), so they do not fail the check."""
        return not self.ungrounded_decimals

    def summary(self) -> str:
        if not self.had_tool_text:
            return (
                f"[numeric-audit:{self.agent}] no tool output to check against "
                f"({self.checked} numbers emitted) — skipped."
            )
        if self.ok and not self.ungrounded_integers:
            return (
                f"[numeric-audit:{self.agent}] OK — all {self.checked} numbers "
                "trace to tool output."
            )
        parts = [f"[numeric-audit:{self.agent}] {self.checked} numbers checked"]
        if self.ungrounded_decimals:
            vals = ", ".join(
                f"{c.raw}({c.where})" for c in self.ungrounded_decimals[:8]
            )
            parts.append(f"UNGROUNDED DECIMALS: {vals}")
        if self.ungrounded_integers:
            vals = ", ".join(
                f"{c.raw}" for c in self.ungrounded_integers[:8]
            )
            parts.append(f"advisory integers not in tool text: {vals}")
        return " | ".join(parts)


def _iter_payload_numbers(
    obj: Any, path: str = ""
) -> Iterable[Tuple[float, str, str]]:
    """Yield (value, raw_str, json-ish path) for every number under `obj`.

    Numbers that are stored as JSON numbers (not strings) are yielded with
    their repr; numbers embedded in strings are extracted via the regex so a
    payload like ``"RRCS 4.21"`` is still checked."""
    if isinstance(obj, bool):
        return  # bool is an int subclass — never a measurement
    if isinstance(obj, (int, float)):
        yield (float(obj), _trim_float(obj), path or "<root>")
        return
    if isinstance(obj, str):
        for m in _NUMBER.finditer(obj):
            raw = m.group(1)
            try:
                yield (float(raw), raw, path or "<root>")
            except ValueError:
                continue
        return
    if isinstance(obj, dict):
        for k, v in obj.items():
            child = f"{path}.{k}" if path else str(k)
            yield from _iter_payload_numbers(v, child)
        return
    if isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            child = f"{path}[{i}]"
            yield from _iter_payload_numbers(v, child)
        return


def _trim_float(x: Any) -> str:
    if isinstance(x, int):
        return str(x)
    s = f"{x}"
    return s


def _looks_like_identifier(raw: str, value: float) -> bool:
    """Heuristic: is this integer an id / index rather than a measurement?

    We treat as identifiers: years (1900-2100), 4+ digit integers (PMIDs,
    resids occasionally), and the integers 0/1/2 that pepper prose. These are
    pushed to the advisory bucket, never the hard-fail bucket."""
    if "." in raw or "e" in raw.lower():
        return False
    iv = int(value)
    if 1900 <= iv <= 2100:
        return True
    if abs(iv) >= 1000:
        return True
    if abs(iv) <= 2:
        return True
    return False


def _build_haystack(tool_texts: Sequence[str]) -> Tuple[str, List[float]]:
    """Concatenate tool text and pre-extract its numeric values for matching."""
    joined = "\n".join(t for t in tool_texts if t)
    nums: List[float] = []
    for m in _NUMBER.finditer(joined):
        try:
            nums.append(float(m.group(1)))
        except ValueError:
            continue
    return joined, nums


def _matches(value: float, haystack_nums: Sequence[float], haystack_text: str) -> bool:
    """Is `value` present in the tool output (with rounding / % tolerance)?"""
    # 1) direct numeric proximity to any number the tools emitted
    for hv in haystack_nums:
        if _close(value, hv):
            return True
    # 2) occupancy rendered as a percentage: payload 0.87 vs tool "87%"
    if 0.0 < value < 1.0:
        pct = value * 100.0
        for hv in haystack_nums:
            if _close(pct, hv):
                return True
    # 3) integer percentage the other way: payload 87 vs tool 0.87
    if value > 1.0:
        frac = value / 100.0
        for hv in haystack_nums:
            if _close(frac, hv):
                return True
    return False


def _close(a: float, b: float) -> bool:
    if a == b:
        return True
    diff = abs(a - b)
    if diff <= _ABS_TOL:
        return True
    scale = max(abs(a), abs(b), 1e-9)
    return diff / scale <= _REL_TOL


def audit_numbers(
    agent: str,
    narrative: str,
    structured_payload: Optional[Dict[str, Any]],
    tool_texts: Sequence[str],
) -> NumericAuditReport:
    """Check every number in `narrative` + `structured_payload` against the
    text of the tool results this agent received.

    `tool_texts` is the list of `SkillCallRecord.raw_result_text` for the
    agent's successful (non-error) calls this run.
    """
    report = NumericAuditReport(agent=agent)
    haystack_text, haystack_nums = _build_haystack(tool_texts)
    report.had_tool_text = bool(haystack_text.strip())

    candidates: List[NumberCheck] = []

    # narrative numbers
    for m in _NUMBER.finditer(narrative or ""):
        raw = m.group(1)
        try:
            val = float(raw)
        except ValueError:
            continue
        candidates.append(
            NumberCheck(value=val, raw=raw, where="narrative",
                        grounded=False, is_integer="." not in raw)
        )

    # payload numbers (including numbers embedded in payload strings)
    if structured_payload:
        for val, raw, path in _iter_payload_numbers(structured_payload):
            candidates.append(
                NumberCheck(value=val, raw=raw, where=f"payload:{path}",
                            grounded=False, is_integer="." not in raw)
            )

    report.checked = len(candidates)
    if not report.had_tool_text:
        return report

    for c in candidates:
        c.grounded = _matches(c.value, haystack_nums, haystack_text)
        if c.grounded:
            continue
        if c.is_integer:
            # identifiers (years/PMIDs) and plain integers are advisory only —
            # hard to distinguish ids from measurements.
            report.ungrounded_integers.append(c)
        elif _is_derived_path(c.where):
            # a decimal the agent legitimately computed (a sum/aggregate) will
            # not appear verbatim in tool text — advisory, not a hard fail.
            report.ungrounded_integers.append(c)
        else:
            report.ungrounded_decimals.append(c)

    return report


def audit_agent_message(
    agent: str,
    narrative: str,
    structured_payload: Optional[Dict[str, Any]],
    skill_calls: Sequence[Any],
) -> NumericAuditReport:
    """Convenience wrapper: pull this agent's tool text out of `skill_calls`
    (a list of `SkillCallRecord`) and run `audit_numbers`."""
    texts: List[str] = []
    for rec in skill_calls:
        if getattr(rec, "agent", None) != agent:
            continue
        if getattr(rec, "is_error", False):
            continue
        txt = getattr(rec, "raw_result_text", None)
        if txt:
            texts.append(txt)
    return audit_numbers(agent, narrative, structured_payload, texts)


__all__ = [
    "NumberCheck",
    "NumericAuditReport",
    "audit_agent_message",
    "audit_numbers",
]
