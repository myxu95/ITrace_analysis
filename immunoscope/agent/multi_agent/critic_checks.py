"""Deterministic, model-free critic checks (2026-06-02).

The LLM `DesignCriticAgent` and the `RecommendationAgent` are the same model
(DeepSeek), so a flaw both share — a number both misread, a constraint both
forget — can slip past even an adversarial review. This module removes the
*mechanical* checks from the LLM entirely: they are pure Python over the
structured payloads, so no model blind spot can reach them.

Two families:

  1. Hard-constraint violations (severity ``blocking``)
     - a germline-conserved framework residue (Bio conservation score ≥ 0.8)
       recommended for mutation (design_hint_class ≠ ``preserve``).

  2. Coverage gaps (severity ``major``)
     - a top-ranked contact hotspot (top-3 by RRCS in the Interaction brief)
       that is neither recommended nor explicitly skipped.
     - a *strong* cavity candidate (Interface brief, severity ``strong``)
       that is neither recommended nor skipped.

Residues are matched across briefs and the recommendation by their
canonical ``(chain, resid)`` key (see ``residue_id``), so the same physical
residue written in different label dialects still lines up.

A ``blocking`` finding here is meant to FORCE a revision regardless of what
the LLM critic concluded — that is the part the shared-model blind spot can
no longer defeat.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set, Tuple

from .context import MutationDesignContext
from .residue_id import ChainAliasMap, canonical_key

_CONSERVED_MIN = 0.80
_TOP_HOTSPOTS = 3

Key = Tuple[str, int]


@dataclass(frozen=True)
class DeterministicFinding:
    issue_type: str   # "constraint_violation" | "missed_signal"
    severity: str     # "blocking" | "major"
    target: str       # residue label / rec_id
    detail: str

    def line(self) -> str:
        return f"[{self.severity}] {self.issue_type} @ {self.target}: {self.detail}"


def _latest_payload(ctx: MutationDesignContext, sender: str) -> Dict[str, Any]:
    msgs = [m for m in ctx.message_history if m.sender == sender]
    if not msgs:
        return {}
    return msgs[-1].structured_payload or {}


def _keys_from(items, alias: ChainAliasMap, field: str = "residue") -> Dict[Key, str]:
    """Map canonical key -> original label for a list of dicts with `field`."""
    out: Dict[Key, str] = {}
    for it in items or []:
        if not isinstance(it, dict):
            continue
        label = it.get(field)
        if not label:
            continue
        k = canonical_key(str(label), alias)
        if k is not None and k not in out:
            out[k] = str(label)
    return out


def _recommended_keys(rec: Dict[str, Any], alias: ChainAliasMap) -> Dict[Key, dict]:
    out: Dict[Key, dict] = {}
    for r in rec.get("recommendations", []) or []:
        if not isinstance(r, dict):
            continue
        k = canonical_key(str(r.get("residue", "")), alias)
        if k is not None:
            out[k] = r
    return out


def _skipped_keys(rec: Dict[str, Any], alias: ChainAliasMap) -> Set[Key]:
    out: Set[Key] = set()
    for s in rec.get("skipped_candidates", []) or []:
        if isinstance(s, dict):
            k = canonical_key(str(s.get("residue", "")), alias)
            if k is not None:
                out.add(k)
    return out


def _is_mutation(rec_entry: dict) -> bool:
    """A recommendation targets a *mutation* unless its direction is preserve."""
    hint = ((rec_entry.get("handoff_to_physchem") or {}).get("design_hint_class")
            or rec_entry.get("design_hint_class") or "").lower()
    role = (rec_entry.get("role") or "").lower()
    return "preserve" not in hint and role != "preserve"


def check_conserved_framework(
    rec: Dict[str, Any], bio: Dict[str, Any], alias: ChainAliasMap
) -> List[DeterministicFinding]:
    conserved = _keys_from(
        ((bio.get("conservation") or {}).get("highly_conserved")), alias
    )
    # also fold in any score>=cutoff entries that carry an explicit score
    findings: List[DeterministicFinding] = []
    for key, entry in _recommended_keys(rec, alias).items():
        if key in conserved and _is_mutation(entry):
            findings.append(
                DeterministicFinding(
                    issue_type="constraint_violation",
                    severity="blocking",
                    target=str(entry.get("residue", key)),
                    detail=(
                        f"recommended for mutation but Bio marks {conserved[key]} "
                        f"germline-conserved (score ≥ {_CONSERVED_MIN}) — a "
                        "fold-load-bearing framework position."
                    ),
                )
            )
    return findings


def check_missed_top_hotspots(
    rec: Dict[str, Any], interaction: Dict[str, Any], alias: ChainAliasMap
) -> List[DeterministicFinding]:
    hotspots = [h for h in (interaction.get("hotspots") or []) if isinstance(h, dict)]
    hotspots.sort(key=lambda h: _as_float(h.get("rrcs_mean")), reverse=True)
    covered = set(_recommended_keys(rec, alias)) | _skipped_keys(rec, alias)
    findings: List[DeterministicFinding] = []
    for h in hotspots[:_TOP_HOTSPOTS]:
        k = canonical_key(str(h.get("residue", "")), alias)
        if k is None or k in covered:
            continue
        findings.append(
            DeterministicFinding(
                issue_type="missed_signal",
                severity="major",
                target=str(h.get("residue")),
                detail=(
                    f"top-{_TOP_HOTSPOTS} hotspot (RRCS "
                    f"{h.get('rrcs_mean')}) neither recommended nor skipped."
                ),
            )
        )
    return findings


def check_missed_strong_cavities(
    rec: Dict[str, Any], interface: Dict[str, Any], alias: ChainAliasMap
) -> List[DeterministicFinding]:
    covered = set(_recommended_keys(rec, alias)) | _skipped_keys(rec, alias)
    findings: List[DeterministicFinding] = []
    for c in interface.get("cavity_candidates", []) or []:
        if not isinstance(c, dict) or str(c.get("severity", "")).lower() != "strong":
            continue
        k = canonical_key(str(c.get("residue", "")), alias)
        if k is None or k in covered:
            continue
        findings.append(
            DeterministicFinding(
                issue_type="missed_signal",
                severity="major",
                target=str(c.get("residue")),
                detail="strong cavity candidate neither recommended nor skipped.",
            )
        )
    return findings


def _as_float(x: Any) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return -1.0


def check_empty_recommendation(rec: Dict[str, Any]) -> List[DeterministicFinding]:
    """A recommendation that carries zero sites is a blocking emptiness.

    The agent occasionally emits a well-formed message whose ``recommendations``
    list is empty (or emits an empty payload entirely). That is not a valid
    design — not even a ``preserve`` call — yet nothing downstream flags it, so
    the LLM critic can rubber-stamp it. Flag it ``blocking`` so the orchestrator
    forces a revision (and, if it stays empty, fails the run) instead of
    silently approving an empty plan.
    """
    if not (rec or {}).get("recommendations"):
        return [
            DeterministicFinding(
                issue_type="empty_plan",
                severity="blocking",
                target="(plan)",
                detail="recommendation produced zero actionable sites "
                       "(no mutations and no preserve calls).",
            )
        ]
    return []


def run_deterministic_checks(
    ctx: MutationDesignContext, alias: Optional[ChainAliasMap] = None
) -> List[DeterministicFinding]:
    """Run all model-free checks against the latest recommendation in `ctx`."""
    alias = alias or ChainAliasMap()
    rec = _latest_payload(ctx, "recommendation_agent")
    # Emptiness is checked first and unconditionally — an empty/missing payload
    # must not short-circuit past it (that was how empty plans slipped through).
    findings: List[DeterministicFinding] = list(check_empty_recommendation(rec))
    if not rec:
        return findings
    bio = _latest_payload(ctx, "bio_agent")
    interaction = _latest_payload(ctx, "interaction_reader")
    interface = _latest_payload(ctx, "interface_exposure_reader")
    findings += check_conserved_framework(rec, bio, alias)
    findings += check_missed_top_hotspots(rec, interaction, alias)
    findings += check_missed_strong_cavities(rec, interface, alias)
    return findings


def summarize(findings: List[DeterministicFinding]) -> Optional[str]:
    if not findings:
        return None
    return "\n".join(f.line() for f in findings)


def has_blocking(findings: List[DeterministicFinding]) -> bool:
    return any(f.severity == "blocking" for f in findings)


__all__ = [
    "DeterministicFinding",
    "run_deterministic_checks",
    "check_empty_recommendation",
    "check_conserved_framework",
    "check_missed_top_hotspots",
    "check_missed_strong_cavities",
    "summarize",
    "has_blocking",
]
