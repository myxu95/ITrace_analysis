"""A2A message envelope and evidence reference types (v1, 2026-05-28).

Agents communicate by emitting `AgentMessage` instances into a shared
`MutationDesignContext.message_history`. The protocol is intentionally
lightweight:

  - `narrative`     : natural language (REQUIRED). The primary semantic
                      carrier between LLM agents.
  - `evidence_refs` : tagged-union references (REQUIRED, ≥1). Every
                      non-common-knowledge claim in `narrative` must be
                      backed by an inline `[ref:<id>]` token whose id
                      appears in this list. The consistency check is
                      enforced by `validate_message`.
  - `structured_payload` : optional dict, only filled when downstream
                      consumers do numerical computation (e.g. Writer
                      emits the 5-field recommendation JSON here).

Why both `[ref:<id>]` inline + side list?
    Inline tokens keep the prose readable and let downstream LLMs see
    "this claim is sourced from x" without parsing JSON. The side list
    carries the resolved evidence (PMID + quote, or a context path +
    excerpt) so an auditor can reach the source without re-fetching it.

The four `EvidenceRef` flavors:

    ContextFieldRef — points at a field path inside MutationDesignContext
                      (e.g. "feature_views.hotspots.top_pairs[3].rrcs_mean").
                      Carries a `quote` excerpt so audit doesn't have to
                      walk the context again.
    LiteratureRef   — PMID + title + relevant quote from the paper.
    SkillCallRef    — external skill / tool invocation, identified by
                      `call_id` matching the audit log; `result_summary`
                      is a one-line digest.
    MessageRef      — points at a prior AgentMessage by msg_id; used when
                      one agent's claim leans on another agent's prior
                      narrative (e.g. Interaction Reader cites Bio Agent).
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Tuple, Union


# ---------------------------------------------------------------------------
# Evidence reference tagged union
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ContextFieldRef:
    """Reference to a field in MutationDesignContext."""

    ref_id: str
    path: str
    quote: str
    kind: Literal["context_field"] = "context_field"


@dataclass(frozen=True)
class LiteratureRef:
    """Reference to a paper in the literature RAG."""

    ref_id: str
    pmid: str
    title: str
    relevant_quote: str
    kind: Literal["literature"] = "literature"


@dataclass(frozen=True)
class SkillCallRef:
    """Reference to an external skill / tool invocation."""

    ref_id: str
    skill_name: str
    call_id: str
    result_summary: str
    kind: Literal["skill_call"] = "skill_call"


@dataclass(frozen=True)
class MessageRef:
    """Reference to a prior AgentMessage by msg_id."""

    ref_id: str
    msg_id: str
    excerpt: str
    kind: Literal["message"] = "message"


EvidenceRef = Union[ContextFieldRef, LiteratureRef, SkillCallRef, MessageRef]


# ---------------------------------------------------------------------------
# Agent message envelope
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AgentMessage:
    """One A2A message from a producing agent to one or more consumers.

    Frozen so that downstream agents can hold a reference without worrying
    about upstream mutation. To "edit" a message, emit a new one.
    """

    msg_id: str
    sender: str
    recipients: Union[List[str], Literal["broadcast"]]
    narrative: str
    evidence_refs: List[EvidenceRef]
    structured_payload: Dict[str, Any] | None = None
    timestamp: float = field(default_factory=time.time)


class MessageValidationError(ValueError):
    """Raised when narrative <-> evidence_refs consistency check fails."""


_INLINE_REF_PATTERN = re.compile(r"\[ref:([^\]\s]+)\]")


def extract_inline_refs(narrative: str) -> List[str]:
    """Return ordered list of ref_ids that appear inline in `narrative`.

    Duplicate ids are kept (one inline occurrence per token) but the side
    list only needs to declare each id once.
    """
    return _INLINE_REF_PATTERN.findall(narrative)


def validate_message(msg: AgentMessage) -> None:
    """Enforce the narrative <-> evidence_refs contract.

    Rules:
      1. evidence_refs is non-empty.
      2. Every ref_id appearing inline as `[ref:<id>]` resolves in
         evidence_refs.
      3. Every ref_id in evidence_refs is unique.

    The reverse direction (every evidence_ref must appear inline) is
    NOT enforced — an agent may attach background refs that supported its
    reasoning even if no single sentence cites them directly.
    """
    if not msg.evidence_refs:
        raise MessageValidationError(
            f"Message {msg.msg_id} from {msg.sender} has empty evidence_refs; "
            "at least one ref is required."
        )

    declared = [r.ref_id for r in msg.evidence_refs]
    if len(declared) != len(set(declared)):
        seen: Dict[str, int] = {}
        for rid in declared:
            seen[rid] = seen.get(rid, 0) + 1
        dups = sorted(rid for rid, n in seen.items() if n > 1)
        raise MessageValidationError(
            f"Message {msg.msg_id} declares duplicate ref_ids: {dups}"
        )

    declared_set = set(declared)
    inline_ids = extract_inline_refs(msg.narrative)
    unresolved = sorted({rid for rid in inline_ids if rid not in declared_set})
    if unresolved:
        raise MessageValidationError(
            f"Message {msg.msg_id} has inline ref(s) {unresolved} that are not "
            f"declared in evidence_refs (declared: {sorted(declared_set)})"
        )


# ---------------------------------------------------------------------------
# Serialization helpers (JSON-friendly)
# ---------------------------------------------------------------------------


def evidence_ref_to_dict(ref: EvidenceRef) -> Dict[str, Any]:
    """Convert any EvidenceRef variant to a plain dict for JSON output."""
    if isinstance(ref, ContextFieldRef):
        return {
            "ref_id": ref.ref_id,
            "kind": ref.kind,
            "path": ref.path,
            "quote": ref.quote,
        }
    if isinstance(ref, LiteratureRef):
        return {
            "ref_id": ref.ref_id,
            "kind": ref.kind,
            "pmid": ref.pmid,
            "title": ref.title,
            "relevant_quote": ref.relevant_quote,
        }
    if isinstance(ref, SkillCallRef):
        return {
            "ref_id": ref.ref_id,
            "kind": ref.kind,
            "skill_name": ref.skill_name,
            "call_id": ref.call_id,
            "result_summary": ref.result_summary,
        }
    if isinstance(ref, MessageRef):
        return {
            "ref_id": ref.ref_id,
            "kind": ref.kind,
            "msg_id": ref.msg_id,
            "excerpt": ref.excerpt,
        }
    raise TypeError(f"Unknown EvidenceRef type: {type(ref).__name__}")


def evidence_ref_from_dict(data: Dict[str, Any]) -> EvidenceRef:
    """Reverse of `evidence_ref_to_dict`."""
    kind = data.get("kind")
    if kind == "context_field":
        return ContextFieldRef(
            ref_id=data["ref_id"], path=data["path"], quote=data["quote"]
        )
    if kind == "literature":
        return LiteratureRef(
            ref_id=data["ref_id"],
            pmid=data["pmid"],
            title=data["title"],
            relevant_quote=data["relevant_quote"],
        )
    if kind == "skill_call":
        return SkillCallRef(
            ref_id=data["ref_id"],
            skill_name=data["skill_name"],
            call_id=data["call_id"],
            result_summary=data["result_summary"],
        )
    if kind == "message":
        return MessageRef(
            ref_id=data["ref_id"], msg_id=data["msg_id"], excerpt=data["excerpt"]
        )
    raise ValueError(f"Unknown EvidenceRef kind: {kind!r}")


def message_to_dict(msg: AgentMessage) -> Dict[str, Any]:
    return {
        "msg_id": msg.msg_id,
        "sender": msg.sender,
        "recipients": msg.recipients,
        "narrative": msg.narrative,
        "evidence_refs": [evidence_ref_to_dict(r) for r in msg.evidence_refs],
        "structured_payload": msg.structured_payload,
        "timestamp": msg.timestamp,
    }


def message_from_dict(data: Dict[str, Any]) -> AgentMessage:
    return AgentMessage(
        msg_id=data["msg_id"],
        sender=data["sender"],
        recipients=data["recipients"],
        narrative=data["narrative"],
        evidence_refs=[evidence_ref_from_dict(r) for r in data["evidence_refs"]],
        structured_payload=data.get("structured_payload"),
        timestamp=data.get("timestamp", time.time()),
    )


__all__ = [
    "AgentMessage",
    "ContextFieldRef",
    "EvidenceRef",
    "LiteratureRef",
    "MessageRef",
    "MessageValidationError",
    "SkillCallRef",
    "evidence_ref_from_dict",
    "evidence_ref_to_dict",
    "extract_inline_refs",
    "message_from_dict",
    "message_to_dict",
    "validate_message",
]
