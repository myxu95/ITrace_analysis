"""Shared `MutationDesignContext` for the multi-agent pipeline (v1, 2026-05-28).

The four round-1 agents (Bio, Interaction, Conformation, Interface &
Exposure) read and write a single `MutationDesignContext` instance
during a run. The contract:

  - The orchestrator constructs it once at the top of a run with the
    user-facing inputs (`case_dir`, `case_id`, `design_intent`).
  - Each agent appends its produced `AgentMessage` via `add_message(msg)`
    and stashes the returned msg_id in the per-agent slot
    (`bio_msg_id`, `interaction_msg_id`, ...). The slots are an
    indexing convenience; the canonical record is `message_history`.
  - Any tool / skill invocation is logged into `skill_calls` so the audit
    trail is reachable from a single object.

The four round-1 agents run in parallel and do **not** read each other's
messages — the slots are populated as agents finish, but a later
integration stage (deferred) is what consumes all four briefs together.

The context is intentionally mutable (unlike `AgentMessage`) because the
orchestrator threads it through every stage and each stage extends it.
Agents must NOT mutate prior `AgentMessage` entries — to "edit" a message
they emit a new one and (optionally) reference the old one via MessageRef.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from .messaging import AgentMessage


# ---------------------------------------------------------------------------
# Skill / tool invocation audit record
# ---------------------------------------------------------------------------


@dataclass
class SkillCallRecord:
    """One external tool / skill invocation made by an agent.

    Kept separate from `AgentMessage.evidence_refs` because a single message
    may cite multiple skill calls, and the same call may be referenced by
    multiple downstream messages. `call_id` is the join key used by
    `SkillCallRef.call_id`.
    """

    call_id: str
    agent: str
    skill_name: str
    arguments: Dict[str, Any]
    result_summary: str
    raw_result_path: Optional[str] = None
    # Full stringified tool output, kept in-memory only (not persisted). The
    # numeric-fidelity audit (`numeric_audit.py`) re-reads this to confirm that
    # numbers an agent put in its narrative / structured_payload actually came
    # back from a tool it called. `result_summary` stays a short audit digest.
    raw_result_text: Optional[str] = None
    is_error: bool = False


# ---------------------------------------------------------------------------
# Shared mutation-design context
# ---------------------------------------------------------------------------


@dataclass
class MutationDesignContext:
    """Single-system mutation design state, shared across all 4 round-1 agents.

    Fields fall into three groups:
      1. Inputs (set by orchestrator at construction): `case_dir`,
         `case_id`, `design_intent`.
      2. Per-agent message pointers (set as each agent finishes):
         `bio_msg_id`, `interaction_msg_id`, `conformation_msg_id`,
         `interface_exposure_msg_id` (round-1) and `recommendation_msg_id`
         (round-2). These are convenience handles into `message_history`;
         the canonical record is still the list.
      3. Cumulative state: `message_history` (every emitted
         `AgentMessage`, in order) and `skill_calls` (audit log).
      4. Preambles: `system_preamble` holds the pre-rendered
         `overview` + `quality` view markdown injected into every round-1
         agent. `recommendation_preamble` holds the orchestrator's
         deterministically pre-aggregated candidate table (built from the
         four round-1 `structured_payload`s) that is injected into the
         round-2 RecommendationAgent.
    """

    case_dir: str
    case_id: str
    design_intent: str

    bio_msg_id: Optional[str] = None
    interaction_msg_id: Optional[str] = None
    conformation_msg_id: Optional[str] = None
    interface_exposure_msg_id: Optional[str] = None
    recommendation_msg_id: Optional[str] = None
    critique_msg_id: Optional[str] = None

    message_history: List[AgentMessage] = field(default_factory=list)
    skill_calls: List[SkillCallRecord] = field(default_factory=list)
    system_preamble: Optional[str] = None
    recommendation_preamble: Optional[str] = None
    # Deterministic numeric-fidelity flags from round-1 (ungrounded decimals
    # per reader). Injected into the DesignCritic so it can escalate a
    # fabricated number the recommendation relies on to a blocking finding.
    numeric_audit_summary: Optional[str] = None
    # Deterministic (model-free) constraint/coverage findings on the current
    # recommendation (see `critic_checks`). A blocking entry here FORCES a
    # revision regardless of the LLM critic's verdict — the shared-model blind
    # spot cannot reach a check the model never makes. Injected into both the
    # critic and the revision payloads.
    deterministic_findings_summary: Optional[str] = None

    # ------------------------------------------------------------------
    # Message helpers
    # ------------------------------------------------------------------

    def add_message(self, msg: AgentMessage) -> str:
        """Append `msg` to history and return its msg_id.

        Raises `ValueError` if a message with the same msg_id already
        exists — msg_ids are the audit join key and must stay unique.
        """
        for existing in self.message_history:
            if existing.msg_id == msg.msg_id:
                raise ValueError(
                    f"msg_id collision: {msg.msg_id!r} already present "
                    f"(sender={existing.sender!r})"
                )
        self.message_history.append(msg)
        return msg.msg_id

    def get_message(self, msg_id: str) -> AgentMessage:
        """Return the message with the given msg_id, or raise KeyError."""
        for msg in self.message_history:
            if msg.msg_id == msg_id:
                return msg
        raise KeyError(f"No message with msg_id={msg_id!r} in history")

    def messages_from(self, sender: str) -> List[AgentMessage]:
        """Return all messages emitted by `sender`, in original order."""
        return [m for m in self.message_history if m.sender == sender]

    # ------------------------------------------------------------------
    # Skill-call helpers
    # ------------------------------------------------------------------

    def add_skill_call(self, record: SkillCallRecord) -> str:
        """Append `record` to the skill-call audit log and return its call_id."""
        for existing in self.skill_calls:
            if existing.call_id == record.call_id:
                raise ValueError(
                    f"call_id collision: {record.call_id!r} already present "
                    f"(agent={existing.agent!r}, skill={existing.skill_name!r})"
                )
        self.skill_calls.append(record)
        return record.call_id

    def get_skill_call(self, call_id: str) -> SkillCallRecord:
        """Return the skill-call record with the given call_id, or raise KeyError."""
        for rec in self.skill_calls:
            if rec.call_id == call_id:
                return rec
        raise KeyError(f"No skill call with call_id={call_id!r}")


__all__ = [
    "MutationDesignContext",
    "SkillCallRecord",
]
