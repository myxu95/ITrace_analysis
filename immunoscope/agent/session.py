from __future__ import annotations
import asyncio
import time
import uuid
from dataclasses import dataclass, field

import anyio

from immunoscope.agent.cost_tracker import CostTracker


@dataclass
class Session:
    id: str
    messages: list[dict] = field(default_factory=list)
    created_at: float = field(default_factory=time.time)
    abort_event: anyio.Event = field(default_factory=anyio.Event)
    pending_permissions: dict[str, asyncio.Future] = field(default_factory=dict)
    turn_idx: int = 0
    # Per-session cost tracker; created lazily by engine.run_turn on first turn
    # because we don't know the LLM model until then. Lives until session is dropped.
    cost_tracker: CostTracker | None = None
    # Design Copilot mode: when set, system prompt switches to design-focused
    # behavior and includes pre-loaded MD evidence summary.
    design_context: dict | None = None


class InProcessSessionStore:
    """Holds active sessions for the lifetime of the FastAPI process. A session
    only exists while a WebSocket is open; on disconnect routers/chat.py drops
    the session. Audit data persists separately to SQLite."""

    def __init__(self) -> None:
        self._sessions: dict[str, Session] = {}

    def create(self) -> Session:
        sid = uuid.uuid4().hex[:16]
        sess = Session(id=sid)
        self._sessions[sid] = sess
        return sess

    def get(self, session_id: str) -> Session | None:
        return self._sessions.get(session_id)

    def drop(self, session_id: str) -> None:
        sess = self._sessions.pop(session_id, None)
        if sess is None:
            return
        for fut in sess.pending_permissions.values():
            if not fut.done():
                fut.set_exception(asyncio.CancelledError())
        if not sess.abort_event.is_set():
            sess.abort_event.set()

    def __len__(self) -> int:
        return len(self._sessions)


session_store = InProcessSessionStore()
