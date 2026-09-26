"""Agent session store for managing WebSocket agent sessions."""

from __future__ import annotations
import asyncio
import time
from typing import Dict, Optional
from dataclasses import dataclass, field

import anyio

from immunoscope.agent.session import Session
from immunoscope.agent.tool import ToolContext
from immunoscope.agent.config import get_settings, AGENT_AUDIT_DB_PATH


@dataclass
class AgentSessionState:
    """State for a single agent session."""

    session: Session
    context: ToolContext
    abort_event: anyio.Event
    created_at: float = field(default_factory=time.time)
    last_activity: float = field(default_factory=time.time)
    active_connections: int = 0

    def touch(self):
        """Update last activity timestamp."""
        self.last_activity = time.time()

    def is_expired(self, timeout: float = 3600.0) -> bool:
        """Check if session has expired."""
        return (time.time() - self.last_activity) > timeout


class AgentSessionStore:
    """Store for managing agent sessions across WebSocket connections."""

    def __init__(self):
        self._sessions: Dict[str, AgentSessionState] = {}
        self._lock = asyncio.Lock()
        self._cleanup_task: Optional[asyncio.Task] = None

    async def start(self):
        """Start background cleanup task."""
        if self._cleanup_task is None:
            self._cleanup_task = asyncio.create_task(self._cleanup_loop())

    async def stop(self):
        """Stop background cleanup task."""
        if self._cleanup_task:
            self._cleanup_task.cancel()
            try:
                await self._cleanup_task
            except asyncio.CancelledError:
                pass
            self._cleanup_task = None

    async def create_session(
        self,
        session_id: str,
        db_path: str = AGENT_AUDIT_DB_PATH,
        design_context: dict | None = None,
    ) -> AgentSessionState:
        """Create a new agent session.

        Args:
            session_id: Unique session identifier
            db_path: Audit DB path
            design_context: When set, session runs in Design Copilot mode
                            with MD evidence pre-loaded into system prompt
        """
        async with self._lock:
            if session_id in self._sessions:
                # Return existing session
                state = self._sessions[session_id]
                state.touch()
                state.active_connections += 1
                # Update design context if provided (allows fresh injection)
                if design_context is not None:
                    state.session.design_context = design_context
                return state

            # Create new session
            settings = get_settings()
            session = Session(id=session_id)
            if design_context:
                session.design_context = design_context
            abort_event = anyio.Event()

            # Event queue for WebSocket communication
            event_queue: asyncio.Queue = asyncio.Queue()

            async def on_event(event: dict):
                """Queue events for WebSocket transmission."""
                await event_queue.put(event)

            context = ToolContext(
                session_id=session_id,
                abort_event=abort_event,
                on_event=on_event,
                settings=settings,
                db_path=db_path,
            )

            # Store event queue in context for retrieval
            context.event_queue = event_queue

            state = AgentSessionState(
                session=session,
                context=context,
                abort_event=abort_event,
                active_connections=1,
            )

            self._sessions[session_id] = state
            return state

    async def get_session(self, session_id: str) -> Optional[AgentSessionState]:
        """Get an existing session."""
        async with self._lock:
            state = self._sessions.get(session_id)
            if state:
                state.touch()
            return state

    async def release_session(self, session_id: str):
        """Release a connection from a session."""
        async with self._lock:
            state = self._sessions.get(session_id)
            if state:
                state.active_connections -= 1
                state.touch()

    async def delete_session(self, session_id: str):
        """Delete a session."""
        async with self._lock:
            if session_id in self._sessions:
                state = self._sessions[session_id]
                # Abort any running operations
                state.abort_event.set()
                del self._sessions[session_id]

    async def list_sessions(self) -> list[dict]:
        """List all active sessions."""
        async with self._lock:
            return [
                {
                    "session_id": session_id,
                    "created_at": state.created_at,
                    "last_activity": state.last_activity,
                    "active_connections": state.active_connections,
                    "turn_count": len(state.session.messages) // 2,
                }
                for session_id, state in self._sessions.items()
            ]

    async def _cleanup_loop(self):
        """Background task to clean up expired sessions."""
        while True:
            try:
                await asyncio.sleep(300)  # Check every 5 minutes
                await self._cleanup_expired()
            except asyncio.CancelledError:
                break
            except Exception as e:
                print(f"Error in cleanup loop: {e}")

    async def _cleanup_expired(self):
        """Remove expired sessions."""
        async with self._lock:
            expired = [
                session_id
                for session_id, state in self._sessions.items()
                if state.active_connections == 0 and state.is_expired()
            ]

            for session_id in expired:
                state = self._sessions[session_id]
                state.abort_event.set()
                del self._sessions[session_id]
                print(f"Cleaned up expired session: {session_id}")


# Global session store instance
_store: Optional[AgentSessionStore] = None


def get_agent_store() -> AgentSessionStore:
    """Get the global agent session store."""
    global _store
    if _store is None:
        _store = AgentSessionStore()
    return _store
