"""WebSocket router for agent chat."""

from __future__ import annotations
import asyncio
import json
import uuid
from typing import Optional

from fastapi import APIRouter, WebSocket, WebSocketDisconnect, Query
from fastapi.responses import JSONResponse

from immunoscope.agent.engine import run_turn
from immunoscope.agent.exceptions import (
    AgentError,
    MaxTurnsExceeded,
    WallClockExceeded,
)
from immunoscope.web.services.agent_store import get_agent_store


router = APIRouter(prefix="/agent", tags=["agent"])


@router.websocket("/ws")
async def agent_websocket(
    websocket: WebSocket,
    session_id: Optional[str] = Query(None),
    design_task_id: Optional[str] = Query(None),
):
    """WebSocket endpoint for agent chat.

    When `design_task_id` is provided, the session is created with Design
    Copilot context (MD evidence + design goals injected into system prompt).

    Protocol:
    - Client sends: {"type": "message", "content": "user message"}
    - Client sends: {"type": "abort"} to cancel current operation
    - Server sends: {"type": "text_delta", "text": "..."}
    - Server sends: {"type": "tool_use_start", "name": "tool_name", "id": "..."}
    - Server sends: {"type": "tool_use_end", "name": "tool_name", "ok": true}
    - Server sends: {"type": "error", "message": "..."}
    - Server sends: {"type": "session_info", "session_id": "...", "turn_count": 0}
    - Server sends: {"type": "draft_saved", "task_id": "...", "n_drafts": N} (design mode)
    """
    await websocket.accept()

    # Generate or use provided session ID
    if not session_id:
        session_id = f"ws_{uuid.uuid4().hex[:16]}"

    # Resolve design context if this is a Design Copilot session
    design_context = None
    if design_task_id:
        try:
            from immunoscope.web.routers.design import _load_task, _extract_design_context
            from pathlib import Path
            task = _load_task(design_task_id)
            if task and task.get("source_run_dir"):
                design_context = _extract_design_context(
                    Path(task["source_run_dir"]),
                    task.get("design_goals", []),
                    task.get("constraints", {}),
                    design_task_id,
                )
        except Exception as e:
            await websocket.send_json({
                "type": "error",
                "message": f"Failed to load design context: {e}",
            })
            await websocket.close()
            return

    store = get_agent_store()
    state = None

    try:
        # Create or get session (with design context if applicable)
        state = await store.create_session(session_id, design_context=design_context)

        # Send session info
        await websocket.send_json({
            "type": "session_info",
            "session_id": session_id,
            "turn_count": len(state.session.messages) // 2,
        })

        # Start event forwarding task
        event_task = asyncio.create_task(
            _forward_events(websocket, state.context.event_queue)
        )

        # Main message loop
        while True:
            try:
                # Receive message from client
                data = await websocket.receive_text()
                message = json.loads(data)

                msg_type = message.get("type")

                if msg_type == "message":
                    # User message - run agent turn
                    user_content = message.get("content", "").strip()
                    if not user_content:
                        continue

                    try:
                        # Run agent turn (events sent via event_queue)
                        await run_turn(state.session, user_content, state.context)

                        # Send turn complete
                        await websocket.send_json({
                            "type": "turn_complete",
                            "turn_count": len(state.session.messages) // 2,
                        })

                    except MaxTurnsExceeded:
                        await websocket.send_json({
                            "type": "error",
                            "message": "Maximum turns exceeded. Session ended.",
                        })
                        break

                    except WallClockExceeded:
                        await websocket.send_json({
                            "type": "error",
                            "message": "Session time limit exceeded.",
                        })
                        break

                    except AgentError as e:
                        await websocket.send_json({
                            "type": "error",
                            "message": str(e),
                        })

                elif msg_type == "abort":
                    # Abort current operation
                    state.abort_event.set()
                    await websocket.send_json({
                        "type": "aborted",
                        "message": "Operation aborted by user.",
                    })
                    # Reset abort event for next turn
                    state.abort_event = asyncio.Event()
                    state.context.abort_event = state.abort_event

                elif msg_type == "ping":
                    # Heartbeat
                    await websocket.send_json({"type": "pong"})

                else:
                    await websocket.send_json({
                        "type": "error",
                        "message": f"Unknown message type: {msg_type}",
                    })

            except json.JSONDecodeError:
                await websocket.send_json({
                    "type": "error",
                    "message": "Invalid JSON message",
                })

    except WebSocketDisconnect:
        pass

    except Exception as e:
        try:
            await websocket.send_json({
                "type": "error",
                "message": f"Server error: {str(e)}",
            })
        except:
            pass

    finally:
        # Cleanup
        if state:
            await store.release_session(session_id)

        # Cancel event forwarding task
        if 'event_task' in locals():
            event_task.cancel()
            try:
                await event_task
            except asyncio.CancelledError:
                pass

        try:
            await websocket.close()
        except:
            pass


async def _forward_events(websocket: WebSocket, event_queue: asyncio.Queue):
    """Forward events from queue to WebSocket."""
    while True:
        try:
            event = await event_queue.get()
            await websocket.send_json(event)
        except asyncio.CancelledError:
            break
        except Exception as e:
            print(f"Error forwarding event: {e}")
            break


@router.get("/sessions")
async def list_sessions():
    """List all active agent sessions."""
    store = get_agent_store()
    sessions = await store.list_sessions()
    return {"sessions": sessions}


@router.delete("/sessions/{session_id}")
async def delete_session(session_id: str):
    """Delete an agent session."""
    store = get_agent_store()
    await store.delete_session(session_id)
    return {"message": f"Session {session_id} deleted"}


@router.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "ok", "service": "agent"}
