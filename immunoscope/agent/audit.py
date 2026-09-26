from __future__ import annotations
import json
import logging
import time
from typing import Any

import aiosqlite

from immunoscope.agent.tool import ToolResult

log = logging.getLogger("prism.agent.audit")


CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS agent_audit (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    session_id TEXT NOT NULL,
    turn_idx INTEGER,
    tool_use_id TEXT,
    tool_name TEXT NOT NULL,
    args_json TEXT,
    result_summary TEXT,
    is_error INTEGER NOT NULL DEFAULT 0,
    error_class TEXT,
    duration_ms INTEGER,
    created_at REAL NOT NULL,
    -- Sprint A: per-turn cost telemetry (denormalized for fast aggregation)
    input_tokens INTEGER,
    output_tokens INTEGER,
    cache_read_tokens INTEGER,
    cache_creation_tokens INTEGER,
    cost_usd REAL
)
"""

CREATE_INDEX_SQL = (
    "CREATE INDEX IF NOT EXISTS idx_agent_audit_session ON agent_audit(session_id)"
)

# Idempotent migrations for existing DBs that pre-date the cost columns.
_MIGRATIONS = [
    "ALTER TABLE agent_audit ADD COLUMN input_tokens INTEGER",
    "ALTER TABLE agent_audit ADD COLUMN output_tokens INTEGER",
    "ALTER TABLE agent_audit ADD COLUMN cache_read_tokens INTEGER",
    "ALTER TABLE agent_audit ADD COLUMN cache_creation_tokens INTEGER",
    "ALTER TABLE agent_audit ADD COLUMN cost_usd REAL",
]


async def init_audit_db(db_path: str) -> None:
    async with aiosqlite.connect(db_path) as db:
        await db.execute(CREATE_TABLE_SQL)
        await db.execute(CREATE_INDEX_SQL)
        for stmt in _MIGRATIONS:
            try:
                await db.execute(stmt)
            except Exception as e:
                # SQLite raises OperationalError for "duplicate column name".
                # Any other failure also gets swallowed because schema additions
                # must never block the audit table from being usable.
                msg = str(e).lower()
                if "duplicate column" not in msg:
                    log.warning("audit migration failed: %s — %s", stmt, e)
        await db.commit()


def _summarize_result(result: ToolResult, max_chars: int = 500) -> str:
    if isinstance(result.content, str):
        s = result.content
    elif result.content is None:
        s = ""
    else:
        try:
            s = json.dumps(result.content, ensure_ascii=False, default=str)
        except Exception:
            s = repr(result.content)
    if len(s) > max_chars:
        return s[:max_chars] + f"...[+{len(s) - max_chars} chars]"
    return s


async def write_tool_call_audit(
    db_path: str,
    *,
    session_id: str,
    turn_idx: int,
    tool_use_id: str,
    tool_name: str,
    args: Any,
    result: ToolResult,
    duration_seconds: float,
    error: BaseException | None,
    # Optional per-turn telemetry. Only the FIRST tool of a turn carries
    # these; subsequent tool calls in the same turn pass None to avoid
    # double-counting in SUM aggregations.
    input_tokens: int | None = None,
    output_tokens: int | None = None,
    cache_read_tokens: int | None = None,
    cache_creation_tokens: int | None = None,
    cost_usd: float | None = None,
) -> None:
    """Persist one tool invocation. Failures are caller-handled (caller wraps
    this in asyncio.wait_for and logs warnings without aborting the agent loop)."""
    try:
        args_json = json.dumps(args, ensure_ascii=False, default=str)
    except Exception:
        args_json = repr(args)

    async with aiosqlite.connect(db_path) as db:
        await db.execute(
            """INSERT INTO agent_audit (
                session_id, turn_idx, tool_use_id, tool_name,
                args_json, result_summary, is_error, error_class,
                duration_ms, created_at,
                input_tokens, output_tokens, cache_read_tokens,
                cache_creation_tokens, cost_usd
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                session_id,
                turn_idx,
                tool_use_id,
                tool_name,
                args_json,
                _summarize_result(result),
                1 if result.is_error else 0,
                type(error).__name__ if error else None,
                int(duration_seconds * 1000),
                time.time(),
                input_tokens,
                output_tokens,
                cache_read_tokens,
                cache_creation_tokens,
                cost_usd,
            ),
        )
        await db.commit()
