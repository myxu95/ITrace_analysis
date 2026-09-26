"""User settings API for the ImmunoScope web app.

P0 scope: LLM configuration only (provider / model / base_url / api_key /
max_tokens / temperature) + test connection. Other settings live elsewhere
or will be added in P1.
"""

from __future__ import annotations

import time
import logging
from typing import Any, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from immunoscope.agent.config import get_settings
from immunoscope.agent.config_store import ALLOWED_PROVIDERS, llm_status, save_config
from immunoscope.agent.llm import make_llm_client
from immunoscope.agent.stream_events import MessageStop, TextDelta

log = logging.getLogger("immunoscope.web.settings")

router = APIRouter(prefix="/settings", tags=["settings"])


class LLMSettingsUpdate(BaseModel):
    provider: Optional[str] = None
    model: Optional[str] = None
    base_url: Optional[str] = None
    api_key: Optional[str] = None  # "" clears (fallback to env); None leaves unchanged
    max_tokens: Optional[int] = Field(None, ge=128, le=32768)
    temperature: Optional[float] = Field(None, ge=0.0, le=2.0)


@router.get("/llm")
def get_llm_settings() -> dict[str, Any]:
    return llm_status()


@router.put("/llm")
def update_llm_settings(payload: LLMSettingsUpdate) -> dict[str, Any]:
    if payload.provider is not None and payload.provider not in ALLOWED_PROVIDERS:
        raise HTTPException(
            status_code=400,
            detail=f"provider must be one of {ALLOWED_PROVIDERS}",
        )
    updates = {k: v for k, v in payload.model_dump().items() if v is not None}
    if updates:
        save_config({"llm": updates})
    return llm_status()


@router.post("/llm/test")
async def test_llm_connection() -> dict[str, Any]:
    settings = get_settings()
    if not settings.LLM_API_KEY:
        raise HTTPException(
            status_code=400,
            detail="No API key configured (neither in ~/.immunoscope/config.json nor in environment).",
        )

    client = make_llm_client(settings)
    t0 = time.monotonic()
    try:
        collected = ""
        async for event in client.stream(
            system="You are a connectivity test.",
            messages=[{"role": "user", "content": "Reply with exactly: ok"}],
            tools=None,
            max_tokens=8,
        ):
            if isinstance(event, TextDelta):
                collected += event.text
            elif isinstance(event, MessageStop):
                break
        latency_ms = int((time.monotonic() - t0) * 1000)
        return {
            "ok": True,
            "latency_ms": latency_ms,
            "provider": settings.LLM_PROVIDER,
            "model": settings.LLM_MODEL,
            "base_url": settings.LLM_BASE_URL,
            "response": collected.strip()[:80],
        }
    except Exception as e:
        log.warning("LLM connectivity test failed: %s", e)
        return {
            "ok": False,
            "error": str(e)[:300],
            "provider": settings.LLM_PROVIDER,
            "model": settings.LLM_MODEL,
            "base_url": settings.LLM_BASE_URL,
        }
    finally:
        if hasattr(client, "aclose"):
            try:
                await client.aclose()
            except Exception:
                pass
