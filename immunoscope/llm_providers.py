"""共享的 LLM provider 元数据。

recommendation/engine.py（同步批处理）和 agent/llm.py（异步流式）两套客户端
各自维护过一份"provider → env 变量名 / 默认模型 / 默认 base_url"的映射，
本模块把它们收敛到一张表，避免新增 provider 时漏改一边。

两套客户端的能力差异（同步 vs 异步、是否支持 tool_use）仍各自实现；
这里只做配置层的去重。
"""
from __future__ import annotations

import os
from typing import Optional, TypedDict


class ProviderConfig(TypedDict):
    env_var: str
    default_model: str
    default_base_url: Optional[str]


PROVIDER_DEFAULTS: dict[str, ProviderConfig] = {
    "deepseek": {
        "env_var": "DEEPSEEK_API_KEY",
        "default_model": "deepseek-chat",
        "default_base_url": "https://api.deepseek.com/v1",
    },
    "anthropic": {
        "env_var": "ANTHROPIC_API_KEY",
        "default_model": "claude-3-5-sonnet-20241022",
        "default_base_url": None,
    },
    "openai": {
        "env_var": "OPENAI_API_KEY",
        "default_model": "gpt-4-turbo-preview",
        "default_base_url": None,
    },
    "gemini": {
        "env_var": "GOOGLE_API_KEY",
        "default_model": "gemini-1.5-pro",
        "default_base_url": (
            "https://generativelanguage.googleapis.com/v1beta/openai/"
        ),
    },
}


def supported_providers() -> list[str]:
    """List known provider names, lowercase."""
    return sorted(PROVIDER_DEFAULTS.keys())


def _normalize(provider: str) -> str:
    name = (provider or "").lower()
    if name not in PROVIDER_DEFAULTS:
        raise ValueError(
            f"Unsupported LLM provider: {provider!r}. "
            f"Known: {', '.join(supported_providers())}."
        )
    return name


def resolve_api_key(provider: str, api_key: Optional[str] = None) -> str:
    """Return ``api_key`` if given, else read from the provider's env var.

    Raises ValueError if neither is available. Callers that pre-resolve via
    a Settings object (agent/llm.py) can skip this and just call the
    underlying client; this helper exists for the env-driven sync path.
    """
    name = _normalize(provider)
    if api_key:
        return api_key
    env_var = PROVIDER_DEFAULTS[name]["env_var"]
    key = os.environ.get(env_var)
    if not key:
        raise ValueError(
            f"API key not found for {name}: set environment variable {env_var} "
            f"or pass api_key explicitly."
        )
    return key


def default_model(provider: str) -> str:
    return PROVIDER_DEFAULTS[_normalize(provider)]["default_model"]


def default_base_url(provider: str) -> Optional[str]:
    return PROVIDER_DEFAULTS[_normalize(provider)]["default_base_url"]
