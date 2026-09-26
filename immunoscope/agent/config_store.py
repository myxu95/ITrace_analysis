"""Persistent LLM config: ~/.immunoscope/config.json.

Priority: config file > env > built-in default. UI writes take effect
immediately on the next call to `get_settings()`.
"""

from __future__ import annotations

import copy
import json
import logging
import os
from pathlib import Path
from typing import Any

log = logging.getLogger("immunoscope.agent.config_store")

CONFIG_DIR = Path.home() / ".immunoscope"
CONFIG_FILE = CONFIG_DIR / "config.json"

ALLOWED_PROVIDERS = ("deepseek", "anthropic", "openai", "gemini")

DEFAULTS: dict[str, Any] = {
    "llm": {
        "provider": "deepseek",
        "model": "",
        "base_url": "",
        "api_key": "",
        "max_tokens": 4096,
        "temperature": 0.7,
    }
}


def _provider_default_model(provider: str) -> str:
    if provider == "anthropic":
        return "claude-opus-4-7"
    if provider == "openai":
        return "gpt-5.5"
    if provider == "gemini":
        return "gemini-3.1-pro"
    return "deepseek-chat"


def _provider_default_base_url(provider: str) -> str:
    if provider == "deepseek":
        return "https://api.deepseek.com/v1"
    if provider == "gemini":
        return "https://generativelanguage.googleapis.com/v1beta/openai/"
    return ""


def _provider_env_key_name(provider: str) -> str:
    return {
        "anthropic": "ANTHROPIC_API_KEY",
        "openai": "OPENAI_API_KEY",
        "deepseek": "DEEPSEEK_API_KEY",
        "gemini": "GEMINI_API_KEY",
    }.get(provider, "DEEPSEEK_API_KEY")


def load_config() -> dict[str, Any]:
    """Read config.json, merging into DEFAULTS. Missing file returns DEFAULTS copy."""
    cfg = copy.deepcopy(DEFAULTS)
    if not CONFIG_FILE.exists():
        return cfg
    try:
        raw = json.loads(CONFIG_FILE.read_text())
    except (json.JSONDecodeError, OSError) as e:
        log.warning("Failed to parse %s, falling back to defaults: %s", CONFIG_FILE, e)
        return cfg
    if isinstance(raw, dict):
        for section, values in raw.items():
            if section in cfg and isinstance(values, dict):
                cfg[section].update({k: v for k, v in values.items() if k in cfg[section]})
    return cfg


def save_config(updates: dict[str, Any]) -> None:
    """Merge `updates` into existing config and persist. Sets mode 0600."""
    cfg = load_config()
    for section, values in updates.items():
        if section not in cfg:
            cfg[section] = {}
        if isinstance(values, dict):
            cfg[section].update(values)
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    CONFIG_FILE.write_text(json.dumps(cfg, indent=2))
    try:
        os.chmod(CONFIG_FILE, 0o600)
    except OSError:
        log.warning("Failed to chmod 600 on %s", CONFIG_FILE)


def _resolve_api_key(file_key: str, provider: str) -> tuple[str, str]:
    """Return (api_key, source). Source ∈ file|env_explicit|env_provider|none."""
    if file_key:
        return file_key, "file"
    explicit = os.getenv("IMMUNOSCOPE_LLM_API_KEY", "")
    if explicit:
        return explicit, "env_explicit"
    env_name = _provider_env_key_name(provider)
    provider_key = os.getenv(env_name, "")
    if provider_key:
        return provider_key, "env_provider"
    return "", "none"


def effective_llm_config() -> dict[str, Any]:
    """Merge priority: file > env > built-in default.

    Returns final {provider, model, base_url, api_key, max_tokens, temperature,
    key_source}. Empty model/base_url fall back to provider default.
    """
    cfg = load_config()["llm"]

    provider = cfg["provider"] or os.getenv("IMMUNOSCOPE_LLM_PROVIDER", "") or DEFAULTS["llm"]["provider"]
    if provider not in ALLOWED_PROVIDERS:
        provider = DEFAULTS["llm"]["provider"]

    model = cfg["model"] or os.getenv("IMMUNOSCOPE_LLM_MODEL", "") or _provider_default_model(provider)
    base_url = cfg["base_url"] or os.getenv("IMMUNOSCOPE_LLM_BASE_URL", "") or _provider_default_base_url(provider)
    api_key, key_source = _resolve_api_key(cfg.get("api_key", ""), provider)

    max_tokens = cfg.get("max_tokens") or DEFAULTS["llm"]["max_tokens"]
    temperature = cfg.get("temperature")
    if temperature is None:
        temperature = DEFAULTS["llm"]["temperature"]

    return {
        "provider": provider,
        "model": model,
        "base_url": base_url,
        "api_key": api_key,
        "max_tokens": int(max_tokens),
        "temperature": float(temperature),
        "key_source": key_source,
    }


def _mask_last4(key: str) -> str:
    if not key:
        return ""
    if len(key) <= 4:
        return "*" * len(key)
    return key[-4:]


def llm_status() -> dict[str, Any]:
    """Snapshot for GET /api/settings/llm."""
    eff = effective_llm_config()
    file_cfg = load_config()["llm"]
    return {
        "active": {
            "provider": eff["provider"],
            "model": eff["model"],
            "base_url": eff["base_url"],
            "max_tokens": eff["max_tokens"],
            "temperature": eff["temperature"],
            "has_key": bool(eff["api_key"]),
            "key_last4": _mask_last4(eff["api_key"]),
            "key_source": eff["key_source"],
        },
        "file": {
            "provider": file_cfg["provider"],
            "model": file_cfg["model"],
            "base_url": file_cfg["base_url"],
            "max_tokens": file_cfg["max_tokens"],
            "temperature": file_cfg["temperature"],
            "has_file_key": bool(file_cfg.get("api_key")),
        },
        "env": {
            "explicit_set": bool(os.getenv("IMMUNOSCOPE_LLM_API_KEY")),
            "deepseek": bool(os.getenv("DEEPSEEK_API_KEY")),
            "anthropic": bool(os.getenv("ANTHROPIC_API_KEY")),
            "openai": bool(os.getenv("OPENAI_API_KEY")),
            "gemini": bool(os.getenv("GEMINI_API_KEY")),
            "provider_override": os.getenv("IMMUNOSCOPE_LLM_PROVIDER", ""),
            "model_override": os.getenv("IMMUNOSCOPE_LLM_MODEL", ""),
            "base_url_override": os.getenv("IMMUNOSCOPE_LLM_BASE_URL", ""),
        },
        "config_file": str(CONFIG_FILE),
        "config_exists": CONFIG_FILE.exists(),
    }
