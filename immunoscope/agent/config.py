"""Agent configuration for ImmunoScope.

LLM-related fields are resolved through `immunoscope.agent.config_store`
which merges (file > env > built-in default). `get_settings()` returns a
fresh `Settings` instance each call so UI updates take effect without a
server restart.
"""

from __future__ import annotations
import os
from typing import Literal

from immunoscope.agent.config_store import (
    effective_llm_config,
    _provider_default_model as _default_model,
)


def _bootstrap_llm() -> dict:
    """Snapshot effective LLM config once at import time for module-level constants."""
    try:
        return effective_llm_config()
    except Exception:
        return {
            "provider": os.getenv("IMMUNOSCOPE_LLM_PROVIDER", "deepseek"),
            "model": "deepseek-chat",
            "base_url": "https://api.deepseek.com/v1",
            "api_key": os.getenv("DEEPSEEK_API_KEY", ""),
            "max_tokens": 4096,
            "temperature": 0.7,
            "key_source": "none",
        }


_BOOT = _bootstrap_llm()

LLM_PROVIDER: str = _BOOT["provider"]
LLM_MODEL: str = _BOOT["model"]
LLM_BASE_URL: str = _BOOT["base_url"]
LLM_API_KEY: str = _BOOT["api_key"]

# Provider-specific aliases kept for diagnostics and backwards compatibility.
ANTHROPIC_API_KEY: str = os.getenv("ANTHROPIC_API_KEY", "")
ANTHROPIC_MODEL: str = _default_model("anthropic")
OPENAI_API_KEY: str = os.getenv("OPENAI_API_KEY", "")
OPENAI_MODEL: str = _default_model("openai")

# Agent Behavior
# "agent": general analysis chat. "design": Design Copilot session with
# recommendation context. Tools may restrict themselves via Tool.allowed_modes
# (e.g. save_recommendation / visualize_residue / show_comparison_card are
# design-only). Gating is enforced in get_enabled_tools(mode=...) and as a
# defense-in-depth check inside permissions.check_permission.
AGENT_MODE: Literal["agent", "design"] = "agent"
AGENT_MAX_TURNS: int = int(os.getenv("IMMUNOSCOPE_AGENT_MAX_TURNS", "50"))
AGENT_TURN_TIMEOUT: float = float(os.getenv("IMMUNOSCOPE_AGENT_TURN_TIMEOUT", "180.0"))
AGENT_SESSION_WALL_CLOCK: float = float(
    os.getenv("IMMUNOSCOPE_AGENT_SESSION_WALL_CLOCK", "3600.0")
)

# Tool Settings
AGENT_TOOL_TIMEOUT_DEFAULT: float = 300.0  # 5 minutes
AGENT_TOOL_RESULT_MAX_CHARS: int = 50_000
AGENT_PERMISSION_TIMEOUT: float = 60.0

# Budget Management
AGENT_TOKEN_BUDGET: int = int(os.getenv("IMMUNOSCOPE_AGENT_TOKEN_BUDGET", "100000"))
AGENT_MAX_COST_USD: float = float(os.getenv("IMMUNOSCOPE_AGENT_MAX_COST_USD", "5.0"))
AGENT_DIMINISHING_THRESHOLD: int = 500
AGENT_DIMINISHING_TURNS: int = 3

# Audit
AGENT_AUDIT_ENABLED: bool = False
AGENT_AUDIT_DB_PATH: str = os.getenv(
    "IMMUNOSCOPE_AGENT_AUDIT_DB", ".immunoscope/agent_audit.db"
)

# Data Directory
DATA_DIR: str = os.getenv("IMMUNOSCOPE_DATA_DIR", ".immunoscope/data")

# Permissions
AGENT_PERMISSION_DEFAULT: Literal["allow", "ask", "deny"] = "ask"


class Settings:
    """Settings object compatible with PRISM agent framework.

    LLM fields are read freshly from `config_store.effective_llm_config()`
    each time a Settings instance is created. Other AGENT_* fields remain
    env-driven (P0 scope).
    """

    def __init__(self):
        eff = effective_llm_config()
        self.LLM_PROVIDER = eff["provider"]
        self.LLM_MODEL = eff["model"]
        self.LLM_BASE_URL = eff["base_url"]
        self.LLM_API_KEY = eff["api_key"]
        self.LLM_MAX_TOKENS = eff["max_tokens"]
        self.LLM_TEMPERATURE = eff["temperature"]
        self.LLM_KEY_SOURCE = eff["key_source"]

        self.AGENT_MODE = AGENT_MODE
        self.AGENT_MAX_TURNS = AGENT_MAX_TURNS
        self.AGENT_TURN_TIMEOUT = AGENT_TURN_TIMEOUT
        self.AGENT_SESSION_WALL_CLOCK = AGENT_SESSION_WALL_CLOCK

        self.AGENT_TOOL_TIMEOUT_DEFAULT = AGENT_TOOL_TIMEOUT_DEFAULT
        self.AGENT_TOOL_RESULT_MAX_CHARS = AGENT_TOOL_RESULT_MAX_CHARS
        self.AGENT_PERMISSION_TIMEOUT = AGENT_PERMISSION_TIMEOUT

        self.AGENT_TOKEN_BUDGET = AGENT_TOKEN_BUDGET
        self.AGENT_DIMINISHING_THRESHOLD = AGENT_DIMINISHING_THRESHOLD
        self.AGENT_DIMINISHING_TURNS = AGENT_DIMINISHING_TURNS

        self.AGENT_AUDIT_ENABLED = AGENT_AUDIT_ENABLED
        self.AGENT_AUDIT_DB_PATH = AGENT_AUDIT_DB_PATH

        self.DATA_DIR = DATA_DIR


def get_settings() -> Settings:
    """Get agent settings instance (re-reads LLM config each call)."""
    return Settings()
