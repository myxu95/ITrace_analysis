"""Configuration loader for ImmunoScope Reporter."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False


@dataclass
class LLMConfig:
    """LLM configuration."""

    enabled: bool = False
    provider: str = "openai"
    model: str | None = None
    api_key: str | None = None
    base_url: str | None = None
    temperature: float = 0.3
    max_tokens: int = 800
    fallback_to_rules: bool = True


@dataclass
class ServerConfig:
    """Server configuration."""

    host: str = "127.0.0.1"
    port: int = 8770


@dataclass
class ReporterConfig:
    """Complete reporter configuration."""

    llm: LLMConfig = field(default_factory=LLMConfig)
    server: ServerConfig = field(default_factory=ServerConfig)


class ConfigLoader:
    """Load configuration from multiple sources with priority."""

    CONFIG_PATHS = [
        Path.home() / ".immunoscope" / "config.yaml",
        Path.cwd() / "immunoscope_config.yaml",
        Path.cwd() / ".immunoscope_config.yaml",
    ]

    @classmethod
    def load(cls, config_file: str | Path | None = None) -> ReporterConfig:
        """Load configuration with priority: CLI args > config file > env vars > defaults."""
        config = ReporterConfig()

        # 1. Load from config file (lowest priority)
        if config_file:
            cls._load_from_file(Path(config_file), config)
        else:
            # Try default locations
            for path in cls.CONFIG_PATHS:
                if path.exists():
                    cls._load_from_file(path, config)
                    break

        # 2. Override with environment variables (higher priority)
        cls._load_from_env(config)

        return config

    @classmethod
    def _load_from_file(cls, path: Path, config: ReporterConfig) -> None:
        """Load configuration from YAML file."""
        if not YAML_AVAILABLE:
            print("Warning: PyYAML not installed, skipping config file loading")
            return

        try:
            data = yaml.safe_load(path.read_text(encoding="utf-8"))
            if not data:
                return

            # Load LLM config
            if "reporter" in data and "llm" in data["reporter"]:
                llm_data = data["reporter"]["llm"]
                config.llm.enabled = llm_data.get("enabled", config.llm.enabled)
                config.llm.provider = llm_data.get("provider", config.llm.provider)
                config.llm.model = llm_data.get("model", config.llm.model)
                config.llm.api_key = llm_data.get("api_key", config.llm.api_key)
                config.llm.base_url = llm_data.get("base_url", config.llm.base_url)
                config.llm.temperature = llm_data.get("temperature", config.llm.temperature)
                config.llm.max_tokens = llm_data.get("max_tokens", config.llm.max_tokens)
                config.llm.fallback_to_rules = llm_data.get("fallback_to_rules", config.llm.fallback_to_rules)

            # Load server config
            if "server" in data:
                server_data = data["server"]
                config.server.host = server_data.get("host", config.server.host)
                config.server.port = server_data.get("port", config.server.port)

        except Exception as e:
            print(f"Warning: Failed to load config from {path}: {e}")

    @classmethod
    def _load_from_env(cls, config: ReporterConfig) -> None:
        """Load configuration from environment variables."""
        # LLM config
        if os.getenv("IMMUNOSCOPE_LLM_ENABLED"):
            config.llm.enabled = os.getenv("IMMUNOSCOPE_LLM_ENABLED", "").lower() in ("true", "1", "yes")

        if os.getenv("IMMUNOSCOPE_LLM_PROVIDER"):
            config.llm.provider = os.getenv("IMMUNOSCOPE_LLM_PROVIDER")

        if os.getenv("IMMUNOSCOPE_LLM_MODEL"):
            config.llm.model = os.getenv("IMMUNOSCOPE_LLM_MODEL")

        if os.getenv("IMMUNOSCOPE_LLM_BASE_URL"):
            config.llm.base_url = os.getenv("IMMUNOSCOPE_LLM_BASE_URL")

        if os.getenv("IMMUNOSCOPE_LLM_TEMPERATURE"):
            try:
                config.llm.temperature = float(os.getenv("IMMUNOSCOPE_LLM_TEMPERATURE"))
            except ValueError:
                pass

        if os.getenv("IMMUNOSCOPE_LLM_MAX_TOKENS"):
            try:
                config.llm.max_tokens = int(os.getenv("IMMUNOSCOPE_LLM_MAX_TOKENS"))
            except ValueError:
                pass

        # API keys. The unified key has highest priority; otherwise use the
        # provider-specific key matching the active backend.
        if os.getenv("IMMUNOSCOPE_LLM_API_KEY"):
            config.llm.api_key = os.getenv("IMMUNOSCOPE_LLM_API_KEY")
        elif config.llm.provider == "openai" and os.getenv("OPENAI_API_KEY"):
            config.llm.api_key = os.getenv("OPENAI_API_KEY")
        elif config.llm.provider in {"deepseek", "openai_compatible"} and os.getenv("DEEPSEEK_API_KEY"):
            config.llm.api_key = os.getenv("DEEPSEEK_API_KEY")
        elif config.llm.provider == "anthropic" and os.getenv("ANTHROPIC_API_KEY"):
            config.llm.api_key = os.getenv("ANTHROPIC_API_KEY")

        # Server config
        if os.getenv("IMMUNOSCOPE_SERVER_HOST"):
            config.server.host = os.getenv("IMMUNOSCOPE_SERVER_HOST")

        if os.getenv("IMMUNOSCOPE_SERVER_PORT"):
            try:
                config.server.port = int(os.getenv("IMMUNOSCOPE_SERVER_PORT"))
            except ValueError:
                pass


def load_config(config_file: str | Path | None = None) -> ReporterConfig:
    """Convenience function to load configuration."""
    return ConfigLoader.load(config_file)
