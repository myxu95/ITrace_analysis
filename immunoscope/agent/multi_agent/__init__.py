"""Multi-agent mutation design pipeline (v1, 2026-05-28).

Round-1 architecture: four parallel agents that each read a slice of
the evidence and emit an `AgentMessage`. A later integration round
(deferred) consumes all four briefs together; this package implements
round 1 only.

Round 1 agents (run in parallel):
    1. BioAgent              — literature-grounded biological briefing.
    2. InteractionReader     — hotspots / pair / fingerprint views.
    3. ConformationReader    — flexibility / dihedrals views.
    4. InterfaceExposureReader — interface (BSA) + exposure (per-residue SASA,
                                 cavity candidates) views.

Deferred (iteration mode that re-models and re-runs MD on prior
recommendations): FES / clustering, docking-angle change.

See:
    `messaging.py`     — AgentMessage envelope + EvidenceRef tagged union
    `context.py`       — MutationDesignContext shared state
    `base.py`          — BaseAgent skeleton (LLM + tool loop + emit_message)
    `orchestrator.py`  — parallel fan-out + overview/quality preamble
    `agents/*`         — concrete agent implementations
    `knowledge/*`      — per-agent domain knowledge slices (markdown)
    `tools.py`         — multi-agent-only tools (e.g. SearchLiteratureTool)

Quick start:

    from immunoscope.agent.llm import make_llm_client
    from immunoscope.agent.config import get_settings
    from immunoscope.agent.multi_agent import (
        MutationDesignContext, build_default_pipeline,
    )

    llm = make_llm_client(get_settings())
    pipeline = build_default_pipeline(llm)
    ctx = MutationDesignContext(
        case_dir="output/5c0a_run2_full_analysis",
        case_id="5c0a_run2",
        design_intent="improve TCR affinity for MART-1 / HLA-A*02:01",
    )
    result = await pipeline.run(ctx)
    for msg in result.context.message_history:
        print(msg.sender, msg.msg_id)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .context import MutationDesignContext, SkillCallRecord
from .messaging import (
    AgentMessage,
    ContextFieldRef,
    EvidenceRef,
    LiteratureRef,
    MessageRef,
    MessageValidationError,
    SkillCallRef,
    validate_message,
)

# Agent / pipeline / base classes pull in the `LLMClient` import chain
# (anthropic + openai SDKs). Keep them off the eager import path so the
# pure-data modules above can be used without those runtime deps.
if TYPE_CHECKING:  # pragma: no cover
    from immunoscope.agent.llm import LLMClient


def __getattr__(name: str):  # PEP 562 lazy module attribute resolution
    """Lazy-load LLM-bound names so importing this package does not pull
    in the anthropic / openai SDKs unless the caller actually needs an
    agent or pipeline instance."""
    if name == "BaseAgent":
        from .base import BaseAgent
        return BaseAgent
    if name == "AgentRunResult":
        from .base import AgentRunResult
        return AgentRunResult
    if name == "AgentTimeoutError":
        from .base import AgentTimeoutError
        return AgentTimeoutError
    if name == "BioAgent":
        from .agents import BioAgent
        return BioAgent
    if name == "InteractionReader":
        from .agents import InteractionReader
        return InteractionReader
    if name == "ConformationReader":
        from .agents import ConformationReader
        return ConformationReader
    if name == "InterfaceExposureReader":
        from .agents import InterfaceExposureReader
        return InterfaceExposureReader
    if name == "RecommendationAgent":
        from .agents import RecommendationAgent
        return RecommendationAgent
    if name == "DesignCriticAgent":
        from .agents import DesignCriticAgent
        return DesignCriticAgent
    if name == "MutationDesignPipeline":
        from .orchestrator import MutationDesignPipeline
        return MutationDesignPipeline
    if name == "PipelineResult":
        from .orchestrator import PipelineResult
        return PipelineResult
    if name == "make_pipeline":
        from .orchestrator import make_pipeline
        return make_pipeline
    if name == "build_default_pipeline":
        return _build_default_pipeline
    raise AttributeError(f"module 'immunoscope.agent.multi_agent' has no attribute {name!r}")


def _make_independent_critic_llm(default_llm: "LLMClient"):
    """Build a SEPARATE LLM client for the critic when configured, so the
    reviewer is a different model than the generator (mitigating the
    shared-model blind spot). Opt-in via env so we never gamble on an
    unknown endpoint:

        IMMUNOSCOPE_CRITIC_PROVIDER  openai | anthropic | deepseek | gemini
        IMMUNOSCOPE_CRITIC_MODEL     model id
        IMMUNOSCOPE_CRITIC_BASE_URL  (optional) override endpoint
        IMMUNOSCOPE_CRITIC_API_KEY   (optional) else the provider's env key

    Returns `default_llm` (same model) when not configured or on any failure.
    """
    import logging
    import os

    log = logging.getLogger("immunoscope.agent.multi_agent")
    provider = os.getenv("IMMUNOSCOPE_CRITIC_PROVIDER", "").strip().lower()
    if not provider:
        return default_llm
    model = os.getenv("IMMUNOSCOPE_CRITIC_MODEL", "").strip()
    base_url = os.getenv("IMMUNOSCOPE_CRITIC_BASE_URL", "").strip() or None
    key = os.getenv("IMMUNOSCOPE_CRITIC_API_KEY", "").strip() or os.getenv(
        {"openai": "OPENAI_API_KEY", "anthropic": "ANTHROPIC_API_KEY",
         "deepseek": "DEEPSEEK_API_KEY", "gemini": "GEMINI_API_KEY"}.get(provider, ""),
        "",
    )
    if not model or not key:
        log.warning(
            "IMMUNOSCOPE_CRITIC_PROVIDER=%s set but MODEL/API key missing; "
            "critic will use the same model as the generator.", provider,
        )
        return default_llm
    try:
        from immunoscope.agent.llm import (
            AnthropicCompatibleClient, OpenAICompatibleClient,
        )
        if provider == "anthropic":
            client = AnthropicCompatibleClient(api_key=key, base_url=base_url, model=model)
        else:  # openai / deepseek / gemini all speak the OpenAI wire format
            client = OpenAICompatibleClient(api_key=key, base_url=base_url, model=model)
        log.info("DesignCritic using INDEPENDENT model: provider=%s model=%s", provider, model)
        return client
    except Exception as exc:  # noqa: BLE001
        log.warning("failed to build independent critic LLM (%s); same model: %s",
                    provider, exc)
        return default_llm


def _build_default_pipeline(llm: "LLMClient"):
    """Construct a pipeline wired with the v1 default agent classes."""
    from .agents import (
        BioAgent,
        ConformationReader,
        DesignCriticAgent,
        InteractionReader,
        InterfaceExposureReader,
        RecommendationAgent,
    )
    from .orchestrator import MutationDesignPipeline

    critic_llm = _make_independent_critic_llm(llm)
    return MutationDesignPipeline(
        bio_agent=BioAgent(llm),
        interaction_reader=InteractionReader(llm),
        conformation_reader=ConformationReader(llm),
        interface_exposure_reader=InterfaceExposureReader(llm),
        recommendation_agent=RecommendationAgent(llm),
        design_critic=DesignCriticAgent(critic_llm),
    )


__all__ = [
    "AgentMessage",
    "AgentRunResult",
    "AgentTimeoutError",
    "BaseAgent",
    "BioAgent",
    "ConformationReader",
    "ContextFieldRef",
    "DesignCriticAgent",
    "EvidenceRef",
    "InteractionReader",
    "InterfaceExposureReader",
    "LiteratureRef",
    "MessageRef",
    "MessageValidationError",
    "MutationDesignContext",
    "MutationDesignPipeline",
    "PipelineResult",
    "RecommendationAgent",
    "SkillCallRecord",
    "SkillCallRef",
    "build_default_pipeline",
    "make_pipeline",
    "validate_message",
]
