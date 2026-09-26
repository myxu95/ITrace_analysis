from __future__ import annotations


class AgentError(Exception):
    """Base for all agent-internal errors."""


class ToolError(AgentError):
    """Recoverable tool error. Surfaced to the model as tool_result(is_error=True)
    so the model can retry, switch tools, or explain the failure to the user."""


class ToolValidationError(ToolError):
    """Tool input failed pydantic / business validation."""


class ToolTimeout(ToolError):
    """Tool exceeded its per-call timeout."""


class ToolFatal(AgentError):
    """Infrastructure-level failure that should abort the current turn rather
    than be fed back to the model. Use sparingly: most failures should be
    recoverable ToolError so the model can self-correct."""


class PermissionDenied(ToolError):
    """User declined or system policy denied a tool call."""


class PermissionTimeout(ToolError):
    """User did not respond to a permission_request within the timeout."""


class MaxTurnsExceeded(AgentError):
    """The agent loop hit AGENT_MAX_TURNS without converging."""


class WallClockExceeded(AgentError):
    """The session exceeded AGENT_SESSION_WALL_CLOCK."""


class LLMStreamError(AgentError):
    """The LLM stream raised or terminated unexpectedly."""
