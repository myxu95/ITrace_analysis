from __future__ import annotations
import enum
import json
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Awaitable, Callable, ClassVar, Literal, TYPE_CHECKING

import anyio
from pydantic import BaseModel

if TYPE_CHECKING:
    from immunoscope.agent.config import Settings


class Permission(str, enum.Enum):
    ALLOW = "allow"
    ASK = "ask"
    DENY = "deny"


@dataclass
class PermissionDecision:
    behavior: Permission
    reason: str = ""
    updated_input: dict | None = None


@dataclass
class ToolResult:
    """Result from Tool.call(). The engine handles serialization to a tool_result
    block and optional truncation against AGENT_TOOL_RESULT_MAX_CHARS."""
    content: Any                     # str | dict | list | None — engine will str-ify
    is_error: bool = False
    truncated: bool = False
    saved_to: str | None = None      # path if engine spilled to disk


@dataclass
class ToolContext:
    """Per-session context handed to every tool. Tools should treat this as
    read-only except for `on_event` (emit progress) and `pending_permissions`
    (managed by permissions.py)."""
    session_id: str
    abort_event: anyio.Event
    on_event: Callable[[dict], Awaitable[None]]
    settings: "Settings"
    db_path: str
    pending_permissions: dict[str, Any] = field(default_factory=dict)

    # engine.run_turn writes here after each LLM stream end so execute_tools
    # can attach the turn's token + cost deltas to the FIRST tool's audit
    # row. Resetting to None after the first tool of the turn prevents
    # double-counting when several tools fire in the same turn.
    last_turn_usage: dict | None = None
    last_turn_cost: float | None = None

    @property
    def aborted(self) -> bool:
        return self.abort_event.is_set()


def _normalize_for_anthropic(schema: dict) -> dict:
    """Strip pydantic JSON-schema artifacts the Anthropic API does not accept.

    - Drops `title` keys (cosmetic, but Anthropic's strict mode complains)
    - Forces `additionalProperties: False` on object schemas (matches Anthropic's
      strict-tool expectations)
    - Inlines $defs (Anthropic accepts $defs but DeepSeek's compat layer is less
      forgiving; inlining keeps the schema portable)
    """
    def _walk(node):
        if isinstance(node, dict):
            node.pop("title", None)
            if node.get("type") == "object" and "additionalProperties" not in node:
                node["additionalProperties"] = False
            for v in node.values():
                _walk(v)
        elif isinstance(node, list):
            for v in node:
                _walk(v)

    out = json.loads(json.dumps(schema))  # deep copy
    _walk(out)
    return out


class Tool(ABC):
    """Base for all tools. Subclass and set the ClassVars + implement call()."""

    name: ClassVar[str]
    description: ClassVar[str]
    Input: ClassVar[type[BaseModel]]

    is_read_only: ClassVar[bool] = True
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float | None] = None  # None = use settings default

    # Marks tools that mutate persistent state in a non-recoverable way
    # (delete, overwrite, send). Surfaced in permission_request as a UI nudge
    # and recorded in audit so reviewers can find write paths quickly.
    is_destructive: ClassVar[bool] = False

    # Behavior on user abort while this tool is running:
    #   - "cancel": cancel the running tool task immediately. Right for
    #               long-running read-only ops (HTTP fetch, log tail).
    #   - "block":  let the tool finish, then unwind. Right for filesystem
    #               writes / DB transactions where mid-flight cancellation
    #               leaves dirty state.
    # Default is "cancel" because all current tools are read-only and safe
    # to restart. Action tools (create_build, delete_build) added in later
    # sprints must explicitly opt into "block".
    interrupt_behavior: ClassVar[Literal["cancel", "block"]] = "cancel"

    # Which agent modes may register / call this tool. Default is "both
    # known modes" so existing analysis tools stay available everywhere.
    # Design Copilot tools (save_recommendation, visualize_residue,
    # show_comparison_card) override this to {"design"} so they don't leak
    # into the generic agent mode where the recommendation context is
    # absent. Mode gating is enforced in agent/permissions.check_mode().
    allowed_modes: ClassVar[frozenset[str]] = frozenset({"agent", "design"})

    def to_anthropic_schema(self) -> dict:
        """Return the tool description in Anthropic's tool-schema format."""
        schema = self.Input.model_json_schema()
        return {
            "name": self.name,
            "description": self.description.strip(),
            "input_schema": _normalize_for_anthropic(schema),
        }

    def system_prompt_section(self) -> str:
        """Optional extra instructions appended to the main system prompt."""
        return ""

    async def validate(self, args: BaseModel, ctx: ToolContext) -> None:
        """Business-level validation. Raise ToolValidationError to reject."""
        return None

    @abstractmethod
    async def call(self, args: BaseModel, ctx: ToolContext) -> ToolResult: ...

    def inputs_equivalent(self, a: dict, b: dict) -> bool:
        """Two raw inputs map to the same effect.

        Override when key ordering or auxiliary fields don't change semantics.
        Used by permission caching: when the user has already approved a
        tool call with input X, an in-flight call with input X' that this
        method considers equivalent gets the same decision without
        re-prompting. Default == comparison is strict and correct for
        small, schema-validated payloads."""
        return a == b

    def serialize_result(self, result: ToolResult, tool_use_id: str) -> dict:
        """Convert a ToolResult into a tool_result content block for the API."""
        if isinstance(result.content, str):
            content_str = result.content
        elif result.content is None:
            content_str = ""
        else:
            content_str = json.dumps(result.content, ensure_ascii=False, default=str)
        return {
            "type": "tool_result",
            "tool_use_id": tool_use_id,
            "content": content_str,
            "is_error": result.is_error,
        }
