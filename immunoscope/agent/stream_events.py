from __future__ import annotations
from dataclasses import dataclass
from typing import Union


@dataclass
class TextDelta:
    text: str


@dataclass
class ToolUseStart:
    id: str
    name: str


@dataclass
class ToolUseInputDelta:
    """Partial JSON fragment of a tool_use input. For frontend observability
    only — engine.py does NOT accumulate these; it relies on the SDK's final
    message accumulator to deliver the parsed input dict."""
    id: str
    partial_json: str


@dataclass
class ToolUseEnd:
    """Emitted after the SDK has fully accumulated a tool_use block.
    `input` is the parsed dict, ready to validate."""
    id: str
    name: str
    input: dict


@dataclass
class MessageStop:
    stop_reason: str | None
    usage: dict
    # Raw assistant content blocks from SDK final message (already dict-ified).
    # Contains text, tool_use, thinking, redacted_thinking, etc. blocks in
    # arrival order. Engine MUST forward this verbatim as the next assistant
    # message — re-encoding loses provider-specific blocks (e.g. DeepSeek's
    # `thinking` block, which the server requires to be echoed back on the
    # next request).
    content: list[dict] | None = None


@dataclass
class StreamError:
    error: BaseException


StreamEvent = Union[
    TextDelta, ToolUseStart, ToolUseInputDelta, ToolUseEnd, MessageStop, StreamError
]
