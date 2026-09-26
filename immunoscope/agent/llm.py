from __future__ import annotations
import json
import logging
import os
from abc import ABC, abstractmethod
from typing import AsyncIterator, TYPE_CHECKING

import anthropic
import httpx
from openai import AsyncOpenAI

from immunoscope.agent.exceptions import LLMStreamError
from immunoscope.agent.stream_events import (
    StreamEvent,
    TextDelta,
    ToolUseStart,
    ToolUseInputDelta,
    ToolUseEnd,
    MessageStop,
)

if TYPE_CHECKING:
    from immunoscope.agent.config import Settings

log = logging.getLogger("immunoscope.agent.llm")


class LLMClient(ABC):
    """Provider-agnostic streaming interface used by engine.run_turn.
    Both real Anthropic and DeepSeek's Anthropic-compatible endpoint share the
    same wire format, so AnthropicCompatibleClient handles both. The abstraction
    is here so a future HttpxAnthropicClient (or OpenAI-style adapter) can be
    swapped in without touching engine.py."""

    @abstractmethod
    def stream(
        self,
        *,
        system: str,
        messages: list[dict],
        tools: list[dict] | None,
        max_tokens: int = 4096,
    ) -> AsyncIterator[StreamEvent]: ...


class AnthropicCompatibleClient(LLMClient):
    def __init__(self, api_key: str, base_url: str | None, model: str):
        if not api_key:
            raise ValueError("LLM_API_KEY is empty; agent mode cannot start")

        # Defensive: the anthropic SDK falls back to ANTHROPIC_AUTH_TOKEN /
        # ANTHROPIC_API_KEY / ANTHROPIC_BASE_URL when its kwargs are omitted,
        # and in some dev environments those are set for an unrelated purpose.
        # We pass api_key explicitly, base_url explicitly when set, and pop the
        # auth_token env to avoid the SDK silently appending an Authorization
        # bearer header that overrides our x-api-key.
        os.environ.pop("ANTHROPIC_AUTH_TOKEN", None)

        kwargs: dict = {"api_key": api_key}
        if base_url:
            kwargs["base_url"] = base_url
        self._client = anthropic.AsyncAnthropic(**kwargs)
        self._model = model
        log.info("LLM client ready: model=%s base_url=%s", model, base_url or "<default>")

    async def aclose(self) -> None:
        try:
            await self._client.close()
        except Exception:
            log.warning("LLM client close failed", exc_info=True)

    async def stream(
        self,
        *,
        system: str,
        messages: list[dict],
        tools: list[dict] | None,
        max_tokens: int = 4096,
    ) -> AsyncIterator[StreamEvent]:
        """Translate the Anthropic SDK's raw event stream into our normalized
        StreamEvent union. The SDK's `messages.stream` accumulator handles
        input_json_delta concatenation internally; we read the parsed `input`
        dict from `final.content` after the stream ends — that means we do
        NOT try to parse partial JSON ourselves."""

        kwargs: dict = {
            "model": self._model,
            "max_tokens": max_tokens,
            "system": system,
            "messages": messages,
        }
        if tools:
            kwargs["tools"] = tools

        try:
            async with self._client.messages.stream(**kwargs) as stream:
                async for event in stream:
                    et = getattr(event, "type", None)
                    if et == "content_block_start":
                        block = getattr(event, "content_block", None)
                        if block is not None and getattr(block, "type", None) == "tool_use":
                            yield ToolUseStart(id=block.id, name=block.name)
                    elif et == "content_block_delta":
                        delta = getattr(event, "delta", None)
                        dt = getattr(delta, "type", None)
                        if dt == "text_delta":
                            yield TextDelta(text=delta.text)
                        elif dt == "input_json_delta":
                            # Frontend observability only; engine does not accumulate.
                            block_id = _resolve_tool_use_id(stream, event)
                            yield ToolUseInputDelta(
                                id=block_id or "",
                                partial_json=delta.partial_json,
                            )
                        # other delta types (thinking, signature) ignored

                final = await stream.get_final_message()
            for block in final.content:
                if getattr(block, "type", None) == "tool_use":
                    yield ToolUseEnd(id=block.id, name=block.name, input=block.input)
            usage = final.usage.model_dump() if final.usage else {}
            # Preserve full assistant content (including provider-specific blocks
            # like `thinking` that must be echoed back) by serializing each
            # block via its Pydantic model_dump.
            content = [b.model_dump(mode="python", exclude_none=False) for b in final.content]
            yield MessageStop(stop_reason=final.stop_reason, usage=usage, content=content)
        except (httpx.ReadError, httpx.RemoteProtocolError, httpx.ReadTimeout,
                anthropic.APIConnectionError) as e:
            # Transport-level interruption — almost always us cancelling the
            # stream task on abort/disconnect. The engine wraps this generator
            # in a Task it cancels explicitly, so silently ending the stream
            # is the correct behavior; raising would re-route through error
            # handling and pollute the logs with a stack frame for an event
            # that isn't actually an error.
            log.debug("LLM stream interrupted (transport): %s", type(e).__name__)
            return
        except anthropic.APIStatusError as e:
            raise LLMStreamError(f"API status error: {e.status_code} {e.message}") from e
        except anthropic.APIError as e:
            raise LLMStreamError(f"API error: {e}") from e


def _resolve_tool_use_id(stream, event) -> str | None:
    """Find the id of the tool_use content block referenced by a content_block_delta
    event. The SDK's accumulator exposes the in-flight blocks via stream.current_message_snapshot."""
    try:
        idx = getattr(event, "index", None)
        snapshot = getattr(stream, "current_message_snapshot", None)
        if snapshot and idx is not None:
            blocks = getattr(snapshot, "content", None) or []
            if 0 <= idx < len(blocks):
                blk = blocks[idx]
                if getattr(blk, "type", None) == "tool_use":
                    return blk.id
    except Exception:
        pass
    return None


class OpenAICompatibleClient(LLMClient):
    """Client for OpenAI-compatible APIs (OpenAI, DeepSeek, etc.)."""

    def __init__(self, api_key: str, base_url: str | None, model: str):
        if not api_key:
            raise ValueError("LLM_API_KEY is empty; agent mode cannot start")

        kwargs: dict = {"api_key": api_key}
        if base_url:
            kwargs["base_url"] = base_url
        self._client = AsyncOpenAI(**kwargs)
        self._model = model
        log.info("OpenAI-compatible LLM client ready: model=%s base_url=%s", model, base_url or "<default>")

    async def aclose(self) -> None:
        try:
            await self._client.close()
        except Exception:
            log.warning("LLM client close failed", exc_info=True)

    async def stream(
        self,
        *,
        system: str,
        messages: list[dict],
        tools: list[dict] | None,
        max_tokens: int = 4096,
    ) -> AsyncIterator[StreamEvent]:
        """Convert OpenAI streaming format to our normalized StreamEvent union."""

        # Convert Anthropic-style messages to OpenAI format
        openai_messages = []
        if system:
            openai_messages.append({"role": "system", "content": system})

        for msg in messages:
            role = msg.get("role")
            content = msg.get("content")

            # Handle text content
            if isinstance(content, str):
                openai_messages.append({"role": role, "content": content})
            elif isinstance(content, list):
                # Anthropic-style user messages with tool_result blocks must be
                # converted to OpenAI `tool` role messages directly — do NOT
                # emit an empty user message in between, or OpenAI/DeepSeek
                # will reject the sequence ("tool_calls must be followed by
                # tool messages").
                if role == "user" and all(
                    isinstance(b, dict) and b.get("type") == "tool_result"
                    for b in content
                ):
                    for block in content:
                        openai_messages.append({
                            "role": "tool",
                            "tool_call_id": block.get("tool_use_id"),
                            "content": str(block.get("content", "")),
                        })
                    continue

                # Otherwise convert Anthropic content blocks (assistant
                # text/tool_use, or mixed user content) to OpenAI format
                text_parts = []
                tool_calls = []
                for block in content:
                    if isinstance(block, dict):
                        if block.get("type") == "text":
                            text_parts.append(block.get("text", ""))
                        elif block.get("type") == "tool_use":
                            tool_calls.append({
                                "id": block.get("id"),
                                "type": "function",
                                "function": {
                                    "name": block.get("name"),
                                    "arguments": json.dumps(block.get("input", {}))
                                }
                            })

                msg_dict = {"role": role}
                # OpenAI/DeepSeek require `content` to be present and to be a
                # string. Use empty string when the assistant turn is purely
                # tool calls (DeepSeek rejects null).
                msg_dict["content"] = "".join(text_parts) if text_parts else ""
                if tool_calls:
                    msg_dict["tool_calls"] = tool_calls
                openai_messages.append(msg_dict)

                # Handle any remaining tool_result blocks in a mixed user msg
                if role == "user":
                    for block in content:
                        if isinstance(block, dict) and block.get("type") == "tool_result":
                            openai_messages.append({
                                "role": "tool",
                                "tool_call_id": block.get("tool_use_id"),
                                "content": str(block.get("content", "")),
                            })

        # Convert Anthropic tools to OpenAI format
        openai_tools = None
        if tools:
            openai_tools = []
            for tool in tools:
                openai_tools.append({
                    "type": "function",
                    "function": {
                        "name": tool.get("name"),
                        "description": tool.get("description", ""),
                        "parameters": tool.get("input_schema", {})
                    }
                })

        kwargs: dict = {
            "model": self._model,
            "messages": openai_messages,
            "max_tokens": max_tokens,
            "stream": True,
        }
        if openai_tools:
            kwargs["tools"] = openai_tools

        try:
            stream = await self._client.chat.completions.create(**kwargs)

            # Track tool calls being built
            tool_calls_buffer = {}

            async for chunk in stream:
                if not chunk.choices:
                    continue

                choice = chunk.choices[0]
                delta = choice.delta

                # Handle text content
                if delta.content:
                    yield TextDelta(text=delta.content)

                # Handle tool calls
                if delta.tool_calls:
                    for tool_call in delta.tool_calls:
                        idx = tool_call.index

                        # Start new tool call
                        if tool_call.id:
                            tool_calls_buffer[idx] = {
                                "id": tool_call.id,
                                "name": tool_call.function.name if tool_call.function else "",
                                "arguments": ""
                            }
                            yield ToolUseStart(
                                id=tool_call.id,
                                name=tool_call.function.name if tool_call.function else ""
                            )

                        # Accumulate arguments
                        if tool_call.function and tool_call.function.arguments:
                            if idx in tool_calls_buffer:
                                tool_calls_buffer[idx]["arguments"] += tool_call.function.arguments
                                yield ToolUseInputDelta(
                                    id=tool_calls_buffer[idx]["id"],
                                    partial_json=tool_call.function.arguments
                                )

                # Check for stream end
                if choice.finish_reason:
                    # Emit ToolUseEnd for each completed tool call
                    for tool_data in tool_calls_buffer.values():
                        try:
                            parsed_input = json.loads(tool_data["arguments"])
                        except json.JSONDecodeError:
                            parsed_input = {}

                        yield ToolUseEnd(
                            id=tool_data["id"],
                            name=tool_data["name"],
                            input=parsed_input
                        )

                    # Build content for MessageStop
                    content = []
                    for tool_data in tool_calls_buffer.values():
                        try:
                            parsed_input = json.loads(tool_data["arguments"])
                        except json.JSONDecodeError:
                            parsed_input = {}
                        content.append({
                            "type": "tool_use",
                            "id": tool_data["id"],
                            "name": tool_data["name"],
                            "input": parsed_input
                        })

                    # Map OpenAI finish reasons to Anthropic format
                    stop_reason_map = {
                        "stop": "end_turn",
                        "tool_calls": "tool_use",
                        "length": "max_tokens",
                    }
                    stop_reason = stop_reason_map.get(choice.finish_reason, "end_turn")

                    yield MessageStop(
                        stop_reason=stop_reason,
                        usage={},
                        content=content
                    )
                    break

        except httpx.ReadError as e:
            log.debug("LLM stream interrupted (transport): %s", type(e).__name__)
            return
        except Exception as e:
            raise LLMStreamError(f"OpenAI API error: {e}") from e


def make_llm_client(settings: "Settings") -> LLMClient:
    """Create appropriate LLM client based on provider configuration.

    Gemini uses Google's OpenAI-compatible endpoint
    (https://generativelanguage.googleapis.com/v1beta/openai/), so it goes
    through OpenAICompatibleClient with the appropriate base_url.
    """
    provider = settings.LLM_PROVIDER.lower()

    if provider in ("openai", "deepseek", "gemini"):
        return OpenAICompatibleClient(
            api_key=settings.LLM_API_KEY,
            base_url=settings.LLM_BASE_URL or None,
            model=settings.LLM_MODEL,
        )
    elif provider == "anthropic":
        return AnthropicCompatibleClient(
            api_key=settings.LLM_API_KEY,
            base_url=settings.LLM_BASE_URL or None,
            model=settings.LLM_MODEL,
        )
    else:
        raise ValueError(
            f"Unsupported LLM provider: {provider}. "
            "Use 'anthropic', 'openai', 'deepseek', or 'gemini'."
        )
