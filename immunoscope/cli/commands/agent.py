"""Agent command for interactive analysis sessions."""

from __future__ import annotations
import argparse
import asyncio
import os
import sys
from pathlib import Path

import anyio

from immunoscope.agent.config import get_settings
from immunoscope.agent.engine import run_turn
from immunoscope.agent.session import Session
from immunoscope.agent.tool import ToolContext
from immunoscope.agent.exceptions import (
    AgentError,
    MaxTurnsExceeded,
    WallClockExceeded,
)


def main(args):
    """Main entry point for agent command."""
    parser = argparse.ArgumentParser(
        description="Start an interactive agent session for MD analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    ims agent
    ims agent --llm-provider anthropic
    ims agent --llm-model claude-sonnet-4-6
    ims agent --debug

Environment variables:
    IMMUNOSCOPE_LLM_API_KEY - API key for DeepSeek
    ANTHROPIC_API_KEY - API key for Anthropic Claude
    OPENAI_API_KEY - API key for OpenAI
        """,
    )

    parser.add_argument(
        "--llm-provider",
        choices=["deepseek", "anthropic", "openai"],
        default=None,
        help="LLM provider to use (default: from config)",
    )
    parser.add_argument(
        "--llm-model",
        type=str,
        default=None,
        help="LLM model to use (default: from config)",
    )
    parser.add_argument(
        "--max-turns",
        type=int,
        default=None,
        help="Maximum number of turns (default: 50)",
    )
    parser.add_argument(
        "--audit-db",
        type=str,
        default=None,
        help="Path to audit database (default: .immunoscope/agent_audit.db)",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging",
    )

    parsed_args = parser.parse_args(args)

    # Setup logging
    import logging

    if parsed_args.debug:
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        )
    else:
        logging.basicConfig(
            level=logging.WARNING,
            format="%(levelname)s: %(message)s",
        )

    # Run async main
    try:
        asyncio.run(
            async_main(
                parsed_args.llm_provider,
                parsed_args.llm_model,
                parsed_args.max_turns,
                parsed_args.audit_db,
                parsed_args.debug,
            )
        )
        return 0
    except KeyboardInterrupt:
        print("\n\nSession interrupted by user.")
        return 0
    except Exception as e:
        print(f"\nFatal error: {e}", file=sys.stderr)
        if parsed_args.debug:
            import traceback
            traceback.print_exc()
        return 1


async def async_main(
    llm_provider: str | None,
    llm_model: str | None,
    max_turns: int | None,
    audit_db: str | None,
    debug: bool,
):
    """Async main function for agent session."""
    settings = get_settings()

    # Override settings if provided. Provider switching must update the coupled
    # values together, otherwise a CLI override can leave a DeepSeek URL/key
    # paired with an OpenAI or Anthropic provider.
    if llm_provider:
        settings.LLM_PROVIDER = llm_provider
        if not llm_model:
            settings.LLM_MODEL = _default_model_for_provider(llm_provider)
        settings.LLM_BASE_URL = _base_url_for_provider(llm_provider)
        settings.LLM_API_KEY = _api_key_for_provider(llm_provider)
    if llm_model:
        settings.LLM_MODEL = llm_model
    if max_turns:
        settings.AGENT_MAX_TURNS = max_turns

    # Setup audit database
    if audit_db:
        db_path = audit_db
    else:
        db_path = ".immunoscope/agent_audit.db"

    Path(db_path).parent.mkdir(parents=True, exist_ok=True)

    # Create session
    import time
    session = Session(id=f"cli_{int(time.time())}")
    abort_event = anyio.Event()

    # Create tool context with event handler
    async def on_event(event: dict):
        """Handle events from agent engine."""
        event_type = event.get("type")

        if event_type == "text_delta":
            # Stream text output
            print(event.get("text", ""), end="", flush=True)
        elif event_type == "tool_use_start":
            # Tool execution started
            tool_name = event.get("name", "unknown")
            print(f"\n[Calling tool: {tool_name}]", file=sys.stderr)
        elif event_type == "tool_use_end":
            # Tool execution completed
            tool_name = event.get("name", "unknown")
            ok = event.get("ok", False)
            status = "✓" if ok else "✗"
            print(f"[{status} {tool_name}]", file=sys.stderr)
        elif event_type == "progress":
            # Progress message
            print(f"  [{event.get('message')}]", file=sys.stderr)

    ctx = ToolContext(
        session_id=session.id,
        abort_event=abort_event,
        on_event=on_event,
        settings=settings,
        db_path=db_path,
    )

    # Print welcome message
    print("=" * 70)
    print("ImmunoScope Agent - Interactive MD Analysis Assistant")
    print("=" * 70)
    print(f"LLM Provider: {settings.LLM_PROVIDER}")
    print(f"Model: {settings.LLM_MODEL}")
    print(f"Max Turns: {settings.AGENT_MAX_TURNS}")
    print(f"Audit DB: {db_path}")
    print()
    print("Type your questions or requests. The agent will help you analyze")
    print("MD trajectories through natural language conversation.")
    print()
    print("Commands:")
    print("  /quit or /exit - Exit the session")
    print("  /help - Show this help message")
    print("=" * 70)
    print()

    turn_count = 0

    try:
        while turn_count < settings.AGENT_MAX_TURNS:
            # Get user input
            try:
                user_input = input("You> ")
            except (EOFError, KeyboardInterrupt):
                break

            # Handle commands
            if user_input.strip() in ["/quit", "/exit"]:
                break
            elif user_input.strip() == "/help":
                print("\nAvailable commands:")
                print("  /quit, /exit - Exit the session")
                print("  /help - Show this help")
                print()
                continue
            elif not user_input.strip():
                continue

            # Run agent turn
            print()
            print("Agent> ", end="", flush=True)

            try:
                # Run turn (events are handled via on_event callback)
                await run_turn(session, user_input, ctx)

                print()  # Newline after response
                print()

                turn_count += 1

            except MaxTurnsExceeded:
                print("\n\nMaximum turns exceeded. Session ended.", file=sys.stderr)
                break
            except WallClockExceeded:
                print("\n\nSession time limit exceeded. Session ended.", file=sys.stderr)
                break
            except AgentError as e:
                print(f"\n\nAgent error: {e}", file=sys.stderr)
                print("You can continue or type /quit to exit.", file=sys.stderr)
                print()

    except Exception as e:
        print(f"\n\nUnexpected error: {e}", file=sys.stderr)
        if debug:
            import traceback
            traceback.print_exc()
        raise

    # Goodbye message
    print()
    print("=" * 70)
    print(f"Session ended. Total turns: {turn_count}")
    print(f"Audit log: {db_path}")
    print("=" * 70)


def _default_model_for_provider(provider: str) -> str:
    if provider == "anthropic":
        return "claude-sonnet-4-6"
    if provider == "openai":
        return "gpt-4o"
    return "deepseek-chat"


def _base_url_for_provider(provider: str) -> str:
    explicit = os.getenv("IMMUNOSCOPE_LLM_BASE_URL")
    if explicit is not None:
        return explicit
    if provider == "deepseek":
        return "https://api.deepseek.com/v1"
    return ""


def _api_key_for_provider(provider: str) -> str:
    explicit = os.getenv("IMMUNOSCOPE_LLM_API_KEY", "")
    if explicit:
        return explicit
    if provider == "anthropic":
        return os.getenv("ANTHROPIC_API_KEY", "")
    if provider == "openai":
        return os.getenv("OPENAI_API_KEY", "")
    return os.getenv("DEEPSEEK_API_KEY", "")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
