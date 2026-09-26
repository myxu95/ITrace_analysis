#!/usr/bin/env python3
"""IMS Agent Web - browser interface for the ImmunoScope Agent."""

from __future__ import annotations

import argparse
import os
import socket
import sys
import threading
import webbrowser

from immunoscope.web.server import main as run_web_server


DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8891


def _has_llm_key() -> bool:
    return any(
        os.getenv(name)
        for name in (
            "IMMUNOSCOPE_LLM_API_KEY",
            "DEEPSEEK_API_KEY",
            "ANTHROPIC_API_KEY",
            "OPENAI_API_KEY",
            "GEMINI_API_KEY",
        )
    )


def _local_ip() -> str:
    try:
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        sock.connect(("8.8.8.8", 80))
        value = sock.getsockname()[0]
        sock.close()
        return value
    except Exception:
        return "YOUR_SERVER_IP"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the ImmunoScope Agent web interface",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ims agent-web
  ims agent-web --no-open
  ims agent-web --port 8891
  ims agent-web --remote --port 8891

SSH forwarding from your laptop:
  ssh -L 8891:127.0.0.1:8891 user@server
  open http://127.0.0.1:8891/#/agent
""",
    )
    parser.add_argument(
        "--host",
        default=DEFAULT_HOST,
        help="Host to bind to. Default is 127.0.0.1 for safe SSH forwarding.",
    )
    parser.add_argument("--port", type=int, default=DEFAULT_PORT, help="Port to bind to.")
    parser.add_argument(
        "--remote",
        action="store_true",
        help="Bind to 0.0.0.0 for direct LAN access. Prefer SSH forwarding on untrusted networks.",
    )
    parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable uvicorn auto-reload. Default is disabled for a stable user-facing server.",
    )
    parser.add_argument(
        "--open",
        action="store_true",
        help="Force opening the Agent page in the default browser after the server starts.",
    )
    parser.add_argument(
        "--no-open",
        action="store_true",
        help="Do not open a browser automatically.",
    )
    parser.add_argument(
        "--require-api-key",
        action="store_true",
        help="Fail fast if no LLM API key is configured.",
    )
    parser.add_argument(
        "--tunnel-user",
        default="user",
        help="Username shown in the SSH tunnel hint.",
    )
    parser.add_argument(
        "--tunnel-host",
        default="server",
        help="Server hostname shown in the SSH tunnel hint.",
    )
    return parser


def _should_auto_open(args: argparse.Namespace, host: str) -> bool:
    if args.no_open or host == "0.0.0.0":
        return False
    if args.open:
        return True
    if os.getenv("SSH_CONNECTION") or os.getenv("SSH_TTY"):
        return False
    return bool(os.getenv("DISPLAY") or os.getenv("WAYLAND_DISPLAY"))


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    host = "0.0.0.0" if args.remote else args.host

    if args.require_api_key and not _has_llm_key():
        print(
            "No LLM API key detected. Set one of: "
            "IMMUNOSCOPE_LLM_API_KEY, DEEPSEEK_API_KEY, ANTHROPIC_API_KEY, "
            "OPENAI_API_KEY, GEMINI_API_KEY.",
            file=sys.stderr,
        )
        return 1

    print("ImmunoScope Agent Web")
    print("=" * 72)
    if host == "0.0.0.0":
        local_ip = _local_ip()
        print("Mode: direct LAN access")
        print(f"URL:  http://{local_ip}:{args.port}/#/agent")
        print("Security: this exposes the server to your network.")
    else:
        local_url = f"http://127.0.0.1:{args.port}/#/agent"
        print("Mode: local bind, SSH-forward friendly")
        print(f"URL:  {local_url}")
        print("SSH tunnel:")
        print(
            f"  ssh -L {args.port}:127.0.0.1:{args.port} "
            f"{args.tunnel_user}@{args.tunnel_host}"
        )
        if _should_auto_open(args, host):
            print("Browser: opening automatically")
            threading.Timer(1.0, lambda: webbrowser.open(local_url)).start()
        else:
            print("Browser: open the URL above after SSH forwarding")
    print("=" * 72)

    server_args = ["--host", host, "--port", str(args.port)]
    if not args.reload:
        server_args.append("--no-reload")
    return run_web_server(server_args)


if __name__ == "__main__":
    raise SystemExit(main())
