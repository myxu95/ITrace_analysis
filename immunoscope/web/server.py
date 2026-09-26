"""Development server entry point for ImmunoScope Web."""

from __future__ import annotations

import argparse
import os
import socket
import sys

DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8890
DEFAULT_DEBUG = True


def _port_is_free(port: int, host: str = "127.0.0.1") -> bool:
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    try:
        sock.bind((host, port))
        return True
    except OSError:
        return False
    finally:
        sock.close()


def _pick_port(preferred: int, host: str) -> int:
    if _port_is_free(preferred, host):
        return preferred
    for candidate in range(preferred + 1, preferred + 50):
        if _port_is_free(candidate, host):
            return candidate
    return preferred


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run ImmunoScope Web")
    parser.add_argument("--host", default=DEFAULT_HOST, help="Host to bind to (use 0.0.0.0 for remote access)")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT)
    parser.add_argument("--no-reload", action="store_true", help="Disable uvicorn reload")
    args = parser.parse_args(argv)

    try:
        import uvicorn
        from immunoscope.web.config import settings
    except ModuleNotFoundError as exc:
        missing = exc.name or "web runtime dependency"
        print(
            f"Missing dependency for ImmunoScope Web: {missing}\n"
            "Install web dependencies with: pip install -r requirements.txt\n"
            "Or update the conda environment with: conda env update -f environment.yml --prune",
            file=sys.stderr,
        )
        return 1

    host = args.host or settings.HOST
    debug = getattr(settings, "DEBUG", DEFAULT_DEBUG)
    port = _pick_port(args.port, args.host)

    # Display access URLs
    if host == "0.0.0.0":
        # Get local IP for remote access instructions
        import socket
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            s.connect(("8.8.8.8", 80))
            local_ip = s.getsockname()[0]
            s.close()
        except Exception:
            local_ip = "YOUR_SERVER_IP"

        print(f"ImmunoScope Web Server Started")
        print(f"  Local:  http://127.0.0.1:{port}/")
        print(f"  Remote: http://{local_ip}:{port}/")
        print(f"\nAgent Interface:")
        print(f"  Local:  http://127.0.0.1:{port}/#/agent")
        print(f"  Remote: http://{local_ip}:{port}/#/agent")
    else:
        print(f"ImmunoScope Web: http://{host}:{port}/")
        print(f"Agent Interface: http://{host}:{port}/#/agent")

    # Check for agent API key configuration
    has_api_key = any([
        os.getenv("IMMUNOSCOPE_LLM_API_KEY"),
        os.getenv("DEEPSEEK_API_KEY"),
        os.getenv("ANTHROPIC_API_KEY"),
        os.getenv("OPENAI_API_KEY"),
    ])

    if not has_api_key:
        print("\n⚠️  Warning: No LLM API key detected!")
        print("   Agent features require an API key. Set one of:")
        print("   - IMMUNOSCOPE_LLM_API_KEY (recommended)")
        print("   - DEEPSEEK_API_KEY")
        print("   - ANTHROPIC_API_KEY")
        print("   - OPENAI_API_KEY")
        print("   See docs/agent/CONFIGURATION_GUIDE.md for details.\n")

    uvicorn.run(
        "immunoscope.web.app:app",
        host=host,
        port=port,
        reload=not args.no_reload and debug,
        log_level="info",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
