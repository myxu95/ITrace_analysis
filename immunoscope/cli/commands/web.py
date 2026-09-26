#!/usr/bin/env python3
"""IMS Web - lightweight local web interface."""

from __future__ import annotations

from immunoscope.web.server import main as run_web_server


def main(argv: list[str] | None = None) -> int:
    return run_web_server(argv)


if __name__ == "__main__":
    raise SystemExit(main())
