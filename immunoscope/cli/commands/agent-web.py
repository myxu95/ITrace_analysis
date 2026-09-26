"""Compatibility shim for the hyphenated ``ims agent-web`` command.

Python module names cannot be imported with a hyphen through normal import
syntax, but ``importlib.import_module`` can load this file. The implementation
lives in ``agent_web.py`` so it can also be imported normally in tests.
"""

from __future__ import annotations

from immunoscope.cli.commands.agent_web import main


if __name__ == "__main__":
    raise SystemExit(main())
