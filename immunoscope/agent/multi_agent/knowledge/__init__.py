"""Per-agent knowledge slices for the multi-agent mutation-design pipeline.

Each agent receives a small markdown brief that frames its task in
domain terms (which views to consult, what numeric thresholds matter,
what counts as evidence). The briefs live as `.md` files in this package
so they can be read by humans without booting Python; the helper
`load_knowledge(name)` returns the file contents.

The legacy `recommendation/prompts.py` flow uses a single monolithic
system prompt. Multi-agent v1 deliberately splits that knowledge along
reader / writer lines so each agent's prompt stays short and on-topic.
"""

from __future__ import annotations

from importlib import resources
from typing import Final

_PACKAGE: Final[str] = "immunoscope.agent.multi_agent.knowledge"


def load_knowledge(name: str) -> str:
    """Return the markdown body for the named knowledge slice.

    `name` is the base filename without `.md` (e.g. "bio_agent",
    "interaction_reader"). Raises FileNotFoundError when the slice is
    not packaged with this module.
    """
    filename = f"{name}.md"
    pkg = resources.files(_PACKAGE)
    target = pkg / filename
    if not target.is_file():
        raise FileNotFoundError(
            f"Knowledge slice {name!r} not found in {_PACKAGE} "
            f"(looked for {filename})."
        )
    return target.read_text(encoding="utf-8")


__all__ = ["load_knowledge"]
