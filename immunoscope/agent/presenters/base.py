"""Base classes for presenters."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, Any, Optional, Tuple
from pathlib import Path


class PresenterError(Exception):
    """Raised when presentation fails."""
    pass


# D-B1 (2026-05-26) spatial-hierarchy vocabulary. Each Presenter subclass
# must declare which spatial layer it serves; the LLM-facing
# `query_analysis_results` tool uses this to support `spatial_layer=` and
# `sub_flavor=` filters and to drive the §6.5.2 ablation cuts.
SPATIAL_LAYERS = ("complex", "interface", "region", "residue", "pair")
SUB_FLAVORS = ("static", "dynamic")


@dataclass
class RenderContext:
    """Context passed to presenters during rendering."""

    locator: "CaseLocator"  # Forward reference
    filters: Dict[str, Any]

    def get_filter(self, key: str, default: Any = None) -> Any:
        """Get a filter value with default."""
        return self.filters.get(key, default)


class Presenter(ABC):
    """
    Base class for view presenters.

    Each presenter renders one type of view (e.g., hotspots, interface)
    as Markdown text optimized for LLM consumption.

    D-B1 (2026-05-26) hierarchy attributes:
      `spatial_layer`        — one of SPATIAL_LAYERS; the primary layer
                               this view operates at.
      `sub_flavors_served`   — tuple of SUB_FLAVORS values; whether the
                               view surfaces static, dynamic, or both
                               sub-flavors at that layer.
    Subclasses MUST override both (an empty tuple is rejected by the
    registry guard). See `feature_registry.py` for the canonical
    per-metric mapping.
    """

    # Subclasses should set these
    view_name: str = ""
    max_length: int = 5000

    # D-B1 attributes — subclasses must override.
    spatial_layer: str = ""
    sub_flavors_served: Tuple[str, ...] = ()

    @abstractmethod
    def render(self, context: RenderContext) -> str:
        """
        Render the view as Markdown.

        Args:
            context: Rendering context with locator and filters

        Returns:
            Markdown string (≤ max_length chars)

        Raises:
            PresenterError: If required data is missing or rendering fails
        """
        pass

    def _cap_length(self, text: str, note: str = "") -> str:
        """
        Cap text to max_length, appending truncation note if needed.

        Args:
            text: Text to cap
            note: Optional note to append (e.g., "Use filters to narrow")

        Returns:
            Capped text
        """
        if len(text) <= self.max_length:
            return text

        truncation_msg = "\n\n---\n**Truncated** (output exceeded length limit)"
        if note:
            truncation_msg += f": {note}"

        # Reserve space for truncation message
        available = self.max_length - len(truncation_msg)
        return text[:available] + truncation_msg

    def _format_sources(self, *paths: Path) -> str:
        """
        Format source file paths as Markdown list.

        Args:
            *paths: Paths relative to case_dir

        Returns:
            Markdown sources section
        """
        if not paths:
            return ""

        lines = ["\n## Sources\n"]
        for path in paths:
            if path:
                lines.append(f"- `{path}`")

        return "\n".join(lines)


def cap_length(text: str, max_len: int = 5000, note: str = "") -> str:
    """
    Standalone helper to cap text length.

    Args:
        text: Text to cap
        max_len: Maximum length
        note: Optional truncation note

    Returns:
        Capped text
    """
    if len(text) <= max_len:
        return text

    truncation_msg = "\n\n---\n**Truncated**"
    if note:
        truncation_msg += f": {note}"

    available = max_len - len(truncation_msg)
    return text[:available] + truncation_msg
