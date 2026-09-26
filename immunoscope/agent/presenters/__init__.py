"""
Analysis data presentation layer for LLM consumption.

Converts raw analysis outputs (CSV/JSON) into standardized Markdown views
optimized for Agent reasoning. Each view is ≤5000 characters with consistent
formatting, field naming, and provenance tracking.
"""

from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path

from .base import (
    Presenter,
    RenderContext,
    PresenterError,
    SPATIAL_LAYERS,
    SUB_FLAVORS,
)
from .locator import CaseLocator

# View registry populated by view modules
VIEW_REGISTRY: Dict[str, type[Presenter]] = {}


def register_view(name: str):
    """Decorator to register a view presenter.

    D-B1 (2026-05-26) guard: rejects a presenter that does not declare a
    valid `spatial_layer` ∈ SPATIAL_LAYERS and a non-empty
    `sub_flavors_served` ⊆ SUB_FLAVORS. This keeps the registry honest so
    that callers can rely on `views_for_layer(...)` /
    `views_for_sub_flavor(...)` returning a complete answer.
    """

    def decorator(cls: type[Presenter]) -> type[Presenter]:
        layer = getattr(cls, "spatial_layer", "")
        flavors = tuple(getattr(cls, "sub_flavors_served", ()) or ())
        if layer not in SPATIAL_LAYERS:
            raise PresenterError(
                f"Presenter '{name}' has invalid spatial_layer={layer!r}; "
                f"must be one of {SPATIAL_LAYERS} (D-B1)."
            )
        if not flavors or any(f not in SUB_FLAVORS for f in flavors):
            raise PresenterError(
                f"Presenter '{name}' has invalid sub_flavors_served="
                f"{flavors!r}; must be a non-empty subset of {SUB_FLAVORS} "
                f"(D-B1)."
            )
        VIEW_REGISTRY[name] = cls
        return cls

    return decorator


# Import views to trigger registration (after register_view is defined)
from . import views as _views


def views_for_layer(layer: str) -> List[str]:
    """D-B1: return registered view names served by a given spatial layer.

    Args:
        layer: one of SPATIAL_LAYERS (case-insensitive).

    Returns:
        Sorted list of view names whose `spatial_layer` matches.
    """
    key = (layer or "").lower()
    return sorted(
        name for name, cls in VIEW_REGISTRY.items()
        if cls.spatial_layer == key
    )


def views_for_sub_flavor(
    sub_flavor: str,
    layer: Optional[str] = None,
) -> List[str]:
    """D-B1: return views serving a sub-flavor, optionally within one layer."""
    key = (sub_flavor or "").lower()
    layer_key = (layer or "").lower() if layer else None
    return sorted(
        name for name, cls in VIEW_REGISTRY.items()
        if key in cls.sub_flavors_served
        and (layer_key is None or cls.spatial_layer == layer_key)
    )


def hierarchy_summary() -> Dict[str, Dict[str, List[str]]]:
    """D-B1: return the canonical view → (layer, sub_flavors) map.

    Returns: {view_name: {"spatial_layer": str, "sub_flavors": [str, ...]}}.
    Used by `query_analysis_results` to emit help text and by tests to
    audit that every registered view declares its D-B1 attributes.
    """
    return {
        name: {
            "spatial_layer": cls.spatial_layer,
            "sub_flavors": list(cls.sub_flavors_served),
        }
        for name, cls in VIEW_REGISTRY.items()
    }


def render(
    view: str,
    locator: CaseLocator,
    filters: Optional[Dict[str, Any]] = None
) -> str:
    """
    Render a view as Markdown.

    Args:
        view: View name (overview, hotspots, interface, etc.)
        locator: CaseLocator for resolving module paths
        filters: Optional filters (top_n, region, residue, partner, etc.)

    Returns:
        Markdown string (≤5000 chars)

    Raises:
        PresenterError: If view not found or rendering fails
    """
    if view not in VIEW_REGISTRY:
        available = ", ".join(sorted(VIEW_REGISTRY.keys()))
        raise PresenterError(
            f"Unknown view '{view}'. Available: {available}"
        )

    presenter_cls = VIEW_REGISTRY[view]
    presenter = presenter_cls()

    context = RenderContext(
        locator=locator,
        filters=filters or {}
    )

    return presenter.render(context)


__all__ = [
    "Presenter",
    "RenderContext",
    "PresenterError",
    "CaseLocator",
    "VIEW_REGISTRY",
    "SPATIAL_LAYERS",
    "SUB_FLAVORS",
    "register_view",
    "render",
    "views_for_layer",
    "views_for_sub_flavor",
    "hierarchy_summary",
]
