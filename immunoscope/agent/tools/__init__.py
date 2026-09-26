"""Tool registry for ImmunoScope Agent."""

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from immunoscope.agent.tool import Tool

# Global tool registry
TOOL_REGISTRY: dict[str, Tool] = {}


def register_tool(tool_class: type[Tool]) -> type[Tool]:
    """Decorator to register a tool in the global registry.

    Usage:
        @register_tool
        class MyTool(Tool):
            name = "my_tool"
            ...
    """
    if not hasattr(tool_class, "name"):
        raise ValueError(f"Tool {tool_class.__name__} must define 'name' class variable")

    tool_name = tool_class.name
    if tool_name in TOOL_REGISTRY:
        raise ValueError(f"Tool '{tool_name}' is already registered")

    # Instantiate and register
    tool_instance = tool_class()
    TOOL_REGISTRY[tool_name] = tool_instance

    return tool_class


def get_tool(name: str) -> Tool | None:
    """Get a tool by name."""
    return TOOL_REGISTRY.get(name)


def get_enabled_tools(mode: str | None = None) -> list[Tool]:
    """Return tools registered and visible in the given agent mode.

    When `mode` is None, returns the full registry (back-compat for callers
    that don't yet route mode). When `mode` is set, filters out tools whose
    `allowed_modes` does not include it — this is how Design Copilot tools
    stay invisible in generic agent mode."""
    tools = list(TOOL_REGISTRY.values())
    if mode is None:
        return tools
    return [t for t in tools if mode in t.allowed_modes]


def list_tool_names() -> list[str]:
    """List all registered tool names."""
    return list(TOOL_REGISTRY.keys())


# Import all tool modules to trigger registration
# This ensures all @register_tool decorated classes are registered

# Basic tools
from immunoscope.agent.tools import list_files  # noqa: F401, E402
from immunoscope.agent.tools import read_file  # noqa: F401, E402

# Phase 2: Core workflow tools
from immunoscope.agent.tools import preprocess_trajectory  # noqa: F401, E402
from immunoscope.agent.tools import check_quality  # noqa: F401, E402
from immunoscope.agent.tools import calculate_rmsd  # noqa: F401, E402
from immunoscope.agent.tools import calculate_rmsf  # noqa: F401, E402

# Phase 3: Interaction analysis tools
from immunoscope.agent.tools import calculate_rrcs  # noqa: F401, E402
from immunoscope.agent.tools import calculate_bsa  # noqa: F401, E402
from immunoscope.agent.tools import analyze_hbonds  # noqa: F401, E402
from immunoscope.agent.tools import analyze_salt_bridges  # noqa: F401, E402
from immunoscope.agent.tools import analyze_hydrophobic  # noqa: F401, E402
from immunoscope.agent.tools import analyze_pi_interactions  # noqa: F401, E402
from immunoscope.agent.tools import analyze_contacts  # noqa: F401, E402

# Phase 4: Conformation analysis tools
from immunoscope.agent.tools import analyze_landscape  # noqa: F401, E402
from immunoscope.agent.tools import cluster_conformations  # noqa: F401, E402
from immunoscope.agent.tools import analyze_angles  # noqa: F401, E402

# Phase 5: Reporting and comparison tools
from immunoscope.agent.tools import generate_report  # noqa: F401, E402
from immunoscope.agent.tools import compare_systems  # noqa: F401, E402
from immunoscope.agent.tools import query_analysis_results  # noqa: F401, E402

# Upload and web interface tools
from immunoscope.agent.tools import analyze_uploaded  # noqa: F401, E402

# Design Copilot tools (only meaningful in design mode)
from immunoscope.agent.tools import save_recommendation  # noqa: F401, E402
from immunoscope.agent.tools import visualize_residue  # noqa: F401, E402
from immunoscope.agent.tools import show_comparison_card  # noqa: F401, E402
