"""Formatting helpers for Markdown rendering."""

from typing import List, Dict, Any, Optional
import pandas as pd


def pair_display(
    tcr_residue: str,
    partner_residue: str,
    partner_component: Optional[str] = None
) -> str:
    """
    Format a residue pair for display.

    Args:
        tcr_residue: TCR residue label (e.g., "ASP92")
        partner_residue: Partner residue label (e.g., "LYS66")
        partner_component: Optional component (e.g., "HLA_alpha", "peptide")

    Returns:
        Formatted string (e.g., "ASP92 ↔ LYS66 (HLA_alpha)")
    """
    display = f"{tcr_residue} ↔ {partner_residue}"
    if partner_component:
        display += f" ({partner_component})"
    return display


def pct(value: float, decimals: int = 0) -> str:
    """
    Format a fraction as percentage.

    Args:
        value: Fraction (0.0-1.0)
        decimals: Decimal places

    Returns:
        Formatted percentage (e.g., "90%")
    """
    return f"{value * 100:.{decimals}f}%"


def num(value: float, decimals: int = 2, unit: str = "") -> str:
    """
    Format a number with optional unit.

    Args:
        value: Number to format
        decimals: Decimal places
        unit: Optional unit (e.g., "Å", "Å²")

    Returns:
        Formatted string (e.g., "4.39", "627.7 Å²")
    """
    formatted = f"{value:.{decimals}f}"
    if unit:
        formatted += f" {unit}"
    return formatted


def mk_table(
    df: pd.DataFrame,
    columns: List[str],
    headers: Optional[List[str]] = None,
    max_rows: int = 10
) -> str:
    """
    Format DataFrame as Markdown table.

    Args:
        df: DataFrame to format
        columns: Column names to include
        headers: Optional custom headers (defaults to column names)
        max_rows: Maximum rows to include

    Returns:
        Markdown table string
    """
    if df.empty:
        return "_No data_\n"

    # Use provided headers or column names
    if headers is None:
        headers = columns

    # Build header row
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]

    # Add data rows (up to max_rows)
    for idx, row in df.head(max_rows).iterrows():
        cells = []
        for col in columns:
            value = row.get(col, "")
            # Format based on type
            if pd.isna(value):
                cells.append("")
            elif isinstance(value, float):
                cells.append(f"{value:.2f}")
            else:
                cells.append(str(value))

        lines.append("| " + " | ".join(cells) + " |")

    # Add truncation note if needed
    if len(df) > max_rows:
        lines.append(f"\n_Showing {max_rows} of {len(df)} rows_")

    return "\n".join(lines) + "\n"


def classify_fingerprint(
    mean_value: float,
    occupancy: float,
    max_value: Optional[float] = None
) -> str:
    """
    Classify interaction dynamics fingerprint.

    Args:
        mean_value: Mean interaction strength (e.g., RRCS)
        occupancy: Fraction of frames with nonzero interaction
        max_value: Optional max value for fluctuation assessment

    Returns:
        Classification string (e.g., "stable", "fluctuating", "weak")
    """
    # High occupancy + high mean = stable
    if occupancy >= 0.8 and mean_value >= 3.0:
        return "stable"

    # High mean but lower occupancy = fluctuating
    if mean_value >= 3.0 and occupancy < 0.8:
        return "fluctuating"

    # High occupancy but lower mean = weak but persistent
    if occupancy >= 0.8 and mean_value < 3.0:
        return "weak-persistent"

    # Check for high fluctuation if max provided
    if max_value and max_value > mean_value * 2:
        return "highly-fluctuating"

    # Default
    if mean_value >= 2.0:
        return "moderate"
    else:
        return "weak"


def format_region(region: str) -> str:
    """
    Format region name for display.

    Args:
        region: Region name (e.g., "CDR3_alpha", "Framework_beta")

    Returns:
        Formatted region (e.g., "CDR3α", "FWβ")
    """
    # Replace underscores
    region = region.replace("_", "")

    # Greek letters
    region = region.replace("alpha", "α").replace("beta", "β")

    # Abbreviations
    region = region.replace("Framework", "FW")

    return region


def truncate_list(
    items: List[str],
    max_items: int = 10,
    note: str = "items"
) -> List[str]:
    """
    Truncate a list with a note if needed.

    Args:
        items: List of items
        max_items: Maximum items to keep
        note: Description for truncation note

    Returns:
        Truncated list with note appended if needed
    """
    if len(items) <= max_items:
        return items

    truncated = items[:max_items]
    truncated.append(f"_...and {len(items) - max_items} more {note}_")
    return truncated
