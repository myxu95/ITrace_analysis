"""Fingerprint view - interface occupancy distribution."""

from pathlib import Path
import pandas as pd

from ..base import Presenter, RenderContext
from ..loaders import read_csv, read_json
from ..formatters import pct, num
from .. import register_view


@register_view("fingerprint")
class FingerprintPresenter(Presenter):
    """
    Fingerprint view - interface occupancy distribution.

    Shows distribution of interaction occupancies and families across
    the interface, revealing overall binding mode characteristics.
    """

    view_name = "fingerprint"
    max_length = 3000
    # D-B1: Region-layer aggregation of Residue-level occupancy / RMSF /
    # interaction family per CDR / MHC helix; carries dynamic (distribution)
    # information by design.
    spatial_layer = "region"
    sub_flavors_served = ("dynamic",)

    def render(self, context: RenderContext) -> str:
        """Render fingerprint view."""
        locator = context.locator
        case_id = locator.get_case_id()

        lines = [f"# Interface Fingerprint: {case_id}\n"]

        # Load RRCS data
        rrcs_root = locator.get_module_root("rrcs")
        if not rrcs_root:
            lines.append("_RRCS analysis not found._")
            return "\n".join(lines)

        csv_path = rrcs_root / "analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv"
        df = read_csv(csv_path)

        if df.empty:
            lines.append("_RRCS data not available._")
            return "\n".join(lines)

        # Render sections
        lines.append(self._render_occupancy_distribution(df))
        lines.append(self._render_interaction_families(df))
        lines.append(self._render_region_breakdown(df))
        lines.append(self._render_interpretation(df))

        # Sources
        lines.append(self._format_sources(
            csv_path.relative_to(locator.case_dir)
        ))

        output = "\n".join(lines)
        return self._cap_length(output)

    def _render_occupancy_distribution(self, df: pd.DataFrame) -> str:
        """Render occupancy distribution."""
        lines = ["## Occupancy distribution\n"]

        if "rrcs_nonzero_fraction" not in df.columns:
            lines.append("_Occupancy data not available_")
            return "\n".join(lines)

        # Bin occupancies
        bins = [
            (0.8, 1.0, "High (≥80%)"),
            (0.5, 0.8, "Medium (50-80%)"),
            (0.2, 0.5, "Low (20-50%)"),
            (0.0, 0.2, "Transient (<20%)"),
        ]

        for low, high, label in bins:
            count = len(df[(df["rrcs_nonzero_fraction"] >= low) & (df["rrcs_nonzero_fraction"] < high)])
            if low == 0.8:  # Include 1.0 in high bin
                count = len(df[df["rrcs_nonzero_fraction"] >= low])

            lines.append(f"- **{label}**: {count} pairs")

        return "\n".join(lines) + "\n"

    def _render_interaction_families(self, df: pd.DataFrame) -> str:
        """Render interaction family breakdown."""
        lines = ["## Interaction families\n"]

        if "interaction_family" not in df.columns:
            lines.append("_Interaction family data not available_")
            return "\n".join(lines)

        # Count by family
        family_counts = df["interaction_family"].value_counts()

        if family_counts.empty:
            lines.append("_No interaction families found_")
            return "\n".join(lines)

        total = len(df)
        for family, count in family_counts.items():
            if pd.notna(family):
                lines.append(f"- **{family}**: {count} pairs ({pct(count/total)})")

        return "\n".join(lines) + "\n"

    def _render_region_breakdown(self, df: pd.DataFrame) -> str:
        """Render breakdown by TCR region."""
        lines = ["## By TCR region\n"]

        region_col = "tcr_region_detailed" if "tcr_region_detailed" in df.columns else "tcr_region"

        if region_col not in df.columns:
            lines.append("_Region data not available_")
            return "\n".join(lines)

        # Count by region
        region_counts = df[region_col].value_counts()

        if region_counts.empty:
            lines.append("_No region data found_")
            return "\n".join(lines)

        # Show top regions
        for region, count in region_counts.head(8).items():
            if pd.notna(region):
                # Calculate mean RRCS for this region
                region_df = df[df[region_col] == region]
                mean_rrcs = region_df["mean_rrcs"].mean() if "mean_rrcs" in region_df.columns else 0

                lines.append(f"- **{region}**: {count} pairs, mean RRCS {num(mean_rrcs, 2)}")

        return "\n".join(lines) + "\n"

    def _render_interpretation(self, df: pd.DataFrame) -> str:
        """Render interpretation of binding mode."""
        lines = ["## Binding mode interpretation\n"]

        if "rrcs_nonzero_fraction" not in df.columns or "mean_rrcs" not in df.columns:
            lines.append("_Insufficient data for interpretation_")
            return "\n".join(lines)

        # Calculate metrics
        high_occupancy = len(df[df["rrcs_nonzero_fraction"] >= 0.8])
        total_pairs = len(df)
        mean_rrcs = df["mean_rrcs"].mean()

        # Classify binding mode
        if high_occupancy / total_pairs >= 0.5:
            lines.append("- **Mode**: Stable binding — many persistent contacts")
        elif high_occupancy / total_pairs >= 0.3:
            lines.append("- **Mode**: Mixed — combination of stable and dynamic contacts")
        else:
            lines.append("- **Mode**: Dynamic binding — mostly transient contacts")

        # Strength assessment
        if mean_rrcs >= 3.5:
            lines.append("- **Strength**: Strong interface (mean RRCS ≥ 3.5)")
        elif mean_rrcs >= 2.5:
            lines.append("- **Strength**: Moderate interface (mean RRCS 2.5-3.5)")
        else:
            lines.append("- **Strength**: Weak interface (mean RRCS < 2.5)")

        return "\n".join(lines)
