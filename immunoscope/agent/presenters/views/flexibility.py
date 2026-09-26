"""Flexibility view - RMSF by region."""

from pathlib import Path
import pandas as pd

from ..base import Presenter, RenderContext
from ..loaders import read_json, read_csv
from ..formatters import num, mk_table
from .. import register_view


@register_view("flexibility")
class FlexibilityPresenter(Presenter):
    """
    Flexibility view - RMSF analysis.

    Shows RMSF statistics by region and highlights highly flexible
    residues that may tolerate mutations.
    """

    view_name = "flexibility"
    max_length = 3000
    # D-B1: per-residue Cα / backbone / sidechain RMSF — Residue-level
    # dynamic (per-frame fluctuation magnitude). Per-region RMSF aggregates
    # surfaced here are computed on-the-fly from the Residue rows.
    spatial_layer = "residue"
    sub_flavors_served = ("dynamic",)

    def render(self, context: RenderContext) -> str:
        """Render flexibility view."""
        locator = context.locator
        case_id = locator.get_case_id()

        # Get filters
        min_rmsf = context.get_filter("min_rmsf", 3.0)
        region_filter = context.get_filter("region")
        top_n = context.get_filter("top_n", 15)

        lines = [f"# Flexibility: {case_id}\n"]

        # Load RMSF data
        rmsf_root = locator.get_module_root("rmsf")
        if not rmsf_root:
            lines.append("_RMSF analysis not found. Run RMSF analysis first._")
            return "\n".join(lines)

        # Load summary
        summary_path = rmsf_root / "analysis/rmsf/rmsf_summary.json"
        summary = read_json(summary_path)

        if summary:
            lines.append(self._render_summary_stats(summary))

        # Load region summary
        region_path = rmsf_root / "analysis/rmsf/region_rmsf_summary.csv"
        region_df = read_csv(region_path)

        if not region_df.empty:
            lines.append(self._render_by_region(region_df, region_filter))

        # Load per-residue data for flexible residues
        residue_path = rmsf_root / "analysis/rmsf/residue_rmsf.csv"
        residue_df = read_csv(residue_path)

        if not residue_df.empty:
            lines.append(self._render_flexible_residues(
                residue_df, min_rmsf, region_filter, top_n
            ))

        # Interpretation guide
        lines.append(self._render_interpretation())

        # Sources
        sources = [p.relative_to(locator.case_dir) for p in [summary_path, region_path, residue_path] if p.exists()]
        if sources:
            lines.append(self._format_sources(*sources))

        output = "\n".join(lines)
        return self._cap_length(output, "Use filters: region, min_rmsf, top_n")

    def _render_summary_stats(self, summary: dict) -> str:
        """Render overall RMSF statistics."""
        lines = ["## Overall statistics\n"]

        # Try angstrom-suffixed names first (actual format), then bare names
        tcr_mean = summary.get("tcr_mean_rmsf_angstrom",
                               summary.get("tcr_mean_rmsf"))
        phla_mean = summary.get("phla_mean_rmsf_angstrom",
                                summary.get("phla_mean_rmsf"))
        overall_mean = summary.get("mean_rmsf_angstrom",
                                   summary.get("overall_mean_rmsf",
                                               summary.get("mean_rmsf")))

        if tcr_mean is not None:
            lines.append(f"- **TCR mean RMSF**: {num(tcr_mean, 2, 'Å')}")
        if phla_mean is not None:
            lines.append(f"- **pMHC mean RMSF**: {num(phla_mean, 2, 'Å')}")
        if overall_mean is not None:
            lines.append(f"- **Overall mean**: {num(overall_mean, 2, 'Å')}")

        return "\n".join(lines) + "\n"

    def _render_by_region(self, df: pd.DataFrame, region_filter: str) -> str:
        """Render RMSF by region."""
        lines = ["## By region\n"]

        # Filter if requested
        if region_filter and "region" in df.columns:
            df = df[df["region"].str.contains(region_filter, case=False, na=False)]

        if df.empty:
            lines.append("_No regions match filter_")
            return "\n".join(lines)

        # Sort by mean RMSF
        if "mean_rmsf" in df.columns:
            df = df.sort_values("mean_rmsf", ascending=False)

        # Build table
        table_data = []
        for _, row in df.iterrows():
            table_data.append({
                "Region": row.get("region", ""),
                "Mean RMSF": f"{row.get('mean_rmsf', 0):.2f} Å",
                "Max RMSF": f"{row.get('max_rmsf', 0):.2f} Å",
                "N residues": int(row.get("n_residues", 0)),
            })

        table_df = pd.DataFrame(table_data)
        lines.append(mk_table(
            table_df,
            columns=["Region", "Mean RMSF", "Max RMSF", "N residues"],
            max_rows=15
        ))

        return "\n".join(lines)

    def _render_flexible_residues(
        self,
        df: pd.DataFrame,
        min_rmsf: float,
        region_filter: str,
        top_n: int
    ) -> str:
        """Render highly flexible residues."""
        lines = [f"## Flexible residues (RMSF ≥ {min_rmsf} Å)\n"]

        # Determine RMSF column name (actual data uses rmsf_angstrom)
        rmsf_col = None
        for candidate in ["rmsf_angstrom", "rmsf"]:
            if candidate in df.columns:
                rmsf_col = candidate
                break

        if not rmsf_col:
            lines.append("_RMSF column not found in data_")
            return "\n".join(lines)

        # Determine region column
        region_col = None
        for candidate in ["tcr_region_detailed", "tcr_region", "region_group", "region"]:
            if candidate in df.columns:
                region_col = candidate
                break

        # Filter by RMSF
        df = df[df[rmsf_col] >= min_rmsf]

        # Filter by region if requested
        if region_filter and region_col:
            df = df[df[region_col].astype(str).str.contains(region_filter, case=False, na=False)]

        if df.empty:
            lines.append(f"_No residues with RMSF ≥ {min_rmsf} Å_")
            return "\n".join(lines)

        # Sort and limit
        df = df.sort_values(rmsf_col, ascending=False).head(top_n)

        lines.append(f"Found {len(df)} flexible residues:\n")

        # List residues
        for _, row in df.iterrows():
            # Build residue label from resname+resid (actual data format)
            if "resname" in row and "resid" in row:
                residue = f"{row.get('resname', '')}{row.get('resid', '')}"
            else:
                residue = row.get("residue", "")

            rmsf = row.get(rmsf_col, 0)
            region = row.get(region_col, "") if region_col else ""

            lines.append(f"- **{residue}** ({region}): {num(rmsf, 2, 'Å')}")

        return "\n".join(lines) + "\n"

    def _render_interpretation(self) -> str:
        """Render interpretation guide."""
        lines = [
            "## Interpretation\n",
            "- **RMSF > 3 Å**: High flexibility, induced fit, tolerates mutations",
            "- **RMSF 1-3 Å**: Moderate flexibility, typical for CDR loops",
            "- **RMSF < 1 Å**: Rigid, structurally constrained, risky to mutate",
        ]
        return "\n".join(lines)
