"""Residue view - single-residue drill-down."""

from pathlib import Path
import pandas as pd

from ..base import Presenter, RenderContext, PresenterError
from ..loaders import read_csv
from ..formatters import pair_display, pct, num
from .. import register_view
from immunoscope.analysis.features import (
    CaseLocator as FeatureCaseLocator,
    FeatureSet,
    compute_feature,
)


@register_view("residue")
class ResiduePresenter(Presenter):
    """
    Residue view - single-residue drill-down.

    Shows all interactions for a specific residue, plus its RMSF
    and context within its region.
    """

    view_name = "residue"
    # D-B1: per-residue drill-down — SASA / interaction counts (static) +
    # persistence profile / chi entropy / RMSF (dynamic).
    spatial_layer = "residue"
    sub_flavors_served = ("static", "dynamic")
    max_length = 2500

    def render(self, context: RenderContext) -> str:
        """Render residue view."""
        locator = context.locator

        # Get required filter
        residue = context.get_filter("residue")
        if not residue:
            raise PresenterError(
                "Residue view requires filter: residue (e.g., 'ASP92', 'TRP97')"
            )

        lines = [f"# Residue: {residue}\n"]

        # Load RRCS data for interactions
        rrcs_root = locator.get_module_root("rrcs")
        if rrcs_root:
            csv_path = rrcs_root / "analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv"
            df = read_csv(csv_path)

            if not df.empty:
                lines.append(self._render_interactions(df, residue))
            else:
                lines.append("_Interaction data not available_\n")
        else:
            lines.append("_RRCS analysis not found_\n")

        # Load RMSF data
        rmsf_root = locator.get_module_root("rmsf")
        if rmsf_root:
            rmsf_path = rmsf_root / "analysis/rmsf/residue_rmsf.csv"
            rmsf_df = read_csv(rmsf_path)

            if not rmsf_df.empty:
                lines.append(self._render_flexibility(rmsf_df, residue))
            else:
                lines.append("_RMSF data not available_\n")
        else:
            lines.append("_RMSF analysis not found_\n")

        lines.append(self._render_design_features(locator.case_dir, residue))

        # Sources
        sources = []
        if rrcs_root:
            sources.append(csv_path.relative_to(locator.case_dir))
        if rmsf_root:
            sources.append(rmsf_path.relative_to(locator.case_dir))

        if sources:
            lines.append(self._format_sources(*sources))

        output = "\n".join(lines)
        return self._cap_length(output)

    def _render_design_features(self, case_dir: Path, residue: str) -> str:
        """Render design-layer features for this residue."""
        lines = ["## Design features\n"]
        try:
            feature_locator = FeatureCaseLocator(case_dir)
            features = FeatureSet(case_id=feature_locator.get_case_id())
            # D-B5 (2026-05-26): `mutability_score` no longer requested on the
            # LLM-facing residue view. Only evidence features + deterministic
            # categorical annotations (chemistry_tags, risk_flag) are pulled.
            for name in (
                "contact_count",
                "rrcs_contribution",
                "chemistry_tags",
                "risk_flag",
                "persistence_profile",
            ):
                compute_feature(name, feature_locator, features)
        except Exception as exc:
            return f"## Design features\n\n_Feature layer unavailable: {exc}_\n"

        matches = [
            item for item in features
            if residue.lower() in item.residue.label().lower()
            or residue.lower() == f"{item.residue.chain}{item.residue.resid}".lower()
        ]
        if not matches:
            return "## Design features\n\n_No design features found for this residue._\n"

        item = matches[0]
        if item.rrcs_contribution is not None:
            lines.append(f"- **RRCS contribution**: {num(item.rrcs_contribution, 2)}")
        if item.contact_count is not None:
            lines.append(f"- **Contact partners**: {item.contact_count}")
        if item.persistence_profile:
            lines.append(f"- **Persistence profile**: {item.persistence_profile}")
        if item.persistence_mean is not None:
            lines.append(f"- **Mean persistence**: {pct(item.persistence_mean)}")
        if item.region:
            lines.append(f"- **Region**: {item.region}")
        if item.chemistry_tags:
            lines.append(f"- **Chemistry tags**: {', '.join(item.chemistry_tags[:8])}")
        if item.risk_flags:
            lines.append(f"- **Risk flags**: {', '.join(item.risk_flags[:8])}")

        return "\n".join(lines) + "\n"

    def _render_interactions(self, df: pd.DataFrame, residue: str) -> str:
        """Render interactions for this residue."""
        lines = ["## Interactions\n"]

        # Filter for this residue (could be TCR or partner)
        if "tcr_residue_label" in df.columns:
            residue_df = df[
                df["tcr_residue_label"].str.contains(residue, case=False, na=False)
            ]
        else:
            lines.append("_Residue label column not found_")
            return "\n".join(lines)

        if residue_df.empty:
            lines.append(f"_No interactions found for {residue}_")
            return "\n".join(lines)

        # Sort by RRCS
        if "mean_rrcs" in residue_df.columns:
            residue_df = residue_df.sort_values("mean_rrcs", ascending=False)

        lines.append(f"Found {len(residue_df)} interaction partners:\n")

        # List interactions
        for _, row in residue_df.iterrows():
            partner = row.get("partner_residue_label", "")
            partner_comp = row.get("partner_component", "")
            mean_rrcs = row.get("mean_rrcs", 0)
            occupancy = row.get("rrcs_nonzero_fraction", 0)
            interaction_family = row.get("interaction_family", "")

            lines.append(
                f"- **{partner}** ({partner_comp}) — "
                f"RRCS {num(mean_rrcs, 2)}, {pct(occupancy)} occupancy"
            )

            if interaction_family:
                lines.append(f"  - Type: {interaction_family}")

            lines.append("")

        return "\n".join(lines)

    def _render_flexibility(self, df: pd.DataFrame, residue: str) -> str:
        """Render RMSF for this residue."""
        lines = ["## Flexibility\n"]

        # Find this residue
        if "residue" in df.columns:
            residue_row = df[
                df["residue"].str.contains(residue, case=False, na=False)
            ]
        else:
            lines.append("_Residue column not found_")
            return "\n".join(lines)

        if residue_row.empty:
            lines.append(f"_RMSF data not found for {residue}_")
            return "\n".join(lines)

        row = residue_row.iloc[0]
        rmsf = row.get("rmsf", 0)
        region = row.get("region", "")

        lines.append(f"- **RMSF**: {num(rmsf, 2, 'Å')}")
        if region:
            lines.append(f"- **Region**: {region}")

        # Interpretation
        if rmsf >= 3.0:
            lines.append("- **Interpretation**: High flexibility, tolerates mutations")
        elif rmsf >= 1.0:
            lines.append("- **Interpretation**: Moderate flexibility")
        else:
            lines.append("- **Interpretation**: Rigid, structurally constrained")

        return "\n".join(lines) + "\n"
