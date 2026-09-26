"""Hotspots view - top RRCS pairs."""

from pathlib import Path
import pandas as pd

from ..base import Presenter, RenderContext
from ..loaders import read_json, read_csv
from ..formatters import pair_display, pct, num, mk_table
from .. import register_view
from immunoscope.analysis.features import (
    CaseLocator as FeatureCaseLocator,
    FeatureSet,
    ResidueKey,
    compute_feature,
)


@register_view("hotspots")
class HotspotsPresenter(Presenter):
    """
    RRCS hotspots view.

    Shows top contact pairs ranked by mean RRCS, with breakdown by
    partner component (peptide, HLA_alpha, HLA_beta).
    """

    view_name = "hotspots"
    # D-B1: per-residue interaction counts (static) + persistence /
    # occupancy fractions (dynamic) — Residue layer, both sub-flavors.
    spatial_layer = "residue"
    sub_flavors_served = ("static", "dynamic")
    max_length = 5000

    def render(self, context: RenderContext) -> str:
        """Render hotspots view."""
        locator = context.locator
        case_id = locator.get_case_id()

        # Get filters
        top_n = min(int(context.get_filter("top_n", 10)), 10)
        min_rrcs = context.get_filter("min_rrcs", 2.0)
        min_occupancy = context.get_filter("min_occupancy", 0.5)
        region_filter = context.get_filter("region")
        partner_filter = context.get_filter("partner")

        lines = [f"# Hotspots: {case_id}\n"]

        # Load data
        rrcs_root = locator.get_module_root("rrcs")
        if not rrcs_root:
            lines.append("_RRCS analysis not found. Run RRCS analysis first._")
            return "\n".join(lines)

        # Try annotated CSV first (has all fields)
        csv_path = self._find_rrcs_pair_csv(locator, rrcs_root)
        df = read_csv(csv_path)

        if df.empty:
            # Fallback to rrcs_summary.json
            summary_path = rrcs_root / "analysis/interactions/rrcs/rrcs_summary.json"
            summary = read_json(summary_path)
            if summary:
                lines.append(self._render_from_summary(summary, top_n))
            else:
                lines.append("_RRCS data not available._")
            return "\n".join(lines)

        # Filter
        df = self._apply_filters(
            df, min_rrcs, min_occupancy, region_filter, partner_filter
        )

        if df.empty:
            lines.append(
                f"_No hotspots found with RRCS ≥ {min_rrcs}, "
                f"occupancy ≥ {min_occupancy*100:.0f}%_"
            )
            return "\n".join(lines)

        design_features = self._load_design_features(locator.case_dir)

        filtered_df = df.copy()
        top_df = filtered_df.sort_values("mean_rrcs", ascending=False).head(top_n)

        # Render sections
        lines.append(self._render_top_table(top_df, design_features))
        lines.append(self._render_partner_breakdown(filtered_df, top_df, design_features, top_n=10))
        lines.append(self._render_by_partner(top_df))

        # Sources
        lines.append(self._format_sources(
            csv_path.relative_to(locator.case_dir) if csv_path and csv_path.is_relative_to(locator.case_dir) else csv_path
        ))

        output = "\n".join(lines)
        return self._cap_length(output, "Use filters to narrow (region, partner, top_n)")

    def _find_rrcs_pair_csv(self, locator, rrcs_root: Path) -> Path | None:
        """Find annotated RRCS pair table across supported output layouts."""
        candidates = [
            rrcs_root / "analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv",
            rrcs_root / "annotated_rrcs_pair_summary.csv",
            locator.case_dir / "analysis/rrcs/analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv",
            locator.case_dir / "analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv",
            locator.case_dir / "overview/rrcs/annotated_rrcs_pair_summary.csv",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return candidates[0]

    def _apply_filters(
        self,
        df: pd.DataFrame,
        min_rrcs: float,
        min_occupancy: float,
        region_filter: str,
        partner_filter: str
    ) -> pd.DataFrame:
        """Apply filters to DataFrame."""
        # RRCS threshold
        if "mean_rrcs" in df.columns:
            df = df[df["mean_rrcs"] >= min_rrcs]

        # Occupancy threshold
        if "rrcs_nonzero_fraction" in df.columns:
            df = df[df["rrcs_nonzero_fraction"] >= min_occupancy]

        # Region filter
        if region_filter:
            if "tcr_region_detailed" in df.columns:
                df = df[
                    df["tcr_region_detailed"].str.contains(
                        region_filter, case=False, na=False
                    )
                ]
            elif "tcr_region" in df.columns:
                df = df[
                    df["tcr_region"].str.contains(
                        region_filter, case=False, na=False
                    )
                ]

        # Partner filter
        if partner_filter:
            if "partner_component" in df.columns:
                df = df[
                    df["partner_component"].str.contains(
                        partner_filter, case=False, na=False
                    )
                ]

        return df

    def _load_design_features(self, case_dir: Path) -> dict[ResidueKey, object]:
        """Load feature-layer records used as optional table annotations.

        D-B5 (2026-05-26): the `mutability_score` (design_priority_score)
        computer is no longer invoked here. Only the deterministic categorical
        annotations (`chemistry_tags`, `risk_flag`) are pulled in, per the
        D-B5 policy that the LLM-facing surface exposes evidence and
        categorical risk flags but **no** composite scalar ranking.
        """
        try:
            feature_locator = FeatureCaseLocator(case_dir)
            features = FeatureSet(case_id=feature_locator.get_case_id())
            for name in ("chemistry_tags", "risk_flag"):
                compute_feature(name, feature_locator, features)
            return dict(features.residues)
        except Exception:
            return {}

    def _render_top_table(self, df: pd.DataFrame, design_features: dict[ResidueKey, object]) -> str:
        """Render top pairs table."""
        lines = ["## Top pairs (by mean RRCS)\n"]

        # Build table
        table_data = []
        for _, row in df.iterrows():
            feature = design_features.get(self._tcr_key_from_row(row))
            table_data.append({
                "TCR residue": row.get("tcr_residue_label", ""),
                "Region": self._format_region(row),
                "Partner": row.get("partner_residue_label", ""),
                "Partner region": row.get("partner_component", ""),
                "Mean RRCS": f"{row.get('mean_rrcs', 0):.2f}",
                "Occupancy": pct(row.get("rrcs_nonzero_fraction", 0)),
                "Tags": self._format_feature_tags(feature),
            })

        # D-B5 (2026-05-26): the `Mutability` column (design_priority_score)
        # was removed. Ranking is RRCS-only; chemistry/risk tags remain.
        table_df = pd.DataFrame(table_data)
        lines.append(mk_table(
            table_df,
            columns=[
                "TCR residue", "Region", "Partner", "Partner region",
                "Mean RRCS", "Occupancy", "Tags",
            ],
            max_rows=10
        ))

        return "\n".join(lines)

    def _render_by_partner(self, df: pd.DataFrame) -> str:
        """Render breakdown by partner component."""
        lines = ["## By partner type\n"]

        if "partner_component" not in df.columns:
            lines.append("_Partner breakdown not available_")
            return "\n".join(lines)

        # Group by partner
        for partner in ["peptide", "HLA_alpha", "HLA_beta"]:
            partner_df = df[
                df["partner_component"].str.contains(partner, case=False, na=False)
            ]

            if partner_df.empty:
                continue

            n_pairs = len(partner_df)
            max_rrcs = partner_df["mean_rrcs"].max()

            lines.append(f"**{partner} contacts** ({n_pairs} pairs, max RRCS {max_rrcs:.2f}):")

            # List top 5 for this partner
            for _, row in partner_df.head(5).iterrows():
                tcr_res = row.get("tcr_residue_label", "")
                tcr_region = self._format_region(row)
                partner_res = row.get("partner_residue_label", "")
                mean_rrcs = row.get("mean_rrcs", 0)
                occupancy = row.get("rrcs_nonzero_fraction", 0)

                lines.append(
                    f"- {tcr_region} {tcr_res} ↔ {partner_res} ({partner}) — "
                    f"RRCS {mean_rrcs:.2f}, {pct(occupancy)} occupancy"
                )

            lines.append("")

        return "\n".join(lines)

    def _render_partner_breakdown(
        self,
        filtered_df: pd.DataFrame,
        top_df: pd.DataFrame,
        design_features: dict[ResidueKey, object],
        top_n: int = 5,
    ) -> str:
        """Render per-TCR-residue RRCS distribution across partner components."""
        lines = ["## By TCR residue\n"]

        required = {"tcr_residue_label", "partner_component", "mean_rrcs", "chain_id_1", "resid_1"}
        if not required.issubset(filtered_df.columns):
            lines.append("_Partner breakdown not available_")
            return "\n".join(lines)

        top_residues = list(dict.fromkeys(top_df["tcr_residue_label"].dropna().astype(str)))[:top_n]
        if not top_residues:
            lines.append("_No TCR residues available for breakdown_")
            return "\n".join(lines)

        work = filtered_df[filtered_df["tcr_residue_label"].astype(str).isin(top_residues)].copy()
        if work.empty:
            lines.append("_No TCR residues available for breakdown_")
            return "\n".join(lines)

        shown = 0
        for residue in top_residues:
            residue_df = work[work["tcr_residue_label"].astype(str) == residue].copy()
            if residue_df.empty:
                continue
            residue_df = residue_df.sort_values("mean_rrcs", ascending=False)
            feature = design_features.get(self._tcr_key_from_row(residue_df.iloc[0]))
            region = self._format_region(residue_df.iloc[0])
            flags = getattr(feature, "risk_flags", []) if feature is not None else []
            # D-B5 (2026-05-26): no composite mutability score; status reflects
            # presence/absence of categorical risk flags only.
            status = "risk-flagged" if flags else "no risk flags"
            lines.append(f"**{region} {residue}** ({status})")

            for partner, partner_df in residue_df.groupby("partner_component", sort=False):
                partner_df = partner_df.sort_values("mean_rrcs", ascending=False)
                best = partner_df.iloc[0]
                lines.append(
                    f"- {partner}: {len(partner_df)} pair"
                    f"{'' if len(partner_df) == 1 else 's'}, "
                    f"max RRCS {float(best.get('mean_rrcs', 0.0)):.2f} "
                    f"(↔ {best.get('partner_residue_label', '-')})"
                )
            lines.append("")
            shown += 1

        if len(top_residues) > shown:
            lines.append(f"_(showing top {shown} of {len(top_residues)})_")
        return "\n".join(lines)

    @staticmethod
    def _tcr_key_from_row(row) -> ResidueKey | None:
        try:
            chain = str(row.get("chain_id_1", "")).strip()
            resid = int(row.get("resid_1"))
            resname = str(row.get("resname_1", "")).strip().upper()
            if chain:
                return ResidueKey(chain=chain, resid=resid, resname=resname)
        except Exception:
            return None
        return None

    def _format_feature_tags(self, feature: object) -> str:
        if feature is None:
            return "-"
        tags = [self._short_tag(tag) for tag in getattr(feature, "chemistry_tags", [])[:3]]
        text = ", ".join(tag for tag in tags if tag)
        if getattr(feature, "risk_flags", []):
            text = f"⚠ {text}" if text else "⚠ risk"
        return text or "-"

    @staticmethod
    def _short_tag(tag: str) -> str:
        mapping = {
            "aromatic_residue": "aromatic",
            "positive_residue": "positive",
            "negative_residue": "negative",
            "polar_residue": "polar",
            "hydrophobic_residue": "hydrophobic",
            "hbond_network_hub": "hbond",
            "salt_bridge_participant": "saltbridge",
            "hydrophobic_contact": "hydrophobic-contact",
            "pi_interaction": "pi",
        }
        return mapping.get(str(tag), str(tag).replace("_", "-"))

    def _format_region(self, row) -> str:
        """Format region from row."""
        region = row.get("tcr_region_detailed") or row.get("tcr_region", "")
        # Abbreviate
        region = region.replace("_alpha", "α").replace("_beta", "β")
        region = region.replace("CDR", "CDR").replace("Framework", "FW")
        return region

    def _render_from_summary(self, summary: dict, top_n: int) -> str:
        """Fallback: render from rrcs_summary.json."""
        lines = ["## Top pairs\n"]

        top_pairs = summary.get("top_pairs", [])[:top_n]

        if not top_pairs:
            lines.append("_No hotspots found_")
            return "\n".join(lines)

        for pair in top_pairs:
            residues = pair.get("residues", "")
            mean_rrcs = pair.get("mean_rrcs", 0)
            occupancy = pair.get("nonzero_fraction", 0)

            lines.append(
                f"- **{residues}** — RRCS {mean_rrcs:.2f}, {pct(occupancy)} occupancy"
            )

        lines.append(
            "\n_Note: Limited detail available. "
            "Use annotated CSV for full breakdown._"
        )

        return "\n".join(lines)
