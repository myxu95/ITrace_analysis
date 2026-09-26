"""Pair view - single-pair dynamics fingerprint."""

from pathlib import Path
import pandas as pd

from ..base import Presenter, RenderContext, PresenterError
from ..loaders import read_csv
from ..formatters import pair_display, pct, num, classify_fingerprint
from .. import register_view
from immunoscope.analysis.features import (
    CaseLocator as FeatureCaseLocator,
    FeatureSet,
    ResidueKey,
    compute_feature,
)


@register_view("pair")
class PairPresenter(Presenter):
    """
    Pair view - single-pair dynamics fingerprint.

    Shows detailed dynamics for a specific residue pair: RRCS statistics,
    occupancy, interaction type, and dynamics classification.
    """

    view_name = "pair"
    # D-B1: per-pair RRCS mean/median/max (static) + occupancy fractions /
    # H-bond / salt-bridge / hydrophobic / π-occupancy time-series summaries
    # (dynamic). Pair is a first-class layer per D-B1 (not folded into
    # Residue) because per-pair contact persistence is an independent MD
    # signal that drives "drop pair" ablations.
    spatial_layer = "pair"
    sub_flavors_served = ("static", "dynamic")
    max_length = 5000

    def render(self, context: RenderContext) -> str:
        """Render pair view."""
        locator = context.locator

        # Get required filters
        residue1 = context.get_filter("residue1")
        residue2 = context.get_filter("residue2")

        if not residue1 or not residue2:
            raise PresenterError(
                "Pair view requires filters: residue1, residue2 "
                "(e.g., residue1='ASP92', residue2='LYS66')"
            )

        lines = [f"# Pair: {residue1} ↔ {residue2}\n"]

        # Load RRCS data
        rrcs_root = locator.get_module_root("rrcs")
        if not rrcs_root:
            lines.append("_RRCS analysis not found._")
            return "\n".join(lines)

        csv_path = self._find_rrcs_pair_csv(locator, rrcs_root)
        if csv_path is None:
            lines.append("_RRCS data not available._")
            return "\n".join(lines)
        df = read_csv(csv_path)

        if df.empty:
            lines.append("_RRCS data not available._")
            return "\n".join(lines)

        # Find the pair
        pair_row = self._find_pair(df, residue1, residue2)

        if pair_row is None:
            lines.append(f"_Pair not found in RRCS data. Check residue labels._")
            return "\n".join(lines)

        # Render sections
        feature_records = self._load_design_features(locator.case_dir)
        lines.append(self._render_statistics(pair_row))
        lines.append(self._render_interaction_type(pair_row))
        decomposition, decomposition_sources = self._load_interaction_decomposition(
            locator, pair_row
        )
        lines.append(self._render_interaction_breakdown(decomposition))
        lines.append(self._render_residue_context(pair_row, feature_records))
        lines.append(self._render_dynamics(pair_row))
        lines.append(self._render_context(pair_row))

        # Sources
        sources = [csv_path.relative_to(locator.case_dir)]
        sources.extend(decomposition_sources)
        lines.append(self._format_sources(*sources))

        output = "\n".join(lines)
        return self._cap_length(output)

    @staticmethod
    def _find_rrcs_pair_csv(locator, rrcs_root: Path) -> Path | None:
        candidates = [
            rrcs_root / "analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv",
            rrcs_root / "annotated_rrcs_pair_summary.csv",
            locator.case_dir / "analysis/rrcs/analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv",
            locator.case_dir / "analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv",
            locator.case_dir / "overview/rrcs/annotated_rrcs_pair_summary.csv",
        ]
        for path in candidates:
            if path.exists():
                return path
        return None

    def _find_pair(self, df: pd.DataFrame, res1: str, res2: str) -> pd.Series:
        """Find the pair in the DataFrame."""
        if "tcr_residue_label" not in df.columns or "partner_residue_label" not in df.columns:
            return None

        # Try both directions
        match = df[
            (df["tcr_residue_label"].str.contains(res1, case=False, na=False)) &
            (df["partner_residue_label"].str.contains(res2, case=False, na=False))
        ]

        if match.empty:
            # Try reverse
            match = df[
                (df["tcr_residue_label"].str.contains(res2, case=False, na=False)) &
                (df["partner_residue_label"].str.contains(res1, case=False, na=False))
            ]

        if match.empty:
            return None

        return match.iloc[0]

    def _render_statistics(self, row: pd.Series) -> str:
        """Render RRCS statistics."""
        lines = ["## Statistics\n"]

        mean_rrcs = row.get("mean_rrcs", 0)
        median_rrcs = row.get("median_rrcs", 0)
        max_rrcs = row.get("max_rrcs", 0)
        occupancy = row.get("rrcs_nonzero_fraction", 0)

        lines.append(f"- **Mean RRCS**: {num(mean_rrcs, 2)}")
        if median_rrcs:
            lines.append(f"- **Median RRCS**: {num(median_rrcs, 2)}")
        lines.append(f"- **Max RRCS**: {num(max_rrcs, 2)}")
        lines.append(f"- **Occupancy**: {pct(occupancy)} of frames")

        nonzero_frames = row.get("rrcs_nonzero_frames", 0)
        total_frames = row.get("total_frames", 0)
        if total_frames:
            lines.append(f"- **Frames**: {int(nonzero_frames)} / {int(total_frames)}")

        return "\n".join(lines) + "\n"

    def _render_interaction_type(self, row: pd.Series) -> str:
        """Render interaction type."""
        lines = ["## Interaction type\n"]

        interaction_family = row.get("interaction_family", "")
        if interaction_family:
            lines.append(f"- **Family**: {interaction_family}")

        interaction_class = row.get("interaction_class", "")
        if interaction_class:
            lines.append(f"- **Class**: {interaction_class}")

        return "\n".join(lines) + "\n"

    def _render_interaction_breakdown(self, rows: list[dict]) -> str:
        """Render fixed typed-interaction breakdown for this residue pair."""
        lines = ["## Interaction breakdown\n"]
        by_key = {row["key"]: row for row in rows}
        lines.append("| Interaction | Occupancy | Mean geom |")
        lines.append("| --- | --- | --- |")
        for spec in self._interaction_specs():
            row = by_key.get(spec["key"])
            if row is None:
                lines.append(f"| {spec['label']} | — | — |")
                continue
            occupancy = row.get("occupancy")
            geom = row.get("geom") or "—"
            occupancy_text = pct(float(occupancy)) if occupancy is not None else "—"
            lines.append(f"| {spec['label']} | {occupancy_text} | {geom} |")
        return "\n".join(lines) + "\n"

    def _load_interaction_decomposition(
        self,
        locator,
        pair_row: pd.Series,
    ) -> tuple[list[dict], list[Path]]:
        rows_by_key: dict[str, dict] = {}
        sources: list[Path] = []
        for spec in self._interaction_specs():
            path = self._find_interaction_pair_table(locator, spec)
            if path is None:
                continue
            df = read_csv(path)
            if df.empty:
                continue
            matches = df[df.apply(lambda row: self._same_residue_pair(pair_row, row), axis=1)]
            if matches.empty:
                continue
            occupancy = self._best_numeric(matches, "contact_frequency", "occupancy", "rrcs_nonzero_fraction")
            min_distance = self._best_numeric(
                matches, "min_distance_observed", "min_distance", "min_distance_angstrom"
            )
            mean_angle = self._best_numeric(
                matches, "mean_angle_observed", "min_geometry_angle_observed", "angle"
            )
            geom_parts = []
            if min_distance is not None:
                geom_parts.append(num(min_distance, 2, "Å"))
            if mean_angle is not None:
                geom_parts.append(num(mean_angle, 1, "deg"))
            rows_by_key[spec["key"]] = {
                "key": spec["key"],
                "label": spec["label"],
                "occupancy": occupancy if occupancy is not None else 0.0,
                "geom": ", ".join(geom_parts),
            }
            try:
                sources.append(path.relative_to(locator.case_dir))
            except ValueError:
                sources.append(path)

        return list(rows_by_key.values()), sources

    @staticmethod
    def _interaction_specs() -> list[dict]:
        return [
            {
                "key": "hbond",
                "label": "Hydrogen bond",
                "paths": [
                    "analysis/interactions/hydrogen_bonds/residue_pair_hbonds.csv",
                    "analysis/interactions/hydrogen_bonds/hbond_pair_summary.csv",
                    "analysis/interactions/hydrogen_bonds/residue_pair_hydrogen_bonds.csv",
                    "analysis/**/residue_pair_hbonds.csv",
                    "analysis/**/hbond_pair_summary.csv",
                    "analysis/**/residue_pair_hydrogen_bonds.csv",
                    "analysis/**/*hbond*pair*.csv",
                    "analysis/**/*hydrogen*pair*.csv",
                ],
            },
            {
                "key": "salt_bridge",
                "label": "Salt bridge",
                "paths": [
                    "analysis/interactions/salt_bridges/salt_bridge_pair_summary.csv",
                    "analysis/interactions/salt_bridges/residue_pair_salt_bridges.csv",
                    "analysis/**/salt_bridge_pair_summary.csv",
                    "analysis/**/residue_pair_salt_bridges.csv",
                    "analysis/**/*salt*pair*.csv",
                ],
            },
            {
                "key": "hydrophobic",
                "label": "Hydrophobic",
                "paths": [
                    "analysis/interactions/hydrophobic_contacts/hydrophobic_pair_summary.csv",
                    "analysis/interactions/hydrophobic_contacts/residue_pair_hydrophobic_contacts.csv",
                    "analysis/**/hydrophobic_pair_summary.csv",
                    "analysis/**/residue_pair_hydrophobic_contacts.csv",
                    "analysis/**/*hydrophobic*pair*.csv",
                ],
            },
            {
                "key": "pi_pi",
                "label": "π-π",
                "paths": [
                    "analysis/interactions/pi_interactions/pi_pi_pair_summary.csv",
                    "analysis/interactions/pi_interactions/residue_pair_pi_pi.csv",
                    "analysis/**/pi_pi_pair_summary.csv",
                    "analysis/**/residue_pair_pi_pi.csv",
                    "analysis/**/*pi_pi*pair*.csv",
                ],
            },
            {
                "key": "cation_pi",
                "label": "Cation-π",
                "paths": [
                    "analysis/interactions/cation_pi_interactions/cation_pi_pair_summary.csv",
                    "analysis/interactions/cation_pi_interactions/residue_pair_cation_pi.csv",
                    "analysis/**/cation_pi_pair_summary.csv",
                    "analysis/**/residue_pair_cation_pi.csv",
                    "analysis/**/*cation*pi*pair*.csv",
                ],
            },
        ]

    @staticmethod
    def _find_interaction_pair_table(locator, spec: dict) -> Path | None:
        seen: set[Path] = set()
        for pattern in spec["paths"]:
            matches = sorted(locator.case_dir.glob(pattern))
            for path in matches:
                resolved = path.resolve()
                if resolved in seen or not path.is_file():
                    continue
                seen.add(resolved)
                return path
        return None

    @staticmethod
    def _same_residue_pair(rrcs_row: pd.Series, interaction_row: pd.Series) -> bool:
        rrcs_pairs = {
            (
                str(rrcs_row.get("chain_id_1", "")),
                str(rrcs_row.get("resid_1", "")),
            ),
            (
                str(rrcs_row.get("chain_id_2", "")),
                str(rrcs_row.get("resid_2", "")),
            ),
        }
        interaction_pairs = {
            (
                str(interaction_row.get("chain_id_1", "")),
                str(interaction_row.get("resid_1", "")),
            ),
            (
                str(interaction_row.get("chain_id_2", "")),
                str(interaction_row.get("resid_2", "")),
            ),
        }
        if rrcs_pairs == interaction_pairs and all(chain and resid for chain, resid in interaction_pairs):
            return True

        rrcs_labels = {
            PairPresenter._normalize_label(rrcs_row.get("tcr_residue_label")),
            PairPresenter._normalize_label(rrcs_row.get("partner_residue_label")),
        }
        interaction_labels = {
            PairPresenter._normalize_label(interaction_row.get("residue_label_1")),
            PairPresenter._normalize_label(interaction_row.get("residue_label_2")),
            PairPresenter._normalize_label(interaction_row.get("phla_residue")),
            PairPresenter._normalize_label(interaction_row.get("tcr_residue")),
        }
        interaction_labels.discard("")
        rrcs_labels.discard("")
        return bool(rrcs_labels) and rrcs_labels.issubset(interaction_labels)

    @staticmethod
    def _normalize_label(value: object) -> str:
        if value is None or pd.isna(value):
            return ""
        return "".join(ch.lower() for ch in str(value) if ch.isalnum())

    @staticmethod
    def _best_numeric(df: pd.DataFrame, *columns: str) -> float | None:
        for column in columns:
            if column not in df.columns:
                continue
            values = pd.to_numeric(df[column], errors="coerce").dropna()
            if not values.empty:
                return float(values.max())
        return None

    def _render_residue_context(
        self,
        row: pd.Series,
        feature_records: dict[ResidueKey, object],
    ) -> str:
        """Render feature-layer context for both residues in the pair."""
        lines = ["## Residue context\n"]
        tcr_key = self._residue_key_from_row(row, side=1)
        partner_key = self._residue_key_from_row(row, side=2)
        tcr_feature = feature_records.get(tcr_key) if tcr_key is not None else None
        partner_feature = feature_records.get(partner_key) if partner_key is not None else None

        lines.extend(self._render_one_residue_context("TCR side", row.get("tcr_residue_label", ""), tcr_feature))
        lines.append("")
        lines.extend(self._render_one_residue_context("Partner side", row.get("partner_residue_label", ""), partner_feature))
        lines.append("")
        lines.append(self._pair_verdict(tcr_feature, partner_feature, row))
        return "\n".join(lines) + "\n"

    def _render_one_residue_context(self, side: str, label: object, feature: object) -> list[str]:
        # D-B5 (2026-05-26): no composite mutability score; descriptor lists
        # only deterministic categorical evidence (pocket chemistry, chemistry
        # tags, risk flags) — the LLM does its own integration.
        pocket = getattr(feature, "pocket_chemistry", None) if feature is not None else None
        tags = getattr(feature, "chemistry_tags", []) if feature is not None else []
        risks = getattr(feature, "risk_flags", []) if feature is not None else []
        lines = [f"**{side} — {label or '-'}**"]
        lines.append(f"- Pocket: {pocket or '-'}")
        lines.append(f"- Tags: {', '.join(tags[:3]) if tags else '-'}")
        if risks:
            lines.append(f"- Risk: {', '.join(risks)}")
        return lines

    @staticmethod
    def _load_design_features(case_dir: Path) -> dict[ResidueKey, object]:
        # D-B5 (2026-05-26): drop `mutability_score`; keep evidence + flags.
        try:
            feature_locator = FeatureCaseLocator(case_dir)
            features = FeatureSet(case_id=feature_locator.get_case_id())
            for name in ("chemistry_tags", "risk_flag"):
                compute_feature(name, feature_locator, features)
            return dict(features.residues)
        except Exception:
            return {}

    @staticmethod
    def _residue_key_from_row(row: pd.Series, side: int) -> ResidueKey | None:
        try:
            chain = str(row.get(f"chain_id_{side}", "")).strip()
            resid = int(row.get(f"resid_{side}"))
            resname = str(row.get(f"resname_{side}", "")).strip().upper()
            if chain:
                return ResidueKey(chain=chain, resid=resid, resname=resname)
        except Exception:
            return None
        return None

    @staticmethod
    def _risk_class(feature: object) -> str:
        """Categorical class derived solely from `risk_flags` presence.

        D-B5 (2026-05-26): replaces the old `_mutability_class` which keyed on
        `design_priority_score` bins. The class is now binary (risk-flagged /
        no-risk) and used only to phrase the pair verdict.
        """
        if feature is None:
            return "unknown"
        if getattr(feature, "risk_flags", []):
            return "risk-flagged"
        return "no-risk"

    def _pair_verdict(self, tcr_feature: object, partner_feature: object, row: pd.Series) -> str:
        tcr_class = self._risk_class(tcr_feature)
        partner_class = self._risk_class(partner_feature)
        partner_component = str(row.get("partner_component", "")).lower()

        if tcr_class == "risk-flagged" and partner_class == "risk-flagged":
            return "_Both sides carry risk flags — this is a high-stakes pair; review flags before mutating._"
        if "hla" in partner_component or "mhc" in partner_component:
            return "_Partner is HLA/MHC; prefer TCR-side mutations over editing the partner._"
        if "peptide" in partner_component:
            return "_Partner is peptide; peptide-side and TCR-side mutations are both admissible — inspect risk flags._"
        if tcr_class == "no-risk" and partner_class == "no-risk":
            return "_Neither side is risk-flagged — both sides are admissible mutation targets._"
        return "_Inspect risk flags on the flagged side before committing a mutation._"

    def _render_dynamics(self, row: pd.Series) -> str:
        """Render dynamics fingerprint."""
        lines = ["## Dynamics fingerprint\n"]

        mean_rrcs = row.get("mean_rrcs", 0)
        occupancy = row.get("rrcs_nonzero_fraction", 0)
        max_rrcs = row.get("max_rrcs", 0)

        fingerprint = classify_fingerprint(mean_rrcs, occupancy, max_rrcs)

        fingerprint_desc = {
            "stable": "Stable contact — high occupancy, consistent strength",
            "fluctuating": "Fluctuating — strong when present, but intermittent",
            "weak-persistent": "Weak but persistent — always present, low strength",
            "highly-fluctuating": "Highly fluctuating — large strength variations",
            "moderate": "Moderate contact — typical interface interaction",
            "weak": "Weak contact — low strength, low occupancy",
        }

        desc = fingerprint_desc.get(fingerprint, "Unknown pattern")
        lines.append(f"- **Classification**: {fingerprint}")
        lines.append(f"- **Description**: {desc}")

        return "\n".join(lines) + "\n"

    def _render_context(self, row: pd.Series) -> str:
        """Render context information."""
        lines = ["## Context\n"]

        tcr_region = row.get("tcr_region_detailed") or row.get("tcr_region", "")
        partner_comp = row.get("partner_component", "")

        if tcr_region:
            lines.append(f"- **TCR region**: {tcr_region}")
        if partner_comp:
            lines.append(f"- **Partner component**: {partner_comp}")

        return "\n".join(lines) + "\n"
