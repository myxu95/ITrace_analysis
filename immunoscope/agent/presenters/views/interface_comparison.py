"""Interface comparison view - residue-level diff between two cases."""

from pathlib import Path

import pandas as pd

from .. import register_view
from ..base import Presenter, RenderContext, PresenterError
from ..formatters import mk_table, num
from ..loaders import read_csv
from immunoscope.analysis.features import CaseLocator


@register_view("interface_comparison")
class InterfaceComparisonPresenter(Presenter):
    """
    Residue-level interface comparison view.

    Compares two case directories by aggregating RRCS pair tables into TCR
    residue-level contribution profiles, then renders the largest gains and
    losses in LLM-friendly Markdown.
    """

    view_name = "interface_comparison"
    max_length = 5000
    # D-B1: residue-level diff is built by aggregating Pair-layer RRCS
    # contributions across two cases; the reported quantity is the residue
    # delta, so this view lives at the Residue layer (static — reports the
    # per-residue scalar delta, not a distribution shape).
    spatial_layer = "residue"
    sub_flavors_served = ("static",)

    def render(self, context: RenderContext) -> str:
        """Render residue-level interface differences between two cases."""
        case_b_dir = context.get_filter("case_b_dir") or context.get_filter("case_b")
        if not case_b_dir:
            raise PresenterError(
                "interface_comparison view requires filter: case_b_dir"
            )

        top_n = int(context.get_filter("top_n", 10))
        min_abs_delta = float(context.get_filter("min_abs_delta", 0.0))
        case_a_label = context.get_filter("case_a_label", context.locator.get_case_id())

        try:
            case_b_locator = CaseLocator(Path(case_b_dir))
        except Exception as exc:
            raise PresenterError(f"Cannot open comparison case: {case_b_dir}") from exc

        case_b_label = context.get_filter("case_b_label", case_b_locator.get_case_id())

        case_a_profile, case_a_source = self._load_residue_profile(context.locator)
        case_b_profile, case_b_source = self._load_residue_profile(case_b_locator)

        lines = [f"# Interface comparison: {case_a_label} vs {case_b_label}\n"]

        if case_a_profile.empty:
            lines.append(f"_RRCS residue profile not available for {case_a_label}._")
            return "\n".join(lines)
        if case_b_profile.empty:
            lines.append(f"_RRCS residue profile not available for {case_b_label}._")
            return "\n".join(lines)

        merged = self._merge_profiles(case_a_profile, case_b_profile, case_a_label, case_b_label)
        if min_abs_delta > 0:
            merged = merged[merged["abs_delta"] >= min_abs_delta]

        if merged.empty:
            lines.append("_No residue-level RRCS differences passed the filters._")
            return "\n".join(lines)

        lines.append(self._render_summary(merged, case_a_label, case_b_label))
        lines.append(self._render_delta_table(
            merged.sort_values("delta", ascending=False),
            title=f"## Increased in {case_b_label}",
            case_a_label=case_a_label,
            case_b_label=case_b_label,
            top_n=top_n,
        ))
        lines.append(self._render_delta_table(
            merged.sort_values("delta", ascending=True),
            title=f"## Decreased in {case_b_label}",
            case_a_label=case_a_label,
            case_b_label=case_b_label,
            top_n=top_n,
        ))

        sources = []
        for locator, source in ((context.locator, case_a_source), (case_b_locator, case_b_source)):
            if source:
                try:
                    sources.append(source.relative_to(locator.case_dir))
                except ValueError:
                    sources.append(source)
        lines.append(self._format_sources(*sources))

        return self._cap_length("\n".join(lines), "Use top_n or min_abs_delta to narrow")

    def _load_residue_profile(self, locator: CaseLocator) -> tuple[pd.DataFrame, Path | None]:
        """Load and aggregate RRCS pair table into residue-level profile."""
        path = self._find_rrcs_pair_csv(locator)
        df = read_csv(path)
        if df.empty or "mean_rrcs" not in df.columns:
            return pd.DataFrame(), path

        residue_col = self._first_column(df, "tcr_residue_label", "residue_1", "residue")
        if residue_col is None:
            return pd.DataFrame(), path

        region_col = self._first_column(df, "tcr_region_detailed", "tcr_region", "region")
        partner_col = self._first_column(df, "partner_component", "component_2", "partner_region")
        occupancy_col = self._first_column(df, "rrcs_nonzero_fraction", "occupancy")

        work = df.copy()
        work["residue"] = work[residue_col].astype(str)
        work["mean_rrcs"] = pd.to_numeric(work["mean_rrcs"], errors="coerce").fillna(0.0)
        if occupancy_col:
            work["occupancy"] = pd.to_numeric(work[occupancy_col], errors="coerce").fillna(0.0)
        else:
            work["occupancy"] = 0.0

        grouped = work.groupby("residue", dropna=False).agg(
            total_rrcs=("mean_rrcs", "sum"),
            max_pair_rrcs=("mean_rrcs", "max"),
            mean_occupancy=("occupancy", "mean"),
            n_pairs=("mean_rrcs", "size"),
        )

        if region_col:
            region = work.groupby("residue")[region_col].agg(self._mode_text)
            grouped["region"] = region
        else:
            grouped["region"] = ""

        if partner_col:
            partner = work.groupby("residue")[partner_col].agg(self._partner_breakdown_text)
            grouped["partner_breakdown"] = partner
        else:
            grouped["partner_breakdown"] = ""

        return grouped.reset_index(), path

    def _find_rrcs_pair_csv(self, locator: CaseLocator) -> Path | None:
        root = locator.get_module_root("rrcs")
        candidates = []
        if root:
            candidates.extend([
                root / "analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv",
                root / "annotated_rrcs_pair_summary.csv",
            ])
        candidates.extend([
            locator.case_dir / "analysis/rrcs/analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv",
            locator.case_dir / "analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv",
            locator.case_dir / "overview/rrcs/annotated_rrcs_pair_summary.csv",
        ])
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return candidates[0] if candidates else None

    @staticmethod
    def _first_column(df: pd.DataFrame, *names: str) -> str | None:
        for name in names:
            if name in df.columns:
                return name
        return None

    @staticmethod
    def _mode_text(values: pd.Series) -> str:
        clean = [str(value) for value in values.dropna() if str(value)]
        if not clean:
            return ""
        return pd.Series(clean).mode().iloc[0]

    @staticmethod
    def _partner_breakdown_text(values: pd.Series) -> str:
        clean = [str(value) for value in values.dropna() if str(value)]
        if not clean:
            return ""
        counts = pd.Series(clean).value_counts()
        return ", ".join(f"{key}:{value}" for key, value in counts.head(4).items())

    @staticmethod
    def _merge_profiles(
        case_a: pd.DataFrame,
        case_b: pd.DataFrame,
        case_a_label: str,
        case_b_label: str,
    ) -> pd.DataFrame:
        merged = case_a.merge(
            case_b,
            on="residue",
            how="outer",
            suffixes=(f"_{case_a_label}", f"_{case_b_label}"),
        )
        for label in (case_a_label, case_b_label):
            for column in ("total_rrcs", "max_pair_rrcs", "mean_occupancy", "n_pairs"):
                name = f"{column}_{label}"
                if name in merged.columns:
                    merged[name] = pd.to_numeric(merged[name], errors="coerce").fillna(0.0)

        merged["delta"] = merged[f"total_rrcs_{case_b_label}"] - merged[f"total_rrcs_{case_a_label}"]
        merged["abs_delta"] = merged["delta"].abs()
        merged["region"] = merged.get(f"region_{case_b_label}", "").fillna("")
        if f"region_{case_a_label}" in merged.columns:
            merged["region"] = merged["region"].where(
                merged["region"].astype(bool),
                merged[f"region_{case_a_label}"].fillna(""),
            )
        merged["partner_breakdown"] = merged.get(f"partner_breakdown_{case_b_label}", "").fillna("")
        return merged.sort_values("abs_delta", ascending=False)

    def _render_summary(self, merged: pd.DataFrame, case_a_label: str, case_b_label: str) -> str:
        """Render compact comparison summary."""
        total_a = merged[f"total_rrcs_{case_a_label}"].sum()
        total_b = merged[f"total_rrcs_{case_b_label}"].sum()
        increased = int((merged["delta"] > 0).sum())
        decreased = int((merged["delta"] < 0).sum())
        return "\n".join([
            "## Summary\n",
            f"- **Total residue RRCS**: {case_a_label} {num(total_a, 2)} -> {case_b_label} {num(total_b, 2)}",
            f"- **Residues increased/decreased**: {increased} / {decreased}",
            f"- **Largest absolute shift**: {num(float(merged['abs_delta'].max()), 2)}",
        ])

    def _render_delta_table(
        self,
        df: pd.DataFrame,
        title: str,
        case_a_label: str,
        case_b_label: str,
        top_n: int,
    ) -> str:
        """Render delta table."""
        rows = []
        for _, row in df.head(top_n).iterrows():
            if row["delta"] == 0:
                continue
            rows.append({
                "Residue": row.get("residue", ""),
                "Region": row.get("region", ""),
                case_a_label: f"{row.get(f'total_rrcs_{case_a_label}', 0):.2f}",
                case_b_label: f"{row.get(f'total_rrcs_{case_b_label}', 0):.2f}",
                "Delta": f"{row.get('delta', 0):+.2f}",
                "Partners": row.get("partner_breakdown", ""),
            })
        if not rows:
            return f"{title}\n\n_No non-zero deltas._\n"
        return "\n".join([
            title,
            "",
            mk_table(
                pd.DataFrame(rows),
                columns=["Residue", "Region", case_a_label, case_b_label, "Delta", "Partners"],
                max_rows=top_n,
            ),
        ])

