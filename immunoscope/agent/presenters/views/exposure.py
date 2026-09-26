"""Exposure view — per-residue SASA and cavity-candidate ranking.

Reads `analysis/interface/residue_sasa.csv` produced by `ims residue_sasa`
and surfaces the per-residue SASA fields that the chain-level BSA in
the `interface` view cannot answer:

  - **Cavity candidates**: residues classified as `interface_core` whose
    `relative_exposure` is anomalously high suggest a packing defect at
    the buried face of the interface — a site where solvent has infiltrated
    or where the sidechain leaves a void. These are the residues a mutation
    designer reaches for when proposing bulkier / more complementary
    substitutions to fill the cavity.
  - **Burial-state breakdown**: counts per burial_state (core /
    interface_core / interface_rim / exposed) so the agent can see how
    much of the interface is core-buried vs. rim-exposed.

D-B1 (locked 2026-05-26): spatial_layer='residue', sub_flavors_served=
('static',). The SASA fields are registered at the residue × static cell
in `analysis/feature_registry.py:123-128`.
"""

from __future__ import annotations

import pandas as pd

from ..base import Presenter, RenderContext
from ..formatters import num
from ..loaders import read_csv
from .. import register_view


# Cavity-anomaly threshold: an interface_core residue whose relative
# exposure exceeds this is flagged as a candidate filler site. The cutoff
# is intentionally loose so the view surfaces signals rather than
# pre-filtering them; the agent decides which are real cavities.
_CAVITY_RELEXP_CUTOFF = 0.20


@register_view("exposure")
class ExposurePresenter(Presenter):
    """Per-residue SASA / burial-state with cavity-candidate ranking."""

    view_name = "exposure"
    max_length = 4000
    spatial_layer = "residue"
    sub_flavors_served = ("static",)

    def render(self, context: RenderContext) -> str:
        locator = context.locator
        case_id = locator.get_case_id()
        top_n = int(context.get_filter("top_n", 15))

        lines = [f"# Exposure / SASA: {case_id}\n"]

        interface_root = locator.get_module_root("interface")
        if not interface_root:
            lines.append("_Interface module not found. Run `ims residue_sasa` first._")
            return "\n".join(lines)

        csv_path = interface_root / "residue_sasa.csv"
        df = read_csv(csv_path)
        if df.empty:
            lines.append("_Per-residue SASA data not available._")
            return "\n".join(lines)

        df = self._normalize(df)

        lines.append(self._render_burial_breakdown(df))
        lines.append(self._render_cavity_candidates(df, top_n))
        lines.append(self._render_most_exposed(df, top_n))
        lines.append(self._render_interpretation())
        lines.append(self._format_sources(csv_path.relative_to(locator.case_dir)))

        return self._cap_length(
            "\n".join(lines),
            "Use filters: top_n",
        )

    @staticmethod
    def _normalize(df: pd.DataFrame) -> pd.DataFrame:
        """Build a uniform `residue` label column and drop unusable rows."""
        if {"chain_id", "resname", "resid"}.issubset(df.columns):
            df = df.copy()
            df["residue"] = (
                df["resname"].astype(str)
                + df["resid"].astype("Int64").astype(str)
                + " ("
                + df["chain_id"].astype(str)
                + ")"
            )
        elif "residue" not in df.columns:
            df = df.copy()
            df["residue"] = df.index.astype(str)
        return df

    def _render_burial_breakdown(self, df: pd.DataFrame) -> str:
        lines = ["## Burial-state breakdown\n"]
        if "burial_state" not in df.columns:
            lines.append("_burial_state column missing._\n")
            return "\n".join(lines)
        counts = df["burial_state"].fillna("(unknown)").value_counts()
        lines.append("| Burial state | Count |")
        lines.append("|--------------|-------|")
        for state, n in counts.items():
            lines.append(f"| {state} | {int(n)} |")
        return "\n".join(lines) + "\n"

    def _render_cavity_candidates(self, df: pd.DataFrame, top_n: int) -> str:
        lines = [
            "## Cavity candidates "
            f"(interface_core residues with relative_exposure ≥ {_CAVITY_RELEXP_CUTOFF})\n",
            "_High exposure at the buried face of the interface implies a packing "
            "defect — solvent infiltration or a sidechain void. Filling it with a "
            "bulkier / more complementary residue can tighten the interface._\n",
        ]
        if not {"burial_state", "relative_exposure"}.issubset(df.columns):
            lines.append("_Insufficient columns to compute cavity candidates._")
            return "\n".join(lines) + "\n"

        cavities = df[
            (df["burial_state"] == "interface_core")
            & (df["relative_exposure"] >= _CAVITY_RELEXP_CUTOFF)
        ].sort_values("relative_exposure", ascending=False)

        if cavities.empty:
            lines.append("_No interface_core residue exceeds the cavity threshold._")
            return "\n".join(lines) + "\n"

        lines.append(
            "| Residue | relative_exposure | sasa_bound (Å²) | "
            "delta_sasa (Å²) | sidechain_sasa_bound (Å²) |"
        )
        lines.append("|---------|-------------------|-----------------|"
                     "-----------------|---------------------------|")
        for _, row in cavities.head(top_n).iterrows():
            lines.append(
                f"| {row['residue']} "
                f"| {num(row.get('relative_exposure', 0), 2)} "
                f"| {num(row.get('sasa_bound', 0), 1)} "
                f"| {num(row.get('delta_sasa', 0), 1)} "
                f"| {num(row.get('sasa_sidechain_bound', 0), 1)} |"
            )
        if len(cavities) > top_n:
            lines.append(f"\n_Showing {top_n} of {len(cavities)} cavity candidates._")
        return "\n".join(lines) + "\n"

    def _render_most_exposed(self, df: pd.DataFrame, top_n: int) -> str:
        """Top residues by absolute delta_sasa — the ones that lose the most
        SASA on binding. These are not cavities but the strongest 'committed'
        burial events; useful context for the agent to anchor its picture
        of which residues actually carry the interface."""
        lines = [
            "## Top burial events (highest |delta_sasa|)\n",
            "_Residues that lose the most SASA on binding — the load-bearing "
            "burial sites._\n",
        ]
        if "delta_sasa" not in df.columns:
            lines.append("_delta_sasa column missing._")
            return "\n".join(lines) + "\n"

        burial = df.assign(_abs_delta=df["delta_sasa"].abs()).sort_values(
            "_abs_delta", ascending=False
        )
        burial = burial[burial["_abs_delta"] > 0]
        if burial.empty:
            lines.append("_No nonzero delta_sasa rows._")
            return "\n".join(lines) + "\n"

        lines.append("| Residue | delta_sasa (Å²) | sasa_bound (Å²) | burial_state |")
        lines.append("|---------|-----------------|-----------------|--------------|")
        for _, row in burial.head(top_n).iterrows():
            lines.append(
                f"| {row['residue']} "
                f"| {num(row.get('delta_sasa', 0), 1)} "
                f"| {num(row.get('sasa_bound', 0), 1)} "
                f"| {row.get('burial_state', '-')} |"
            )
        return "\n".join(lines) + "\n"

    @staticmethod
    def _render_interpretation() -> str:
        return (
            "## Interpretation\n\n"
            "- **interface_core + high relative_exposure**: cavity / packing "
            "defect; mutation candidates → bulkier or chemically complementary "
            "residues to fill the void.\n"
            "- **interface_rim + high relative_exposure**: rim residue with "
            "solvent access; mutations here primarily tune specificity rather "
            "than affinity.\n"
            "- **large |delta_sasa|**: load-bearing burial site; mutating away "
            "from the wild-type sidechain risks losing affinity unless the "
            "replacement preserves the burial.\n"
        )
