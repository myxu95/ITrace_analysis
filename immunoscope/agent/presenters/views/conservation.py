"""Conservation view — per-residue TCR germline conservation.

Reads `analysis/conservation/conservation_tcr.csv` (produced by the
`TcrConservationNode` / `write_conservation_csv`) and surfaces how
conserved each TCR residue is across the IMGT germline V-gene reference:

  - **Per-region summary**: framework (FR) regions are evolutionarily
    load-bearing and score high; CDR1/CDR2 are germline-variable and score
    lower. This is the headline signal a mutation designer wants.
  - **Most-conserved residues**: high-confidence "do-not-touch" sites —
    mutating them risks destabilising the Ig fold.
  - **Most-variable residues**: naturally tolerant positions — a wider
    substitution menu is admissible (affinity-tuning knobs).

Scope caveat surfaced in the view: germline V conservation does NOT cover
the somatic CDR3 junction (those residues are reported `covered=False`);
only the germline-encoded CDR3 anchor carries a score.

D-B1: spatial_layer='residue', sub_flavors_served=('static',). Conservation
is sequence-derived (not an MD trajectory feature), so it is intentionally
NOT part of the locked MD feature_registry grid; the residue/static view
cell is already covered there, so this view is registry-consistent.
"""

from __future__ import annotations

import pandas as pd

from ..base import Presenter, RenderContext
from ..formatters import num
from ..loaders import read_csv
from .. import register_view


# Conservation-class cutoffs (mirror tcr_germline_conservation thresholds).
_CONSERVED_MIN = 0.80
_VARIABLE_MAX = 0.50


@register_view("conservation")
class ConservationPresenter(Presenter):
    """Per-residue TCR germline conservation (FR / CDR1 / CDR2)."""

    view_name = "conservation"
    max_length = 4000
    spatial_layer = "residue"
    sub_flavors_served = ("static",)

    def render(self, context: RenderContext) -> str:
        locator = context.locator
        case_id = locator.get_case_id()
        top_n = int(context.get_filter("top_n", 15))

        lines = [f"# TCR germline conservation: {case_id}\n"]

        root = locator.get_module_root("conservation")
        if not root:
            lines.append(
                "_Conservation module not found. Run the TcrConservationNode "
                "(or `write_conservation_csv`) first._"
            )
            return "\n".join(lines)

        csv_path = root / "conservation_tcr.csv"
        df = read_csv(csv_path)
        if df.empty:
            lines.append(
                "_TCR conservation data not available (empty CSV — ANARCI "
                "unavailable or TCR chains unresolved)._"
            )
            return "\n".join(lines)

        try:
            df = self._normalize(df)
            covered = df[df["covered"] == True]  # noqa: E712 — pandas mask
            lines.append(self._render_region_summary(covered))
            lines.append(self._render_most_conserved(covered, top_n))
            lines.append(self._render_most_variable(covered, top_n))
            lines.append(self._render_coverage_note(df))
        except Exception as exc:  # noqa: BLE001 — malformed/partial CSV
            lines.append(
                f"_Conservation CSV present but could not be rendered "
                f"({type(exc).__name__}: {exc}). The file may be malformed._"
            )
            return "\n".join(lines)
        lines.append(self._render_interpretation())
        lines.append(self._format_sources(csv_path.relative_to(locator.case_dir)))

        return self._cap_length("\n".join(lines), "Use filters: top_n")

    @staticmethod
    def _normalize(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        if {"chain_id", "resname", "resid"}.issubset(df.columns):
            df["residue"] = (
                df["resname"].astype(str)
                + df["resid"].astype("Int64").astype(str)
                + " ("
                + df["chain_id"].astype(str)
                + ")"
            )
        elif "residue" not in df.columns:
            df["residue"] = df.index.astype(str)
        # coerce score to numeric (covered rows have it; uncovered are blank)
        if "conservation_score" in df.columns:
            df["conservation_score"] = pd.to_numeric(
                df["conservation_score"], errors="coerce"
            )
        if "covered" not in df.columns:
            df["covered"] = df.get("conservation_score").notna()
        return df

    def _render_region_summary(self, df: pd.DataFrame) -> str:
        lines = ["## Per-region conservation (mean over covered residues)\n"]
        if df.empty or "region" not in df.columns:
            lines.append("_No covered residues to summarise._\n")
            return "\n".join(lines)
        lines.append("| Region | n | mean conservation |")
        lines.append("|--------|---|-------------------|")
        order = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
        present = [r for r in order if r in set(df["region"])]
        present += [r for r in sorted(set(df["region"])) if r not in order]
        for reg in present:
            sub = df[df["region"] == reg]["conservation_score"].dropna()
            if len(sub):
                lines.append(f"| {reg} | {len(sub)} | {num(sub.mean(), 3)} |")
        return "\n".join(lines) + "\n"

    def _render_most_conserved(self, df: pd.DataFrame, top_n: int) -> str:
        lines = [
            f"## Most-conserved residues (score ≥ {_CONSERVED_MIN} → avoid)\n",
            "_Evolutionarily load-bearing — mutating these risks destabilising "
            "the fold. Treat as do-not-touch unless evidence is strong._\n",
        ]
        return self._residue_table(
            df[df["conservation_score"] >= _CONSERVED_MIN].sort_values(
                "conservation_score", ascending=False
            ),
            top_n,
            lines,
            "_No residue reaches the conserved threshold._",
        )

    def _render_most_variable(self, df: pd.DataFrame, top_n: int) -> str:
        lines = [
            f"## Most-variable residues (score < {_VARIABLE_MAX} → tolerant)\n",
            "_Germline-variable positions — a wider substitution menu is "
            "admissible; candidate affinity-tuning knobs._\n",
        ]
        return self._residue_table(
            df[df["conservation_score"] < _VARIABLE_MAX].sort_values(
                "conservation_score", ascending=True
            ),
            top_n,
            lines,
            "_No residue falls below the variable threshold._",
        )

    @staticmethod
    def _residue_table(
        sub: pd.DataFrame, top_n: int, lines: list, empty_msg: str
    ) -> str:
        if sub.empty:
            lines.append(empty_msg)
            return "\n".join(lines) + "\n"
        lines.append("| Residue | locus | region | IMGT | conservation | class |")
        lines.append("|---------|-------|--------|------|--------------|-------|")
        for _, row in sub.head(top_n).iterrows():
            lines.append(
                f"| {row['residue']} "
                f"| {row.get('locus', '-')} "
                f"| {row.get('region', '-')} "
                f"| {row.get('imgt_position', '-')} "
                f"| {num(row.get('conservation_score', 0), 3)} "
                f"| {row.get('conservation_class', '-')} |"
            )
        if len(sub) > top_n:
            lines.append(f"\n_Showing {top_n} of {len(sub)}._")
        return "\n".join(lines) + "\n"

    @staticmethod
    def _render_coverage_note(df: pd.DataFrame) -> str:
        total = len(df)
        covered = int((df["covered"] == True).sum())  # noqa: E712
        uncovered = total - covered
        return (
            "## Coverage\n\n"
            f"- {covered}/{total} residues covered by the germline V reference.\n"
            f"- {uncovered} residues uncovered — chiefly the somatic **CDR3 "
            "junction**, which is not encoded in the germline V gene and so has "
            "no germline conservation score (only the germline-encoded CDR3 "
            "anchor is scored).\n"
        )

    @staticmethod
    def _render_interpretation() -> str:
        return (
            "## Interpretation\n\n"
            "- **FR (framework) high conservation**: structural scaffold; "
            "mutations risk destabilising the Ig fold → avoid.\n"
            "- **CDR1 / CDR2 lower conservation**: germline-variable contact "
            "loops; more substitution freedom.\n"
            "- **CDR3 (covered=False)**: somatic junction, not assessed by "
            "germline conservation — defer to repertoire / specificity-group "
            "analysis.\n"
            "- Cross-check with MD evidence: **conserved + contact hotspot** = "
            "load-bearing (high risk); **variable + hotspot** = affinity-tuning "
            "knob (high value).\n"
        )
