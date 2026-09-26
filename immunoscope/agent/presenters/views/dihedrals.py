"""Dihedrals view — Ramachandran (φ/ψ) backbone conformation per residue.

Reports the time-averaged φ/ψ statistics, region populations
(α-helix / β-sheet / left-handed / other), and the residues with the
greatest backbone flexibility (large angular spread). Useful for
mutation design: residues with rigid backbone tolerate fewer
substitutions; residues sampling multiple Ramachandran basins are
candidate flexibility hotspots.
"""

from __future__ import annotations

from ..base import Presenter, RenderContext
from ..loaders import read_json
from .. import register_view


@register_view("dihedrals")
class DihedralsPresenter(Presenter):
    view_name = "dihedrals"
    max_length = 3500
    # D-B1: per-residue φ/ψ statistics, region populations, χ entropies —
    # Residue-level dynamic (distribution / ensemble shape over rotamer
    # bins).
    spatial_layer = "residue"
    sub_flavors_served = ("dynamic",)

    def render(self, context: RenderContext) -> str:
        locator = context.locator
        case_id = locator.get_case_id()
        lines = [f"# Ramachandran (φ/ψ): {case_id}", ""]

        dih_root = locator.get_module_root("dihedrals")
        if not dih_root:
            lines.append("_Dihedral analysis not run for this case._")
            return "\n".join(lines)

        summary_path = locator.first_existing(
            dih_root / "dihedrals_summary.json",
        )
        if not summary_path:
            lines.append("_Dihedrals summary not generated yet._")
            return "\n".join(lines)

        summary = read_json(summary_path) or {}
        n_frames = summary.get("n_frames", 0)
        n_res = summary.get("n_residues_analyzed", 0)
        stride = summary.get("stride", 1)
        lines.append(f"**Frames analyzed**: {n_frames} (stride {stride}) — **Residues**: {n_res}")
        mean_spread = summary.get("mean_angular_spread_deg")
        if mean_spread is not None:
            lines.append(f"**Average backbone angular spread**: {mean_spread:.2f}°")
        lines.append("")

        # Region populations
        regions = summary.get("region_summary", {})
        if regions:
            lines.append("## Ramachandran region populations")
            lines.append("")
            lines.append("| Region | Mean fraction | Residues with this as dominant |")
            lines.append("|--------|---------------|--------------------------------|")
            for r in ("alpha", "beta", "left_alpha", "other"):
                stats = regions.get(r) or {}
                mf = stats.get("mean_fraction", 0)
                nr = stats.get("n_residues_dominant", 0)
                label = {
                    "alpha": "α-helix",
                    "beta": "β-sheet",
                    "left_alpha": "Left-handed α",
                    "other": "Other (loops, turns)",
                }[r]
                lines.append(f"| {label} | {mf*100:.1f}% | {nr} |")
            lines.append("")

        # Most flexible residues — candidates for mutation
        flex = summary.get("most_flexible_residues") or []
        if flex:
            lines.append("## Top 10 most backbone-flexible residues")
            lines.append("_(High angular spread = backbone samples multiple basins; often tolerates mutations.)_")
            lines.append("")
            lines.append("| Chain | Resid | Residue | Angular spread (°) | Dominant region |")
            lines.append("|-------|-------|---------|--------------------|------------------|")
            for row in flex[:10]:
                lines.append(
                    f"| {row.get('chain_id', '-')} "
                    f"| {row.get('resid', '-')} "
                    f"| {row.get('resname', '-')} "
                    f"| {row.get('angular_spread_deg', 0):.1f} "
                    f"| {row.get('dominant_region', '-')} |"
                )
            lines.append("")

        # Most rigid residues — risky for mutation
        rigid = summary.get("most_rigid_residues") or []
        if rigid:
            lines.append("## Top 10 most backbone-rigid residues")
            lines.append("_(Low angular spread = backbone locked in one conformation; "
                         "mutations may disrupt local structure.)_")
            lines.append("")
            lines.append("| Chain | Resid | Residue | Angular spread (°) | Dominant region |")
            lines.append("|-------|-------|---------|--------------------|------------------|")
            for row in rigid[:10]:
                lines.append(
                    f"| {row.get('chain_id', '-')} "
                    f"| {row.get('resid', '-')} "
                    f"| {row.get('resname', '-')} "
                    f"| {row.get('angular_spread_deg', 0):.1f} "
                    f"| {row.get('dominant_region', '-')} |"
                )
            lines.append("")

        lines.append("## Sources")
        lines.append(f"- `{summary_path}`")
        return "\n".join(lines)
