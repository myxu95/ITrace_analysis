"""Interface view - BSA statistics and composition."""

from pathlib import Path
import pandas as pd

from ..base import Presenter, RenderContext
from ..loaders import read_json, read_csv
from ..formatters import num, pct
from .. import register_view


@register_view("interface")
class InterfacePresenter(Presenter):
    """
    Interface view - BSA statistics and composition.

    Shows buried surface area statistics and interface composition
    (which regions contribute to the interface).
    """

    view_name = "interface"
    max_length = 2000
    # D-B1: BSA mean/std + interface-ratio — Interface-level static scalars.
    # H-bond/salt-bridge/etc. interface totals are on-the-fly aggregates of
    # the Pair layer; they are surfaced here for convenience but not
    # persisted as separate metrics (see Methods Table A1.2).
    spatial_layer = "interface"
    sub_flavors_served = ("static",)

    def render(self, context: RenderContext) -> str:
        """Render interface view."""
        locator = context.locator
        case_id = locator.get_case_id()

        lines = [f"# Interface: {case_id}\n"]

        # Load BSA data
        bsa_root = locator.get_module_root("bsa")
        if not bsa_root:
            lines.append("_BSA analysis not found. Run BSA analysis first._")
            return "\n".join(lines)

        summary_path = bsa_root / "analysis/interface/interface_summary.json"
        summary = read_json(summary_path)

        if not summary:
            lines.append("_BSA data not available._")
            return "\n".join(lines)

        # Render sections
        lines.append(self._render_bsa_stats(summary))
        lines.append(self._render_composition(summary))

        # Try to add contact overview if available
        contact_overview = self._render_contact_overview(locator)
        if contact_overview:
            lines.append(contact_overview)

        # Sources
        lines.append(self._format_sources(
            summary_path.relative_to(locator.case_dir)
        ))

        output = "\n".join(lines)
        return self._cap_length(output)

    def _render_bsa_stats(self, summary: dict) -> str:
        """Render BSA statistics."""
        lines = ["## BSA statistics\n"]

        # Try nested structure (buried_surface_area.mean/std/min/max) first
        bsa_data = summary.get("buried_surface_area", {})
        if isinstance(bsa_data, dict) and bsa_data:
            mean_bsa = bsa_data.get("mean", 0)
            std_bsa = bsa_data.get("std", 0)
            min_bsa = bsa_data.get("min", 0)
            max_bsa = bsa_data.get("max", 0)
        else:
            # Flat structure fallback
            mean_bsa = summary.get("mean_bsa", 0)
            std_bsa = summary.get("std_bsa", 0)
            min_bsa = summary.get("min_bsa", 0)
            max_bsa = summary.get("max_bsa", 0)

        lines.append(f"- **Mean BSA**: {num(mean_bsa, 1, 'Å²')}")
        lines.append(f"- **Std dev**: {num(std_bsa, 1, 'Å²')}")
        lines.append(f"- **Range**: {num(min_bsa, 1, 'Å²')} – {num(max_bsa, 1, 'Å²')}")

        # Convergence info if available
        if "is_converged" in summary:
            converged = summary["is_converged"]
            icon = "✅" if converged else "⚠️"
            lines.append(f"- **Convergence**: {icon} {'Yes' if converged else 'No'}")

        return "\n".join(lines) + "\n"

    def _render_composition(self, summary: dict) -> str:
        """Render interface composition."""
        lines = ["## Interface composition\n"]

        # Try to extract component contributions
        components = summary.get("components", {})
        if not components:
            # Try alternative structure
            tcr_bsa = summary.get("tcr_bsa")
            phla_bsa = summary.get("phla_bsa")

            if tcr_bsa is not None and phla_bsa is not None:
                total = tcr_bsa + phla_bsa
                if total > 0:
                    lines.append(f"- **TCR contribution**: {num(tcr_bsa, 1, 'Å²')} ({pct(tcr_bsa/total)})")
                    lines.append(f"- **pMHC contribution**: {num(phla_bsa, 1, 'Å²')} ({pct(phla_bsa/total)})")
                return "\n".join(lines) + "\n"

        if components:
            for component, bsa in components.items():
                lines.append(f"- **{component}**: {num(bsa, 1, 'Å²')}")
            return "\n".join(lines) + "\n"

        lines.append("_Composition breakdown not available_")
        return "\n".join(lines) + "\n"

    def _render_contact_overview(self, locator) -> str:
        """Render contact overview if available."""
        # Try to load interaction_overview.csv from landscape or contact module
        landscape_root = locator.get_module_root("landscape")
        contact_root = locator.get_module_root("contact")

        overview_path = None
        if landscape_root:
            candidate = landscape_root / "analysis/contacts/interaction_overview.csv"
            if candidate.exists():
                overview_path = candidate

        if not overview_path and contact_root:
            candidate = contact_root / "analysis/contacts/interaction_overview.csv"
            if candidate.exists():
                overview_path = candidate

        if not overview_path:
            return ""

        df = read_csv(overview_path)
        if df.empty:
            return ""

        lines = ["## Contact overview\n"]

        # Show top interaction families
        if "interaction_family" in df.columns and "count" in df.columns:
            family_counts = df.groupby("interaction_family")["count"].sum().sort_values(ascending=False)

            lines.append("**Interaction families**:")
            for family, count in family_counts.head(5).items():
                lines.append(f"- {family}: {int(count)} contacts")

        return "\n".join(lines) + "\n"
