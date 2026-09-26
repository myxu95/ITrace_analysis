"""Docking angles view - TCR-pMHC binding geometry (crossing/incident/tilt)."""

from __future__ import annotations

from ..base import Presenter, RenderContext
from ..loaders import read_json
from .. import register_view


@register_view("angles")
class DockingAnglesPresenter(Presenter):
    """Docking angles view.

    Reports the time-averaged TCR-pMHC docking geometry: crossing angle
    (TCR vs peptide main axis), incident angle (TCR vs groove plane
    normal), and tilt. These define the binding orientation and are
    closely correlated with TCR specificity and contact geometry.
    """

    view_name = "angles"
    max_length = 1800
    # D-B1: crossing / incident / tilt angles + CDR3 loop geometric span —
    # Complex-level static scalars (means + std reported, but reported value
    # is the scalar, not the distribution shape).
    spatial_layer = "complex"
    sub_flavors_served = ("static",)

    def render(self, context: RenderContext) -> str:
        locator = context.locator
        case_id = locator.get_case_id()
        lines = [f"# Docking angles: {case_id}", ""]

        angles_root = locator.get_module_root("angles")
        if not angles_root:
            lines.append("_Docking angle analysis not run for this case._")
            lines.append("")
            lines.append("To compute: include `angles` in the run modules.")
            return "\n".join(lines)

        summary_path = locator.first_existing(
            angles_root / "docking_angles_summary.json",
        )
        if not summary_path:
            lines.append("_Angle summary not generated yet._")
            return "\n".join(lines)

        summary = read_json(summary_path) or {}
        n_frames = summary.get("n_frames", 0)
        stride = summary.get("stride", 1)
        lines.append(f"**Frames analyzed**: {n_frames} (stride {stride})")
        lines.append("")

        # Render each angle type as a row
        lines.append("| Angle | Mean (°) | Std (°) | Range (°) | Interpretation |")
        lines.append("|-------|----------|---------|-----------|----------------|")

        for key, label in (
            ("crossing_angle", "Crossing"),
            ("incident_angle", "Incident"),
            ("tilt_angle", "Tilt"),
        ):
            stats = summary.get(key)
            if not stats:
                continue
            mean = stats.get("mean")
            std = stats.get("std")
            rmin = stats.get("min")
            rmax = stats.get("max")
            interp = self._interpret(key, mean, std)
            lines.append(
                f"| {label} | {mean:.2f} | {std:.2f} | {rmin:.1f}–{rmax:.1f} | {interp} |"
            )

        lines.append("")
        # Stability summary — angle std < 5° usually means very stable binding
        crossings_std = (summary.get("crossing_angle") or {}).get("std", 0)
        if crossings_std:
            if crossings_std < 5:
                lines.append("_Binding orientation is **highly stable** "
                             "(crossing std < 5°). Strong specificity signal._")
            elif crossings_std < 10:
                lines.append("_Binding orientation is **moderately stable** "
                             "(crossing std 5–10°)._")
            else:
                lines.append("_Binding orientation is **flexible** "
                             "(crossing std > 10°). May indicate weak or "
                             "transient binding._")

        lines.append("")
        lines.append("## Sources")
        lines.append(f"- `{summary_path}`")

        return "\n".join(lines)

    def _interpret(self, key: str, mean: float, std: float) -> str:
        """Quick interpretation hint for each angle."""
        if mean is None:
            return ""
        if key == "crossing_angle":
            # Typical TCR-pMHC crossing angle: 20°-70°, canonical ~40°
            if mean < 20:
                return "Low crossing (parallel-like)"
            if mean < 50:
                return "Canonical TCR orientation"
            if mean < 70:
                return "Reverse/diagonal binding"
            return "Atypical orientation"
        if key == "incident_angle":
            # Tilt of TCR axis vs groove plane normal
            if mean < 30:
                return "Steep dive into groove"
            if mean < 60:
                return "Typical incident"
            return "Shallow approach"
        if key == "tilt_angle":
            if std < 5:
                return "Stable tilt"
            return "Variable tilt"
        return ""
