"""Quality view - trajectory quality metrics."""

from pathlib import Path

from ..base import Presenter, RenderContext
from ..loaders import read_json
from ..formatters import num
from .. import register_view


@register_view("quality")
class QualityPresenter(Presenter):
    """
    Quality view - trajectory quality assessment.

    Shows RMSD convergence, trajectory length, and quality warnings.
    """

    view_name = "quality"
    max_length = 1500
    # D-B1: convergence grade, energy drift, T/P stability — Complex-level
    # static scalars only.
    spatial_layer = "complex"
    sub_flavors_served = ("static",)

    def render(self, context: RenderContext) -> str:
        """Render quality view."""
        locator = context.locator
        case_id = locator.get_case_id()

        lines = [f"# Quality: {case_id}\n"]

        # Load quality data
        quality_root = locator.get_module_root("quality")
        if not quality_root:
            lines.append("_Quality analysis not found._")
            return "\n".join(lines)

        # Try multiple possible paths
        quality_path = locator.first_existing(
            quality_root / "preprocess_quality_report.json",
            quality_root / "quality_report.json",
        )

        if not quality_path:
            lines.append("_Quality report not available._")
            return "\n".join(lines)

        quality = read_json(quality_path)

        if not quality:
            lines.append("_Quality data not available._")
            return "\n".join(lines)

        # Render sections
        lines.append(self._render_convergence(quality))
        lines.append(self._render_trajectory_info(quality))
        lines.append(self._render_warnings(quality))

        # Sources
        lines.append(self._format_sources(
            quality_path.relative_to(locator.case_dir)
        ))

        output = "\n".join(lines)
        return self._cap_length(output)

    def _render_convergence(self, quality: dict) -> str:
        """Render RMSD convergence info."""
        lines = ["## RMSD convergence\n"]

        # Check convergence status
        is_converged = quality.get("is_converged", quality.get("rmsd_converged"))
        if is_converged is not None:
            icon = "✅" if is_converged else "⚠️"
            status = "Converged" if is_converged else "Not converged"
            lines.append(f"- **Status**: {icon} {status}")

        # RMSD statistics
        rmsd_mean = quality.get("rmsd_mean", quality.get("mean_rmsd"))
        rmsd_std = quality.get("rmsd_std", quality.get("std_rmsd"))

        if rmsd_mean is not None:
            lines.append(f"- **Mean RMSD**: {num(rmsd_mean, 2, 'nm')}")
        if rmsd_std is not None:
            lines.append(f"- **Std dev**: {num(rmsd_std, 2, 'nm')}")

        # Convergence time if available
        convergence_time = quality.get("convergence_time_ns")
        if convergence_time is not None:
            lines.append(f"- **Converged at**: {num(convergence_time, 1, 'ns')}")

        return "\n".join(lines) + "\n"

    def _render_trajectory_info(self, quality: dict) -> str:
        """Render trajectory information."""
        lines = ["## Trajectory\n"]

        n_frames = quality.get("n_frames", quality.get("total_frames"))
        if n_frames is not None:
            lines.append(f"- **Frames**: {int(n_frames)}")

        # Time info
        total_time = quality.get("total_time_ns")
        if total_time is not None:
            lines.append(f"- **Total time**: {num(total_time, 1, 'ns')}")

        dt = quality.get("dt_ps")
        if dt is not None:
            lines.append(f"- **Time step**: {num(dt, 1, 'ps')}")

        # Grade if available
        grade = quality.get("overall_grade", quality.get("grade"))
        if grade:
            lines.append(f"- **Quality grade**: {grade}")

        return "\n".join(lines) + "\n"

    def _render_warnings(self, quality: dict) -> str:
        """Render quality warnings."""
        lines = ["## Warnings\n"]

        warnings = quality.get("warnings", [])
        issues = quality.get("issues", [])

        all_warnings = warnings + issues

        if not all_warnings:
            lines.append("✅ No quality issues detected")
            return "\n".join(lines) + "\n"

        for warning in all_warnings:
            lines.append(f"- ⚠️ {warning}")

        return "\n".join(lines) + "\n"
