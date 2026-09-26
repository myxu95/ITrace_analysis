"""Compatibility facade for the legacy PBC plus RMSD quality pipeline."""

from __future__ import annotations

from typing import Dict

from md_analysis.analysis.trajectory import PBCProcessor

from .quality_assessment_pipeline import QualityAssessmentPipeline


class PBCRMSDPipeline(QualityAssessmentPipeline):
    """
    Thin compatibility facade over current preprocessing and quality modules.

    New workflow code should prefer ``PreprocessQualityPipeline`` for node-based
    orchestration. This class preserves the old helper methods used by existing
    tests and notebooks.
    """

    def __init__(
        self,
        gmx_executable: str = "gmx",
        max_com_drift: float = 1.0,
        max_rg_std_ratio: float = 0.15,
        max_frame_jump_rmsd: float = 0.5,
    ):
        self.pbc_processor = PBCProcessor(gmx_executable=gmx_executable)
        super().__init__(
            max_com_drift=max_com_drift,
            max_rg_std_ratio=max_rg_std_ratio,
            max_frame_jump_rmsd=max_frame_jump_rmsd,
        )

    def _determine_overall_grade(self, validation_results: Dict | None, rmsd_metrics: Dict) -> str:
        """Determine the overall grade, allowing skipped post-PBC validation."""
        rmsd_grade = self._assign_legacy_rmsd_grade(rmsd_metrics)
        if validation_results is None:
            return rmsd_grade

        pbc_grade = validation_results.get("overall_grade", "D")
        if pbc_grade == "D" or rmsd_grade == "D":
            return "D"
        grade_order = {"A": 1, "B": 2, "C": 3}
        worst = max(grade_order.get(pbc_grade, 4), grade_order.get(rmsd_grade, 4))
        return {1: "A", 2: "B", 3: "C"}.get(worst, "D")

    def _generate_markdown_report(self, results: Dict, trajectory_name: str | None = None) -> str:
        """Generate a quality report while accepting the legacy trajectory argument."""
        results = self._normalize_legacy_report_results(results)
        report = super()._generate_markdown_report(results)
        if trajectory_name:
            return report.replace(
                "# MD Trajectory Quality Assessment Report",
                f"# MD Trajectory Quality Assessment Report\n\n**Trajectory**: {trajectory_name}",
                1,
            )
        return report

    @staticmethod
    def _assign_legacy_rmsd_grade(rmsd_metrics: Dict) -> str:
        """Assign the historical RMSD grade used by the old PBC-RMSD tests."""
        if not rmsd_metrics.get("is_converged", False):
            return "D"
        mean = float(rmsd_metrics.get("mean_rmsd", 0.0))
        std = float(rmsd_metrics.get("std_rmsd", 0.0))
        fraction = float(rmsd_metrics.get("convergence_fraction", 1.0))
        if mean <= 0.30 and std <= 0.05 and fraction <= 0.20:
            return "A"
        if mean <= 0.50 and std <= 0.10 and fraction <= 0.30:
            return "B"
        if mean <= 0.75 and std <= 0.15 and fraction <= 0.50:
            return "C"
        return "D"

    @staticmethod
    def _normalize_legacy_report_results(results: Dict) -> Dict:
        """Fill optional fields that old tests omitted but current reports expect."""
        normalized = dict(results)
        metrics = dict(normalized.get("rmsd_metrics") or {})
        mean = float(metrics.get("mean_rmsd", 0.0))
        std = float(metrics.get("std_rmsd", 0.0))
        metrics.setdefault("min_rmsd", max(0.0, mean - std))
        metrics.setdefault("max_rmsd", mean + std)

        block = metrics.get("block_analysis")
        if isinstance(block, dict):
            block = dict(block)
            block.setdefault("block_stability", 0.0)
            for info in block.get("block_info", []):
                info.setdefault("std_rmsd", 0.0)
            metrics["block_analysis"] = block

        normalized["rmsd_metrics"] = metrics
        return normalized
