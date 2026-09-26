"""
Quality Assessment Pipeline

This module provides a focused pipeline for comprehensive trajectory quality assessment
without PBC processing. It integrates Post-PBC validation and RMSD convergence analysis.

Classes:
    QualityAssessmentPipeline: Quality assessment pipeline for processed trajectories

Author: ImmunoScope Development Team
Date: 2026-03-15
"""

import logging
from pathlib import Path
from typing import Dict, Optional
import json

from ..analysis.quality.post_pbc_validator import PostPBCValidator
from ..analysis.trajectory.rmsd_refactored import RMSDCalculator, RMSDInput
from ..analysis.trajectory.rmsd_convergence import RMSDConvergenceAnalyzer


class QualityAssessmentPipeline:
    """
    Quality Assessment Pipeline

    This pipeline provides comprehensive quality assessment for processed trajectories.
    It combines Post-PBC validation and RMSD convergence analysis.

    Use this pipeline when:
    - Trajectory is already PBC-corrected
    - You need quality assessment only (no PBC processing)
    - You want detailed convergence and stability analysis

    Parameters
    ----------
    max_com_drift : float, default=1.0
        Maximum allowed center-of-mass drift in nm
    max_rg_std_ratio : float, default=0.15
        Maximum allowed Rg std/mean ratio
    max_frame_jump_rmsd : float, default=0.5
        Maximum allowed frame-to-frame RMSD in nm

    Examples
    --------
    >>> pipeline = QualityAssessmentPipeline()
    >>> results = pipeline.run_comprehensive_assessment(
    ...     trajectory="md_pbc.xtc",
    ...     topology="md.tpr",
    ...     rmsd_selection="backbone"
    ... )
    >>> print(f"Quality Grade: {results['overall_grade']}")
    """

    def __init__(self,
                 max_com_drift: float = 1.0,
                 max_rg_std_ratio: float = 0.15,
                 max_frame_jump_rmsd: float = 0.5):
        """Initialize quality assessment pipeline"""
        self.logger = logging.getLogger(__name__)

        # Initialize functional modules
        self.post_pbc_validator = PostPBCValidator(
            max_com_drift=max_com_drift,
            max_rg_std_ratio=max_rg_std_ratio,
            max_frame_jump_rmsd=max_frame_jump_rmsd
        )
        self.convergence_analyzer = RMSDConvergenceAnalyzer()

    def run_comprehensive_assessment(self,
                                     trajectory: str,
                                     topology: str,
                                     rmsd_selection: str = "backbone",
                                     validation_selection: str = "protein",
                                     validation_stride: int = 1,
                                     output_dir: Optional[str] = None,
                                     enable_post_pbc_validation: bool = True) -> Dict:
        """
        Run comprehensive quality assessment

        Workflow:
        1. Post-PBC validation (trajectory integrity, PBC quality, stability)
        2. RMSD calculation
        3. Convergence analysis
        4. Overall quality grading

        Parameters
        ----------
        trajectory : str
            Path to PBC-corrected trajectory file
        topology : str
            Path to topology file
        rmsd_selection : str, default="backbone"
            MDAnalysis selection string for RMSD calculation
        validation_selection : str, default="protein"
            MDAnalysis selection string for validation
        validation_stride : int, default=1
            Stride for validation sampling
        output_dir : str, optional
            Output directory for saving results (if None, no files saved)

        Returns
        -------
        dict
            Assessment results with keys:
            - post_pbc_validation: validation results
            - rmsd_metrics: convergence metrics
            - overall_grade: quality grade ('A', 'B', 'C', 'D')
            - is_qualified: whether trajectory passes quality checks
        """
        self.logger.info(f"Starting comprehensive quality assessment: {trajectory}")

        results = {
            'post_pbc_validation': None,
            'rmsd_metrics': None,
            'overall_grade': None,
            'is_qualified': False
        }

        try:
            # Step 1: Post-PBC validation
            if enable_post_pbc_validation:
                self.logger.info("Step 1/3: Post-PBC validation")
                validation_results = self.post_pbc_validator.comprehensive_validation(
                    trajectory=trajectory,
                    topology=topology,
                    selection=validation_selection,
                    stride=validation_stride
                )
                results['post_pbc_validation'] = validation_results

                if validation_results['overall_status'] == 'fail':
                    self.logger.warning(
                        f"Post-PBC validation failed: {validation_results['all_issues']}"
                    )
            else:
                self.logger.info("Step 1/3: Post-PBC validation skipped")
                validation_results = {
                    'overall_status': 'skip',
                    'overall_grade': 'A',
                    'all_issues': [],
                    'skipped': True,
                }
                results['post_pbc_validation'] = validation_results

            # Step 2: RMSD calculation
            self.logger.info("Step 2/3: RMSD calculation")
            rmsd_calc = RMSDCalculator()

            rmsd_output = None
            if output_dir:
                output_path = Path(output_dir)
                output_path.mkdir(parents=True, exist_ok=True)
                rmsd_output = str(output_path / "rmsd.csv")

            rmsd_result = rmsd_calc.calculate(
                RMSDInput(
                    topology=topology,
                    trajectory=trajectory,
                    selection=rmsd_selection,
                    output_file=rmsd_output,
                )
            )
            if not rmsd_result.success or rmsd_result.times is None or rmsd_result.rmsd_values is None:
                raise RuntimeError(rmsd_result.error_message or "RMSD calculation failed")

            times = rmsd_result.times
            rmsd_values = rmsd_result.rmsd_values

            # Step 3: Convergence analysis
            self.logger.info("Step 3/3: Convergence analysis")
            rmsd_metrics = self.convergence_analyzer.calculate_convergence_metrics(
                times=times,
                rmsd_values=rmsd_values
            )
            results['rmsd_metrics'] = rmsd_metrics

            # Overall quality grading
            overall_grade = self._determine_overall_grade(
                validation_results, rmsd_metrics
            )
            results['overall_grade'] = overall_grade
            results['is_qualified'] = overall_grade in ['A', 'B', 'C']

            self.logger.info(
                f"Assessment complete: Grade={overall_grade}, "
                f"Qualified={results['is_qualified']}"
            )

        except Exception as e:
            self.logger.error(f"Quality assessment failed: {str(e)}")
            results['error'] = str(e)
            raise

        return results

    def _determine_overall_grade(self,
                                validation_results: Dict,
                                rmsd_metrics: Dict) -> str:
        """
        Determine overall quality grade

        Uses the worse grade from Post-PBC validation and RMSD analysis.
        """
        pbc_grade = validation_results.get('overall_grade', 'D')
        rmsd_grade = self.convergence_analyzer.assign_quality_grade(rmsd_metrics)

        # Grade D has highest priority
        if pbc_grade == 'D' or rmsd_grade == 'D':
            return 'D'

        # Take the worse grade
        grade_order = {'A': 1, 'B': 2, 'C': 3}
        pbc_score = grade_order.get(pbc_grade, 4)
        rmsd_score = grade_order.get(rmsd_grade, 4)

        worst_score = max(pbc_score, rmsd_score)
        grade_map = {1: 'A', 2: 'B', 3: 'C'}

        return grade_map.get(worst_score, 'D')

    def generate_quality_report(self,
                                results: Dict,
                                output_file: str,
                                format: str = "markdown") -> str:
        """
        Generate quality assessment report

        Parameters
        ----------
        results : dict
            Assessment results from run_comprehensive_assessment()
        output_file : str
            Output file path
        format : str, default="markdown"
            Report format: "markdown" or "json"

        Returns
        -------
        str
            Path to generated report file
        """
        if format == "markdown":
            content = self._generate_markdown_report(results)
        elif format == "json":
            content = json.dumps(self._serialize_results(results), indent=2)
        else:
            raise ValueError(f"Unknown format: {format}")

        with open(output_file, 'w') as f:
            f.write(content)

        self.logger.info(f"Generated {format} report: {output_file}")
        return output_file

    def _generate_markdown_report(self, results: Dict) -> str:
        """Generate Markdown quality report"""
        lines = [
            "# MD Trajectory Quality Assessment Report",
            "",
            f"**Overall Grade**: {results['overall_grade']}",
            f"**Status**: {'✅ Qualified' if results['is_qualified'] else '❌ Not Qualified'}",
            "",
            "## Post-PBC Validation Results",
            ""
        ]

        # Post-PBC validation section
        if results.get('post_pbc_validation'):
            pbc_val = results['post_pbc_validation']

            lines.extend([
                f"**Status**: {pbc_val['overall_status'].upper()}",
                f"**Grade**: {pbc_val['overall_grade']}",
                ""
            ])

            # Trajectory integrity
            if pbc_val.get('integrity'):
                integrity = pbc_val['integrity']
                lines.extend([
                    "### Trajectory Integrity",
                    "",
                    f"- Frames: {integrity.get('n_frames', 0)}",
                    f"- File Size: {integrity.get('file_size_mb', 0):.2f} MB",
                    f"- Status: {integrity.get('status', 'unknown').upper()}",
                    ""
                ])

            # PBC quality
            if pbc_val.get('pbc_quality'):
                pbc_q = pbc_val['pbc_quality']
                coord_range = pbc_q.get('coordinate_range_nm', (0, 0))
                lines.extend([
                    "### PBC Correction Quality",
                    "",
                    f"- COM Drift: {pbc_q.get('com_drift_nm', 0):.3f} nm "
                    f"{'✅' if pbc_q.get('status') == 'pass' else '❌'}",
                    f"- Coordinate Range: [{coord_range[0]:.2f}, {coord_range[1]:.2f}] nm",
                    ""
                ])

            # Structural stability
            if pbc_val.get('stability'):
                stab = pbc_val['stability']
                lines.extend([
                    "### Structural Stability",
                    "",
                    f"- Rg Mean: {stab.get('rg_mean_nm', 0):.3f} nm",
                    f"- Rg Std: {stab.get('rg_std_nm', 0):.3f} nm",
                    f"- Rg Std/Mean: {stab.get('rg_std_ratio', 0):.3f} "
                    f"{'✅' if stab.get('rg_std_ratio', 1) < 0.15 else '❌'}",
                    f"- Max Frame Jump: {stab.get('max_frame_jump_rmsd_nm', 0):.3f} nm",
                    ""
                ])

            # Issues
            if pbc_val.get('all_issues'):
                lines.extend([
                    "### Detected Issues",
                    ""
                ])
                for issue in pbc_val['all_issues']:
                    lines.append(f"- {issue}")
                lines.append("")

        # RMSD analysis section
        if results.get('rmsd_metrics'):
            metrics = results['rmsd_metrics']
            rmsd_grade = self.convergence_analyzer.assign_quality_grade(metrics)

            lines.extend([
                "## RMSD Quality Assessment",
                "",
                f"**Grade**: {rmsd_grade}",
                "",
                "### Overall Statistics",
                "",
                f"- Mean RMSD: {metrics['mean_rmsd']:.3f} nm",
                f"- Std RMSD: {metrics['std_rmsd']:.3f} nm",
                f"- Min RMSD: {metrics['min_rmsd']:.3f} nm",
                f"- Max RMSD: {metrics['max_rmsd']:.3f} nm",
                "",
                "### Convergence",
                "",
                f"- Converged: {'Yes ✅' if metrics['is_converged'] else 'No ❌'}",
                f"- Convergence Time: {metrics['convergence_time']:.1f}",
                f"- Convergence Fraction: {metrics['convergence_fraction']:.1%}",
                ""
            ])

            # Block analysis
            if metrics.get('block_analysis'):
                block = metrics['block_analysis']
                lines.extend([
                    f"### Block Analysis ({block['n_blocks']} blocks)",
                    "",
                    f"- Block Stability (CV): {block['block_stability']:.3f}",
                    "",
                    "| Block | Time Range | Mean (nm) | Std (nm) |",
                    "|-------|------------|-----------|----------|"
                ])

                for info in block['block_info']:
                    t_range = f"{info['time_range'][0]:.0f}-{info['time_range'][1]:.0f}"
                    lines.append(
                        f"| {info['block_id']} | {t_range} | "
                        f"{info['mean_rmsd']:.3f} | {info['std_rmsd']:.3f} |"
                    )

                lines.append("")

            # Moving average
            if metrics.get('moving_avg_analysis'):
                ma = metrics['moving_avg_analysis']
                lines.extend([
                    "### Moving Average Analysis",
                    "",
                    f"- Window Size: {ma['window_size']} frames",
                    f"- Trend Slope: {ma['trend_slope']:.6f} nm/unit",
                    f"- Has Drift: {'Yes' if ma['has_drift'] else 'No'}",
                    ""
                ])

        lines.extend([
            "---",
            f"*Report generated by ImmunoScope Quality Assessment Pipeline*"
        ])

        return "\n".join(lines)

    def _serialize_results(self, results: Dict) -> Dict:
        """Convert results to JSON-serializable format"""
        import numpy as np

        def serialize_value(value):
            if isinstance(value, np.ndarray):
                return value.tolist()
            elif isinstance(value, dict):
                return {k: serialize_value(v) for k, v in value.items()}
            elif isinstance(value, list):
                return [serialize_value(v) for v in value]
            elif isinstance(value, (str, int, float, bool, type(None))):
                return value
            else:
                return str(value)

        return {k: serialize_value(v) for k, v in results.items()}
