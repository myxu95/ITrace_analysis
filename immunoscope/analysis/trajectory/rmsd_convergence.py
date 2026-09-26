"""
RMSD Convergence Analyzer

This module provides tools for analyzing RMSD convergence and stability in MD trajectories.
It includes block analysis, moving average analysis, and quality grading.

Classes:
    RMSDConvergenceAnalyzer: Main analyzer for RMSD convergence assessment

Author: ImmunoScope Development Team
Date: 2026-03-15
"""

import logging
from typing import Dict, Tuple, Optional
import numpy as np


class RMSDConvergenceAnalyzer:
    """
    RMSD Convergence and Stability Analyzer

    Analyzes RMSD time series to determine:
    1. Convergence status and convergence time
    2. Stability metrics (mean, std, block-to-block variation)
    3. Quality grade assignment (A/B/C/D)

    Parameters
    ----------
    n_blocks : int, default=5
        Number of blocks for block analysis
    moving_avg_window : int, default=100
        Window size for moving average analysis
    convergence_threshold : float, default=0.05
        Threshold for convergence detection (nm)

    Quality Grading Criteria (Stability-Based)
    ------------------------------------------
    Uses Coefficient of Variation (CV) = Std / Mean to assess stability:

    Grade A (Excellent): CV < 0.10, Converged, Conv_frac < 0.2
    Grade B (Good): CV < 0.15 (convergence not required)
    Grade C (Acceptable): CV < 0.20
    Grade D (Poor): CV >= 0.20

    Note: Absolute RMSD values are NOT used for grading. A trajectory with
    high RMSD but low CV (e.g., mean=4.0nm, std=0.6nm, CV=0.15) indicates
    stable dynamics and receives Grade B.

    Examples
    --------
    >>> analyzer = RMSDConvergenceAnalyzer()
    >>> metrics = analyzer.calculate_convergence_metrics(times, rmsd_values)
    >>> grade = analyzer.assign_quality_grade(metrics)
    >>> print(f"Quality Grade: {grade}")
    """

    def __init__(self,
                 n_blocks: int = 5,
                 moving_avg_window: int = 100,
                 convergence_threshold: float = 0.05):
        """Initialize RMSDConvergenceAnalyzer"""
        self.n_blocks = n_blocks
        self.moving_avg_window = moving_avg_window
        self.convergence_threshold = convergence_threshold

        self.logger = logging.getLogger(__name__)

    def block_analysis(self,
                      times: np.ndarray,
                      rmsd_values: np.ndarray,
                      n_blocks: Optional[int] = None) -> Dict:
        """
        Perform block analysis to assess convergence

        Divides trajectory into blocks and calculates statistics for each block.
        Convergence is indicated by stable block-to-block statistics.

        Parameters
        ----------
        times : np.ndarray
            Time points (ps or ns)
        rmsd_values : np.ndarray
            RMSD values (nm)
        n_blocks : int, optional
            Number of blocks (uses self.n_blocks if None)

        Returns
        -------
        dict
            Block analysis results with keys:
            - n_blocks: number of blocks
            - block_info: list of dicts with block statistics
            - block_means: array of block mean RMSD values
            - block_stds: array of block standard deviations
            - block_stability: coefficient of variation of block means
        """
        if n_blocks is None:
            n_blocks = self.n_blocks

        n_frames = len(rmsd_values)
        block_size = n_frames // n_blocks

        if block_size < 10:
            self.logger.warning(
                f"Block size {block_size} too small, reducing n_blocks"
            )
            n_blocks = max(2, n_frames // 10)
            block_size = n_frames // n_blocks

        block_info = []
        block_means = []
        block_stds = []

        for i in range(n_blocks):
            start_idx = i * block_size
            if i == n_blocks - 1:
                end_idx = n_frames  # Include remaining frames in last block
            else:
                end_idx = (i + 1) * block_size

            block_times = times[start_idx:end_idx]
            block_rmsd = rmsd_values[start_idx:end_idx]

            block_mean = np.mean(block_rmsd)
            block_std = np.std(block_rmsd)

            block_info.append({
                'block_id': i + 1,
                'time_range': (block_times[0], block_times[-1]),
                'n_frames': len(block_rmsd),
                'mean_rmsd': block_mean,
                'std_rmsd': block_std
            })

            block_means.append(block_mean)
            block_stds.append(block_std)

        block_means = np.array(block_means)
        block_stds = np.array(block_stds)

        # Calculate block-to-block stability (coefficient of variation)
        mean_of_means = np.mean(block_means)
        std_of_means = np.std(block_means)
        block_stability = std_of_means / mean_of_means if mean_of_means > 0 else 0.0

        results = {
            'n_blocks': n_blocks,
            'block_info': block_info,
            'block_means': block_means,
            'block_stds': block_stds,
            'block_stability': block_stability
        }

        return results

    def moving_average_analysis(self,
                                times: np.ndarray,
                                rmsd_values: np.ndarray,
                                window_size: Optional[int] = None) -> Dict:
        """
        Calculate moving average to analyze RMSD trends

        Parameters
        ----------
        times : np.ndarray
            Time points (ps or ns)
        rmsd_values : np.ndarray
            RMSD values (nm)
        window_size : int, optional
            Window size for moving average (uses self.moving_avg_window if None)

        Returns
        -------
        dict
            Moving average results with keys:
            - window_size: size of moving average window
            - moving_avg: moving average values
            - trend_slope: linear trend slope (nm per time unit)
            - has_drift: whether significant drift is detected
        """
        if window_size is None:
            window_size = self.moving_avg_window

        n_frames = len(rmsd_values)
        if window_size > n_frames:
            window_size = max(10, n_frames // 10)

        # Calculate moving average using convolution
        kernel = np.ones(window_size) / window_size
        moving_avg = np.convolve(rmsd_values, kernel, mode='valid')

        # Fit linear trend to moving average
        valid_times = times[window_size - 1:]
        if len(moving_avg) > 10:
            trend_coeffs = np.polyfit(valid_times, moving_avg, 1)
            trend_slope = trend_coeffs[0]
        else:
            trend_slope = 0.0

        # Detect significant drift (slope > threshold)
        drift_threshold = self.convergence_threshold / (times[-1] - times[0])
        has_drift = abs(trend_slope) > drift_threshold

        results = {
            'window_size': window_size,
            'moving_avg': moving_avg,
            'moving_avg_times': valid_times,
            'trend_slope': trend_slope,
            'has_drift': has_drift
        }

        return results

    def calculate_convergence_metrics(self,
                                     times: np.ndarray,
                                     rmsd_values: np.ndarray) -> Dict:
        """
        Calculate comprehensive convergence metrics

        Parameters
        ----------
        times : np.ndarray
            Time points (ps or ns)
        rmsd_values : np.ndarray
            RMSD values (nm)

        Returns
        -------
        dict
            Convergence metrics with keys:
            - mean_rmsd: overall mean RMSD (nm)
            - std_rmsd: overall standard deviation (nm)
            - cv_rmsd: coefficient of variation (std/mean, dimensionless)
            - min_rmsd: minimum RMSD value (nm)
            - max_rmsd: maximum RMSD value (nm)
            - is_converged: whether trajectory is converged
            - convergence_time: estimated convergence time (same units as input times)
            - convergence_fraction: fraction of trajectory before convergence
            - block_analysis: results from block_analysis()
            - moving_avg_analysis: results from moving_average_analysis()
        """
        self.logger.info("Calculating RMSD convergence metrics")

        # Basic statistics
        mean_rmsd = np.mean(rmsd_values)
        std_rmsd = np.std(rmsd_values)
        min_rmsd = np.min(rmsd_values)
        max_rmsd = np.max(rmsd_values)

        # Calculate coefficient of variation (stability metric)
        cv_rmsd = std_rmsd / mean_rmsd if mean_rmsd > 0 else 0.0

        # Run block analysis
        block_results = self.block_analysis(times, rmsd_values)

        # Run moving average analysis
        moving_avg_results = self.moving_average_analysis(times, rmsd_values)

        # Estimate convergence time
        # Method: Find when moving average stabilizes (variance < threshold)
        convergence_time, convergence_fraction, is_converged = self._estimate_convergence_time(
            times, rmsd_values, moving_avg_results
        )

        metrics = {
            'mean_rmsd': mean_rmsd,
            'std_rmsd': std_rmsd,
            'cv_rmsd': cv_rmsd,
            'min_rmsd': min_rmsd,
            'max_rmsd': max_rmsd,
            'is_converged': is_converged,
            'convergence_time': convergence_time,
            'convergence_fraction': convergence_fraction,
            'block_analysis': block_results,
            'moving_avg_analysis': moving_avg_results
        }

        self.logger.info(
            f"Metrics: Mean={mean_rmsd:.3f} nm, Std={std_rmsd:.3f} nm, CV={cv_rmsd:.3f}, "
            f"Converged={is_converged}, Conv_frac={convergence_fraction:.2%}"
        )

        return metrics

    def _estimate_convergence_time(self,
                                   times: np.ndarray,
                                   rmsd_values: np.ndarray,
                                   moving_avg_results: Dict) -> Tuple[float, float, bool]:
        """
        Estimate convergence time from moving average

        Convergence is detected when the rolling standard deviation
        of the moving average falls below threshold.

        Returns
        -------
        tuple
            (convergence_time, convergence_fraction, is_converged)
        """
        moving_avg = moving_avg_results['moving_avg']
        moving_avg_times = moving_avg_results['moving_avg_times']

        if len(moving_avg) < 20:
            # Too few points to determine convergence
            return times[-1], 1.0, False

        # Calculate rolling std of moving average
        window = max(10, len(moving_avg) // 10)
        rolling_std = np.array([
            np.std(moving_avg[max(0, i - window):i + 1])
            for i in range(len(moving_avg))
        ])

        # Find first point where rolling std stays below threshold
        stable_indices = np.where(rolling_std < self.convergence_threshold)[0]

        if len(stable_indices) > 0:
            # Check if stability is maintained for at least 20% of remaining trajectory
            for idx in stable_indices:
                remaining_length = len(rolling_std) - idx
                if remaining_length > len(rolling_std) * 0.2:
                    # Check if mostly stable after this point
                    stable_fraction = np.sum(
                        rolling_std[idx:] < self.convergence_threshold
                    ) / remaining_length

                    if stable_fraction > 0.8:
                        convergence_time = moving_avg_times[idx]
                        convergence_fraction = convergence_time / times[-1]
                        return convergence_time, convergence_fraction, True

        # Not converged
        return times[-1], 1.0, False

    def assign_quality_grade(self, metrics: Dict) -> str:
        """
        Assign quality grade (A/B/C/D) based on convergence metrics

        Grading criteria (Stability-Based):
        - Grade A: CV < 0.10, Converged, Conv_frac < 0.2
        - Grade B: CV < 0.15 (regardless of convergence status)
        - Grade C: CV < 0.20
        - Grade D: CV >= 0.20

        Uses Coefficient of Variation (CV = Std/Mean) to assess stability.
        A trajectory with low CV indicates stable dynamics and deserves a
        good grade, even if formal "convergence" is not detected.

        Parameters
        ----------
        metrics : dict
            Convergence metrics from calculate_convergence_metrics()

        Returns
        -------
        str
            Quality grade: 'A', 'B', 'C', or 'D'
        """
        mean_rmsd = metrics['mean_rmsd']
        std_rmsd = metrics['std_rmsd']
        conv_frac = metrics['convergence_fraction']
        is_converged = metrics['is_converged']

        if mean_rmsd > 0.8 or std_rmsd > 0.15:
            return 'D'

        # Calculate coefficient of variation (CV)
        cv = std_rmsd / mean_rmsd if mean_rmsd > 0 else float('inf')

        # Grade A: Excellent stability AND early convergence
        if cv < 0.10 and is_converged and conv_frac < 0.2:
            return 'A'

        # Grade B: Good stability (convergence not required)
        if cv < 0.15:
            return 'B'

        # Grade C: Acceptable stability
        if cv < 0.20:
            return 'C'

        # Grade D: Poor stability
        return 'D'

    def generate_convergence_report(self, metrics: Dict) -> str:
        """
        Generate human-readable convergence report

        Parameters
        ----------
        metrics : dict
            Convergence metrics from calculate_convergence_metrics()

        Returns
        -------
        str
            Formatted convergence report
        """
        grade = self.assign_quality_grade(metrics)

        report_lines = [
            "=" * 60,
            "RMSD Convergence Analysis Report",
            "=" * 60,
            "",
            f"Quality Grade: {grade}",
            "",
            "Overall Statistics:",
            f"  Mean RMSD: {metrics['mean_rmsd']:.3f} nm",
            f"  Std RMSD: {metrics['std_rmsd']:.3f} nm",
            f"  CV (Std/Mean): {metrics['cv_rmsd']:.3f}",
            f"  Min RMSD: {metrics['min_rmsd']:.3f} nm",
            f"  Max RMSD: {metrics['max_rmsd']:.3f} nm",
            "",
            "Convergence:",
            f"  Converged: {'Yes' if metrics['is_converged'] else 'No'}",
            f"  Convergence Time: {metrics['convergence_time']:.1f}",
            f"  Convergence Fraction: {metrics['convergence_fraction']:.1%}",
            "",
            f"Block Analysis ({metrics['block_analysis']['n_blocks']} blocks):",
            f"  Block Stability (CV): {metrics['block_analysis']['block_stability']:.3f}",
            ""
        ]

        # Add block details
        for block in metrics['block_analysis']['block_info']:
            report_lines.append(
                f"  Block {block['block_id']}: "
                f"Mean={block['mean_rmsd']:.3f} nm, "
                f"Std={block['std_rmsd']:.3f} nm"
            )

        report_lines.extend([
            "",
            "Moving Average Analysis:",
            f"  Window Size: {metrics['moving_avg_analysis']['window_size']}",
            f"  Trend Slope: {metrics['moving_avg_analysis']['trend_slope']:.6f} nm/time_unit",
            f"  Has Drift: {'Yes' if metrics['moving_avg_analysis']['has_drift'] else 'No'}",
            "",
            "=" * 60
        ])

        return "\n".join(report_lines)
