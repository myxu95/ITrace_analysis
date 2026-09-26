"""
RMSD Analysis and Quality Control Module

This module provides tools to analyze RMSD trajectories, detect conformational
transitions, and identify potential artifacts.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


class RMSDAnalyzer:
    """Advanced RMSD analysis and quality control."""

    def __init__(self, rmsd_file: str):
        """
        Initialize RMSD analyzer.

        Args:
            rmsd_file: Path to RMSD XVG file
        """
        self.rmsd_file = rmsd_file
        self.time = None
        self.rmsd = None
        self._load_data()

    def _load_data(self):
        """Load RMSD data from XVG file."""
        data = []
        with open(self.rmsd_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    data.append([float(parts[0]), float(parts[1])])

        if not data:
            raise ValueError(f"No data found in {self.rmsd_file}")

        data = np.array(data)
        self.time = data[:, 0]
        self.rmsd = data[:, 1]

        logger.info(f"Loaded {len(self.rmsd)} frames from {self.rmsd_file}")

    def detect_transitions(self,
                          window_size: int = 50,
                          threshold_factor: float = 2.0) -> List[Dict]:
        """
        Detect major conformational transitions in RMSD trajectory.

        A transition is detected when RMSD changes significantly beyond
        normal fluctuations.

        Args:
            window_size: Window for calculating local statistics
            threshold_factor: Multiplier for std to define "significant change"

        Returns:
            List of transition events with details
        """
        transitions = []

        # Calculate rolling statistics
        rmsd_series = pd.Series(self.rmsd)
        rolling_mean = rmsd_series.rolling(window=window_size, center=True).mean()
        rolling_std = rmsd_series.rolling(window=window_size, center=True).std()

        # Detect significant changes
        for i in range(window_size, len(self.rmsd) - window_size):
            # Change from previous window
            delta = abs(rolling_mean.iloc[i] - rolling_mean.iloc[i - window_size])
            threshold = threshold_factor * rolling_std.iloc[i]

            if delta > threshold:
                transitions.append({
                    'frame': i,
                    'time_ps': self.time[i],
                    'time_ns': self.time[i] / 1000,
                    'rmsd_before': rolling_mean.iloc[i - window_size],
                    'rmsd_after': rolling_mean.iloc[i],
                    'delta': delta,
                    'type': 'increase' if rolling_mean.iloc[i] > rolling_mean.iloc[i - window_size] else 'decrease'
                })

        # Remove duplicates (keep only major transitions)
        filtered_transitions = []
        for trans in transitions:
            # Check if this is too close to a previous transition
            if not filtered_transitions or \
               (trans['frame'] - filtered_transitions[-1]['frame']) > window_size * 2:
                filtered_transitions.append(trans)

        logger.info(f"Detected {len(filtered_transitions)} major conformational transitions")
        return filtered_transitions

    def identify_stable_regions(self,
                               max_std: float = 0.1) -> List[Tuple[int, int, float]]:
        """
        Identify stable conformational regions (low RMSD fluctuation).

        Args:
            max_std: Maximum standard deviation for "stable" region (nm)

        Returns:
            List of (start_frame, end_frame, mean_rmsd) tuples
        """
        window_size = 100
        stable_regions = []

        rmsd_series = pd.Series(self.rmsd)
        rolling_std = rmsd_series.rolling(window=window_size, center=True).std()

        in_stable_region = False
        region_start = 0

        for i in range(len(self.rmsd)):
            if pd.notna(rolling_std.iloc[i]):
                is_stable = rolling_std.iloc[i] < max_std

                if is_stable and not in_stable_region:
                    # Start of stable region
                    region_start = i
                    in_stable_region = True
                elif not is_stable and in_stable_region:
                    # End of stable region
                    if i - region_start > window_size:  # Minimum length
                        mean_rmsd = np.mean(self.rmsd[region_start:i])
                        stable_regions.append((region_start, i, mean_rmsd))
                    in_stable_region = False

        # Close last region if still open
        if in_stable_region and len(self.rmsd) - region_start > window_size:
            mean_rmsd = np.mean(self.rmsd[region_start:])
            stable_regions.append((region_start, len(self.rmsd), mean_rmsd))

        logger.info(f"Identified {len(stable_regions)} stable conformational regions")
        return stable_regions

    def flag_outliers(self,
                     z_score_threshold: float = 3.0,
                     min_consecutive: int = 1) -> Dict[str, List]:
        """
        Flag frames with unusually high or low RMSD values.

        Args:
            z_score_threshold: Z-score threshold for outlier detection
            min_consecutive: Minimum consecutive frames to flag (1=isolated spikes)

        Returns:
            Dictionary with 'high_outliers' and 'low_outliers' lists
        """
        mean_rmsd = np.mean(self.rmsd)
        std_rmsd = np.std(self.rmsd)

        z_scores = (self.rmsd - mean_rmsd) / std_rmsd

        high_outliers = []
        low_outliers = []

        # Find consecutive outliers
        consecutive_count = 0
        current_outliers = []

        for i, z in enumerate(z_scores):
            if abs(z) > z_score_threshold:
                consecutive_count += 1
                current_outliers.append(i)
            else:
                if consecutive_count >= min_consecutive:
                    # Determine if high or low
                    if z_scores[current_outliers[0]] > 0:
                        high_outliers.extend(current_outliers)
                    else:
                        low_outliers.extend(current_outliers)

                consecutive_count = 0
                current_outliers = []

        # Close last group
        if consecutive_count >= min_consecutive:
            if z_scores[current_outliers[0]] > 0:
                high_outliers.extend(current_outliers)
            else:
                low_outliers.extend(current_outliers)

        result = {
            'high_outliers': high_outliers,
            'low_outliers': low_outliers,
            'z_scores': z_scores
        }

        logger.info(f"Flagged {len(high_outliers)} high outliers and {len(low_outliers)} low outliers")
        return result

    def classify_outliers(self,
                         outliers: List[int],
                         consecutive_threshold: int = 10) -> Dict:
        """
        Classify outliers as artifacts or conformational changes.

        Args:
            outliers: List of outlier frame indices
            consecutive_threshold: Frames threshold for "sustained" vs "spike"

        Returns:
            Classification results
        """
        if not outliers:
            return {'artifacts': [], 'conformational_changes': []}

        # Group consecutive outliers
        groups = []
        current_group = [outliers[0]]

        for i in range(1, len(outliers)):
            if outliers[i] - outliers[i-1] <= 2:  # Allow 1-frame gap
                current_group.append(outliers[i])
            else:
                groups.append(current_group)
                current_group = [outliers[i]]
        groups.append(current_group)

        # Classify groups
        artifacts = []
        conformational_changes = []

        for group in groups:
            if len(group) < consecutive_threshold:
                # Short spike = likely artifact
                artifacts.append({
                    'frames': group,
                    'duration': len(group),
                    'time_range': (self.time[group[0]], self.time[group[-1]]),
                    'type': 'isolated_spike'
                })
            else:
                # Sustained high RMSD = likely conformational change
                conformational_changes.append({
                    'frames': group,
                    'duration': len(group),
                    'time_range': (self.time[group[0]], self.time[group[-1]]),
                    'mean_rmsd': np.mean(self.rmsd[group]),
                    'type': 'sustained_change'
                })

        return {
            'artifacts': artifacts,
            'conformational_changes': conformational_changes
        }

    def generate_report(self, output_file: Optional[str] = None) -> str:
        """
        Generate comprehensive RMSD analysis report.

        Args:
            output_file: Optional file to save report

        Returns:
            Report as string
        """
        report_lines = []

        report_lines.append("=" * 80)
        report_lines.append("RMSD Analysis Report")
        report_lines.append("=" * 80)
        report_lines.append(f"Input file: {self.rmsd_file}")
        report_lines.append(f"Total frames: {len(self.rmsd)}")
        report_lines.append(f"Time range: {self.time[0]/1000:.1f} - {self.time[-1]/1000:.1f} ns")
        report_lines.append("")

        # Basic statistics
        report_lines.append("Basic Statistics:")
        report_lines.append(f"  Mean RMSD: {np.mean(self.rmsd):.3f} nm")
        report_lines.append(f"  Std RMSD:  {np.std(self.rmsd):.3f} nm")
        report_lines.append(f"  Min RMSD:  {np.min(self.rmsd):.3f} nm")
        report_lines.append(f"  Max RMSD:  {np.max(self.rmsd):.3f} nm")
        report_lines.append("")

        # Detect transitions
        transitions = self.detect_transitions()
        report_lines.append(f"Conformational Transitions: {len(transitions)}")
        for i, trans in enumerate(transitions, 1):
            report_lines.append(f"  {i}. Frame {trans['frame']} ({trans['time_ns']:.1f} ns)")
            report_lines.append(f"     Type: {trans['type']}")
            report_lines.append(f"     RMSD: {trans['rmsd_before']:.2f} -> {trans['rmsd_after']:.2f} nm")
            report_lines.append(f"     Delta: {trans['delta']:.2f} nm")
        report_lines.append("")

        # Stable regions
        stable_regions = self.identify_stable_regions()
        report_lines.append(f"Stable Conformational Regions: {len(stable_regions)}")
        for i, (start, end, mean_rmsd) in enumerate(stable_regions, 1):
            report_lines.append(f"  {i}. Frames {start}-{end} ({self.time[start]/1000:.1f}-{self.time[end-1]/1000:.1f} ns)")
            report_lines.append(f"     Duration: {(self.time[end-1] - self.time[start])/1000:.1f} ns")
            report_lines.append(f"     Mean RMSD: {mean_rmsd:.3f} nm")
        report_lines.append("")

        # Outliers
        outliers_info = self.flag_outliers()
        high_outliers = outliers_info['high_outliers']
        classification = self.classify_outliers(high_outliers)

        report_lines.append(f"Outlier Analysis:")
        report_lines.append(f"  High outliers: {len(high_outliers)} frames")
        report_lines.append(f"  Artifacts (isolated spikes): {len(classification['artifacts'])}")
        report_lines.append(f"  Conformational changes: {len(classification['conformational_changes'])}")

        if classification['artifacts']:
            report_lines.append("\n  Potential Artifacts:")
            for artifact in classification['artifacts']:
                report_lines.append(f"    - Frames {artifact['frames'][0]}-{artifact['frames'][-1]}")
                report_lines.append(f"      Time: {artifact['time_range'][0]/1000:.1f}-{artifact['time_range'][1]/1000:.1f} ns")

        if classification['conformational_changes']:
            report_lines.append("\n  Major Conformational Changes:")
            for change in classification['conformational_changes']:
                report_lines.append(f"    - Frames {change['frames'][0]}-{change['frames'][-1]}")
                report_lines.append(f"      Time: {change['time_range'][0]/1000:.1f}-{change['time_range'][1]/1000:.1f} ns")
                report_lines.append(f"      Mean RMSD: {change['mean_rmsd']:.2f} nm")

        report_lines.append("\n" + "=" * 80)

        report = "\n".join(report_lines)

        if output_file:
            with open(output_file, 'w') as f:
                f.write(report)
            logger.info(f"Report saved to {output_file}")

        return report
