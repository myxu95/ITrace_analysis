#!/usr/bin/env python3
"""
Trajectory Convergence Checker

Lightweight convergence assessment using energy data from md.edr.
Pre-PBC compatible - no trajectory processing required.

Core indicators:
1. Potential energy block analysis
2. Total energy drift detection
3. Temperature stability
4. Convergence time estimation
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from scipy import stats
from scipy.ndimage import uniform_filter1d


class TrajectoryConvergenceChecker:
    """
    Assess trajectory convergence using energy-based analysis.

    Provides quick Pre-PBC quality screening to determine if trajectory
    has reached equilibrium and identify the convergence time point.
    """

    def __init__(self,
                 n_blocks: int = 5,
                 convergence_threshold: float = 0.005,
                 gmx_executable: str = "gmx"):
        """
        Initialize convergence checker.

        Args:
            n_blocks: Number of blocks to divide trajectory (default 5)
            convergence_threshold: Threshold for convergence (0.005 = 0.5%)
            gmx_executable: GROMACS executable name
        """
        self.n_blocks = n_blocks
        self.threshold = convergence_threshold
        self.gmx_exe = gmx_executable

    def check_convergence(self, edr_file: str) -> Dict:
        """
        Perform comprehensive convergence analysis.

        Args:
            edr_file: Path to GROMACS energy file (md.edr)

        Returns:
            Dictionary with convergence assessment results
        """
        edr_path = Path(edr_file)
        if not edr_path.exists():
            return {
                'status': 'error',
                'message': f'Energy file not found: {edr_file}'
            }

        # Extract energy data
        potential = self._extract_energy_term(edr_file, 'Potential')
        total = self._extract_energy_term(edr_file, 'Total-Energy')
        temperature = self._extract_energy_term(edr_file, 'Temperature')

        if potential is None or total is None:
            return {
                'status': 'error',
                'message': 'Failed to extract energy data from EDR file'
            }

        # 1. Block analysis on potential energy
        pot_analysis = self._block_analysis(potential)

        # 2. Drift analysis on total energy
        drift_analysis = self._drift_analysis(total)

        # 3. Temperature stability
        temp_analysis = self._temperature_stability(temperature) if temperature is not None else None

        # 4. Estimate convergence time
        conv_time_ps = self._estimate_convergence_time(potential)

        # Overall convergence decision
        is_converged = self._make_convergence_decision(
            pot_analysis, drift_analysis, temp_analysis
        )

        # Assign grade
        grade = self._assign_grade(is_converged, pot_analysis, drift_analysis)

        return {
            'status': 'success',
            'is_converged': is_converged,
            'convergence_time_ps': conv_time_ps,
            'convergence_time_ns': conv_time_ps / 1000 if conv_time_ps else None,
            'grade': grade,
            'potential_blocks': pot_analysis,
            'energy_drift': drift_analysis,
            'temperature_stability': temp_analysis,
            'recommendation': self._generate_recommendation(
                is_converged, conv_time_ps, grade
            )
        }

    def _extract_energy_term(self, edr_file: str, term: str) -> Optional[np.ndarray]:
        """
        Extract energy term from EDR file using gmx energy.

        Args:
            edr_file: Path to energy file
            term: Energy term name (e.g., 'Potential', 'Temperature')

        Returns:
            Numpy array of energy values, or None if extraction fails
        """
        try:
            # Create temporary XVG file
            import tempfile
            with tempfile.NamedTemporaryFile(suffix='.xvg', delete=False) as tmp:
                tmp_xvg = tmp.name

            # Run gmx energy
            cmd = f'echo "{term}" | {self.gmx_exe} energy -f {edr_file} -o {tmp_xvg}'
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=60
            )

            if result.returncode != 0:
                return None

            # Parse XVG file
            data = []
            with open(tmp_xvg, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.startswith('@'):
                        continue
                    try:
                        values = line.split()
                        if len(values) >= 2:
                            data.append(float(values[1]))
                    except (ValueError, IndexError):
                        continue

            # Clean up
            Path(tmp_xvg).unlink(missing_ok=True)

            return np.array(data) if data else None

        except Exception as e:
            print(f"Warning: Failed to extract {term}: {e}")
            return None

    def _block_analysis(self, data: np.ndarray) -> Dict:
        """
        Divide trajectory into blocks and analyze convergence.

        Args:
            data: Energy data array

        Returns:
            Dictionary with block analysis results
        """
        n = len(data)
        block_size = n // self.n_blocks

        block_means = []
        block_stds = []

        for i in range(self.n_blocks):
            start = i * block_size
            end = start + block_size if i < self.n_blocks - 1 else n
            block = data[start:end]

            block_means.append(np.mean(block))
            block_stds.append(np.std(block))

        block_means = np.array(block_means)
        block_stds = np.array(block_stds)

        # Check if last two blocks are similar (local convergence)
        if len(block_means) >= 2:
            last_two_diff = abs(block_means[-1] - block_means[-2]) / abs(block_means[-1])
            last_two_converged = last_two_diff < self.threshold
        else:
            last_two_diff = 0.0
            last_two_converged = True

        # Check overall drift (first vs last block)
        first_last_diff = abs(block_means[-1] - block_means[0]) / abs(block_means[0])
        no_strong_drift = first_last_diff < 0.02  # 2% threshold

        # Check if blocks show decreasing trend (still equilibrating)
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            range(len(block_means)), block_means
        )

        has_trend = abs(r_value) > 0.7 and p_value < 0.05

        return {
            'block_means': block_means.tolist(),
            'block_stds': block_stds.tolist(),
            'last_two_diff_percent': last_two_diff * 100,
            'last_two_converged': last_two_converged,
            'first_last_diff_percent': first_last_diff * 100,
            'no_strong_drift': no_strong_drift,
            'linear_slope': slope,
            'r_squared': r_value**2,
            'has_significant_trend': has_trend,
            'is_converged': last_two_converged and no_strong_drift and not has_trend
        }

    def _drift_analysis(self, data: np.ndarray) -> Dict:
        """
        Analyze long-term energy drift.

        Args:
            data: Energy data array

        Returns:
            Dictionary with drift analysis results
        """
        # Linear regression to detect drift
        x = np.arange(len(data))
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, data)

        # Calculate drift percentage
        mean_value = np.mean(data)
        drift_per_frame = slope
        total_drift = drift_per_frame * len(data)
        drift_percent = abs(total_drift / mean_value) * 100

        # Check if drift is significant
        has_drift = drift_percent > 2.0 and p_value < 0.05

        return {
            'mean': float(mean_value),
            'std': float(np.std(data)),
            'slope': float(slope),
            'drift_percent': float(drift_percent),
            'r_squared': float(r_value**2),
            'p_value': float(p_value),
            'has_significant_drift': has_drift,
            'drift_direction': 'increasing' if slope > 0 else 'decreasing'
        }

    def _temperature_stability(self, data: np.ndarray) -> Dict:
        """
        Check temperature stability throughout simulation.

        Args:
            data: Temperature data array

        Returns:
            Dictionary with temperature stability results
        """
        mean_temp = np.mean(data)
        std_temp = np.std(data)

        # Temperature should not drift significantly
        slope, _, r_value, p_value, _ = stats.linregress(
            np.arange(len(data)), data
        )

        has_drift = abs(r_value) > 0.5 and p_value < 0.05

        # Check stability (std should be small)
        is_stable = std_temp < 5.0 and not has_drift

        return {
            'mean': float(mean_temp),
            'std': float(std_temp),
            'min': float(np.min(data)),
            'max': float(np.max(data)),
            'has_drift': has_drift,
            'is_stable': is_stable
        }

    def _estimate_convergence_time(self, data: np.ndarray) -> Optional[float]:
        """
        Estimate when the trajectory converged.

        Uses moving average to identify when energy stabilizes.

        Args:
            data: Potential energy data

        Returns:
            Convergence time in ps, or None if not converged
        """
        if len(data) < 100:
            return None

        # Calculate moving average and moving std
        window = max(100, len(data) // 20)
        moving_avg = uniform_filter1d(data, size=window)

        # Calculate moving standard deviation
        moving_std = np.array([
            np.std(data[max(0, i-window):min(len(data), i+window)])
            for i in range(len(data))
        ])

        # Find where std stabilizes (use last 20% as reference)
        ref_start = int(len(data) * 0.8)
        ref_std = np.mean(moving_std[ref_start:])

        # Find first point where std stays within 20% of reference
        threshold = ref_std * 1.2

        for i in range(window, len(data) - window):
            if moving_std[i] <= threshold:
                # Check if it stays stable for at least 10% of trajectory
                check_length = max(100, len(data) // 10)
                if all(moving_std[i:i+check_length] <= threshold):
                    # Convert frame index to time (assume 2 ps per frame)
                    return float(i * 2.0)

        # If no clear convergence point, use 20% mark as estimate
        return float(len(data) * 0.2 * 2.0)

    def _make_convergence_decision(self,
                                   pot_analysis: Dict,
                                   drift_analysis: Dict,
                                   temp_analysis: Optional[Dict]) -> bool:
        """
        Make final convergence decision based on all analyses.

        Args:
            pot_analysis: Potential energy block analysis
            drift_analysis: Energy drift analysis
            temp_analysis: Temperature stability analysis

        Returns:
            True if converged, False otherwise
        """
        # Core criteria
        pot_converged = pot_analysis['is_converged']
        no_drift = not drift_analysis['has_significant_drift']

        # Temperature check (if available)
        temp_ok = True
        if temp_analysis is not None:
            temp_ok = temp_analysis['is_stable']

        return pot_converged and no_drift and temp_ok

    def _assign_grade(self,
                     is_converged: bool,
                     pot_analysis: Dict,
                     drift_analysis: Dict) -> str:
        """
        Assign quality grade based on convergence analysis.

        Grades:
        A: Excellent convergence
        B: Good convergence with minor issues
        C: Acceptable but needs attention
        D: Poor convergence

        Args:
            is_converged: Overall convergence status
            pot_analysis: Potential energy analysis
            drift_analysis: Energy drift analysis

        Returns:
            Grade string (A/B/C/D)
        """
        if not is_converged:
            # Check severity
            if drift_analysis['drift_percent'] > 5.0:
                return 'D'  # Severe drift
            elif pot_analysis['has_significant_trend']:
                return 'D'  # Still equilibrating
            else:
                return 'C'  # Minor convergence issues

        # Converged - check quality
        last_two_diff = pot_analysis['last_two_diff_percent']
        drift_percent = drift_analysis['drift_percent']

        if last_two_diff < 0.1 and drift_percent < 0.5:
            return 'A'  # Excellent
        elif last_two_diff < 0.5 and drift_percent < 1.0:
            return 'B'  # Good
        else:
            return 'C'  # Acceptable

    def _generate_recommendation(self,
                                is_converged: bool,
                                conv_time_ps: Optional[float],
                                grade: str) -> str:
        """
        Generate actionable recommendation based on analysis.

        Args:
            is_converged: Convergence status
            conv_time_ps: Convergence time in ps
            grade: Quality grade

        Returns:
            Recommendation string
        """
        if not is_converged:
            if grade == 'D':
                return (
                    "Trajectory has not converged. System may still be equilibrating. "
                    "Consider extending simulation time or checking for system instabilities."
                )
            else:
                return (
                    "Trajectory shows minor convergence issues. "
                    "Proceed with caution and validate results carefully."
                )

        conv_time_ns = conv_time_ps / 1000 if conv_time_ps else 0

        if grade == 'A':
            return (
                f"Excellent convergence achieved (estimated ~{conv_time_ns:.1f} ns). "
                f"Recommend discarding first {conv_time_ns:.1f} ns for production analysis."
            )
        elif grade == 'B':
            return (
                f"Good convergence achieved (estimated ~{conv_time_ns:.1f} ns). "
                f"Recommend discarding first {conv_time_ns:.1f} ns for production analysis."
            )
        else:
            return (
                f"Acceptable convergence (estimated ~{conv_time_ns:.1f} ns). "
                f"Recommend careful validation and discarding first {conv_time_ns:.1f} ns."
            )

    def plot_convergence_analysis(self,
                                 edr_file: str,
                                 output_path: Optional[str] = None,
                                 dpi: int = 300) -> None:
        """
        Generate comprehensive convergence analysis plot.

        Creates a 3-panel figure showing:
        1. Potential energy with block divisions
        2. Block means comparison
        3. Temperature stability (if available)

        Args:
            edr_file: Path to energy file
            output_path: Path to save plot (optional)
            dpi: Plot resolution
        """
        # Extract data
        potential = self._extract_energy_term(edr_file, 'Potential')
        temperature = self._extract_energy_term(edr_file, 'Temperature')

        if potential is None:
            print("Error: Could not extract potential energy data")
            return

        # Perform analysis
        result = self.check_convergence(edr_file)
        pot_analysis = result['potential_blocks']

        # Setup plot
        has_temp = temperature is not None
        n_panels = 3 if has_temp else 2

        fig = plt.figure(figsize=(12, 4*n_panels))
        gs = fig.add_gridspec(n_panels, 1, hspace=0.3)

        # Panel 1: Potential energy with blocks
        ax1 = fig.add_subplot(gs[0])

        time_ps = np.arange(len(potential)) * 2.0  # Assume 2 ps per frame
        time_ns = time_ps / 1000

        ax1.plot(time_ns, potential, linewidth=0.5, alpha=0.6, label='Potential Energy')

        # Add smoothed line
        if len(potential) > 200:
            window = len(potential) // 50
            smoothed = uniform_filter1d(potential, size=window)
            ax1.plot(time_ns, smoothed, linewidth=2.0, color='orange',
                    label='Smoothed', alpha=0.9)

        # Mark block boundaries
        block_size = len(potential) // self.n_blocks
        colors = sns.color_palette("husl", self.n_blocks)

        for i in range(self.n_blocks):
            start_idx = i * block_size
            end_idx = start_idx + block_size if i < self.n_blocks - 1 else len(potential)
            start_time = time_ns[start_idx]
            end_time = time_ns[end_idx-1]

            ax1.axvspan(start_time, end_time, alpha=0.1, color=colors[i])

            # Add block label
            mid_time = (start_time + end_time) / 2
            ax1.text(mid_time, ax1.get_ylim()[1], f'Block {i+1}',
                    ha='center', va='top', fontsize=9, color=colors[i])

        # Mark convergence time
        if result['convergence_time_ns']:
            conv_time = result['convergence_time_ns']
            ax1.axvline(conv_time, color='red', linestyle='--', linewidth=2,
                       label=f'Est. Convergence: {conv_time:.1f} ns')

        ax1.set_xlabel('Time (ns)', fontsize=12)
        ax1.set_ylabel('Potential Energy (kJ/mol)', fontsize=12)
        ax1.set_title('Potential Energy Evolution', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Panel 2: Block means comparison
        ax2 = fig.add_subplot(gs[1])

        block_means = pot_analysis['block_means']
        block_stds = pot_analysis['block_stds']
        x_pos = np.arange(len(block_means))

        bars = ax2.bar(x_pos, block_means, yerr=block_stds,
                      color=colors[:len(block_means)], alpha=0.7,
                      capsize=5, edgecolor='black', linewidth=1.5)

        # Add value labels on bars
        for i, (bar, mean, std) in enumerate(zip(bars, block_means, block_stds)):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{mean:.0f}\n±{std:.0f}',
                    ha='center', va='bottom', fontsize=9)

        # Add horizontal line for overall mean
        overall_mean = np.mean(potential)
        ax2.axhline(overall_mean, color='red', linestyle='--', linewidth=2,
                   label=f'Overall Mean: {overall_mean:.1f}')

        ax2.set_xlabel('Block Number', fontsize=12)
        ax2.set_ylabel('Mean Potential Energy (kJ/mol)', fontsize=12)
        ax2.set_title('Block-wise Energy Comparison', fontsize=14, fontweight='bold')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels([f'Block {i+1}' for i in range(len(block_means))])
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3, axis='y')

        # Panel 3: Temperature (if available)
        if has_temp and temperature is not None:
            ax3 = fig.add_subplot(gs[2])

            temp_time_ns = np.arange(len(temperature)) * 2.0 / 1000
            ax3.plot(temp_time_ns, temperature, linewidth=0.5, alpha=0.6,
                    label='Temperature')

            # Add smoothed line
            if len(temperature) > 200:
                window = len(temperature) // 50
                temp_smoothed = uniform_filter1d(temperature, size=window)
                ax3.plot(temp_time_ns, temp_smoothed, linewidth=2.0,
                        color='orange', label='Smoothed', alpha=0.9)

            # Add mean line
            temp_mean = np.mean(temperature)
            temp_std = np.std(temperature)
            ax3.axhline(temp_mean, color='red', linestyle='--', linewidth=2,
                       label=f'Mean: {temp_mean:.1f}±{temp_std:.1f} K')

            # Add ±1 std band
            ax3.fill_between(temp_time_ns,
                           temp_mean - temp_std,
                           temp_mean + temp_std,
                           alpha=0.2, color='red')

            ax3.set_xlabel('Time (ns)', fontsize=12)
            ax3.set_ylabel('Temperature (K)', fontsize=12)
            ax3.set_title('Temperature Stability', fontsize=14, fontweight='bold')
            ax3.legend(fontsize=10)
            ax3.grid(True, alpha=0.3)

        # Add overall summary
        summary_text = (
            f"Convergence Status: {'CONVERGED' if result['is_converged'] else 'NOT CONVERGED'}\n"
            f"Grade: {result['grade']}\n"
            f"Convergence Time: ~{result['convergence_time_ns']:.1f} ns\n"
            f"Last 2 Blocks Diff: {pot_analysis['last_two_diff_percent']:.2f}%"
        )

        fig.text(0.02, 0.98, summary_text, transform=fig.transFigure,
                fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if output_path:
            plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
            print(f"Convergence plot saved to: {output_path}")

        plt.show()
