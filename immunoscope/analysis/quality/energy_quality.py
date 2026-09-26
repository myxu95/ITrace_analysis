"""
Energy Quality Checker

Validates MD simulation quality based on energy file (md.edr) analysis.
This checker is independent of PBC processing and can be used for rapid
pre-screening of MD trajectories.
"""

import subprocess
import re
import tempfile
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import logging

logger = logging.getLogger(__name__)


class EnergyQualityChecker:
    """
    Energy-based quality checker for MD simulations.

    Advantages:
    - Independent of PBC processing
    - Fast screening based on md.edr file
    - Reliable thermodynamic indicators

    Quality checks:
    - Temperature stability and accuracy
    - Pressure stability (for NPT ensemble)
    - Total energy conservation
    - Potential/kinetic energy balance
    """

    def __init__(self,
                 target_temperature: float = 300.0,
                 temp_tolerance: float = 5.0,
                 target_pressure: float = 1.0,
                 pressure_tolerance: float = 50.0,
                 energy_drift_threshold: float = 2.0,
                 gmx_executable: str = "gmx"):
        """
        Initialize energy quality checker.

        Args:
            target_temperature: Target simulation temperature in K
            temp_tolerance: Temperature deviation tolerance in K
            target_pressure: Target pressure in bar (for NPT)
            pressure_tolerance: Pressure deviation tolerance in bar
            energy_drift_threshold: Maximum acceptable energy drift percentage
            gmx_executable: GROMACS executable command
        """
        self.target_temp = target_temperature
        self.temp_tolerance = temp_tolerance
        self.target_pressure = target_pressure
        self.pressure_tolerance = pressure_tolerance
        self.energy_drift_threshold = energy_drift_threshold
        self.gmx = gmx_executable

        self._check_gromacs()

    def _check_gromacs(self):
        """Check if GROMACS is available."""
        try:
            result = subprocess.run(
                [self.gmx, "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            logger.info("GROMACS found and accessible for energy analysis")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("GROMACS not found - energy analysis will not work")

    def _extract_energy_data(self,
                            edr_file: str,
                            energy_term: str,
                            output_file: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract energy data from .edr file using gmx energy.

        Args:
            edr_file: Path to .edr file
            energy_term: Energy term name (e.g., "Temperature", "Pressure", "Total-Energy")
            output_file: Optional output XVG file path

        Returns:
            Tuple of (time_array, energy_array)
        """
        if not Path(edr_file).exists():
            raise FileNotFoundError(f"Energy file not found: {edr_file}")

        # Create temporary output file if not specified
        if output_file is None:
            temp_fd, output_file = tempfile.mkstemp(suffix='.xvg')
            os.close(temp_fd)
            temp_output = True
        else:
            temp_output = False

        try:
            # Run gmx energy
            cmd = [
                self.gmx, "energy",
                "-f", edr_file,
                "-o", output_file
            ]

            stdin_input = f"{energy_term}\n"

            result = subprocess.run(
                cmd,
                input=stdin_input,
                text=True,
                capture_output=True,
                check=True
            )

            # Parse XVG file
            time_data, energy_data = self._parse_xvg_file(output_file)

            logger.info(f"Extracted {len(time_data)} data points for {energy_term}")
            return time_data, energy_data

        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to extract {energy_term} from {edr_file}")
            logger.error(f"stderr: {e.stderr}")
            raise

        finally:
            # Clean up temporary file
            if temp_output and Path(output_file).exists():
                try:
                    os.remove(output_file)
                except Exception as e:
                    logger.warning(f"Failed to remove temporary file {output_file}: {e}")

    def _parse_xvg_file(self, xvg_file: str) -> Tuple[np.ndarray, np.ndarray]:
        """Parse GROMACS XVG file and extract data."""
        time_data = []
        value_data = []

        with open(xvg_file, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if line.startswith('#') or line.startswith('@') or not line:
                    continue

                parts = line.split()
                if len(parts) >= 2:
                    try:
                        time_data.append(float(parts[0]))
                        value_data.append(float(parts[1]))
                    except ValueError:
                        continue

        return np.array(time_data), np.array(value_data)

    def _calculate_statistics(self, data: np.ndarray) -> Dict:
        """Calculate basic statistics for data array."""
        return {
            'mean': float(np.mean(data)),
            'std': float(np.std(data)),
            'min': float(np.min(data)),
            'max': float(np.max(data)),
            'median': float(np.median(data))
        }

    def _check_drift(self, time: np.ndarray, data: np.ndarray) -> Dict:
        """
        Check for systematic drift in data.

        Uses linear regression to detect trends.
        """
        if len(data) < 10:
            return {
                'has_drift': False,
                'drift_rate': 0.0,
                'drift_percentage': 0.0,
                'reason': 'insufficient_data'
            }

        # Linear regression
        coeffs = np.polyfit(time, data, 1)
        slope = coeffs[0]

        # Calculate drift as percentage change over simulation
        if len(data) > 0 and abs(data[0]) > 1e-10:
            time_span = time[-1] - time[0]
            total_drift = slope * time_span
            drift_percentage = abs(total_drift / data[0]) * 100
        else:
            drift_percentage = 0.0

        has_drift = drift_percentage > self.energy_drift_threshold

        return {
            'has_drift': bool(has_drift),
            'drift_rate': float(slope),
            'drift_percentage': float(drift_percentage),
            'assessment': 'significant_drift' if has_drift else 'stable'
        }

    def check_temperature_stability(self, edr_file: str) -> Dict:
        """
        Check temperature stability and accuracy.

        Args:
            edr_file: Path to md.edr file

        Returns:
            Dictionary with temperature quality assessment
        """
        logger.info(f"Checking temperature stability: {edr_file}")

        try:
            time, temp = self._extract_energy_data(edr_file, "Temperature")

            if len(temp) == 0:
                return {
                    'status': 'error',
                    'message': 'No temperature data found',
                    'grade': 'D'
                }

            stats = self._calculate_statistics(temp)
            drift = self._check_drift(time, temp)

            # Temperature deviation from target
            temp_deviation = abs(stats['mean'] - self.target_temp)
            temp_in_range = temp_deviation < self.temp_tolerance

            # Grade assessment
            if temp_deviation < 2.0 and stats['std'] < 5.0:
                grade = 'A'
                assessment = 'Excellent temperature control'
            elif temp_deviation < 5.0 and stats['std'] < 10.0:
                grade = 'B'
                assessment = 'Good temperature control'
            elif temp_deviation < 10.0 and stats['std'] < 15.0:
                grade = 'C'
                assessment = 'Acceptable temperature control'
            else:
                grade = 'D'
                assessment = 'Poor temperature control'

            result = {
                'status': 'success',
                'grade': grade,
                'assessment': assessment,
                'target_temperature': self.target_temp,
                'statistics': stats,
                'deviation': float(temp_deviation),
                'in_tolerance': temp_in_range,
                'drift': drift,
                'data_points': len(temp)
            }

            logger.info(f"Temperature check: Grade {grade}, Mean={stats['mean']:.2f}K, "
                       f"Std={stats['std']:.2f}K")

            return result

        except Exception as e:
            logger.error(f"Temperature check failed: {e}")
            return {
                'status': 'error',
                'message': str(e),
                'grade': 'D'
            }

    def check_pressure_stability(self, edr_file: str) -> Dict:
        """
        Check pressure stability (for NPT ensemble).

        Args:
            edr_file: Path to md.edr file

        Returns:
            Dictionary with pressure quality assessment
        """
        logger.info(f"Checking pressure stability: {edr_file}")

        try:
            time, pressure = self._extract_energy_data(edr_file, "Pressure")

            if len(pressure) == 0:
                return {
                    'status': 'skipped',
                    'message': 'No pressure data found (possibly NVT ensemble)',
                    'grade': 'N/A'
                }

            stats = self._calculate_statistics(pressure)

            # Pressure deviation from target
            pressure_deviation = abs(stats['mean'] - self.target_pressure)
            pressure_in_range = pressure_deviation < self.pressure_tolerance

            # Grade assessment (more lenient for pressure)
            if pressure_deviation < 20.0 and stats['std'] < 100.0:
                grade = 'A'
                assessment = 'Excellent pressure control'
            elif pressure_deviation < 50.0 and stats['std'] < 200.0:
                grade = 'B'
                assessment = 'Good pressure control'
            elif pressure_deviation < 100.0 and stats['std'] < 300.0:
                grade = 'C'
                assessment = 'Acceptable pressure control'
            else:
                grade = 'D'
                assessment = 'Poor pressure control'

            result = {
                'status': 'success',
                'grade': grade,
                'assessment': assessment,
                'target_pressure': self.target_pressure,
                'statistics': stats,
                'deviation': float(pressure_deviation),
                'in_tolerance': pressure_in_range,
                'data_points': len(pressure)
            }

            logger.info(f"Pressure check: Grade {grade}, Mean={stats['mean']:.1f}bar, "
                       f"Std={stats['std']:.1f}bar")

            return result

        except Exception as e:
            logger.error(f"Pressure check failed: {e}")
            return {
                'status': 'error',
                'message': str(e),
                'grade': 'N/A'
            }

    def check_total_energy_conservation(self, edr_file: str) -> Dict:
        """
        Check total energy conservation and drift.

        Args:
            edr_file: Path to md.edr file

        Returns:
            Dictionary with energy conservation assessment
        """
        logger.info(f"Checking energy conservation: {edr_file}")

        try:
            time, total_energy = self._extract_energy_data(edr_file, "Total-Energy")

            if len(total_energy) == 0:
                return {
                    'status': 'error',
                    'message': 'No total energy data found',
                    'grade': 'D'
                }

            stats = self._calculate_statistics(total_energy)
            drift = self._check_drift(time, total_energy)

            # Grade based on drift percentage
            drift_pct = drift['drift_percentage']
            if drift_pct < 0.5:
                grade = 'A'
                assessment = 'Excellent energy conservation'
            elif drift_pct < 1.0:
                grade = 'B'
                assessment = 'Good energy conservation'
            elif drift_pct < 2.0:
                grade = 'C'
                assessment = 'Acceptable energy conservation'
            else:
                grade = 'D'
                assessment = 'Poor energy conservation (significant drift)'

            result = {
                'status': 'success',
                'grade': grade,
                'assessment': assessment,
                'statistics': stats,
                'drift': drift,
                'data_points': len(total_energy)
            }

            logger.info(f"Energy conservation: Grade {grade}, Drift={drift_pct:.2f}%")

            return result

        except Exception as e:
            logger.error(f"Energy conservation check failed: {e}")
            return {
                'status': 'error',
                'message': str(e),
                'grade': 'D'
            }

    def check_potential_kinetic_balance(self, edr_file: str) -> Dict:
        """
        Check potential and kinetic energy balance.

        Args:
            edr_file: Path to md.edr file

        Returns:
            Dictionary with energy balance assessment
        """
        logger.info(f"Checking potential/kinetic energy balance: {edr_file}")

        try:
            # Extract both energies
            time_pot, potential = self._extract_energy_data(edr_file, "Potential")
            time_kin, kinetic = self._extract_energy_data(edr_file, "Kinetic-En.")

            if len(potential) == 0 or len(kinetic) == 0:
                return {
                    'status': 'error',
                    'message': 'Missing potential or kinetic energy data',
                    'grade': 'D'
                }

            pot_stats = self._calculate_statistics(potential)
            kin_stats = self._calculate_statistics(kinetic)

            # Calculate ratio and its stability
            if len(potential) == len(kinetic):
                ratio = potential / kinetic
                ratio_mean = float(np.mean(ratio))
                ratio_std = float(np.std(ratio))
                ratio_cv = float(ratio_std / abs(ratio_mean)) if abs(ratio_mean) > 1e-10 else 0.0
            else:
                ratio_mean = float(pot_stats['mean'] / kin_stats['mean'])
                ratio_std = 0.0
                ratio_cv = 0.0

            # Grade based on ratio stability
            if ratio_cv < 0.05:
                grade = 'A'
                assessment = 'Excellent energy balance'
            elif ratio_cv < 0.10:
                grade = 'B'
                assessment = 'Good energy balance'
            elif ratio_cv < 0.15:
                grade = 'C'
                assessment = 'Acceptable energy balance'
            else:
                grade = 'D'
                assessment = 'Poor energy balance'

            result = {
                'status': 'success',
                'grade': grade,
                'assessment': assessment,
                'potential_energy': pot_stats,
                'kinetic_energy': kin_stats,
                'ratio': {
                    'mean': ratio_mean,
                    'std': ratio_std,
                    'coefficient_of_variation': ratio_cv
                }
            }

            logger.info(f"Energy balance: Grade {grade}, Ratio={ratio_mean:.2f}, CV={ratio_cv:.3f}")

            return result

        except Exception as e:
            logger.error(f"Energy balance check failed: {e}")
            return {
                'status': 'error',
                'message': str(e),
                'grade': 'D'
            }

    def comprehensive_energy_check(self, edr_file: str) -> Dict:
        """
        Perform comprehensive energy quality check.

        Args:
            edr_file: Path to md.edr file

        Returns:
            Dictionary with comprehensive energy quality assessment
        """
        logger.info(f"Starting comprehensive energy quality check: {edr_file}")

        if not Path(edr_file).exists():
            return {
                'status': 'error',
                'message': f'Energy file not found: {edr_file}',
                'overall_status': 'fail',
                'energy_grade': 'D',
                'score': 0
            }

        # Run all checks
        temp_result = self.check_temperature_stability(edr_file)
        pressure_result = self.check_pressure_stability(edr_file)
        energy_result = self.check_total_energy_conservation(edr_file)
        balance_result = self.check_potential_kinetic_balance(edr_file)

        # Collect results
        results = {
            'edr_file': edr_file,
            'temperature': temp_result,
            'pressure': pressure_result,
            'energy_conservation': energy_result,
            'energy_balance': balance_result
        }

        # Calculate overall grade and score
        grades = [
            temp_result.get('grade', 'D'),
            energy_result.get('grade', 'D'),
            balance_result.get('grade', 'D')
        ]

        # Add pressure grade only if it was checked (NPT ensemble)
        if pressure_result.get('status') == 'success':
            grades.append(pressure_result.get('grade', 'D'))

        # Convert grades to scores
        grade_to_score = {'A': 95, 'B': 85, 'C': 70, 'D': 50, 'N/A': 0}
        scores = [grade_to_score.get(g, 50) for g in grades if g != 'N/A']

        if scores:
            overall_score = float(np.mean(scores))
        else:
            overall_score = 0.0

        # Determine overall grade
        if overall_score >= 90:
            overall_grade = 'A'
            overall_status = 'excellent'
        elif overall_score >= 75:
            overall_grade = 'B'
            overall_status = 'good'
        elif overall_score >= 60:
            overall_grade = 'C'
            overall_status = 'acceptable'
        else:
            overall_grade = 'D'
            overall_status = 'poor'

        # Collect issues
        issues = []
        if temp_result.get('grade') in ['C', 'D']:
            issues.append('temperature_instability')
        if pressure_result.get('grade') in ['C', 'D'] and pressure_result.get('status') == 'success':
            issues.append('pressure_instability')
        if energy_result.get('grade') in ['C', 'D']:
            issues.append('energy_drift')
        if balance_result.get('grade') in ['C', 'D']:
            issues.append('energy_imbalance')

        comprehensive_result = {
            'status': 'success',
            'overall_status': overall_status,
            'energy_grade': overall_grade,
            'score': overall_score,
            'issues': issues,
            'details': results,
            'summary': self._generate_summary(results, overall_grade, overall_score)
        }

        logger.info(f"Comprehensive energy check completed: Grade {overall_grade}, Score {overall_score:.1f}")

        return comprehensive_result

    def _generate_summary(self, results: Dict, grade: str, score: float) -> str:
        """Generate human-readable summary of energy quality check."""
        temp_grade = results['temperature'].get('grade', 'D')
        energy_grade = results['energy_conservation'].get('grade', 'D')

        summary_lines = [
            f"Overall Energy Quality: Grade {grade} (Score: {score:.1f}/100)",
            f"Temperature Control: Grade {temp_grade}",
            f"Energy Conservation: Grade {energy_grade}"
        ]

        if results['pressure'].get('status') == 'success':
            pressure_grade = results['pressure'].get('grade', 'N/A')
            summary_lines.append(f"Pressure Control: Grade {pressure_grade}")

        return " | ".join(summary_lines)
