"""
FEL Calculator Module

Core class for computing free energy landscapes from collective variables.
"""

import numpy as np
from typing import Optional, Union, Tuple, List, Dict, Any
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks


class FELCalculator:
    """
    Free Energy Landscape Calculator

    Calculate 2D/1D free energy landscapes from collective variable time series
    using histogram-based methods.

    Attributes:
        temperature (float): Temperature in Kelvin
        kT (float): Boltzmann constant * temperature in kJ/mol
    """

    def __init__(self, temperature: float = 300.0):
        """
        Initialize FEL calculator.

        Args:
            temperature: Temperature in Kelvin (default: 300.0)
        """
        self.temperature = temperature
        # Boltzmann constant: 8.314 J/(mol*K) = 8.314e-3 kJ/(mol*K)
        self.kT = 8.314e-3 * temperature

    def _freedman_diaconis_bins(self, data: np.ndarray) -> int:
        """
        Calculate optimal bin number using Freedman-Diaconis rule.

        Formula: bin_width = 2 * IQR / N^(1/3)

        Args:
            data: 1D array of CV values

        Returns:
            Number of bins (constrained to [10, 100])
        """
        q75, q25 = np.percentile(data, [75, 25])
        iqr = q75 - q25

        if iqr == 0:
            # Fallback to Sturges' rule if IQR is zero
            n_bins = int(np.ceil(np.log2(len(data))) + 1)
        else:
            bin_width = 2.0 * iqr / (len(data) ** (1.0/3.0))
            data_range = data.max() - data.min()
            n_bins = int(np.ceil(data_range / bin_width))

        # Constrain to reasonable range
        n_bins = max(10, min(n_bins, 100))

        return n_bins

    def compute_2d_fel(
        self,
        cv1_data: np.ndarray,
        cv2_data: np.ndarray,
        bins: Optional[Union[int, Tuple[int, int], str]] = 'auto',
        smooth: bool = False,
        sigma: float = 1.0,
        vmax: Optional[float] = None
    ) -> Dict[str, Any]:
        """
        Compute 2D free energy landscape.

        Formula: ΔG(cv1, cv2) = -kT ln[P(cv1, cv2)]
        where P(cv1, cv2) is normalized probability from histogram.

        Args:
            cv1_data: First CV time series (1D array)
            cv2_data: Second CV time series (1D array)
            bins: Bin specification:
                  - 'auto': Freedman-Diaconis rule (default)
                  - int: Same number of bins for both CVs
                  - tuple: (n_bins_cv1, n_bins_cv2)
            smooth: Apply Gaussian smoothing to FEL
            sigma: Gaussian kernel standard deviation (if smooth=True)
            vmax: Maximum free energy value (kJ/mol), higher values set to vmax

        Returns:
            Dictionary containing:
                - 'free_energy': 2D array of relative free energy (kJ/mol)
                - 'cv1_edges': CV1 bin edges
                - 'cv2_edges': CV2 bin edges
                - 'cv1_centers': CV1 bin centers
                - 'cv2_centers': CV2 bin centers
                - 'probability': Normalized probability distribution
                - 'histogram': Raw count histogram
                - 'n_bins': Actual bins used (n_cv1, n_cv2)
        """
        # Input validation
        if len(cv1_data) != len(cv2_data):
            raise ValueError(f"CV data length mismatch: {len(cv1_data)} vs {len(cv2_data)}")

        if len(cv1_data) == 0:
            raise ValueError("Empty CV data provided")

        # Determine bin numbers
        if bins == 'auto':
            n_bins_cv1 = self._freedman_diaconis_bins(cv1_data)
            n_bins_cv2 = self._freedman_diaconis_bins(cv2_data)
        elif isinstance(bins, int):
            n_bins_cv1 = n_bins_cv2 = bins
        elif isinstance(bins, (tuple, list)) and len(bins) == 2:
            n_bins_cv1, n_bins_cv2 = bins
        else:
            raise ValueError(f"Invalid bins parameter: {bins}")

        # Compute 2D histogram
        H, cv1_edges, cv2_edges = np.histogram2d(
            cv1_data, cv2_data,
            bins=[n_bins_cv1, n_bins_cv2],
            density=False
        )

        # Compute probability distribution (normalize)
        total_count = np.sum(H)
        if total_count == 0:
            raise ValueError("Histogram is empty - check CV data ranges")

        P = H / total_count

        # Handle zero probability (avoid log(0))
        # Strategy: Set to 1% of minimum non-zero probability
        min_nonzero_prob = np.min(P[P > 0]) if np.any(P > 0) else 1e-10
        P_safe = np.where(P > 0, P, min_nonzero_prob * 0.01)

        # Compute free energy: ΔG = -kT ln(P)
        G = -self.kT * np.log(P_safe)

        # Relative free energy (set minimum to zero)
        G_relative = G - np.min(G)

        # Apply maximum cutoff if specified
        if vmax is not None:
            G_relative = np.where(G_relative > vmax, vmax, G_relative)

        # Optional smoothing
        if smooth:
            G_relative = gaussian_filter(G_relative, sigma=sigma)
            # Re-zero after smoothing
            G_relative = G_relative - np.min(G_relative)

        # Compute bin centers
        cv1_centers = 0.5 * (cv1_edges[:-1] + cv1_edges[1:])
        cv2_centers = 0.5 * (cv2_edges[:-1] + cv2_edges[1:])

        return {
            'free_energy': G_relative,
            'cv1_edges': cv1_edges,
            'cv2_edges': cv2_edges,
            'cv1_centers': cv1_centers,
            'cv2_centers': cv2_centers,
            'probability': P,
            'histogram': H,
            'n_bins': (n_bins_cv1, n_bins_cv2)
        }

    def compute_1d_projection(
        self,
        cv_data: np.ndarray,
        bins: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute 1D free energy projection.

        Args:
            cv_data: CV time series
            bins: Number of bins (auto if None)

        Returns:
            Tuple of (cv_centers, free_energy_1d)
        """
        if bins is None:
            bins = self._freedman_diaconis_bins(cv_data)

        # Compute histogram
        H, edges = np.histogram(cv_data, bins=bins, density=False)

        # Probability
        P = H / np.sum(H)

        # Handle zeros
        min_nonzero_prob = np.min(P[P > 0]) if np.any(P > 0) else 1e-10
        P_safe = np.where(P > 0, P, min_nonzero_prob * 0.01)

        # Free energy
        G = -self.kT * np.log(P_safe)
        G_relative = G - np.min(G)

        # Bin centers
        centers = 0.5 * (edges[:-1] + edges[1:])

        return centers, G_relative

    def find_minima(
        self,
        free_energy: np.ndarray,
        cv1_centers: np.ndarray,
        cv2_centers: np.ndarray,
        threshold: float = 5.0,
        method: str = 'local'
    ) -> List[Tuple[float, float, float]]:
        """
        Identify energy minima (basins) in FEL.

        Args:
            free_energy: 2D free energy array
            cv1_centers: CV1 bin centers
            cv2_centers: CV2 bin centers
            threshold: Energy threshold above global minimum (kJ/mol)
            method: Detection method:
                    - 'local': Find all local minima
                    - 'global': Only global minimum

        Returns:
            List of (cv1, cv2, free_energy) tuples, sorted by energy (ascending)
        """
        minima = []

        # Global minimum
        global_min_idx = np.unravel_index(np.argmin(free_energy), free_energy.shape)
        global_min_energy = free_energy[global_min_idx]
        global_min_cv1 = cv1_centers[global_min_idx[0]]
        global_min_cv2 = cv2_centers[global_min_idx[1]]

        minima.append((global_min_cv1, global_min_cv2, global_min_energy))

        if method == 'global':
            return minima

        # Local minima detection
        # Strategy: Find points that are lower than all 8 neighbors
        rows, cols = free_energy.shape

        for i in range(1, rows - 1):
            for j in range(1, cols - 1):
                center_value = free_energy[i, j]

                # Check if it's a local minimum (lower than all neighbors)
                neighborhood = free_energy[i-1:i+2, j-1:j+2]
                neighbors = neighborhood[neighborhood != center_value]

                if len(neighbors) > 0 and center_value < np.min(neighbors):
                    # Check if within threshold
                    if center_value - global_min_energy <= threshold:
                        cv1_val = cv1_centers[i]
                        cv2_val = cv2_centers[j]

                        # Avoid duplicate of global minimum
                        if not (i == global_min_idx[0] and j == global_min_idx[1]):
                            minima.append((cv1_val, cv2_val, center_value))

        # Sort by energy
        minima.sort(key=lambda x: x[2])

        return minima

    def extract_representative_structures(
        self,
        cv1_data: np.ndarray,
        cv2_data: np.ndarray,
        minima: List[Tuple[float, float, float]],
        distance_threshold: Optional[float] = None
    ) -> Dict[int, List[int]]:
        """
        Extract representative structure frame indices for each basin.

        Args:
            cv1_data: CV1 time series
            cv2_data: CV2 time series
            minima: List of (cv1, cv2, energy) minima
            distance_threshold: Max distance from minimum in CV space
                               (auto-calculated if None)

        Returns:
            Dictionary {basin_id: [frame_idx1, frame_idx2, ...]}
        """
        # Auto-calculate threshold if not provided
        if distance_threshold is None:
            # Use 5% of CV range as threshold
            cv1_range = cv1_data.max() - cv1_data.min()
            cv2_range = cv2_data.max() - cv2_data.min()
            distance_threshold = 0.05 * np.sqrt(cv1_range**2 + cv2_range**2)

        representatives = {}

        for basin_id, (cv1_min, cv2_min, _) in enumerate(minima):
            # Calculate distances from minimum
            distances = np.sqrt(
                (cv1_data - cv1_min)**2 + (cv2_data - cv2_min)**2
            )

            # Find frames within threshold
            in_basin = distances < distance_threshold
            basin_frames = np.where(in_basin)[0]

            if len(basin_frames) > 0:
                # Select frame closest to minimum
                closest_idx = basin_frames[np.argmin(distances[basin_frames])]
                representatives[basin_id] = [int(closest_idx)]
            else:
                # No frames in basin (rare), use closest frame anyway
                closest_idx = np.argmin(distances)
                representatives[basin_id] = [int(closest_idx)]

        return representatives
