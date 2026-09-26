"""
Free Energy Landscape Comparison

Compares FEL between two conditions (e.g., Standard MD vs Enhanced Sampling)
and computes overlap metrics, KL divergence, and difference maps.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, Tuple, Optional
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import entropy

from ..landscape.landscape_analyzer import LandscapeResult


class FELComparator:
    """
    Free Energy Landscape Comparator

    Compares two FEL surfaces and computes:
    - Difference maps (ΔG)
    - Overlap coefficients
    - KL divergence
    - Basin differences
    """

    def __init__(self, energy_cutoff: float = 10.0):
        """
        Initialize FEL comparator

        Args:
            energy_cutoff: Maximum energy to display (kcal/mol)
        """
        self.energy_cutoff = energy_cutoff

    def align_grids(
        self,
        fel_a: LandscapeResult,
        fel_b: LandscapeResult,
        bins: int = 50
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Align two FEL to a common grid

        Args:
            fel_a: First FEL result
            fel_b: Second FEL result
            bins: Number of bins for common grid

        Returns:
            (cv1_common, cv2_common, fel_a_aligned, fel_b_aligned, valid_mask)
        """
        # Determine common CV ranges
        cv1_min = max(fel_a.cv1_edges[0], fel_b.cv1_edges[0])
        cv1_max = min(fel_a.cv1_edges[-1], fel_b.cv1_edges[-1])
        cv2_min = max(fel_a.cv2_edges[0], fel_b.cv2_edges[0])
        cv2_max = min(fel_a.cv2_edges[-1], fel_b.cv2_edges[-1])

        # Create common grid
        cv1_common = np.linspace(cv1_min, cv1_max, bins)
        cv2_common = np.linspace(cv2_min, cv2_max, bins)

        # Interpolate FEL A onto common grid
        interp_a = RegularGridInterpolator(
            (fel_a.cv1_values, fel_a.cv2_values),
            fel_a.free_energy,
            method='linear',
            bounds_error=False,
            fill_value=np.nan
        )

        # Interpolate FEL B onto common grid
        interp_b = RegularGridInterpolator(
            (fel_b.cv1_values, fel_b.cv2_values),
            fel_b.free_energy,
            method='linear',
            bounds_error=False,
            fill_value=np.nan
        )

        # Evaluate on common grid
        cv1_grid, cv2_grid = np.meshgrid(cv1_common, cv2_common, indexing='ij')
        points = np.column_stack([cv1_grid.ravel(), cv2_grid.ravel()])

        fel_a_aligned = interp_a(points).reshape(cv1_grid.shape)
        fel_b_aligned = interp_b(points).reshape(cv2_grid.shape)

        # Create valid mask (both FEL have data)
        valid_mask = ~(np.isnan(fel_a_aligned) | np.isnan(fel_b_aligned))

        return cv1_common, cv2_common, fel_a_aligned, fel_b_aligned, valid_mask

    def compute_overlap_coefficient(
        self,
        fel_a: np.ndarray,
        fel_b: np.ndarray,
        valid_mask: np.ndarray
    ) -> float:
        """
        Compute overlap coefficient between two FEL

        OC = ∫ min(P_A, P_B) dCV1 dCV2

        Args:
            fel_a: First FEL (aligned)
            fel_b: Second FEL (aligned)
            valid_mask: Valid data mask

        Returns:
            Overlap coefficient [0, 1]
        """
        # Convert free energy to probability
        # P = exp(-F/kT), normalized
        kT = 0.001987204 * 300.0  # kcal/mol at 300K

        prob_a = np.exp(-fel_a / kT)
        prob_b = np.exp(-fel_b / kT)

        # Normalize probabilities
        prob_a_norm = prob_a / np.nansum(prob_a[valid_mask])
        prob_b_norm = prob_b / np.nansum(prob_b[valid_mask])

        # Compute overlap
        overlap = np.nansum(np.minimum(prob_a_norm, prob_b_norm)[valid_mask])

        return overlap

    def compute_kl_divergence(
        self,
        fel_a: np.ndarray,
        fel_b: np.ndarray,
        valid_mask: np.ndarray
    ) -> Tuple[float, float]:
        """
        Compute KL divergence between two FEL

        KL(A||B) = ∫ P_A log(P_A / P_B) dCV

        Args:
            fel_a: First FEL (aligned)
            fel_b: Second FEL (aligned)
            valid_mask: Valid data mask

        Returns:
            (KL(A||B), KL(B||A))
        """
        kT = 0.001987204 * 300.0  # kcal/mol at 300K

        # Convert to probability
        prob_a = np.exp(-fel_a / kT)
        prob_b = np.exp(-fel_b / kT)

        # Normalize
        prob_a_norm = prob_a / np.nansum(prob_a[valid_mask])
        prob_b_norm = prob_b / np.nansum(prob_b[valid_mask])

        # Flatten and filter valid points
        prob_a_flat = prob_a_norm[valid_mask].ravel()
        prob_b_flat = prob_b_norm[valid_mask].ravel()

        # Compute KL divergence
        kl_ab = entropy(prob_a_flat, prob_b_flat)
        kl_ba = entropy(prob_b_flat, prob_a_flat)

        return kl_ab, kl_ba

    def compute_difference(
        self,
        fel_a: LandscapeResult,
        fel_b: LandscapeResult
    ) -> Dict:
        """
        Compute FEL difference and metrics

        Args:
            fel_a: First FEL result
            fel_b: Second FEL result

        Returns:
            Dictionary with difference map and metrics
        """
        # Align to common grid
        cv1_common, cv2_common, fel_a_aligned, fel_b_aligned, valid_mask = \
            self.align_grids(fel_a, fel_b)

        # Compute difference (B - A)
        difference = fel_b_aligned - fel_a_aligned

        # Compute metrics
        overlap = self.compute_overlap_coefficient(fel_a_aligned, fel_b_aligned, valid_mask)
        kl_ab, kl_ba = self.compute_kl_divergence(fel_a_aligned, fel_b_aligned, valid_mask)

        # Identify new basins in B
        new_basins_in_b = []
        for basin in fel_b.basins:
            # Check if basin center is in low-energy region of A
            cv1_idx = np.argmin(np.abs(cv1_common - basin['cv1_center']))
            cv2_idx = np.argmin(np.abs(cv2_common - basin['cv2_center']))

            if valid_mask[cv1_idx, cv2_idx]:
                energy_in_a = fel_a_aligned[cv1_idx, cv2_idx]
                if energy_in_a > 3.0:  # High energy in A, low in B
                    new_basins_in_b.append(basin)

        # Identify lost basins from A
        lost_basins_from_a = []
        for basin in fel_a.basins:
            cv1_idx = np.argmin(np.abs(cv1_common - basin['cv1_center']))
            cv2_idx = np.argmin(np.abs(cv2_common - basin['cv2_center']))

            if valid_mask[cv1_idx, cv2_idx]:
                energy_in_b = fel_b_aligned[cv1_idx, cv2_idx]
                if energy_in_b > 3.0:  # Low energy in A, high in B
                    lost_basins_from_a.append(basin)

        return {
            'cv1_common': cv1_common,
            'cv2_common': cv2_common,
            'fel_a_aligned': fel_a_aligned,
            'fel_b_aligned': fel_b_aligned,
            'difference': difference,
            'valid_mask': valid_mask,
            'overlap_coefficient': overlap,
            'kl_divergence_ab': kl_ab,
            'kl_divergence_ba': kl_ba,
            'new_basins_in_b': new_basins_in_b,
            'lost_basins_from_a': lost_basins_from_a
        }

    def plot_fel_comparison(
        self,
        fel_a: LandscapeResult,
        fel_b: LandscapeResult,
        difference: Dict,
        output_path: Path,
        case_a_label: str = "Condition A",
        case_b_label: str = "Condition B",
        dpi: int = 300
    ) -> None:
        """
        Generate 3-panel FEL comparison figure

        Layout: [FEL A] [FEL B] [Difference]

        Args:
            fel_a: First FEL result
            fel_b: Second FEL result
            difference: Difference dictionary from compute_difference()
            output_path: Output file path
            case_a_label: Label for condition A
            case_b_label: Label for condition B
            dpi: Figure resolution
        """
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # Extract data
        cv1 = difference['cv1_common']
        cv2 = difference['cv2_common']
        fel_a_data = difference['fel_a_aligned']
        fel_b_data = difference['fel_b_aligned']
        diff_data = difference['difference']
        valid_mask = difference['valid_mask']

        # Apply energy cutoff
        fel_a_plot = np.where(fel_a_data > self.energy_cutoff, np.nan, fel_a_data)
        fel_b_plot = np.where(fel_b_data > self.energy_cutoff, np.nan, fel_b_data)

        # Panel 1: FEL A
        im1 = axes[0].contourf(
            cv1, cv2, fel_a_plot.T,
            levels=20, cmap='viridis', vmin=0, vmax=self.energy_cutoff
        )
        axes[0].contour(
            cv1, cv2, fel_a_plot.T,
            levels=10, colors='white', linewidths=0.5, alpha=0.3
        )
        # Mark basins
        for i, basin in enumerate(fel_a.basins[:3]):  # Top 3 basins
            axes[0].plot(
                basin['cv1_center'], basin['cv2_center'],
                'r*', markersize=15, markeredgecolor='white', markeredgewidth=1
            )
            axes[0].text(
                basin['cv1_center'], basin['cv2_center'],
                f"A{i+1}", color='white', fontsize=10, ha='center', va='bottom'
            )

        axes[0].set_xlabel(f"{fel_a.cv1_name.upper()} (Å)")
        axes[0].set_ylabel(f"{fel_a.cv2_name.replace('_', ' ').title()}")
        axes[0].set_title(f"{case_a_label}\n(n={fel_a.n_frames} frames)")
        plt.colorbar(im1, ax=axes[0], label="Free Energy (kcal/mol)")

        # Panel 2: FEL B
        im2 = axes[1].contourf(
            cv1, cv2, fel_b_plot.T,
            levels=20, cmap='viridis', vmin=0, vmax=self.energy_cutoff
        )
        axes[1].contour(
            cv1, cv2, fel_b_plot.T,
            levels=10, colors='white', linewidths=0.5, alpha=0.3
        )
        # Mark basins
        for i, basin in enumerate(fel_b.basins[:3]):
            axes[1].plot(
                basin['cv1_center'], basin['cv2_center'],
                'r*', markersize=15, markeredgecolor='white', markeredgewidth=1
            )
            axes[1].text(
                basin['cv1_center'], basin['cv2_center'],
                f"B{i+1}", color='white', fontsize=10, ha='center', va='bottom'
            )

        axes[1].set_xlabel(f"{fel_b.cv1_name.upper()} (Å)")
        axes[1].set_ylabel(f"{fel_b.cv2_name.replace('_', ' ').title()}")
        axes[1].set_title(f"{case_b_label}\n(n={fel_b.n_frames} frames)")
        plt.colorbar(im2, ax=axes[1], label="Free Energy (kcal/mol)")

        # Panel 3: Difference (B - A)
        diff_plot = np.where(~valid_mask, np.nan, diff_data)
        vmax_diff = min(5.0, np.nanmax(np.abs(diff_plot)))
        im3 = axes[2].contourf(
            cv1, cv2, diff_plot.T,
            levels=20, cmap='RdBu_r', vmin=-vmax_diff, vmax=vmax_diff
        )
        axes[2].contour(
            cv1, cv2, diff_plot.T,
            levels=10, colors='black', linewidths=0.5, alpha=0.3
        )

        axes[2].set_xlabel(f"{fel_a.cv1_name.upper()} (Å)")
        axes[2].set_ylabel(f"{fel_a.cv2_name.replace('_', ' ').title()}")
        axes[2].set_title(
            f"Difference (B - A)\n"
            f"Overlap: {difference['overlap_coefficient']:.3f}"
        )
        plt.colorbar(im3, ax=axes[2], label="ΔG (kcal/mol)")

        plt.tight_layout()
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        plt.close()

        print(f"FEL comparison plot saved to {output_path}")
