"""Free energy landscape visualization for comparison analysis."""

from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from .fel_comparison import Dict


class FELVisualizer:
    """Visualizer for free energy landscapes and their comparisons."""

    def __init__(self, dpi: int = 300, figsize: Tuple[float, float] = (12, 4)):
        """Initialize visualizer.

        Args:
            dpi: Resolution for saved figures
            figsize: Figure size (width, height) in inches
        """
        self.dpi = dpi
        self.figsize = figsize

    def plot_fel_comparison(
        self,
        result: Dict,
        output_path: Path,
        title: Optional[str] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
    ) -> None:
        """Create side-by-side FEL comparison plot.

        Args:
            result: FEL comparison result
            output_path: Path to save the figure
            title: Optional overall title
            vmin: Minimum energy for color scale (kcal/mol)
            vmax: Maximum energy for color scale (kcal/mol)
        """
        fig, axes = plt.subplots(1, 3, figsize=self.figsize)

        # Determine color scale limits
        if vmin is None:
            vmin = min(result.fel_a.free_energy.min(), result.fel_b.free_energy.min())
        if vmax is None:
            vmax = max(result.fel_a.free_energy.max(), result.fel_b.free_energy.max())

        # Plot FEL A
        self._plot_single_fel(
            axes[0],
            result.fel_a.cv1_edges,
            result.fel_a.cv2_edges,
            result.fel_a.free_energy,
            result.fel_a.cv1_name,
            result.fel_a.cv2_name,
            f"System A: {result.system_a_name}",
            vmin,
            vmax,
        )

        # Plot FEL B
        self._plot_single_fel(
            axes[1],
            result.fel_b.cv1_edges,
            result.fel_b.cv2_edges,
            result.fel_b.free_energy,
            result.fel_b.cv1_name,
            result.fel_b.cv2_name,
            f"System B: {result.system_b_name}",
            vmin,
            vmax,
        )

        # Plot difference map
        if result.difference_map is not None:
            self._plot_difference_map(
                axes[2],
                result.aligned_cv1_edges,
                result.aligned_cv2_edges,
                result.difference_map,
                result.fel_a.cv1_name,
                result.fel_a.cv2_name,
                "Difference (B - A)",
            )
        else:
            axes[2].text(
                0.5, 0.5,
                "Grids not aligned\n(different CV ranges)",
                ha='center', va='center',
                transform=axes[2].transAxes,
            )
            axes[2].set_title("Difference (B - A)")

        if title:
            fig.suptitle(title, fontsize=14, fontweight='bold')

        plt.tight_layout()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()

    def _plot_single_fel(
        self,
        ax: plt.Axes,
        cv1_edges: np.ndarray,
        cv2_edges: np.ndarray,
        free_energy: np.ndarray,
        cv1_name: str,
        cv2_name: str,
        title: str,
        vmin: float,
        vmax: float,
    ) -> None:
        """Plot a single FEL."""
        # Create custom colormap (blue -> white -> red)
        colors = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0',
                  '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
        n_bins = 256
        cmap = LinearSegmentedColormap.from_list('custom', colors, N=n_bins)

        # Plot 2D histogram
        im = ax.pcolormesh(
            cv1_edges,
            cv2_edges,
            free_energy.T,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            shading='auto',
        )

        # Add contour lines
        cv1_centers = (cv1_edges[:-1] + cv1_edges[1:]) / 2
        cv2_centers = (cv2_edges[:-1] + cv2_edges[1:]) / 2
        contour_levels = np.arange(0, vmax, 1.0)  # Contours every 1 kcal/mol
        ax.contour(
            cv1_centers,
            cv2_centers,
            free_energy.T,
            levels=contour_levels,
            colors='black',
            linewidths=0.5,
            alpha=0.3,
        )

        ax.set_xlabel(cv1_name)
        ax.set_ylabel(cv2_name)
        ax.set_title(title)

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Free Energy (kcal/mol)')

    def _plot_difference_map(
        self,
        ax: plt.Axes,
        cv1_edges: np.ndarray,
        cv2_edges: np.ndarray,
        difference: np.ndarray,
        cv1_name: str,
        cv2_name: str,
        title: str,
    ) -> None:
        """Plot FEL difference map."""
        # Symmetric color scale around zero
        vmax = np.abs(difference).max()
        vmin = -vmax

        # Diverging colormap (blue = A lower, red = B lower)
        im = ax.pcolormesh(
            cv1_edges,
            cv2_edges,
            difference.T,
            cmap='RdBu_r',
            vmin=vmin,
            vmax=vmax,
            shading='auto',
        )

        # Add zero contour
        cv1_centers = (cv1_edges[:-1] + cv1_edges[1:]) / 2
        cv2_centers = (cv2_edges[:-1] + cv2_edges[1:]) / 2
        ax.contour(
            cv1_centers,
            cv2_centers,
            difference.T,
            levels=[0],
            colors='black',
            linewidths=1.5,
        )

        ax.set_xlabel(cv1_name)
        ax.set_ylabel(cv2_name)
        ax.set_title(title)

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('ΔG difference (kcal/mol)')

    def plot_basin_comparison(
        self,
        result: Dict,
        output_path: Path,
        title: Optional[str] = None,
    ) -> None:
        """Create basin comparison visualization.

        Args:
            result: FEL comparison result with basin information
            output_path: Path to save the figure
            title: Optional overall title
        """
        if not result.basin_comparison:
            raise ValueError("No basin comparison data available")

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        # Plot basin locations on FEL A
        self._plot_basins_on_fel(
            axes[0],
            result.fel_a.cv1_edges,
            result.fel_a.cv2_edges,
            result.fel_a.free_energy,
            result.fel_a.basins,
            result.fel_a.cv1_name,
            result.fel_a.cv2_name,
            f"System A: {result.system_a_name}",
        )

        # Plot basin locations on FEL B
        self._plot_basins_on_fel(
            axes[1],
            result.fel_b.cv1_edges,
            result.fel_b.cv2_edges,
            result.fel_b.free_energy,
            result.fel_b.basins,
            result.fel_b.cv1_name,
            result.fel_b.cv2_name,
            f"System B: {result.system_b_name}",
        )

        if title:
            fig.suptitle(title, fontsize=14, fontweight='bold')

        plt.tight_layout()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()

    def _plot_basins_on_fel(
        self,
        ax: plt.Axes,
        cv1_edges: np.ndarray,
        cv2_edges: np.ndarray,
        free_energy: np.ndarray,
        basins: list,
        cv1_name: str,
        cv2_name: str,
        title: str,
    ) -> None:
        """Plot basins overlaid on FEL."""
        # Plot FEL
        im = ax.pcolormesh(
            cv1_edges,
            cv2_edges,
            free_energy.T,
            cmap='viridis',
            shading='auto',
        )

        # Overlay basin markers
        for i, basin in enumerate(basins):
            ax.plot(
                basin['cv1_center'],
                basin['cv2_center'],
                'r*',
                markersize=15,
                markeredgecolor='white',
                markeredgewidth=1,
            )
            ax.text(
                basin['cv1_center'],
                basin['cv2_center'],
                f"{i+1}",
                color='white',
                fontsize=10,
                fontweight='bold',
                ha='center',
                va='center',
            )

        ax.set_xlabel(cv1_name)
        ax.set_ylabel(cv2_name)
        ax.set_title(title)

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Free Energy (kcal/mol)')
