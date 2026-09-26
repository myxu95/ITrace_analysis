"""
FEL Visualizer Module

Visualization tools for free energy landscapes.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
from typing import Optional, Tuple, List
import os


class FELVisualizer:
    """
    Free Energy Landscape Visualizer

    Create publication-quality visualizations of FEL data.

    Attributes:
        style (str): Matplotlib style
        figsize (tuple): Default figure size
        dpi (int): Figure resolution
    """

    def __init__(
        self,
        style: str = "seaborn-v0_8",
        figsize: Tuple[int, int] = (10, 8),
        dpi: int = 300
    ):
        """
        Initialize visualizer.

        Args:
            style: Matplotlib style
            figsize: Default figure size (width, height)
            dpi: Resolution for saved figures
        """
        self.style = style
        self.figsize = figsize
        self.dpi = dpi

        try:
            plt.style.use(style)
        except:
            # Fallback to default
            plt.style.use('default')

    def plot_2d_contour(
        self,
        free_energy: np.ndarray,
        cv1_centers: np.ndarray,
        cv2_centers: np.ndarray,
        cv1_label: str = "CV1",
        cv2_label: str = "CV2",
        title: str = "Free Energy Landscape",
        levels: Optional[np.ndarray] = None,
        vmax: float = 20.0,
        cmap: str = 'viridis',
        output_file: Optional[str] = None,
        annotate_minima: bool = True,
        minima: Optional[List[Tuple]] = None,
        figsize: Optional[Tuple[int, int]] = None,
        show_colorbar: bool = True,
        contour_linewidths: float = 0.5,
        filled: bool = True
    ) -> None:
        """
        Plot 2D free energy landscape as contour map.

        Args:
            free_energy: 2D free energy array (n_cv1, n_cv2)
            cv1_centers: CV1 bin centers
            cv2_centers: CV2 bin centers
            cv1_label: Label for CV1 axis
            cv2_label: Label for CV2 axis
            title: Plot title
            levels: Contour levels (kJ/mol), auto if None
            vmax: Maximum energy value for colormap (kJ/mol)
            cmap: Colormap name
            output_file: Path to save figure (None = display only)
            annotate_minima: Whether to mark energy minima
            minima: List of (cv1, cv2, energy) minima
            figsize: Figure size (uses default if None)
            show_colorbar: Show colorbar
            contour_linewidths: Width of contour lines
            filled: Use filled contours (contourf) vs line contours
        """
        if figsize is None:
            figsize = self.figsize

        fig, ax = plt.subplots(figsize=figsize)

        # Prepare mesh grids for contour plotting
        CV1, CV2 = np.meshgrid(cv1_centers, cv2_centers, indexing='ij')

        # Clip free energy to vmax
        fe_clipped = np.clip(free_energy, 0, vmax)

        # Define contour levels if not provided
        if levels is None:
            # Default: [0, 2, 4, ..., 20] kJ/mol
            levels = np.arange(0, vmax + 2, 2)

        # Plot filled contours or line contours
        if filled:
            contour = ax.contourf(
                CV1, CV2, fe_clipped,
                levels=levels,
                cmap=cmap,
                extend='max'
            )
            # Add line contours on top for clarity
            ax.contour(
                CV1, CV2, fe_clipped,
                levels=levels,
                colors='black',
                linewidths=contour_linewidths,
                alpha=0.3
            )
        else:
            contour = ax.contour(
                CV1, CV2, fe_clipped,
                levels=levels,
                cmap=cmap,
                linewidths=contour_linewidths
            )

        # Colorbar
        if show_colorbar:
            cbar = plt.colorbar(contour, ax=ax)
            cbar.set_label('Free Energy (kJ/mol)', fontsize=12)

        # Annotate minima
        if annotate_minima and minima is not None and len(minima) > 0:
            for i, (cv1_min, cv2_min, energy_min) in enumerate(minima):
                # Plot marker
                ax.plot(cv1_min, cv2_min, 'r*', markersize=15,
                       markeredgecolor='white', markeredgewidth=1.5,
                       label=f'Min {i+1}' if i == 0 else '')

                # Annotate with energy value
                ax.annotate(
                    f'{energy_min:.1f}',
                    xy=(cv1_min, cv2_min),
                    xytext=(5, 5),
                    textcoords='offset points',
                    fontsize=10,
                    color='red',
                    fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='red', alpha=0.8)
                )

            if len(minima) > 0:
                ax.legend(loc='best', fontsize=10)

        # Labels and title
        ax.set_xlabel(cv1_label, fontsize=14, fontweight='bold')
        ax.set_ylabel(cv2_label, fontsize=14, fontweight='bold')
        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)

        # Grid
        ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

        # Tight layout
        plt.tight_layout()

        # Save or show
        if output_file:
            # Create output directory if needed
            output_dir = os.path.dirname(output_file)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)

            plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
            print(f"FEL contour plot saved to: {output_file}")
        else:
            plt.show()

        plt.close()

    def plot_1d_projections(
        self,
        cv1_data: np.ndarray,
        cv2_data: np.ndarray,
        fel_calculator: 'FELCalculator',
        cv1_label: str = "CV1",
        cv2_label: str = "CV2",
        output_file: Optional[str] = None,
        figsize: Tuple[int, int] = (12, 5)
    ) -> None:
        """
        Plot 1D free energy projections for both CVs.

        Args:
            cv1_data: CV1 time series
            cv2_data: CV2 time series
            fel_calculator: FELCalculator instance
            cv1_label: CV1 axis label
            cv2_label: CV2 axis label
            output_file: Path to save figure
            figsize: Figure size
        """
        fig, axes = plt.subplots(1, 2, figsize=figsize)

        # CV1 projection
        cv1_centers, fe_cv1 = fel_calculator.compute_1d_projection(cv1_data)
        axes[0].plot(cv1_centers, fe_cv1, 'b-', linewidth=2)
        axes[0].fill_between(cv1_centers, fe_cv1, alpha=0.3)
        axes[0].set_xlabel(cv1_label, fontsize=12, fontweight='bold')
        axes[0].set_ylabel('Free Energy (kJ/mol)', fontsize=12, fontweight='bold')
        axes[0].set_title(f'{cv1_label} Projection', fontsize=14, fontweight='bold')
        axes[0].grid(True, alpha=0.3)

        # CV2 projection
        cv2_centers, fe_cv2 = fel_calculator.compute_1d_projection(cv2_data)
        axes[1].plot(cv2_centers, fe_cv2, 'r-', linewidth=2)
        axes[1].fill_between(cv2_centers, fe_cv2, alpha=0.3, color='red')
        axes[1].set_xlabel(cv2_label, fontsize=12, fontweight='bold')
        axes[1].set_ylabel('Free Energy (kJ/mol)', fontsize=12, fontweight='bold')
        axes[1].set_title(f'{cv2_label} Projection', fontsize=14, fontweight='bold')
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()

        if output_file:
            output_dir = os.path.dirname(output_file)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)

            plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
            print(f"1D projections saved to: {output_file}")
        else:
            plt.show()

        plt.close()

    def plot_combined_figure(
        self,
        free_energy: np.ndarray,
        cv1_centers: np.ndarray,
        cv2_centers: np.ndarray,
        cv1_data: np.ndarray,
        cv2_data: np.ndarray,
        fel_calculator: 'FELCalculator',
        cv1_label: str = "CV1",
        cv2_label: str = "CV2",
        title: str = "Free Energy Landscape Analysis",
        output_file: Optional[str] = None,
        minima: Optional[List[Tuple]] = None
    ) -> None:
        """
        Create combined 4-panel figure:
        - Top left: CV1 1D projection
        - Top right: 2D FEL contour
        - Bottom left: CV2 1D projection
        - Bottom right: CV time series scatter

        Args:
            free_energy: 2D free energy array
            cv1_centers: CV1 bin centers
            cv2_centers: CV2 bin centers
            cv1_data: CV1 time series
            cv2_data: CV2 time series
            fel_calculator: FELCalculator instance
            cv1_label: CV1 label
            cv2_label: CV2 label
            title: Overall title
            output_file: Save path
            minima: Energy minima for annotation
        """
        fig = plt.figure(figsize=(14, 12))
        gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

        # Top left: CV1 projection
        ax1 = fig.add_subplot(gs[0, 0])
        cv1_proj_centers, fe_cv1 = fel_calculator.compute_1d_projection(cv1_data)
        ax1.plot(cv1_proj_centers, fe_cv1, 'b-', linewidth=2)
        ax1.fill_between(cv1_proj_centers, fe_cv1, alpha=0.3)
        ax1.set_ylabel('Free Energy (kJ/mol)', fontsize=11, fontweight='bold')
        ax1.set_title(f'{cv1_label} Projection', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3)

        # Top right: 2D FEL
        ax2 = fig.add_subplot(gs[0, 1])
        CV1, CV2 = np.meshgrid(cv1_centers, cv2_centers, indexing='ij')
        fe_clipped = np.clip(free_energy, 0, 20)
        levels = np.arange(0, 22, 2)
        contour = ax2.contourf(CV1, CV2, fe_clipped, levels=levels, cmap='viridis')
        ax2.contour(CV1, CV2, fe_clipped, levels=levels, colors='black',
                   linewidths=0.5, alpha=0.3)

        if minima is not None and len(minima) > 0:
            for cv1_min, cv2_min, e_min in minima:
                ax2.plot(cv1_min, cv2_min, 'r*', markersize=12,
                        markeredgecolor='white', markeredgewidth=1)

        cbar = plt.colorbar(contour, ax=ax2)
        cbar.set_label('Free Energy (kJ/mol)', fontsize=10)
        ax2.set_ylabel(cv2_label, fontsize=11, fontweight='bold')
        ax2.set_title('2D FEL', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3)

        # Bottom left: CV2 projection
        ax3 = fig.add_subplot(gs[1, 0])
        cv2_proj_centers, fe_cv2 = fel_calculator.compute_1d_projection(cv2_data)
        ax3.plot(cv2_proj_centers, fe_cv2, 'r-', linewidth=2)
        ax3.fill_between(cv2_proj_centers, fe_cv2, alpha=0.3, color='red')
        ax3.set_xlabel(cv2_label, fontsize=11, fontweight='bold')
        ax3.set_ylabel('Free Energy (kJ/mol)', fontsize=11, fontweight='bold')
        ax3.set_title(f'{cv2_label} Projection', fontsize=12, fontweight='bold')
        ax3.grid(True, alpha=0.3)

        # Bottom right: CV scatter
        ax4 = fig.add_subplot(gs[1, 1])
        scatter = ax4.scatter(cv1_data, cv2_data, c=range(len(cv1_data)),
                             cmap='coolwarm', s=1, alpha=0.5)
        ax4.set_xlabel(cv1_label, fontsize=11, fontweight='bold')
        ax4.set_ylabel(cv2_label, fontsize=11, fontweight='bold')
        ax4.set_title('CV Time Series', fontsize=12, fontweight='bold')
        ax4.grid(True, alpha=0.3)

        cbar_scatter = plt.colorbar(scatter, ax=ax4)
        cbar_scatter.set_label('Frame Index', fontsize=10)

        fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)

        if output_file:
            output_dir = os.path.dirname(output_file)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)

            plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
            print(f"Combined figure saved to: {output_file}")
        else:
            plt.show()

        plt.close()
