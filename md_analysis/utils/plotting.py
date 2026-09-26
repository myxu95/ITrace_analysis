import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union

# Optional plotly import
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


class PlotManager:
    """Plotting utility for visualizing MD analysis results."""
    
    def __init__(self, style: str = "seaborn-v0_8", figsize: Tuple[int, int] = (10, 6)):
        """
        Initialize PlotManager.
        
        Args:
            style: Matplotlib style
            figsize: Default figure size
        """
        plt.style.use(style)
        self.figsize = figsize
        self.colors = sns.color_palette("husl", 10)
        
    def read_xvg_file(self, file_path: str) -> pd.DataFrame:
        """
        Read GROMACS XVG file and return as DataFrame.
        
        Args:
            file_path: Path to XVG file
            
        Returns:
            DataFrame with time series data
        """
        data = []
        headers = []
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('@') and 'legend' in line:
                    legend = line.split('"')[1]
                    headers.append(legend)
                elif line.startswith('#') or line.startswith('@'):
                    continue
                else:
                    values = [float(x) for x in line.split()]
                    data.append(values)
        
        df = pd.DataFrame(data)
        if headers and len(headers) == len(df.columns) - 1:
            df.columns = ['Time'] + headers
        elif len(df.columns) == 2:
            df.columns = ['Time', 'Value']
        
        return df
    
    def plot_time_series(self, 
                        data: Union[str, pd.DataFrame],
                        title: str = "Time Series",
                        xlabel: str = "Time (ps)",
                        ylabel: str = "Value",
                        output_path: Optional[str] = None,
                        interactive: bool = False) -> None:
        """
        Plot time series data.
        
        Args:
            data: XVG file path or DataFrame
            title: Plot title
            xlabel: X-axis label
            ylabel: Y-axis label
            output_path: Path to save plot
            interactive: Use plotly for interactive plot
        """
        if isinstance(data, str):
            df = self.read_xvg_file(data)
        else:
            df = data
            
        if interactive and PLOTLY_AVAILABLE:
            self._plot_interactive_time_series(df, title, xlabel, ylabel, output_path)
        elif interactive and not PLOTLY_AVAILABLE:
            print("⚠️ Plotly not available, falling back to static plot")
            self._plot_static_time_series(df, title, xlabel, ylabel, output_path)
        else:
            self._plot_static_time_series(df, title, xlabel, ylabel, output_path)
    
    def _plot_static_time_series(self, df: pd.DataFrame, title: str, 
                               xlabel: str, ylabel: str, output_path: Optional[str]):
        """Create static time series plot with matplotlib."""
        fig, ax = plt.subplots(figsize=self.figsize)
        
        time_col = df.columns[0]
        for i, col in enumerate(df.columns[1:]):
            ax.plot(df[time_col], df[col], label=col, color=self.colors[i % len(self.colors)])
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        
        if len(df.columns) > 2:
            ax.legend()
            
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            
        plt.show()
    
    def _plot_interactive_time_series(self, df: pd.DataFrame, title: str,
                                    xlabel: str, ylabel: str, output_path: Optional[str]):
        """Create interactive time series plot with plotly."""
        if not PLOTLY_AVAILABLE:
            raise ImportError("Plotly is required for interactive plots. Install with: pip install plotly")
            
        fig = go.Figure()
        
        time_col = df.columns[0]
        for col in df.columns[1:]:
            fig.add_trace(go.Scatter(
                x=df[time_col],
                y=df[col],
                mode='lines',
                name=col
            ))
        
        fig.update_layout(
            title=title,
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            hovermode='x unified'
        )
        
        if output_path:
            if output_path.endswith('.html'):
                fig.write_html(output_path)
            else:
                fig.write_image(output_path)
        
        fig.show()
    
    def plot_distribution(self, 
                         data: Union[List, np.ndarray, pd.Series],
                         title: str = "Distribution",
                         xlabel: str = "Value",
                         bins: int = 50,
                         output_path: Optional[str] = None) -> None:
        """
        Plot distribution histogram.
        
        Args:
            data: Data to plot
            title: Plot title
            xlabel: X-axis label
            bins: Number of histogram bins
            output_path: Path to save plot
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        ax.hist(data, bins=bins, alpha=0.7, density=True, color=self.colors[0])
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Density")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            
        plt.show()
    
    def plot_heatmap(self, 
                    data: Union[np.ndarray, pd.DataFrame],
                    title: str = "Heatmap",
                    xlabel: str = "X",
                    ylabel: str = "Y",
                    output_path: Optional[str] = None) -> None:
        """
        Plot 2D heatmap.
        
        Args:
            data: 2D data array or DataFrame
            title: Plot title
            xlabel: X-axis label
            ylabel: Y-axis label
            output_path: Path to save plot
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        if isinstance(data, pd.DataFrame):
            sns.heatmap(data, ax=ax, cmap='viridis')
        else:
            im = ax.imshow(data, cmap='viridis', aspect='auto')
            plt.colorbar(im, ax=ax)
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            
        plt.show()
    
    def plot_multiple_xvg(self, 
                         file_paths: List[str],
                         labels: Optional[List[str]] = None,
                         title: str = "Multiple Time Series",
                         output_path: Optional[str] = None) -> None:
        """
        Plot multiple XVG files on the same plot.
        
        Args:
            file_paths: List of XVG file paths
            labels: Optional labels for each file
            title: Plot title
            output_path: Path to save plot
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        if labels is None:
            labels = [Path(fp).stem for fp in file_paths]
        
        for i, (file_path, label) in enumerate(zip(file_paths, labels)):
            df = self.read_xvg_file(file_path)
            time_col = df.columns[0]
            value_col = df.columns[1]
            
            ax.plot(df[time_col], df[value_col], 
                   label=label, color=self.colors[i % len(self.colors)])
        
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel("Value")
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')

        plt.show()

    def plot_rmsd(self,
                  rmsd_file: Union[str, List[str]],
                  title: Optional[str] = None,
                  xlabel: str = "Time (ps)",
                  ylabel: str = "RMSD (nm)",
                  output_path: Optional[str] = None,
                  show_stats: bool = True,
                  moving_average: Optional[int] = None,
                  interactive: bool = False) -> None:
        """
        Plot RMSD data from XVG file(s) with enhanced visualization.

        Args:
            rmsd_file: Path to RMSD XVG file or list of paths for comparison
            title: Plot title (auto-generated if None)
            xlabel: X-axis label
            ylabel: Y-axis label
            output_path: Output file path for saving plot
            show_stats: Show statistics (mean, std) on plot
            moving_average: Window size for moving average (optional)
            interactive: Use interactive plotly plot
        """
        if isinstance(rmsd_file, str):
            rmsd_files = [rmsd_file]
        else:
            rmsd_files = rmsd_file

        if interactive and PLOTLY_AVAILABLE:
            self._plot_rmsd_interactive(rmsd_files, title, xlabel, ylabel,
                                      output_path, show_stats, moving_average)
        else:
            self._plot_rmsd_static(rmsd_files, title, xlabel, ylabel,
                                 output_path, show_stats, moving_average)

    def _plot_rmsd_static(self, rmsd_files: List[str], title: Optional[str],
                         xlabel: str, ylabel: str, output_path: Optional[str],
                         show_stats: bool, moving_average: Optional[int]) -> None:
        """Static RMSD plot using matplotlib."""
        fig, ax = plt.subplots(figsize=self.figsize)

        for i, file_path in enumerate(rmsd_files):
            # Read RMSD data
            df = self.read_xvg_file(file_path)
            if df.empty:
                continue

            time_col = df.columns[0]
            rmsd_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]

            time = df[time_col].values
            rmsd = df[rmsd_col].values

            # Convert time from ps to ns for better readability
            time_ns = time / 1000

            file_name = Path(file_path).stem
            color = self.colors[i % len(self.colors)]

            # Plot main RMSD curve
            ax.plot(time_ns, rmsd, color=color, linewidth=1.5,
                   label=file_name, alpha=0.8)

            # Add moving average if requested
            if moving_average and len(rmsd) >= moving_average:
                rmsd_ma = pd.Series(rmsd).rolling(window=moving_average, center=True).mean()
                ax.plot(time_ns, rmsd_ma, color=color, linewidth=2.5,
                       linestyle='--', alpha=0.9,
                       label=f'{file_name} (MA {moving_average})')

            # Add statistics if requested
            if show_stats:
                mean_rmsd = np.mean(rmsd)
                std_rmsd = np.std(rmsd)

                # Add horizontal line for mean
                ax.axhline(y=mean_rmsd, color=color, linestyle=':',
                          alpha=0.6, linewidth=1)

                # Add text with statistics
                stats_text = f'{file_name}: μ={mean_rmsd:.2f}±{std_rmsd:.2f} nm'
                ax.text(0.02, 0.98 - i*0.05, stats_text,
                       transform=ax.transAxes, fontsize=9,
                       bbox=dict(boxstyle="round,pad=0.3", facecolor=color, alpha=0.2))

        # Customize plot
        ax.set_xlabel(f"{xlabel.replace('(ps)', '(ns)')}")
        ax.set_ylabel(ylabel)

        if title is None:
            if len(rmsd_files) == 1:
                title = f"RMSD Analysis - {Path(rmsd_files[0]).stem}"
            else:
                title = f"RMSD Comparison ({len(rmsd_files)} trajectories)"
        ax.set_title(title, fontsize=14, fontweight='bold')

        # Add grid and legend
        ax.grid(True, alpha=0.3)
        if len(rmsd_files) > 1 or moving_average:
            ax.legend(loc='best', framealpha=0.9)

        # Set y-axis to start from 0
        ax.set_ylim(bottom=0)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"RMSD plot saved to: {output_path}")

        plt.show()

    def _plot_rmsd_interactive(self, rmsd_files: List[str], title: Optional[str],
                              xlabel: str, ylabel: str, output_path: Optional[str],
                              show_stats: bool, moving_average: Optional[int]) -> None:
        """Interactive RMSD plot using plotly."""
        fig = go.Figure()

        for i, file_path in enumerate(rmsd_files):
            # Read RMSD data
            df = self.read_xvg_file(file_path)
            if df.empty:
                continue

            time_col = df.columns[0]
            rmsd_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]

            time = df[time_col].values
            rmsd = df[rmsd_col].values

            # Convert time from ps to ns
            time_ns = time / 1000

            file_name = Path(file_path).stem
            color = px.colors.qualitative.Set1[i % len(px.colors.qualitative.Set1)]

            # Add main RMSD trace
            fig.add_trace(go.Scatter(
                x=time_ns, y=rmsd,
                mode='lines',
                name=file_name,
                line=dict(color=color, width=2),
                hovertemplate=f'<b>{file_name}</b><br>Time: %{{x:.1f}} ns<br>RMSD: %{{y:.3f}} nm<extra></extra>'
            ))

            # Add moving average if requested
            if moving_average and len(rmsd) >= moving_average:
                rmsd_ma = pd.Series(rmsd).rolling(window=moving_average, center=True).mean()
                fig.add_trace(go.Scatter(
                    x=time_ns, y=rmsd_ma,
                    mode='lines',
                    name=f'{file_name} (MA {moving_average})',
                    line=dict(color=color, width=3, dash='dash'),
                    hovertemplate=f'<b>{file_name} MA</b><br>Time: %{{x:.1f}} ns<br>RMSD: %{{y:.3f}} nm<extra></extra>'
                ))

            # Add mean line if requested
            if show_stats:
                mean_rmsd = np.mean(rmsd)
                fig.add_hline(y=mean_rmsd, line_dash="dot", line_color=color,
                             annotation_text=f'{file_name} mean: {mean_rmsd:.2f} nm',
                             annotation_position="top left")

        # Update layout
        if title is None:
            if len(rmsd_files) == 1:
                title = f"RMSD Analysis - {Path(rmsd_files[0]).stem}"
            else:
                title = f"RMSD Comparison ({len(rmsd_files)} trajectories)"

        fig.update_layout(
            title=title,
            xaxis_title=xlabel.replace('(ps)', '(ns)'),
            yaxis_title=ylabel,
            hovermode='x unified',
            template='plotly_white',
            width=1000,
            height=600,
            yaxis=dict(rangemode='tozero')  # Start y-axis from 0
        )

        if output_path:
            if output_path.endswith('.html'):
                fig.write_html(output_path)
            else:
                fig.write_image(output_path, width=1000, height=600)
            print(f"Interactive RMSD plot saved to: {output_path}")

        fig.show()

    def plot_rmsd_distribution(self,
                              rmsd_file: Union[str, List[str]],
                              title: Optional[str] = None,
                              output_path: Optional[str] = None,
                              bins: int = 50,
                              show_kde: bool = True) -> None:
        """
        Plot RMSD distribution histogram with KDE.

        Args:
            rmsd_file: Path to RMSD XVG file or list of paths
            title: Plot title
            output_path: Output file path
            bins: Number of histogram bins
            show_kde: Show kernel density estimation
        """
        if isinstance(rmsd_file, str):
            rmsd_files = [rmsd_file]
        else:
            rmsd_files = rmsd_file

        fig, ax = plt.subplots(figsize=self.figsize)

        for i, file_path in enumerate(rmsd_files):
            # Read RMSD data
            df = self.read_xvg_file(file_path)
            if df.empty:
                continue

            rmsd_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]
            rmsd = df[rmsd_col].values

            file_name = Path(file_path).stem
            color = self.colors[i % len(self.colors)]

            # Plot histogram
            ax.hist(rmsd, bins=bins, alpha=0.6, color=color,
                   label=f'{file_name} (n={len(rmsd)})', density=True)

            # Add KDE if requested
            if show_kde:
                try:
                    from scipy.stats import gaussian_kde
                    kde = gaussian_kde(rmsd)
                    x_range = np.linspace(rmsd.min(), rmsd.max(), 200)
                    ax.plot(x_range, kde(x_range), color=color, linewidth=2)
                except ImportError:
                    print("Warning: scipy not available, skipping KDE")

            # Add statistics
            mean_rmsd = np.mean(rmsd)
            ax.axvline(mean_rmsd, color=color, linestyle='--', linewidth=2,
                      label=f'{file_name} mean: {mean_rmsd:.2f} nm')

        # Customize plot
        ax.set_xlabel("RMSD (nm)")
        ax.set_ylabel("Density")

        if title is None:
            if len(rmsd_files) == 1:
                title = f"RMSD Distribution - {Path(rmsd_files[0]).stem}"
            else:
                title = f"RMSD Distribution Comparison"
        ax.set_title(title, fontsize=14, fontweight='bold')

        ax.grid(True, alpha=0.3)
        ax.legend()

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"RMSD distribution plot saved to: {output_path}")

        plt.show()

    def plot_cdr3_rmsd_summary(self,
                              summary_csv: str,
                              output_dir: Optional[str] = None,
                              show_plots: bool = True) -> Dict[str, str]:
        """
        Generate comprehensive CDR3 RMSD summary visualizations.

        Args:
            summary_csv: Path to CDR3 RMSD summary CSV file
            output_dir: Directory to save plots (defaults to same dir as CSV)
            show_plots: Whether to display plots

        Returns:
            Dictionary mapping plot types to saved file paths
        """
        df = pd.read_csv(summary_csv)

        if output_dir is None:
            output_dir = Path(summary_csv).parent
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        output_paths = {}

        # 4-panel summary figure
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))

        # Panel 1: RMSD Distribution Histogram
        ax = axes[0, 0]
        ax.hist(df['rmsd_mean_nm'], bins=30, alpha=0.7, color='steelblue', edgecolor='black')
        mean_val = df['rmsd_mean_nm'].mean()
        median_val = df['rmsd_mean_nm'].median()
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {mean_val:.4f} nm')
        ax.axvline(median_val, color='green', linestyle='--', linewidth=2,
                  label=f'Median: {median_val:.4f} nm')
        ax.set_xlabel('Mean RMSD (nm)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
        ax.set_title('CDR3β RMSD Distribution', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Panel 2: Box plot comparison
        ax = axes[0, 1]
        box_data = [df['rmsd_mean_nm'], df['rmsd_std_nm'], df['rmsd_max_nm']]
        bp = ax.boxplot(box_data, labels=['Mean', 'Std', 'Max'], patch_artist=True)
        for patch, color in zip(bp['boxes'], ['lightblue', 'lightgreen', 'lightcoral']):
            patch.set_facecolor(color)
        ax.set_ylabel('RMSD (nm)', fontsize=12, fontweight='bold')
        ax.set_title('RMSD Statistics Comparison', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

        # Panel 3: Mean vs Std scatter plot
        ax = axes[1, 0]
        scatter = ax.scatter(df['rmsd_mean_nm'], df['rmsd_std_nm'],
                           c=df['rmsd_max_nm'], cmap='viridis',
                           s=100, alpha=0.6, edgecolors='black', linewidth=0.5)
        ax.set_xlabel('Mean RMSD (nm)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Std RMSD (nm)', fontsize=12, fontweight='bold')
        ax.set_title('Mean vs Std RMSD (colored by Max)', fontsize=14, fontweight='bold')
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Max RMSD (nm)', fontsize=10)
        ax.grid(True, alpha=0.3)

        # Panel 4: Top 20 highest RMSD tasks
        ax = axes[1, 1]
        top_20 = df.nlargest(20, 'rmsd_mean_nm')
        colors_gradient = plt.cm.coolwarm(np.linspace(0, 1, 20))
        ax.barh(range(20), top_20['rmsd_mean_nm'].values, color=colors_gradient)
        ax.set_yticks(range(20))
        ax.set_yticklabels(top_20['task_name'].values, fontsize=8)
        ax.set_xlabel('Mean RMSD (nm)', fontsize=12, fontweight='bold')
        ax.set_title('Top 20 Highest RMSD Tasks', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='x')
        ax.invert_yaxis()

        plt.tight_layout()

        # Save 4-panel figure
        summary_path = output_dir / "cdr3_beta_rmsd_summary.png"
        plt.savefig(summary_path, dpi=300, bbox_inches='tight')
        output_paths['summary'] = str(summary_path)

        if show_plots:
            plt.show()
        else:
            plt.close()

        print(f"CDR3 RMSD summary plot saved: {summary_path}")

        return output_paths

    def plot_cdr3_stability_distribution(self,
                                        summary_csv: str,
                                        output_dir: Optional[str] = None,
                                        show_plots: bool = True) -> Dict[str, str]:
        """
        Generate CDR3 stability distribution visualizations.

        Args:
            summary_csv: Path to CDR3 RMSD summary CSV file
            output_dir: Directory to save plots
            show_plots: Whether to display plots

        Returns:
            Dictionary mapping plot types to saved file paths
        """
        df = pd.read_csv(summary_csv)

        if output_dir is None:
            output_dir = Path(summary_csv).parent
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        # Categorize by RMSD
        def categorize_rmsd(rmsd):
            if rmsd < 0.10:
                return 'Very Stable (<0.10 nm)'
            elif rmsd < 0.15:
                return 'Stable (0.10-0.15 nm)'
            elif rmsd < 0.20:
                return 'Moderate (0.15-0.20 nm)'
            elif rmsd < 0.25:
                return 'Flexible (0.20-0.25 nm)'
            else:
                return 'Highly Flexible (>0.25 nm)'

        df['stability_category'] = df['rmsd_mean_nm'].apply(categorize_rmsd)
        category_counts = df['stability_category'].value_counts()

        # Create 2-panel figure
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        # Panel 1: Pie chart
        ax = axes[0]
        colors = ['#2ecc71', '#3498db', '#f39c12', '#e74c3c', '#9b59b6']
        wedges, texts, autotexts = ax.pie(category_counts, labels=category_counts.index,
                                          autopct='%1.1f%%', colors=colors, startangle=90,
                                          textprops={'fontsize': 10, 'fontweight': 'bold'})
        ax.set_title('CDR3β Stability Distribution', fontsize=14, fontweight='bold')

        # Panel 2: Bar chart
        ax = axes[1]
        category_order = ['Very Stable (<0.10 nm)', 'Stable (0.10-0.15 nm)',
                         'Moderate (0.15-0.20 nm)', 'Flexible (0.20-0.25 nm)',
                         'Highly Flexible (>0.25 nm)']
        counts = [category_counts.get(cat, 0) for cat in category_order]
        bars = ax.bar(range(len(category_order)), counts, color=colors)
        ax.set_xticks(range(len(category_order)))
        ax.set_xticklabels(category_order, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel('Number of Complexes', fontsize=12, fontweight='bold')
        ax.set_title('Stability Category Counts', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

        # Add count labels on bars
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(count)}', ha='center', va='bottom', fontweight='bold')

        plt.tight_layout()

        # Save figure
        stability_path = output_dir / "cdr3_beta_stability_categories.png"
        plt.savefig(stability_path, dpi=300, bbox_inches='tight')

        if show_plots:
            plt.show()
        else:
            plt.close()

        print(f"CDR3 stability distribution plot saved: {stability_path}")

        return {'stability': str(stability_path)}

    def plot_cdr3_rmsd_comparison(self,
                                 summary_csv: str,
                                 task_names: List[str],
                                 trajectory_dir: str,
                                 output_path: Optional[str] = None,
                                 show_stats: bool = True) -> None:
        """
        Compare CDR3 RMSD trajectories for multiple tasks.

        Args:
            summary_csv: Path to CDR3 RMSD summary CSV file
            task_names: List of task names to compare
            trajectory_dir: Base directory containing trajectory subdirectories
            output_path: Output file path
            show_stats: Show statistics on plot
        """
        df = pd.read_csv(summary_csv)
        trajectory_dir = Path(trajectory_dir)

        fig, ax = plt.subplots(figsize=(14, 8))

        for i, task_name in enumerate(task_names):
            # Get RMSD file path
            rmsd_file = trajectory_dir / task_name / "cdr3_beta_rmsd.xvg"

            if not rmsd_file.exists():
                print(f"Warning: RMSD file not found for {task_name}")
                continue

            # Read RMSD data
            rmsd_df = self.read_xvg_file(str(rmsd_file))
            if rmsd_df.empty:
                continue

            time = rmsd_df[rmsd_df.columns[0]].values / 1000  # Convert to ns
            rmsd = rmsd_df[rmsd_df.columns[1]].values

            color = self.colors[i % len(self.colors)]

            # Plot trajectory
            ax.plot(time, rmsd, color=color, linewidth=1.5, alpha=0.8, label=task_name)

            # Add statistics if requested
            if show_stats:
                task_stats = df[df['task_name'] == task_name]
                if not task_stats.empty:
                    mean_rmsd = task_stats['rmsd_mean_nm'].values[0]
                    ax.axhline(y=mean_rmsd, color=color, linestyle=':', alpha=0.5, linewidth=1)

        ax.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
        ax.set_ylabel('CDR3β RMSD (nm)', fontsize=12, fontweight='bold')
        ax.set_title(f'CDR3β RMSD Comparison ({len(task_names)} tasks)',
                    fontsize=14, fontweight='bold')
        ax.legend(loc='best', framealpha=0.9, fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"CDR3 RMSD comparison plot saved: {output_path}")

        plt.show()

    def plot_cdr3_rmsd_heatmap(self,
                              summary_csv: str,
                              output_path: Optional[str] = None,
                              top_n: Optional[int] = None,
                              sort_by: str = 'rmsd_mean_nm') -> None:
        """
        Generate heatmap visualization of CDR3 RMSD statistics across tasks.

        Args:
            summary_csv: Path to CDR3 RMSD summary CSV file
            output_path: Output file path
            top_n: Show only top N tasks (by sort_by metric)
            sort_by: Column to sort by ('rmsd_mean_nm', 'rmsd_std_nm', etc.)
        """
        df = pd.read_csv(summary_csv)

        # Select columns for heatmap
        heatmap_cols = ['rmsd_mean_nm', 'rmsd_std_nm', 'rmsd_min_nm',
                       'rmsd_max_nm', 'rmsd_median_nm', 'rmsd_range_nm']

        # Filter columns that exist
        heatmap_cols = [col for col in heatmap_cols if col in df.columns]

        # Sort and optionally limit to top N
        df_sorted = df.sort_values(by=sort_by, ascending=False)
        if top_n:
            df_sorted = df_sorted.head(top_n)

        # Prepare data for heatmap
        heatmap_data = df_sorted[heatmap_cols].values
        task_labels = df_sorted['task_name'].values

        # Create figure
        fig, ax = plt.subplots(figsize=(12, max(8, len(task_labels) * 0.3)))

        # Create heatmap
        im = ax.imshow(heatmap_data, cmap='YlOrRd', aspect='auto')

        # Set ticks and labels
        ax.set_xticks(range(len(heatmap_cols)))
        ax.set_xticklabels([col.replace('rmsd_', '').replace('_nm', '')
                           for col in heatmap_cols], rotation=45, ha='right')
        ax.set_yticks(range(len(task_labels)))
        ax.set_yticklabels(task_labels, fontsize=8)

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('RMSD (nm)', fontsize=10, fontweight='bold')

        # Add title
        title = f'CDR3β RMSD Heatmap'
        if top_n:
            title += f' (Top {top_n} by {sort_by})'
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

        # Add values as text
        for i in range(len(task_labels)):
            for j in range(len(heatmap_cols)):
                text = ax.text(j, i, f'{heatmap_data[i, j]:.3f}',
                             ha="center", va="center", color="black", fontsize=7)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"CDR3 RMSD heatmap saved: {output_path}")

        plt.show()

    def plot_rmsd_convergence(self,
                             rmsd_file: str,
                             window_sizes: List[int] = [100, 500, 1000, 2000],
                             title: Optional[str] = None,
                             output_path: Optional[str] = None) -> None:
        """
        Plot RMSD convergence analysis with different time windows.

        Args:
            rmsd_file: Path to RMSD XVG file
            window_sizes: List of window sizes for moving average
            title: Plot title
            output_path: Output file path
        """
        # Read RMSD data
        df = self.read_xvg_file(rmsd_file)
        if df.empty:
            print(f"Warning: Could not read data from {rmsd_file}")
            return

        time_col = df.columns[0]
        rmsd_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]

        time = df[time_col].values / 1000  # Convert to ns
        rmsd = df[rmsd_col].values

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(self.figsize[0], self.figsize[1]*1.5))

        # Plot 1: Raw RMSD with different moving averages
        ax1.plot(time, rmsd, color='lightgray', alpha=0.6, linewidth=0.8, label='Raw RMSD')

        for i, window in enumerate(window_sizes):
            if len(rmsd) >= window:
                rmsd_ma = pd.Series(rmsd).rolling(window=window, center=True).mean()
                color = self.colors[i % len(self.colors)]
                ax1.plot(time, rmsd_ma, color=color, linewidth=2,
                        label=f'MA {window} frames')

        ax1.set_xlabel("Time (ns)")
        ax1.set_ylabel("RMSD (nm)")
        ax1.set_title("RMSD with Moving Averages", fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        ax1.set_ylim(bottom=0)

        # Plot 2: Cumulative mean convergence
        cumulative_mean = np.cumsum(rmsd) / np.arange(1, len(rmsd) + 1)
        ax2.plot(time, cumulative_mean, color='darkblue', linewidth=2, label='Cumulative Mean')

        # Add horizontal line for final mean
        final_mean = np.mean(rmsd)
        ax2.axhline(final_mean, color='red', linestyle='--', linewidth=2,
                   label=f'Final Mean: {final_mean:.2f} nm')

        # Add convergence threshold (±5% of final mean)
        threshold = 0.05 * final_mean
        ax2.fill_between(time, final_mean - threshold, final_mean + threshold,
                        alpha=0.2, color='red', label='±5% threshold')

        ax2.set_xlabel("Time (ns)")
        ax2.set_ylabel("Cumulative Mean RMSD (nm)")
        ax2.set_title("RMSD Convergence Analysis", fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        ax2.set_ylim(bottom=0)

        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        else:
            file_name = Path(rmsd_file).stem
            fig.suptitle(f"RMSD Convergence Analysis - {file_name}",
                        fontsize=16, fontweight='bold')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"RMSD convergence plot saved to: {output_path}")

        plt.show()

    @staticmethod
    def plot_cdr_rmsf_curves(rmsf_output_dir: str,
                             output_path: Optional[str] = None,
                             show_plot: bool = True,
                             n_points: int = 100,
                             figsize: Tuple[int, int] = (14, 8)) -> None:
        """
        Generate CDR RMSF comparison using normalized position curves.

        Creates overlaid RMSF curves for 6 CDR loops (Alpha/Beta CDR1/2/3) with
        normalized relative positions to handle variable CDR lengths. Shows mean
        curve with standard deviation shaded region.

        Args:
            rmsf_output_dir: Directory containing task subdirectories with cdr_rmsf_extracted.json
            output_path: Output file path for saving plot (optional)
            show_plot: Whether to display the plot interactively
            n_points: Number of interpolation points for normalized position (default: 100)
            figsize: Figure size (width, height) in inches

        Returns:
            None

        Example:
            >>> PlotManager.plot_cdr_rmsf_curves(
            ...     rmsf_output_dir="output/rmsf_cdr",
            ...     output_path="output/rmsf_cdr/cdr_rmsf_curves_comparison.png"
            ... )
        """
        import json
        from scipy.interpolate import interp1d

        rmsf_dir = Path(rmsf_output_dir)
        if not rmsf_dir.exists():
            raise FileNotFoundError(f"RMSF output directory not found: {rmsf_output_dir}")

        # CDR mapping: JSON keys to display names
        cdr_mapping = {
            'alpha_cdr1': 'Alpha CDR1',
            'alpha_cdr2': 'Alpha CDR2',
            'alpha_cdr3': 'Alpha CDR3',
            'beta_cdr1': 'Beta CDR1',
            'beta_cdr2': 'Beta CDR2',
            'beta_cdr3': 'Beta CDR3'
        }

        # Collect RMSF curves for each CDR
        cdr_curves = {key: [] for key in cdr_mapping.keys()}

        # Scan all task directories
        task_dirs = [d for d in rmsf_dir.iterdir() if d.is_dir()]
        print(f"Scanning {len(task_dirs)} task directories...")

        for task_dir in task_dirs:
            json_file = task_dir / "cdr_rmsf_extracted.json"
            if not json_file.exists():
                continue

            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)

                # Extract RMSF values for each CDR
                for cdr_key in cdr_mapping.keys():
                    if cdr_key in data and 'rmsf_values' in data[cdr_key]:
                        rmsf_values = data[cdr_key]['rmsf_values']
                        if len(rmsf_values) > 0:
                            cdr_curves[cdr_key].append(rmsf_values)
            except Exception as e:
                print(f"Warning: Failed to read {json_file}: {e}")
                continue

        # Check if we have data
        total_curves = sum(len(curves) for curves in cdr_curves.values())
        if total_curves == 0:
            raise ValueError("No valid RMSF curves found in any JSON files")

        print(f"Loaded RMSF curves: {total_curves} total")
        for cdr_key, curves in cdr_curves.items():
            if len(curves) > 0:
                print(f"  {cdr_mapping[cdr_key]}: {len(curves)} samples")

        # Normalize and interpolate curves to common grid
        normalized_grid = np.linspace(0, 1, n_points)
        interpolated_curves = {key: [] for key in cdr_mapping.keys()}

        for cdr_key, curves in cdr_curves.items():
            for curve in curves:
                if len(curve) < 2:
                    continue  # Skip too short curves

                # Create normalized position for original curve
                original_positions = np.linspace(0, 1, len(curve))

                # Interpolate to common grid
                try:
                    interp_func = interp1d(original_positions, curve,
                                          kind='linear', bounds_error=False,
                                          fill_value='extrapolate')
                    interpolated_curve = interp_func(normalized_grid)
                    interpolated_curves[cdr_key].append(interpolated_curve)
                except Exception as e:
                    print(f"Warning: Interpolation failed for {cdr_key}: {e}")
                    continue

        # Calculate mean and std for each CDR
        cdr_stats = {}
        for cdr_key in cdr_mapping.keys():
            if len(interpolated_curves[cdr_key]) > 0:
                curves_array = np.array(interpolated_curves[cdr_key])
                cdr_stats[cdr_key] = {
                    'mean': np.mean(curves_array, axis=0),
                    'std': np.std(curves_array, axis=0),
                    'n_samples': len(interpolated_curves[cdr_key])
                }

        # Set style - academic journal style
        plt.style.use('seaborn-v0_8-white')
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.size'] = 11
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['axes.edgecolor'] = '#333333'

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Academic journal color scheme (high contrast, print-friendly)
        # Using distinct colors from ColorBrewer Set1 palette
        colors = {
            'alpha_cdr1': '#E41A1C',  # Red
            'alpha_cdr2': '#377EB8',  # Blue
            'alpha_cdr3': '#4DAF4A',  # Green
            'beta_cdr1': '#984EA3',   # Purple
            'beta_cdr2': '#FF7F00',   # Orange
            'beta_cdr3': '#A65628'    # Brown
        }

        # Line styles for better distinction in grayscale printing
        line_styles = {
            'alpha_cdr1': '-',      # Solid
            'alpha_cdr2': '-',      # Solid
            'alpha_cdr3': '-',      # Solid
            'beta_cdr1': '--',      # Dashed
            'beta_cdr2': '--',      # Dashed
            'beta_cdr3': '--'       # Dashed
        }

        # Markers for additional distinction
        markers = {
            'alpha_cdr1': 'o', 'alpha_cdr2': 's', 'alpha_cdr3': '^',
            'beta_cdr1': 'o', 'beta_cdr2': 's', 'beta_cdr3': '^'
        }

        # Plot each CDR curve
        for cdr_key in cdr_mapping.keys():
            if cdr_key not in cdr_stats:
                continue

            stats = cdr_stats[cdr_key]
            label = f"{cdr_mapping[cdr_key]} (n={stats['n_samples']})"

            # Plot mean curve with markers every 10th point
            markevery = max(1, n_points // 10)
            ax.plot(normalized_grid, stats['mean'],
                   color=colors[cdr_key],
                   linestyle=line_styles[cdr_key],
                   linewidth=2.0,
                   label=label,
                   alpha=0.95,
                   marker=markers[cdr_key],
                   markevery=markevery,
                   markersize=5,
                   markerfacecolor='white',
                   markeredgewidth=1.5,
                   markeredgecolor=colors[cdr_key])

            # Plot very subtle std region (minimal, for reference only)
            ax.fill_between(normalized_grid,
                           stats['mean'] - stats['std'],
                           stats['mean'] + stats['std'],
                           color=colors[cdr_key], alpha=0.08, linewidth=0)

        # Customize axes - clean academic style
        ax.set_xlabel('Normalized Relative Position', fontsize=12, fontweight='normal')
        ax.set_ylabel('RMSF (nm)', fontsize=12, fontweight='normal')
        ax.set_title('CDR Loop Flexibility Patterns',
                    fontsize=14, fontweight='bold', pad=15)

        # Set axis limits
        ax.set_xlim(0, 1)
        ax.set_ylim(bottom=0)

        # Clean grid style
        ax.grid(True, which='major', axis='both', alpha=0.25,
               linestyle='-', linewidth=0.6, color='#CCCCCC')
        ax.set_axisbelow(True)

        # Spine styling
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.2)
        ax.spines['bottom'].set_linewidth(1.2)

        # Tick parameters
        ax.tick_params(axis='both', which='major', labelsize=10,
                      width=1.2, length=5, direction='out')

        # Legend - compact and clean
        legend = ax.legend(loc='upper left', frameon=True, framealpha=1.0,
                          fontsize=9, ncol=2, columnspacing=1.0,
                          handlelength=2.5, handletextpad=0.5,
                          edgecolor='#333333', fancybox=False)
        legend.get_frame().set_linewidth(1.0)

        plt.tight_layout()

        # Save figure
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"CDR RMSF curves comparison plot saved: {output_path}")

        # Display plot
        if show_plot:
            plt.show()
        else:
            plt.close()