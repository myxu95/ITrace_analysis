"""
Free Energy Landscape Analysis Tool

Performs dimensionality reduction and free energy landscape calculation.
"""

import asyncio
from pathlib import Path
from typing import ClassVar

import numpy as np
from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class AnalyzeLandscapeInput(BaseModel):
    """Input schema for landscape analysis."""

    trajectory: str = Field(
        description="Path to trajectory file (.xtc, .trr)"
    )
    topology: str = Field(
        description="Path to topology file (.tpr, .gro, .pdb)"
    )
    structure_pdb: str = Field(
        description="Path to reference PDB for residue annotation"
    )
    chain_mapping: dict[str, str] = Field(
        description="Mapping of biological roles to chain IDs"
    )
    output_csv: str | None = Field(
        default=None,
        description="Optional path to save PC coordinates as CSV"
    )
    output_plot: str | None = Field(
        default=None,
        description="Optional path to save FEL plot"
    )
    reducer: str = Field(
        default="pca",
        description="Dimensionality reduction method: 'pca', 'tica', or 'umap' (default: pca)"
    )
    n_components: int = Field(
        default=2,
        description="Number of components to extract (default: 2)"
    )
    tica_lag: int = Field(
        default=10,
        description="Lag time for TICA in frames (default: 10)"
    )
    temperature: float = Field(
        default=300.0,
        description="Temperature in Kelvin for free energy calculation (default: 300.0)"
    )
    bins: int = Field(
        default=50,
        description="Number of bins for FEL histogram (default: 50)"
    )
    stride: int = Field(
        default=1,
        description="Frame stride for analysis (default: 1)"
    )


@register_tool
class AnalyzeLandscapeTool(Tool):
    """Tool for free energy landscape analysis with dimensionality reduction."""

    name: ClassVar[str] = "analyze_landscape"
    description: ClassVar[str] = (
        "Perform free energy landscape (FEL) analysis using dimensionality reduction. "
        "Reduces high-dimensional conformational space to 2D using PCA/TICA/UMAP, "
        "then calculates free energy surface. Identifies conformational states and transition barriers."
    )
    Input: ClassVar[type[BaseModel]] = AnalyzeLandscapeInput

    # Tool metadata
    is_read_only: ClassVar[bool] = False
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 900.0  # 15 minutes for large trajectories

    async def call(self, args: AnalyzeLandscapeInput, ctx: ToolContext) -> ToolResult:
        """Execute landscape analysis."""
        # Validate inputs
        traj_path = Path(args.trajectory)
        top_path = Path(args.topology)
        pdb_path = Path(args.structure_pdb)

        if not traj_path.exists():
            return ToolResult(
                is_error=True,
                error=f"Trajectory file not found: {args.trajectory}"
            )

        if not top_path.exists():
            return ToolResult(
                is_error=True,
                error=f"Topology file not found: {args.topology}"
            )

        if not pdb_path.exists():
            return ToolResult(
                is_error=True,
                error=f"Structure PDB not found: {args.structure_pdb}"
            )

        if args.reducer not in ['pca', 'tica', 'umap']:
            return ToolResult(
                is_error=True,
                error=f"Invalid reducer: {args.reducer}. Must be 'pca', 'tica', or 'umap'"
            )

        if args.n_components < 2:
            return ToolResult(
                is_error=True,
                error=f"n_components must be >= 2, got {args.n_components}"
            )

        if args.temperature <= 0:
            return ToolResult(
                is_error=True,
                error=f"Temperature must be positive, got {args.temperature}"
            )

        await ctx.on_event({
            "type": "progress",
            "message": f"Starting landscape analysis with {args.reducer.upper()}..."
        })

        loop = asyncio.get_event_loop()
        try:
            result = await loop.run_in_executor(
                None,
                self._run_landscape_analysis,
                args,
                ctx
            )
            return result
        except Exception as e:
            return ToolResult(
                is_error=True,
                error=f"Landscape analysis failed: {str(e)}"
            )

    def _run_landscape_analysis(
        self,
        args: AnalyzeLandscapeInput,
        ctx: ToolContext
    ) -> ToolResult:
        """Run landscape analysis (blocking operation)."""
        from immunoscope.analysis.landscape.feature_matrix import FeatureMatrixBuilder
        from immunoscope.analysis.landscape.landscape_analyzer import LandscapeAnalyzer, LandscapeInput

        # Build feature matrix
        builder = FeatureMatrixBuilder(
            trajectory_path=args.trajectory,
            topology_path=args.topology,
            structure_pdb=args.structure_pdb,
            chain_mapping=args.chain_mapping
        )

        feature_result = builder.build(stride=args.stride)
        feature_matrix = feature_result.feature_matrix
        feature_names = feature_result.feature_names
        times = feature_result.times

        # Prepare landscape input
        landscape_input = LandscapeInput(
            feature_matrix=feature_matrix,
            feature_names=feature_names,
            times=times,
            reducer=args.reducer,
            n_components=args.n_components,
            tica_lag=args.tica_lag,
            temperature=args.temperature,
            bins=args.bins
        )

        # Run landscape analysis
        analyzer = LandscapeAnalyzer()
        landscape_result = analyzer.analyze(landscape_input)

        # Extract statistics
        n_frames = len(landscape_result.pc_coordinates)
        explained_var = landscape_result.explained_variance

        # Get top contributors for PC1 and PC2
        top_pc1 = landscape_result.get_top_contributors(0, top_n=5)
        top_pc2 = landscape_result.get_top_contributors(1, top_n=5) if args.n_components >= 2 else []

        # Free energy statistics
        fel = landscape_result.free_energy
        min_energy = float(np.nanmin(fel))
        max_energy = float(np.nanmax(fel))
        energy_range = max_energy - min_energy

        # Save outputs
        output_files = {}
        if args.output_csv:
            pc_df = landscape_result.to_pca_dataframe()
            pc_df.to_csv(args.output_csv, index=False)
            output_files['csv'] = args.output_csv

        if args.output_plot:
            from immunoscope.analysis.landscape.visualizer import LandscapeVisualizer
            visualizer = LandscapeVisualizer()
            visualizer.plot_fel(
                landscape_result=landscape_result,
                output_path=args.output_plot
            )
            output_files['plot'] = args.output_plot

        # Format top contributors
        top_pc1_formatted = [
            {'feature': name, 'loading': round(loading, 3)}
            for name, loading in top_pc1
        ]
        top_pc2_formatted = [
            {'feature': name, 'loading': round(loading, 3)}
            for name, loading in top_pc2
        ]

        return ToolResult(
            is_error=False,
            data={
                'n_frames': n_frames,
                'n_features': len(feature_names),
                'reducer': args.reducer,
                'n_components': args.n_components,
                'explained_variance': [round(v, 4) for v in explained_var.tolist()] if len(explained_var) > 0 else None,
                'total_variance_explained': round(float(explained_var.sum()), 4) if len(explained_var) > 0 else None,
                'top_pc1_contributors': top_pc1_formatted,
                'top_pc2_contributors': top_pc2_formatted,
                'free_energy_range': round(energy_range, 2),
                'min_free_energy': round(min_energy, 2),
                'max_free_energy': round(max_energy, 2),
                'output_files': output_files
            },
            message=(
                f"Landscape analysis complete using {args.reducer.upper()}. "
                f"Analyzed {n_frames} frames with {len(feature_names)} features. "
                f"Free energy range: {energy_range:.2f} kJ/mol."
            )
        )

    def system_prompt_section(self) -> str:
        """Return system prompt section for this tool."""
        return """
## analyze_landscape

Perform free energy landscape (FEL) analysis with dimensionality reduction.

**When to use:**
- To visualize conformational space explored during simulation
- To identify stable conformational states (energy minima)
- To detect transition barriers between states
- To understand dominant motions driving conformational changes
- After RMSD/RMSF analysis to understand conformational heterogeneity

**Input requirements:**
- Preprocessed trajectory
- Topology and reference PDB
- Chain mapping for feature extraction

**Dimensionality reduction methods:**
1. **PCA (Principal Component Analysis)**: Default, fast, linear
   - Best for: Initial exploration, linear motions
   - Explained variance: quantifies information retained
   - Loadings: show which features contribute to each PC

2. **TICA (Time-lagged Independent Component Analysis)**: Slow dynamics
   - Best for: Identifying slow conformational transitions
   - Requires: tica_lag parameter (typically 10-50 frames)
   - Captures: Kinetically relevant motions

3. **UMAP (Uniform Manifold Approximation and Projection)**: Nonlinear
   - Best for: Complex, nonlinear conformational changes
   - Slower than PCA, no explained variance
   - Better separation of distinct states

**Output interpretation:**
- **Free energy (kJ/mol)**: Lower = more stable/populated
- **Energy minima**: Stable conformational states (basins)
- **Energy barriers**: Transition states between basins
- **Barrier height**: Activation energy for transitions
- **PC loadings**: Features driving conformational changes

**Typical FEL patterns:**
- Single basin: One dominant conformation
- Two basins: Two-state system (e.g., open/closed)
- Multiple basins: Multi-state system
- Flat landscape: High conformational flexibility
- Deep basin: Very stable conformation

**Best practices:**
- Start with PCA for initial exploration
- Use TICA for kinetic analysis (set lag = 10-50 frames)
- Use stride=10 for large trajectories to speed up
- n_components=2 for visualization, 3-5 for detailed analysis
- Compare with clustering results for validation
- Top PC contributors reveal which structural features drive motion
- In pHLA-TCR: CDR loop motions often dominate PC1/PC2
"""
