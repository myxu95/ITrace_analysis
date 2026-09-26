"""
Conformational Clustering Tool

Performs RMSD-based hierarchical clustering of interface conformations.
"""

import asyncio
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class ClusterConformationsInput(BaseModel):
    """Input schema for conformational clustering."""

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
        description="Mapping of biological roles to chain IDs (mhc_alpha, b2m, peptide, tcr_alpha, tcr_beta)"
    )
    cdr_detection: dict = Field(
        description="CDR detection results from TCR analysis"
    )
    output_csv: str | None = Field(
        default=None,
        description="Optional path to save cluster assignments as CSV"
    )
    output_plot: str | None = Field(
        default=None,
        description="Optional path to save cluster visualization"
    )
    stride: int = Field(
        default=10,
        description="Frame stride for analysis (default: 10)"
    )
    contact_cutoff: float = Field(
        default=4.5,
        description="Contact distance cutoff in Angstroms (default: 4.5)"
    )
    distance_cutoff: float = Field(
        default=0.35,
        description="Hierarchical clustering distance cutoff (default: 0.35)"
    )
    linkage_method: str = Field(
        default="average",
        description="Linkage method: 'average', 'single', 'complete', 'ward' (default: average)"
    )
    geometry_weight: float = Field(
        default=0.45,
        description="Weight for geometric features (default: 0.45)"
    )
    sidechain_weight: float = Field(
        default=0.25,
        description="Weight for sidechain orientation (default: 0.25)"
    )
    interaction_weight: float = Field(
        default=0.30,
        description="Weight for interaction patterns (default: 0.30)"
    )


@register_tool
class ClusterConformationsTool(Tool):
    """Tool for RMSD-based conformational clustering of interface states."""

    name: ClassVar[str] = "cluster_conformations"
    description: ClassVar[str] = (
        "Perform hierarchical clustering of interface conformations based on RMSD. "
        "Combines geometric features (backbone RMSD), sidechain orientations, and interaction patterns "
        "to identify distinct conformational states. Useful for understanding conformational heterogeneity."
    )
    Input: ClassVar[type[BaseModel]] = ClusterConformationsInput

    # Tool metadata
    is_read_only: ClassVar[bool] = False
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 600.0  # 10 minutes

    async def call(self, args: ClusterConformationsInput, ctx: ToolContext) -> ToolResult:
        """Execute conformational clustering."""
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

        # Validate chain mapping
        required_chains = ['mhc_alpha', 'peptide', 'tcr_alpha', 'tcr_beta']
        missing = [c for c in required_chains if c not in args.chain_mapping]
        if missing:
            return ToolResult(
                is_error=True,
                error=f"Missing required chains in chain_mapping: {missing}"
            )

        # Validate linkage method
        valid_methods = ['average', 'single', 'complete', 'ward']
        if args.linkage_method not in valid_methods:
            return ToolResult(
                is_error=True,
                error=f"Invalid linkage_method: {args.linkage_method}. Must be one of {valid_methods}"
            )

        # Validate weights sum to 1.0
        weight_sum = args.geometry_weight + args.sidechain_weight + args.interaction_weight
        if abs(weight_sum - 1.0) > 0.01:
            return ToolResult(
                is_error=True,
                error=f"Feature weights must sum to 1.0, got {weight_sum:.3f}"
            )

        await ctx.on_event({
            "type": "progress",
            "message": "Starting conformational clustering analysis..."
        })

        loop = asyncio.get_event_loop()
        try:
            result = await loop.run_in_executor(
                None,
                self._run_clustering,
                args,
                ctx
            )
            return result
        except Exception as e:
            return ToolResult(
                is_error=True,
                error=f"Clustering analysis failed: {str(e)}"
            )

    def _run_clustering(
        self,
        args: ClusterConformationsInput,
        ctx: ToolContext
    ) -> ToolResult:
        """Run clustering analysis (blocking operation)."""
        from immunoscope.analysis.conformation.interface_clustering import InterfaceClusteringAnalyzer

        # Initialize analyzer
        analyzer = InterfaceClusteringAnalyzer(
            topology_file=args.topology,
            trajectory_file=args.trajectory,
            structure_pdb=args.structure_pdb
        )

        # Run clustering
        result = analyzer.calculate(
            chain_mapping=args.chain_mapping,
            cdr_detection=args.cdr_detection,
            stride=args.stride,
            contact_cutoff_angstrom=args.contact_cutoff,
            distance_cutoff=args.distance_cutoff,
            linkage_method=args.linkage_method,
            geometry_weight=args.geometry_weight,
            sidechain_weight=args.sidechain_weight,
            interaction_weight=args.interaction_weight
        )

        # Extract key statistics
        n_clusters = result.summary['n_clusters']
        n_frames = result.summary['n_frames']

        # Get cluster populations
        cluster_pops = result.cluster_summary[['cluster_id', 'n_frames', 'population_percent']].to_dict('records')

        # Get top cluster features
        top_features = []
        if not result.cluster_feature_digest.empty:
            top_features = result.cluster_feature_digest.head(10).to_dict('records')

        # Save outputs
        output_files = {}
        if args.output_csv:
            result.frame_assignments.to_csv(args.output_csv, index=False)
            output_files['assignments'] = args.output_csv

            # Also save cluster summary
            summary_path = str(Path(args.output_csv).with_suffix('')) + '_summary.csv'
            result.cluster_summary.to_csv(summary_path, index=False)
            output_files['summary'] = summary_path

        if args.output_plot:
            # Create visualization if requested
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(1, 2, figsize=(14, 5))

            # Plot 1: Cluster populations
            ax1 = axes[0]
            cluster_ids = result.cluster_summary['cluster_id'].values
            populations = result.cluster_summary['population_percent'].values
            ax1.bar(cluster_ids, populations)
            ax1.set_xlabel('Cluster ID')
            ax1.set_ylabel('Population (%)')
            ax1.set_title('Cluster Populations')
            ax1.grid(axis='y', alpha=0.3)

            # Plot 2: Cluster timeline
            ax2 = axes[1]
            frames = result.frame_assignments['frame'].values
            clusters = result.frame_assignments['cluster_id'].values
            ax2.scatter(frames, clusters, alpha=0.5, s=10)
            ax2.set_xlabel('Frame')
            ax2.set_ylabel('Cluster ID')
            ax2.set_title('Cluster Assignment Timeline')
            ax2.grid(alpha=0.3)

            plt.tight_layout()
            plt.savefig(args.output_plot, dpi=300, bbox_inches='tight')
            plt.close()
            output_files['plot'] = args.output_plot

        # Format cluster descriptions
        cluster_descriptions = []
        for _, row in result.cluster_summary.iterrows():
            desc = {
                'cluster_id': int(row['cluster_id']),
                'n_frames': int(row['n_frames']),
                'population_percent': round(float(row['population_percent']), 2),
                'avg_intra_distance': round(float(row['avg_intra_distance']), 3)
            }
            if 'structural_descriptor' in row:
                desc['descriptor'] = row['structural_descriptor']
            cluster_descriptions.append(desc)

        return ToolResult(
            is_error=False,
            data={
                'n_clusters': n_clusters,
                'n_frames': n_frames,
                'stride': args.stride,
                'distance_cutoff': args.distance_cutoff,
                'linkage_method': args.linkage_method,
                'clusters': cluster_descriptions,
                'top_discriminating_features': top_features[:5] if top_features else [],
                'output_files': output_files
            },
            message=(
                f"Clustering complete. Identified {n_clusters} distinct conformational states "
                f"from {n_frames} frames. Largest cluster: {cluster_descriptions[0]['population_percent']:.1f}%."
            )
        )

    def system_prompt_section(self) -> str:
        """Return system prompt section for this tool."""
        return """
## cluster_conformations

Perform hierarchical clustering of interface conformations based on RMSD and interaction patterns.

**When to use:**
- To identify distinct conformational states in the trajectory
- To understand conformational heterogeneity
- To find representative structures for each state
- After landscape analysis to validate energy basins
- To quantify conformational transitions

**Input requirements:**
- Preprocessed trajectory
- Topology and reference PDB
- Chain mapping (mhc_alpha, b2m, peptide, tcr_alpha, tcr_beta)
- CDR detection results (from TCR analysis)

**Feature types (weighted combination):**
1. **Geometry (default 0.45)**: Backbone RMSD of CDR3 loops and peptide
   - Captures overall conformational changes
   - Most important for structural clustering

2. **Sidechain (default 0.25)**: Sidechain orientation relative to CA
   - Captures rotamer changes
   - Important for interaction specificity

3. **Interaction (default 0.30)**: Contact patterns between CDR3 and peptide/MHC
   - Captures functional interaction changes
   - Important for binding mode classification

**Clustering parameters:**
- **distance_cutoff**: Lower = more clusters (more granular)
  - 0.2-0.3: Fine-grained (many small clusters)
  - 0.35-0.45: Moderate (balanced)
  - 0.5+: Coarse (few large clusters)

- **linkage_method**:
  - 'average': Balanced, most common (default)
  - 'complete': Compact clusters, sensitive to outliers
  - 'single': Can create elongated clusters
  - 'ward': Minimizes variance, good for equal-sized clusters

**Output interpretation:**
- **n_clusters**: Number of distinct conformational states
- **population_percent**: Time spent in each state
- **avg_intra_distance**: Cluster compactness (lower = tighter)
- **structural_descriptor**: Human-readable cluster characterization
- **discriminating_features**: Key differences between clusters

**Typical patterns:**
- 1-2 clusters: Stable, homogeneous binding
- 3-5 clusters: Moderate flexibility, multiple substates
- 6+ clusters: High flexibility, many transient states
- Dominant cluster (>70%): One preferred conformation
- Balanced clusters: Multiple equally stable states

**Best practices:**
- Use stride=10 for large trajectories
- Start with default weights, adjust if needed
- Compare with FEL analysis (clusters should match energy basins)
- Examine top discriminating features to understand differences
- Representative frames from each cluster are good for visualization
- In pHLA-TCR: clusters often correspond to CDR3 loop conformations
"""
