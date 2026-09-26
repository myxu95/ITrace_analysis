"""
Salt Bridge Analysis Tool

Analyzes salt bridges (ionic interactions) at molecular interfaces.
"""

import asyncio
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class AnalyzeSaltBridgesInput(BaseModel):
    """Input schema for salt bridge analysis."""

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
        description="Optional path to save salt bridge results as CSV"
    )
    output_plot: str | None = Field(
        default=None,
        description="Optional path to save salt bridge occupancy plot"
    )
    distance_cutoff: float = Field(
        default=4.0,
        description="Distance cutoff in Angstroms (default: 4.0)"
    )
    stride: int = Field(
        default=1,
        description="Frame stride for analysis (default: 1)"
    )


@register_tool
class AnalyzeSaltBridgesTool(Tool):
    """Tool for analyzing salt bridges at molecular interfaces."""

    name: ClassVar[str] = "analyze_salt_bridges"
    description: ClassVar[str] = (
        "Analyze salt bridges (ionic interactions) between charged residues. "
        "Salt bridges are electrostatic interactions between positively and negatively charged groups. "
        "They contribute significantly to binding specificity and stability."
    )
    Input: ClassVar[type[BaseModel]] = AnalyzeSaltBridgesInput

    # Tool metadata
    is_read_only: ClassVar[bool] = False
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 600.0

    async def call(self, args: AnalyzeSaltBridgesInput, ctx: ToolContext) -> ToolResult:
        """Execute salt bridge analysis."""
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

        if args.distance_cutoff <= 0 or args.distance_cutoff > 6.0:
            return ToolResult(
                is_error=True,
                error=f"Distance cutoff must be between 0 and 6.0 Å, got {args.distance_cutoff}"
            )

        await ctx.on_event({
            "type": "progress",
            "message": f"Starting salt bridge analysis (d<{args.distance_cutoff}Å)..."
        })

        loop = asyncio.get_event_loop()
        try:
            result = await loop.run_in_executor(
                None,
                self._run_saltbridge_analysis,
                args,
                ctx
            )
            return result
        except Exception as e:
            return ToolResult(
                is_error=True,
                error=f"Salt bridge analysis failed: {str(e)}"
            )

    def _run_saltbridge_analysis(
        self,
        args: AnalyzeSaltBridgesInput,
        ctx: ToolContext
    ) -> ToolResult:
        """Run salt bridge analysis (blocking operation)."""
        from immunoscope.analysis.interactions.salt_bridge_pairs import SaltBridgeAnalyzer

        # Initialize analyzer
        analyzer = SaltBridgeAnalyzer(
            trajectory_path=args.trajectory,
            topology_path=args.topology,
            structure_pdb=args.structure_pdb,
            chain_mapping=args.chain_mapping
        )

        # Run analysis
        sb_df = analyzer.analyze(
            distance_cutoff=args.distance_cutoff,
            stride=args.stride
        )

        # Get statistics
        total_saltbridges = len(sb_df)
        persistent_sb = sb_df[sb_df['occupancy'] > 0.5]
        mean_occupancy = sb_df['occupancy'].mean() if len(sb_df) > 0 else 0.0

        # Top salt bridges by occupancy
        top_sb = sb_df.nlargest(10, 'occupancy') if len(sb_df) > 0 else sb_df
        top_sb_list = []
        for _, row in top_sb.iterrows():
            top_sb_list.append({
                'positive': f"{row['pos_chain']}:{row['pos_resname']}{row['pos_resid']}",
                'negative': f"{row['neg_chain']}:{row['neg_resname']}{row['neg_resid']}",
                'occupancy': round(row['occupancy'], 3),
                'mean_distance': round(row['mean_distance'], 2)
            })

        # Save outputs
        output_files = {}
        if args.output_csv:
            sb_df.to_csv(args.output_csv, index=False)
            output_files['csv'] = args.output_csv

        if args.output_plot and len(sb_df) > 0:
            analyzer.plot_occupancy(
                sb_df=sb_df,
                output_path=args.output_plot
            )
            output_files['plot'] = args.output_plot

        return ToolResult(
            is_error=False,
            data={
                'total_salt_bridges': total_saltbridges,
                'persistent_salt_bridges': len(persistent_sb),
                'mean_occupancy': round(mean_occupancy, 3),
                'top_10_salt_bridges': top_sb_list,
                'output_files': output_files,
                'distance_cutoff': args.distance_cutoff
            },
            message=(
                f"Salt bridge analysis complete. Found {total_saltbridges} salt bridges, "
                f"{len(persistent_sb)} persistent (>50% occupancy). "
                f"Mean occupancy: {mean_occupancy:.3f}."
            )
        )

    def system_prompt_section(self) -> str:
        """Return system prompt section for this tool."""
        return """
## analyze_salt_bridges

Analyze salt bridges (ionic interactions) at molecular interfaces.

**When to use:**
- To identify electrostatic interactions stabilizing the interface
- To understand pH-dependent binding
- To predict effects of charge-altering mutations
- After H-bond analysis to get complete interaction picture

**Input requirements:**
- Preprocessed trajectory
- Topology and reference PDB
- Chain mapping for interface definition

**Charged residues detected:**
- Positive: LYS (NZ), ARG (NH1, NH2), HIS (ND1, NE2)
- Negative: ASP (OD1, OD2), GLU (OE1, OE2)
- Distance between charged atoms < 4.0 Å (default)

**Output interpretation:**
- Occupancy: % of frames where salt bridge is present
- Occupancy > 75%: Very stable, critical for binding
- Occupancy 50-75%: Stable, important for specificity
- Occupancy 25-50%: Moderate, may contribute
- Occupancy < 25%: Transient, less important
- Mean distance: Average distance when salt bridge is formed

**Best practices:**
- Use default 4.0 Å cutoff for standard analysis
- Salt bridges are pH-dependent (HIS protonation state matters)
- Compare with H-bonds - salt bridges often accompanied by H-bonds
- Buried salt bridges are more stable than surface-exposed ones
- In pHLA-TCR: peptide charged residues often form salt bridges with TCR
- Charge complementarity is a key specificity determinant
"""
