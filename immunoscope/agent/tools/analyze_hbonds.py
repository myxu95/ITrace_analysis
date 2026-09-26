"""
Hydrogen Bond Analysis Tool

Analyzes hydrogen bonds at molecular interfaces.
"""

import asyncio
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class AnalyzeHBondsInput(BaseModel):
    """Input schema for hydrogen bond analysis."""

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
        description="Optional path to save hydrogen bond results as CSV"
    )
    output_plot: str | None = Field(
        default=None,
        description="Optional path to save hydrogen bond occupancy plot"
    )
    distance_cutoff: float = Field(
        default=3.5,
        description="Distance cutoff in Angstroms (default: 3.5)"
    )
    angle_cutoff: float = Field(
        default=150.0,
        description="Angle cutoff in degrees (default: 150)"
    )
    stride: int = Field(
        default=1,
        description="Frame stride for analysis (default: 1)"
    )


@register_tool
class AnalyzeHBondsTool(Tool):
    """Tool for analyzing hydrogen bonds at molecular interfaces."""

    name: ClassVar[str] = "analyze_hbonds"
    description: ClassVar[str] = (
        "Analyze hydrogen bonds between molecular components. "
        "Identifies persistent H-bonds that stabilize the interface. "
        "Reports occupancy (% of frames where H-bond is present) and geometric properties."
    )
    Input: ClassVar[type[BaseModel]] = AnalyzeHBondsInput

    # Tool metadata
    is_read_only: ClassVar[bool] = False
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 600.0

    async def call(self, args: AnalyzeHBondsInput, ctx: ToolContext) -> ToolResult:
        """Execute hydrogen bond analysis."""
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

        if args.distance_cutoff <= 0 or args.distance_cutoff > 5.0:
            return ToolResult(
                is_error=True,
                error=f"Distance cutoff must be between 0 and 5.0 Å, got {args.distance_cutoff}"
            )

        if args.angle_cutoff < 90 or args.angle_cutoff > 180:
            return ToolResult(
                is_error=True,
                error=f"Angle cutoff must be between 90 and 180 degrees, got {args.angle_cutoff}"
            )

        await ctx.on_event({
            "type": "progress",
            "message": f"Starting H-bond analysis (d<{args.distance_cutoff}Å, angle>{args.angle_cutoff}°)..."
        })

        loop = asyncio.get_event_loop()
        try:
            result = await loop.run_in_executor(
                None,
                self._run_hbond_analysis,
                args,
                ctx
            )
            return result
        except Exception as e:
            return ToolResult(
                is_error=True,
                error=f"H-bond analysis failed: {str(e)}"
            )

    def _run_hbond_analysis(
        self,
        args: AnalyzeHBondsInput,
        ctx: ToolContext
    ) -> ToolResult:
        """Run hydrogen bond analysis (blocking operation)."""
        from immunoscope.analysis.interactions.hydrogen_bond_pairs import HydrogenBondAnalyzer

        # Initialize analyzer
        analyzer = HydrogenBondAnalyzer(
            trajectory_path=args.trajectory,
            topology_path=args.topology,
            structure_pdb=args.structure_pdb,
            chain_mapping=args.chain_mapping
        )

        # Run analysis
        hbond_df = analyzer.analyze(
            distance_cutoff=args.distance_cutoff,
            angle_cutoff=args.angle_cutoff,
            stride=args.stride
        )

        # Get statistics
        total_hbonds = len(hbond_df)
        persistent_hbonds = hbond_df[hbond_df['occupancy'] > 0.5]
        mean_occupancy = hbond_df['occupancy'].mean()

        # Top H-bonds by occupancy
        top_hbonds = hbond_df.nlargest(10, 'occupancy')
        top_hbonds_list = []
        for _, row in top_hbonds.iterrows():
            top_hbonds_list.append({
                'donor': f"{row['donor_chain']}:{row['donor_resname']}{row['donor_resid']}@{row['donor_atom']}",
                'acceptor': f"{row['acceptor_chain']}:{row['acceptor_resname']}{row['acceptor_resid']}@{row['acceptor_atom']}",
                'occupancy': round(row['occupancy'], 3),
                'mean_distance': round(row['mean_distance'], 2),
                'mean_angle': round(row['mean_angle'], 1)
            })

        # Save outputs
        output_files = {}
        if args.output_csv:
            hbond_df.to_csv(args.output_csv, index=False)
            output_files['csv'] = args.output_csv

        if args.output_plot:
            analyzer.plot_occupancy(
                hbond_df=hbond_df,
                output_path=args.output_plot
            )
            output_files['plot'] = args.output_plot

        return ToolResult(
            is_error=False,
            data={
                'total_hbonds': total_hbonds,
                'persistent_hbonds': len(persistent_hbonds),
                'mean_occupancy': round(mean_occupancy, 3),
                'top_10_hbonds': top_hbonds_list,
                'output_files': output_files,
                'distance_cutoff': args.distance_cutoff,
                'angle_cutoff': args.angle_cutoff
            },
            message=(
                f"H-bond analysis complete. Found {total_hbonds} H-bonds, "
                f"{len(persistent_hbonds)} persistent (>50% occupancy). "
                f"Mean occupancy: {mean_occupancy:.3f}."
            )
        )

    def system_prompt_section(self) -> str:
        """Return system prompt section for this tool."""
        return """
## analyze_hbonds

Analyze hydrogen bonds at molecular interfaces.

**When to use:**
- To identify key H-bonds stabilizing the interface
- To understand specificity determinants
- To validate structural models
- After RRCS analysis to understand interaction chemistry

**Input requirements:**
- Preprocessed trajectory
- Topology and reference PDB
- Chain mapping for interface definition

**Geometric criteria:**
- Distance: Donor-Acceptor < 3.5 Å (default)
- Angle: Donor-H-Acceptor > 150° (default)
- Standard criteria from Baker & Hubbard (1984)

**Output interpretation:**
- Occupancy: % of frames where H-bond is present
- Occupancy > 75%: Very persistent, likely critical for binding
- Occupancy 50-75%: Persistent, contributes to stability
- Occupancy 25-50%: Moderate, may be important
- Occupancy < 25%: Transient, less critical
- Mean distance: Average D-A distance when H-bond is formed
- Mean angle: Average D-H-A angle when H-bond is formed

**Best practices:**
- Use default cutoffs (3.5 Å, 150°) for standard analysis
- Focus on H-bonds with occupancy > 50% for key interactions
- Compare with experimental data (mutagenesis, crystallography)
- H-bonds involving backbone atoms are often more stable
- Side chain H-bonds provide specificity
- In pHLA-TCR: peptide-TCR H-bonds are critical for recognition
"""
