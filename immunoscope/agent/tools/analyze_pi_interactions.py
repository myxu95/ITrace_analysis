"""
Pi Interaction Analysis Tool

Analyzes pi-pi stacking and cation-pi interactions at molecular interfaces.
"""

import asyncio
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class AnalyzePiInteractionsInput(BaseModel):
    """Input schema for pi interaction analysis."""

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
        description="Optional path to save pi interaction results as CSV"
    )
    output_plot: str | None = Field(
        default=None,
        description="Optional path to save pi interaction occupancy plot"
    )
    distance_cutoff: float = Field(
        default=6.0,
        description="Distance cutoff in Angstroms (default: 6.0)"
    )
    angle_cutoff: float = Field(
        default=30.0,
        description="Angle cutoff in degrees for pi-pi stacking (default: 30.0)"
    )
    stride: int = Field(
        default=1,
        description="Frame stride for analysis (default: 1)"
    )


@register_tool
class AnalyzePiInteractionsTool(Tool):
    """Tool for analyzing pi-pi stacking and cation-pi interactions."""

    name: ClassVar[str] = "analyze_pi_interactions"
    description: ClassVar[str] = (
        "Analyze pi-pi stacking and cation-pi interactions involving aromatic residues. "
        "These interactions contribute to binding specificity and are important for "
        "recognizing aromatic peptide residues."
    )
    Input: ClassVar[type[BaseModel]] = AnalyzePiInteractionsInput

    # Tool metadata
    is_read_only: ClassVar[bool] = False
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 600.0

    async def call(self, args: AnalyzePiInteractionsInput, ctx: ToolContext) -> ToolResult:
        """Execute pi interaction analysis."""
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

        if args.distance_cutoff <= 0 or args.distance_cutoff > 8.0:
            return ToolResult(
                is_error=True,
                error=f"Distance cutoff must be between 0 and 8.0 Å, got {args.distance_cutoff}"
            )

        if args.angle_cutoff <= 0 or args.angle_cutoff > 90.0:
            return ToolResult(
                is_error=True,
                error=f"Angle cutoff must be between 0 and 90 degrees, got {args.angle_cutoff}"
            )

        await ctx.on_event({
            "type": "progress",
            "message": f"Starting pi interaction analysis (d<{args.distance_cutoff}Å, angle<{args.angle_cutoff}°)..."
        })

        loop = asyncio.get_event_loop()
        try:
            result = await loop.run_in_executor(
                None,
                self._run_pi_analysis,
                args,
                ctx
            )
            return result
        except Exception as e:
            return ToolResult(
                is_error=True,
                error=f"Pi interaction analysis failed: {str(e)}"
            )

    def _run_pi_analysis(
        self,
        args: AnalyzePiInteractionsInput,
        ctx: ToolContext
    ) -> ToolResult:
        """Run pi interaction analysis (blocking operation)."""
        from immunoscope.analysis.interactions.pi_interaction_pairs import PiInteractionAnalyzer

        # Initialize analyzer
        analyzer = PiInteractionAnalyzer(
            trajectory_path=args.trajectory,
            topology_path=args.topology,
            structure_pdb=args.structure_pdb,
            chain_mapping=args.chain_mapping
        )

        # Run analysis
        pi_df = analyzer.analyze(
            distance_cutoff=args.distance_cutoff,
            angle_cutoff=args.angle_cutoff,
            stride=args.stride
        )

        # Separate pi-pi and cation-pi
        pi_pi_df = pi_df[pi_df['interaction_type'] == 'pi-pi']
        cation_pi_df = pi_df[pi_df['interaction_type'] == 'cation-pi']

        # Get statistics
        total_pi_pi = len(pi_pi_df)
        total_cation_pi = len(cation_pi_df)
        persistent_pi_pi = pi_pi_df[pi_pi_df['occupancy'] > 0.5]
        persistent_cation_pi = cation_pi_df[cation_pi_df['occupancy'] > 0.5]

        # Top interactions by occupancy
        top_pi_pi = pi_pi_df.nlargest(5, 'occupancy') if len(pi_pi_df) > 0 else pi_pi_df
        top_cation_pi = cation_pi_df.nlargest(5, 'occupancy') if len(cation_pi_df) > 0 else cation_pi_df

        top_pi_pi_list = []
        for _, row in top_pi_pi.iterrows():
            top_pi_pi_list.append({
                'residue1': f"{row['chain1']}:{row['resname1']}{row['resid1']}",
                'residue2': f"{row['chain2']}:{row['resname2']}{row['resid2']}",
                'occupancy': round(row['occupancy'], 3),
                'mean_distance': round(row['mean_distance'], 2)
            })

        top_cation_pi_list = []
        for _, row in top_cation_pi.iterrows():
            top_cation_pi_list.append({
                'cation': f"{row['chain1']}:{row['resname1']}{row['resid1']}",
                'aromatic': f"{row['chain2']}:{row['resname2']}{row['resid2']}",
                'occupancy': round(row['occupancy'], 3),
                'mean_distance': round(row['mean_distance'], 2)
            })

        # Save outputs
        output_files = {}
        if args.output_csv:
            pi_df.to_csv(args.output_csv, index=False)
            output_files['csv'] = args.output_csv

        if args.output_plot and len(pi_df) > 0:
            analyzer.plot_occupancy(
                pi_df=pi_df,
                output_path=args.output_plot
            )
            output_files['plot'] = args.output_plot

        return ToolResult(
            is_error=False,
            data={
                'total_pi_pi': total_pi_pi,
                'total_cation_pi': total_cation_pi,
                'persistent_pi_pi': len(persistent_pi_pi),
                'persistent_cation_pi': len(persistent_cation_pi),
                'top_5_pi_pi': top_pi_pi_list,
                'top_5_cation_pi': top_cation_pi_list,
                'output_files': output_files,
                'distance_cutoff': args.distance_cutoff,
                'angle_cutoff': args.angle_cutoff
            },
            message=(
                f"Pi interaction analysis complete. Found {total_pi_pi} pi-pi stacking, "
                f"{total_cation_pi} cation-pi interactions. "
                f"Persistent: {len(persistent_pi_pi)} pi-pi, {len(persistent_cation_pi)} cation-pi."
            )
        )

    def system_prompt_section(self) -> str:
        """Return system prompt section for this tool."""
        return """
## analyze_pi_interactions

Analyze pi-pi stacking and cation-pi interactions at molecular interfaces.

**When to use:**
- To identify aromatic interactions
- To understand recognition of aromatic peptide residues (PHE, TRP, TYR)
- To predict effects of aromatic mutations
- After hydrophobic analysis (aromatics contribute to both)

**Input requirements:**
- Preprocessed trajectory
- Topology and reference PDB
- Chain mapping for interface definition

**Interaction types detected:**
1. **Pi-pi stacking**: aromatic-aromatic interactions
   - Aromatics: PHE, TRP, TYR, HIS
   - Distance between ring centers < 6.0 Å
   - Angle between ring planes < 30° (parallel/T-shaped)

2. **Cation-pi**: charged-aromatic interactions
   - Cations: LYS (NZ), ARG (NH1, NH2), HIS (ND1, NE2)
   - Aromatics: PHE, TRP, TYR, HIS
   - Distance between cation and ring center < 6.0 Å

**Output interpretation:**
- Occupancy: % of frames where interaction is present
- Occupancy > 75%: Very stable, critical for recognition
- Occupancy 50-75%: Stable, important for specificity
- Occupancy 25-50%: Moderate, may contribute
- Occupancy < 25%: Transient, less important
- Mean distance: Average distance when interaction is formed

**Best practices:**
- Use default 6.0 Å distance and 30° angle cutoffs
- Pi-pi stacking can be parallel or T-shaped (perpendicular)
- Cation-pi is stronger than typical hydrophobic contacts
- In pHLA-TCR: aromatic peptide residues often recognized via pi interactions
- TRP is the strongest pi donor/acceptor
- Compare with hydrophobic contacts - aromatics contribute to both
"""
