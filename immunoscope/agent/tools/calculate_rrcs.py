"""
RRCS (Residue-Residue Contact Score) Analysis Tool

Calculates interaction hotspots using RRCS methodology.
"""

import asyncio
from pathlib import Path
from typing import ClassVar, Any

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class CalculateRRCSInput(BaseModel):
    """Input schema for RRCS calculation."""

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
        description="Mapping of biological roles to chain IDs, e.g. {'tcr_alpha': 'D', 'tcr_beta': 'E', 'peptide': 'C', 'mhc_alpha': 'A', 'mhc_beta': 'B'}"
    )
    output_csv: str | None = Field(
        default=None,
        description="Optional path to save RRCS results as CSV"
    )
    output_plot: str | None = Field(
        default=None,
        description="Optional path to save RRCS heatmap plot"
    )
    distance_cutoff: float = Field(
        default=4.5,
        description="Distance cutoff in Angstroms for contact detection (default: 4.5)"
    )
    stride: int = Field(
        default=1,
        description="Frame stride for analysis (default: 1, analyze every frame)"
    )


@register_tool
class CalculateRRCSTool(Tool):
    """Tool for calculating RRCS interaction hotspots."""

    name: ClassVar[str] = "calculate_rrcs"
    description: ClassVar[str] = (
        "Calculate RRCS (Residue-Residue Contact Score) to identify interaction hotspots. "
        "RRCS quantifies the strength and persistence of residue-residue contacts across the trajectory. "
        "Higher RRCS values indicate more important interaction sites."
    )
    Input: ClassVar[type[BaseModel]] = CalculateRRCSInput

    # Tool metadata
    is_read_only: ClassVar[bool] = False  # Writes output files
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 600.0  # 10 minutes for large trajectories

    async def call(self, args: CalculateRRCSInput, ctx: ToolContext) -> ToolResult:
        """Execute RRCS calculation."""
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

        if args.distance_cutoff <= 0:
            return ToolResult(
                is_error=True,
                error=f"Distance cutoff must be positive, got {args.distance_cutoff}"
            )

        # Report start
        await ctx.on_event({
            "type": "progress",
            "message": f"Starting RRCS calculation with {args.distance_cutoff}Å cutoff..."
        })

        # Run calculation in executor to avoid blocking
        loop = asyncio.get_event_loop()
        try:
            result = await loop.run_in_executor(
                None,
                self._run_rrcs_calculation,
                args,
                ctx
            )
            return result
        except Exception as e:
            return ToolResult(
                is_error=True,
                error=f"RRCS calculation failed: {str(e)}"
            )

    def _run_rrcs_calculation(
        self,
        args: CalculateRRCSInput,
        ctx: ToolContext
    ) -> ToolResult:
        """Run RRCS calculation (blocking operation)."""
        from immunoscope.analysis.interactions.rrcs import RRCSCalculator

        # Initialize calculator
        calculator = RRCSCalculator(
            trajectory_path=args.trajectory,
            topology_path=args.topology,
            structure_pdb=args.structure_pdb,
            chain_mapping=args.chain_mapping
        )

        # Run calculation
        rrcs_df = calculator.calculate(
            distance_cutoff=args.distance_cutoff,
            stride=args.stride
        )

        # Get statistics
        top_contacts = rrcs_df.nlargest(10, 'rrcs_score')
        mean_rrcs = rrcs_df['rrcs_score'].mean()
        max_rrcs = rrcs_df['rrcs_score'].max()

        # Save outputs if requested
        output_files = {}
        if args.output_csv:
            rrcs_df.to_csv(args.output_csv, index=False)
            output_files['csv'] = args.output_csv

        if args.output_plot:
            calculator.plot_heatmap(
                rrcs_df=rrcs_df,
                output_path=args.output_plot
            )
            output_files['plot'] = args.output_plot

        # Format top contacts for output
        top_contacts_list = []
        for _, row in top_contacts.iterrows():
            top_contacts_list.append({
                'residue1': f"{row['chain1']}:{row['resname1']}{row['resid1']}",
                'residue2': f"{row['chain2']}:{row['resname2']}{row['resid2']}",
                'rrcs_score': round(row['rrcs_score'], 3),
                'occupancy': round(row['occupancy'], 3)
            })

        return ToolResult(
            is_error=False,
            data={
                'total_contacts': len(rrcs_df),
                'mean_rrcs': round(mean_rrcs, 3),
                'max_rrcs': round(max_rrcs, 3),
                'top_10_contacts': top_contacts_list,
                'output_files': output_files,
                'distance_cutoff': args.distance_cutoff,
                'frames_analyzed': len(rrcs_df) // args.stride if args.stride > 1 else 'all'
            },
            message=(
                f"RRCS calculation complete. Found {len(rrcs_df)} contacts. "
                f"Mean RRCS: {mean_rrcs:.3f}, Max RRCS: {max_rrcs:.3f}. "
                f"Top contact: {top_contacts_list[0]['residue1']} - {top_contacts_list[0]['residue2']} "
                f"(RRCS: {top_contacts_list[0]['rrcs_score']:.3f})"
            )
        )

    def system_prompt_section(self) -> str:
        """Return system prompt section for this tool."""
        return """
## calculate_rrcs

Calculate RRCS (Residue-Residue Contact Score) to identify interaction hotspots.

**When to use:**
- After preprocessing to identify key interaction sites
- To find hotspot residues that drive binding
- To compare interaction patterns between simulations
- To validate experimental mutagenesis data

**Input requirements:**
- Preprocessed trajectory (PBC-corrected)
- Topology file
- Reference PDB with correct chain IDs
- Chain mapping dictionary (biological roles → chain IDs)

**Output interpretation:**
- RRCS score: Combined metric of contact occupancy and distance
- Higher RRCS = more important interaction
- RRCS > 0.5: Strong, persistent contact (hotspot candidate)
- RRCS 0.3-0.5: Moderate interaction
- RRCS < 0.3: Weak or transient contact
- Top 10 contacts typically represent the binding interface core

**Best practices:**
- Use default 4.5Å cutoff for most analyses
- For large trajectories, use stride=10 to speed up calculation
- Always save CSV output for downstream analysis
- Generate heatmap plot for visual inspection
- Compare RRCS patterns across different conditions to identify changes

**Common patterns in pHLA-TCR:**
- CDR3 loops typically have highest RRCS with peptide
- CDR1/CDR2 interact more with MHC helices
- Peptide central residues (P4-P6) often show high RRCS
- MHC anchor pockets (P2, P9) show strong contacts with peptide
"""
