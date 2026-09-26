"""
Hydrophobic Contact Analysis Tool

Analyzes hydrophobic interactions at molecular interfaces.
"""

import asyncio
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class AnalyzeHydrophobicInput(BaseModel):
    """Input schema for hydrophobic contact analysis."""

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
        description="Optional path to save hydrophobic contact results as CSV"
    )
    output_plot: str | None = Field(
        default=None,
        description="Optional path to save hydrophobic contact occupancy plot"
    )
    distance_cutoff: float = Field(
        default=4.5,
        description="Distance cutoff in Angstroms (default: 4.5)"
    )
    stride: int = Field(
        default=1,
        description="Frame stride for analysis (default: 1)"
    )


@register_tool
class AnalyzeHydrophobicTool(Tool):
    """Tool for analyzing hydrophobic contacts at molecular interfaces."""

    name: ClassVar[str] = "analyze_hydrophobic"
    description: ClassVar[str] = (
        "Analyze hydrophobic interactions between nonpolar residues. "
        "Hydrophobic contacts are critical for protein-protein binding, "
        "contributing to the hydrophobic effect that drives association."
    )
    Input: ClassVar[type[BaseModel]] = AnalyzeHydrophobicInput

    # Tool metadata
    is_read_only: ClassVar[bool] = False
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 600.0

    async def call(self, args: AnalyzeHydrophobicInput, ctx: ToolContext) -> ToolResult:
        """Execute hydrophobic contact analysis."""
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
            "message": f"Starting hydrophobic contact analysis (d<{args.distance_cutoff}Å)..."
        })

        loop = asyncio.get_event_loop()
        try:
            result = await loop.run_in_executor(
                None,
                self._run_hydrophobic_analysis,
                args,
                ctx
            )
            return result
        except Exception as e:
            return ToolResult(
                is_error=True,
                error=f"Hydrophobic contact analysis failed: {str(e)}"
            )

    def _run_hydrophobic_analysis(
        self,
        args: AnalyzeHydrophobicInput,
        ctx: ToolContext
    ) -> ToolResult:
        """Run hydrophobic contact analysis (blocking operation)."""
        from immunoscope.analysis.interactions.hydrophobic_contact_pairs import HydrophobicContactAnalyzer

        # Initialize analyzer
        analyzer = HydrophobicContactAnalyzer(
            trajectory_path=args.trajectory,
            topology_path=args.topology,
            structure_pdb=args.structure_pdb,
            chain_mapping=args.chain_mapping
        )

        # Run analysis
        hc_df = analyzer.analyze(
            distance_cutoff=args.distance_cutoff,
            stride=args.stride
        )

        # Get statistics
        total_contacts = len(hc_df)
        persistent_contacts = hc_df[hc_df['occupancy'] > 0.5]
        mean_occupancy = hc_df['occupancy'].mean() if len(hc_df) > 0 else 0.0

        # Top contacts by occupancy
        top_contacts = hc_df.nlargest(10, 'occupancy') if len(hc_df) > 0 else hc_df
        top_contacts_list = []
        for _, row in top_contacts.iterrows():
            top_contacts_list.append({
                'residue1': f"{row['chain1']}:{row['resname1']}{row['resid1']}",
                'residue2': f"{row['chain2']}:{row['resname2']}{row['resid2']}",
                'occupancy': round(row['occupancy'], 3),
                'mean_distance': round(row['mean_distance'], 2)
            })

        # Save outputs
        output_files = {}
        if args.output_csv:
            hc_df.to_csv(args.output_csv, index=False)
            output_files['csv'] = args.output_csv

        if args.output_plot and len(hc_df) > 0:
            analyzer.plot_occupancy(
                hc_df=hc_df,
                output_path=args.output_plot
            )
            output_files['plot'] = args.output_plot

        return ToolResult(
            is_error=False,
            data={
                'total_contacts': total_contacts,
                'persistent_contacts': len(persistent_contacts),
                'mean_occupancy': round(mean_occupancy, 3),
                'top_10_contacts': top_contacts_list,
                'output_files': output_files,
                'distance_cutoff': args.distance_cutoff
            },
            message=(
                f"Hydrophobic contact analysis complete. Found {total_contacts} contacts, "
                f"{len(persistent_contacts)} persistent (>50% occupancy). "
                f"Mean occupancy: {mean_occupancy:.3f}."
            )
        )

    def system_prompt_section(self) -> str:
        """Return system prompt section for this tool."""
        return """
## analyze_hydrophobic

Analyze hydrophobic interactions at molecular interfaces.

**When to use:**
- To identify nonpolar contacts driving binding
- To understand the hydrophobic effect contribution
- To predict effects of hydrophobic → polar mutations
- After polar interaction analysis (H-bonds, salt bridges)

**Input requirements:**
- Preprocessed trajectory
- Topology and reference PDB
- Chain mapping for interface definition

**Hydrophobic residues detected:**
- ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO
- Distance between carbon atoms < 4.5 Å (default)
- Uses carbon-carbon distances (not Cα-Cα)

**Output interpretation:**
- Occupancy: % of frames where contact is present
- Occupancy > 75%: Very stable hydrophobic core
- Occupancy 50-75%: Stable, important for binding
- Occupancy 25-50%: Moderate, may contribute
- Occupancy < 25%: Transient, less critical
- Mean distance: Average C-C distance when contact is formed

**Best practices:**
- Use default 4.5 Å cutoff for carbon-carbon contacts
- Hydrophobic contacts are entropy-driven (desolvation)
- Buried hydrophobic contacts are more stable
- Compare with BSA - hydrophobic contacts correlate with buried area
- In pHLA-TCR: peptide hydrophobic anchors are critical
- Aromatic residues (PHE, TRP, TYR) can form both hydrophobic and pi interactions
"""
