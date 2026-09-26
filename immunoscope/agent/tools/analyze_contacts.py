"""
General Contact Analysis Tool

Analyzes all residue-residue contacts at molecular interfaces.
"""

import asyncio
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class AnalyzeContactsInput(BaseModel):
    """Input schema for general contact analysis."""

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
        description="Optional path to save contact results as CSV"
    )
    output_heatmap: str | None = Field(
        default=None,
        description="Optional path to save contact heatmap"
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
class AnalyzeContactsTool(Tool):
    """Tool for analyzing all residue-residue contacts at molecular interfaces."""

    name: ClassVar[str] = "analyze_contacts"
    description: ClassVar[str] = (
        "Analyze all residue-residue contacts at molecular interfaces. "
        "This provides a comprehensive view of interface contacts regardless of interaction type. "
        "Use this for initial interface characterization or when you need all contacts."
    )
    Input: ClassVar[type[BaseModel]] = AnalyzeContactsInput

    # Tool metadata
    is_read_only: ClassVar[bool] = False
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 600.0

    async def call(self, args: AnalyzeContactsInput, ctx: ToolContext) -> ToolResult:
        """Execute general contact analysis."""
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

        if args.distance_cutoff <= 0 or args.distance_cutoff > 10.0:
            return ToolResult(
                is_error=True,
                error=f"Distance cutoff must be between 0 and 10.0 Å, got {args.distance_cutoff}"
            )

        await ctx.on_event({
            "type": "progress",
            "message": f"Starting general contact analysis (d<{args.distance_cutoff}Å)..."
        })

        loop = asyncio.get_event_loop()
        try:
            result = await loop.run_in_executor(
                None,
                self._run_contact_analysis,
                args,
                ctx
            )
            return result
        except Exception as e:
            return ToolResult(
                is_error=True,
                error=f"Contact analysis failed: {str(e)}"
            )

    def _run_contact_analysis(
        self,
        args: AnalyzeContactsInput,
        ctx: ToolContext
    ) -> ToolResult:
        """Run general contact analysis (blocking operation)."""
        from immunoscope.analysis.topology.contact_annotation import InterfaceContactAnnotator

        # Initialize annotator
        annotator = InterfaceContactAnnotator(
            trajectory_path=args.trajectory,
            topology_path=args.topology,
            structure_pdb=args.structure_pdb,
            chain_mapping=args.chain_mapping
        )

        # Run analysis
        contact_df = annotator.analyze(
            distance_cutoff=args.distance_cutoff,
            stride=args.stride
        )

        # Get statistics
        total_contacts = len(contact_df)
        persistent_contacts = contact_df[contact_df['occupancy'] > 0.5]
        mean_occupancy = contact_df['occupancy'].mean() if len(contact_df) > 0 else 0.0

        # Analyze by region (if available)
        region_stats = {}
        if 'region1' in contact_df.columns and 'region2' in contact_df.columns:
            for region_pair in contact_df.groupby(['region1', 'region2']):
                region_name = f"{region_pair[0][0]}-{region_pair[0][1]}"
                region_df = region_pair[1]
                region_stats[region_name] = {
                    'total_contacts': len(region_df),
                    'persistent_contacts': len(region_df[region_df['occupancy'] > 0.5]),
                    'mean_occupancy': round(region_df['occupancy'].mean(), 3)
                }

        # Top contacts by occupancy
        top_contacts = contact_df.nlargest(10, 'occupancy') if len(contact_df) > 0 else contact_df
        top_contacts_list = []
        for _, row in top_contacts.iterrows():
            contact_info = {
                'residue1': f"{row['chain1']}:{row['resname1']}{row['resid1']}",
                'residue2': f"{row['chain2']}:{row['resname2']}{row['resid2']}",
                'occupancy': round(row['occupancy'], 3),
                'mean_distance': round(row['mean_distance'], 2)
            }
            if 'region1' in row and 'region2' in row:
                contact_info['regions'] = f"{row['region1']}-{row['region2']}"
            top_contacts_list.append(contact_info)

        # Save outputs
        output_files = {}
        if args.output_csv:
            contact_df.to_csv(args.output_csv, index=False)
            output_files['csv'] = args.output_csv

        if args.output_heatmap and len(contact_df) > 0:
            from immunoscope.analysis.topology.contact_heatmap import ContactHeatmapPlotter
            plotter = ContactHeatmapPlotter()
            plotter.plot_heatmap(
                contact_df=contact_df,
                output_path=args.output_heatmap
            )
            output_files['heatmap'] = args.output_heatmap

        return ToolResult(
            is_error=False,
            data={
                'total_contacts': total_contacts,
                'persistent_contacts': len(persistent_contacts),
                'mean_occupancy': round(mean_occupancy, 3),
                'region_statistics': region_stats,
                'top_10_contacts': top_contacts_list,
                'output_files': output_files,
                'distance_cutoff': args.distance_cutoff
            },
            message=(
                f"Contact analysis complete. Found {total_contacts} contacts, "
                f"{len(persistent_contacts)} persistent (>50% occupancy). "
                f"Mean occupancy: {mean_occupancy:.3f}."
            )
        )

    def system_prompt_section(self) -> str:
        """Return system prompt section for this tool."""
        return """
## analyze_contacts

Analyze all residue-residue contacts at molecular interfaces.

**When to use:**
- For initial interface characterization
- To get a comprehensive view of all contacts
- To identify contact hotspots before detailed interaction analysis
- When you need region-level contact statistics (CDR-peptide, etc.)

**Input requirements:**
- Preprocessed trajectory
- Topology and reference PDB
- Chain mapping for interface definition

**Contact detection:**
- Any heavy atom distance < 4.5 Å (default)
- All residue types included
- Minimum distance between any atoms of two residues

**Output interpretation:**
- Occupancy: % of frames where contact is present
- Occupancy > 75%: Very stable contact
- Occupancy 50-75%: Stable contact
- Occupancy 25-50%: Moderate contact
- Occupancy < 25%: Transient contact
- Mean distance: Average minimum distance when contact is formed

**Region-aware analysis:**
- If structure has region annotations (CDR1, CDR2, CDR3, peptide, MHC)
- Provides statistics per region pair (e.g., CDR3-peptide contacts)
- Useful for understanding which CDRs dominate binding

**Best practices:**
- Use this tool first for interface overview
- Follow up with specific interaction tools (H-bonds, salt bridges, etc.)
- Compare contact occupancy with interaction occupancy
- High contact occupancy without specific interactions → van der Waals/hydrophobic
- Heatmap visualization shows contact patterns across interface
- In pHLA-TCR: CDR3 typically has highest contact occupancy with peptide
"""
