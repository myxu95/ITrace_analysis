"""Calculate RMSF tool for ImmunoScope Agent."""

from __future__ import annotations
import asyncio
from pathlib import Path
from typing import ClassVar, Optional

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class CalculateRMSFInput(BaseModel):
    """Input schema for calculate_rmsf tool."""

    trajectory: str = Field(
        ...,
        description="Path to trajectory file (.xtc, .trr)",
    )
    topology: str = Field(
        ...,
        description="Path to topology file (.tpr, .gro, .pdb)",
    )
    structure_pdb: str = Field(
        ...,
        description="Path to reference structure PDB file for residue annotation",
    )
    chain_mapping: dict = Field(
        ...,
        description="Chain mapping dict, e.g., {'TCR_alpha': 'D', 'TCR_beta': 'E', 'HLA_alpha': 'A', 'peptide': 'P', 'beta2m': 'B'}",
    )
    cdr_detection: Optional[dict] = Field(
        default=None,
        description="Optional CDR detection dict with 'method' and parameters",
    )
    stride: int = Field(
        default=1,
        ge=1,
        description="Frame stride for analysis (default: 1 = all frames)",
    )
    selection: Optional[str] = Field(
        default=None,
        description="Optional atom selection (default: 'name CA' for mapped chains)",
    )
    output_csv: Optional[str] = Field(
        default=None,
        description="Optional output CSV file path for per-residue RMSF data",
    )
    output_plot: Optional[str] = Field(
        default=None,
        description="Optional output plot file path (.png, .pdf)",
    )


@register_tool
class CalculateRMSFTool(Tool):
    """Calculate Root Mean Square Fluctuation (RMSF) per residue.

    RMSF measures per-residue flexibility over the trajectory. It identifies
    flexible loops, rigid domains, and dynamic regions. For pHLA-TCR complexes,
    it provides region-aware analysis (CDRs, peptide, MHC helices, etc.).

    The tool requires:
    - Trajectory and topology files
    - Reference PDB structure for residue annotation
    - Chain mapping to identify complex components

    Returns per-residue RMSF values with biological region annotations and
    summary statistics for each region (CDRs, peptide, MHC helices, etc.).
    """

    name: ClassVar[str] = "calculate_rmsf"
    description: ClassVar[str] = (
        "Calculate Root Mean Square Fluctuation (RMSF) per residue to identify flexible and rigid regions. "
        "For pHLA-TCR complexes, provides region-aware analysis (CDRs, peptide, MHC helices). "
        "Returns per-residue RMSF values with biological annotations and region summaries. "
        "Use after RMSD to understand which regions are most dynamic."
    )
    Input: ClassVar[type[BaseModel]] = CalculateRMSFInput

    is_read_only: ClassVar[bool] = False  # May write output files
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 300.0  # 5 minutes for large trajectories

    async def call(self, args: CalculateRMSFInput, ctx: ToolContext) -> ToolResult:
        """Execute RMSF calculation."""
        try:
            # Validate inputs
            traj_path = Path(args.trajectory).resolve()
            topo_path = Path(args.topology).resolve()
            struct_path = Path(args.structure_pdb).resolve()

            if not traj_path.exists():
                return ToolResult(
                    content=f"Trajectory file not found: {args.trajectory}",
                    is_error=True,
                )
            if not topo_path.exists():
                return ToolResult(
                    content=f"Topology file not found: {args.topology}",
                    is_error=True,
                )
            if not struct_path.exists():
                return ToolResult(
                    content=f"Structure PDB file not found: {args.structure_pdb}",
                    is_error=True,
                )

            # Validate chain_mapping
            if not args.chain_mapping:
                return ToolResult(
                    content="chain_mapping is required and cannot be empty",
                    is_error=True,
                )

            # Prepare output paths
            output_csv = None
            output_plot = None
            if args.output_csv:
                output_csv = Path(args.output_csv).resolve()
                output_csv.parent.mkdir(parents=True, exist_ok=True)
            if args.output_plot:
                output_plot = Path(args.output_plot).resolve()
                output_plot.parent.mkdir(parents=True, exist_ok=True)

            # Report progress
            await ctx.on_event({
                "type": "progress",
                "message": f"Calculating RMSF for {traj_path.name}",
                "progress": 0.0,
            })

            # Run RMSF calculation in executor
            loop = asyncio.get_event_loop()
            result = await loop.run_in_executor(
                None,
                self._run_rmsf_calculation,
                str(traj_path),
                str(topo_path),
                str(struct_path),
                args.chain_mapping,
                args.cdr_detection,
                args.stride,
                args.selection,
                str(output_csv) if output_csv else None,
                str(output_plot) if output_plot else None,
            )

            if result["success"]:
                rmsf_data = result["data"]

                return ToolResult(content={
                    "success": True,
                    "trajectory": str(traj_path),
                    "topology": str(topo_path),
                    "structure_pdb": str(struct_path),
                    "n_residues": rmsf_data["n_residues"],
                    "n_frames": rmsf_data["n_frames"],
                    "stride": args.stride,
                    "mean_rmsf_angstrom": rmsf_data["mean_rmsf"],
                    "max_rmsf_angstrom": rmsf_data["max_rmsf"],
                    "tcr_mean_rmsf_angstrom": rmsf_data.get("tcr_mean_rmsf"),
                    "phla_mean_rmsf_angstrom": rmsf_data.get("phla_mean_rmsf"),
                    "region_summary": rmsf_data.get("region_summary", {}),
                    "output_csv": str(output_csv) if output_csv else None,
                    "output_plot": str(output_plot) if output_plot else None,
                    "message": self._format_rmsf_message(rmsf_data),
                })
            else:
                return ToolResult(
                    content=f"RMSF calculation failed: {result['error']}",
                    is_error=True,
                )

        except Exception as e:
            return ToolResult(
                content=f"RMSF calculation failed: {str(e)}",
                is_error=True,
            )

    def _run_rmsf_calculation(
        self,
        trajectory: str,
        topology: str,
        structure_pdb: str,
        chain_mapping: dict,
        cdr_detection: Optional[dict],
        stride: int,
        selection: Optional[str],
        output_csv: Optional[str],
        output_plot: Optional[str],
    ) -> dict:
        """Run RMSF calculation (blocking operation)."""
        try:
            from immunoscope.analysis.trajectory.residue_rmsf import ResidueRMSFAnalyzer

            # Create analyzer
            analyzer = ResidueRMSFAnalyzer(
                topology_file=topology,
                trajectory_file=trajectory,
                structure_pdb=structure_pdb,
            )

            # Calculate RMSF
            result = analyzer.calculate(
                chain_mapping=chain_mapping,
                cdr_detection=cdr_detection,
                stride=stride,
                selection=selection,
            )

            # Save CSV if requested
            if output_csv:
                result.residue_frame.to_csv(output_csv, index=False)

            # Generate plot if requested
            if output_plot:
                analyzer.plot(result, output_file=output_plot)

            # Extract region summary
            region_summary = {}
            if not result.region_summary.empty:
                for _, row in result.region_summary.iterrows():
                    region_summary[row['region_group']] = {
                        'mean_rmsf': round(float(row['mean_rmsf_angstrom']), 3),
                        'std_rmsf': round(float(row['std_rmsf_angstrom']), 3),
                        'n_residues': int(row['n_residues']),
                    }

            return {
                "success": True,
                "data": {
                    "n_residues": result.summary["n_residues"],
                    "n_frames": result.summary["n_frames"],
                    "mean_rmsf": round(result.summary["mean_rmsf_angstrom"], 3),
                    "max_rmsf": round(result.summary["max_rmsf_angstrom"], 3),
                    "tcr_mean_rmsf": round(result.summary.get("tcr_mean_rmsf_angstrom", 0.0), 3),
                    "phla_mean_rmsf": round(result.summary.get("phla_mean_rmsf_angstrom", 0.0), 3),
                    "region_summary": region_summary,
                },
            }

        except Exception as e:
            return {
                "success": False,
                "error": str(e),
            }

    def _format_rmsf_message(self, data: dict) -> str:
        """Format RMSF result message."""
        mean = data["mean_rmsf"]
        max_rmsf = data["max_rmsf"]
        n_residues = data["n_residues"]

        msg = f"Analyzed {n_residues} residues. "
        msg += f"Mean RMSF: {mean:.2f} Å, Max RMSF: {max_rmsf:.2f} Å. "

        # Add component-specific info if available
        if data.get("tcr_mean_rmsf") and data.get("phla_mean_rmsf"):
            tcr_mean = data["tcr_mean_rmsf"]
            phla_mean = data["phla_mean_rmsf"]
            msg += f"TCR: {tcr_mean:.2f} Å, pHLA: {phla_mean:.2f} Å. "

        # Highlight most flexible regions
        region_summary = data.get("region_summary", {})
        if region_summary:
            sorted_regions = sorted(
                region_summary.items(),
                key=lambda x: x[1]['mean_rmsf'],
                reverse=True
            )[:3]
            if sorted_regions:
                top_regions = ", ".join([
                    f"{region} ({vals['mean_rmsf']:.2f} Å)"
                    for region, vals in sorted_regions
                ])
                msg += f"Most flexible: {top_regions}."

        return msg

    def system_prompt_section(self) -> str:
        """Tool-specific guidance for the agent."""
        return """
## RMSF Calculation Guidelines

**When to calculate RMSF:**
- After RMSD to understand which regions contribute to structural variation
- To identify flexible loops, CDRs, or binding sites
- When user asks about residue flexibility or dynamic regions

**Interpreting RMSF values (for proteins):**
- < 1.0 Å: Rigid regions (core, secondary structures)
- 1.0-2.0 Å: Moderate flexibility (typical for structured regions)
- 2.0-3.0 Å: High flexibility (loops, termini)
- > 3.0 Å: Very high flexibility (disordered regions, long loops)

**For pHLA-TCR complexes:**
- CDR3 loops typically show highest RMSF (2-4 Å)
- Peptide flexibility varies (1-3 Å depending on binding)
- MHC helices are usually rigid (< 1.5 Å)
- Compare TCR vs pHLA flexibility to assess binding stability

**Required inputs:**
- chain_mapping: Must map all components (TCR_alpha, TCR_beta, HLA_alpha, peptide, beta2m)
- structure_pdb: Reference structure for residue annotation
- cdr_detection: Optional, use if CDR boundaries are known

**Best practices:**
- Use stride > 1 for large trajectories to speed up calculation
- Save CSV output for detailed per-residue analysis
- Generate plot to visualize flexibility patterns
- Compare RMSF across different simulations or conditions
"""
