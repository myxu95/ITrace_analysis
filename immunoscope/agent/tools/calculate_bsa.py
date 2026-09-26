"""
Buried Surface Area (BSA) Analysis Tool

Calculates interface buried surface area between molecular components.
"""

import asyncio
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class CalculateBSAInput(BaseModel):
    """Input schema for BSA calculation."""

    trajectory: str = Field(
        description="Path to trajectory file (.xtc, .trr)"
    )
    topology: str = Field(
        description="Path to topology file (.tpr, .gro, .pdb)"
    )
    structure_pdb: str = Field(
        description="Path to reference PDB for chain annotation"
    )
    chain_mapping: dict[str, str] = Field(
        description="Mapping of biological roles to chain IDs"
    )
    output_csv: str | None = Field(
        default=None,
        description="Optional path to save BSA time series as CSV"
    )
    output_plot: str | None = Field(
        default=None,
        description="Optional path to save BSA plot"
    )
    stride: int = Field(
        default=1,
        description="Frame stride for analysis (default: 1)"
    )


@register_tool
class CalculateBSATool(Tool):
    """Tool for calculating buried surface area at molecular interfaces."""

    name: ClassVar[str] = "calculate_bsa"
    description: ClassVar[str] = (
        "Calculate buried surface area (BSA) at the molecular interface. "
        "BSA quantifies the surface area hidden from solvent upon complex formation. "
        "Larger BSA typically indicates stronger binding and more extensive interfaces."
    )
    Input: ClassVar[type[BaseModel]] = CalculateBSAInput

    # Tool metadata
    is_read_only: ClassVar[bool] = False
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 600.0

    async def call(self, args: CalculateBSAInput, ctx: ToolContext) -> ToolResult:
        """Execute BSA calculation."""
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

        await ctx.on_event({
            "type": "progress",
            "message": "Starting BSA calculation..."
        })

        loop = asyncio.get_event_loop()
        try:
            result = await loop.run_in_executor(
                None,
                self._run_bsa_calculation,
                args,
                ctx
            )
            return result
        except Exception as e:
            return ToolResult(
                is_error=True,
                error=f"BSA calculation failed: {str(e)}"
            )

    def _run_bsa_calculation(
        self,
        args: CalculateBSAInput,
        ctx: ToolContext
    ) -> ToolResult:
        """Run BSA calculation (blocking operation)."""
        from immunoscope.analysis.interface.buried_surface_area import BuriedSurfaceAreaCalculator

        # Initialize calculator
        calculator = BuriedSurfaceAreaCalculator(
            trajectory_path=args.trajectory,
            topology_path=args.topology,
            structure_pdb=args.structure_pdb,
            chain_mapping=args.chain_mapping
        )

        # Run calculation
        bsa_results = calculator.calculate(stride=args.stride)

        # Extract statistics
        mean_bsa = bsa_results['bsa_total'].mean()
        std_bsa = bsa_results['bsa_total'].std()
        min_bsa = bsa_results['bsa_total'].min()
        max_bsa = bsa_results['bsa_total'].max()

        # Component breakdown (if available)
        component_bsa = {}
        for col in bsa_results.columns:
            if col.startswith('bsa_') and col != 'bsa_total':
                component_name = col.replace('bsa_', '')
                component_bsa[component_name] = round(bsa_results[col].mean(), 2)

        # Save outputs
        output_files = {}
        if args.output_csv:
            bsa_results.to_csv(args.output_csv, index=False)
            output_files['csv'] = args.output_csv

        if args.output_plot:
            calculator.plot_timeseries(
                bsa_results=bsa_results,
                output_path=args.output_plot
            )
            output_files['plot'] = args.output_plot

        return ToolResult(
            is_error=False,
            data={
                'mean_bsa': round(mean_bsa, 2),
                'std_bsa': round(std_bsa, 2),
                'min_bsa': round(min_bsa, 2),
                'max_bsa': round(max_bsa, 2),
                'component_bsa': component_bsa,
                'frames_analyzed': len(bsa_results),
                'output_files': output_files
            },
            message=(
                f"BSA calculation complete. Mean BSA: {mean_bsa:.2f} ± {std_bsa:.2f} Ų. "
                f"Range: {min_bsa:.2f} - {max_bsa:.2f} Ų."
            )
        )

    def system_prompt_section(self) -> str:
        """Return system prompt section for this tool."""
        return """
## calculate_bsa

Calculate buried surface area (BSA) at molecular interfaces.

**When to use:**
- To quantify interface size and binding strength
- To assess interface stability over time
- To compare binding modes between different complexes
- After preprocessing, often alongside RRCS analysis

**Input requirements:**
- Preprocessed trajectory
- Topology and reference PDB
- Chain mapping for interface definition

**Output interpretation:**
- BSA in Ų (square Angstroms)
- Typical pHLA-TCR interfaces: 1200-2000 Ų
- Typical antibody-antigen: 1400-2200 Ų
- Typical protein-protein: 1000-3000 Ų
- BSA < 800 Ų: Weak or transient interaction
- BSA > 2500 Ų: Very extensive interface
- Stable BSA over time indicates stable binding
- Large BSA fluctuations suggest conformational changes

**Best practices:**
- Use stride=10 for quick assessment of large trajectories
- Compare BSA with RRCS to understand interface composition
- Monitor BSA stability to assess equilibration
- Component breakdown shows which chains contribute most to interface
"""
