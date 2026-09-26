"""
Docking Angle Analysis Tool

Calculates TCR-pMHC docking angles (crossing angle and tilt angle).
"""

import asyncio
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class AnalyzeAnglesInput(BaseModel):
    """Input schema for docking angle analysis."""

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
        description="Mapping of biological roles to chain IDs (mhc_alpha, tcr_alpha, tcr_beta)"
    )
    output_csv: str | None = Field(
        default=None,
        description="Optional path to save angle time series as CSV"
    )
    output_plot: str | None = Field(
        default=None,
        description="Optional path to save angle distribution plot"
    )
    stride: int = Field(
        default=10,
        description="Frame stride for analysis (default: 10)"
    )


@register_tool
class AnalyzeAnglesTool(Tool):
    """Tool for analyzing TCR-pMHC docking angles."""

    name: ClassVar[str] = "analyze_angles"
    description: ClassVar[str] = (
        "Calculate TCR-pMHC docking angles (crossing angle and tilt angle). "
        "The crossing angle measures the rotation of TCR relative to the MHC groove axis. "
        "The tilt angle measures the inclination of TCR relative to the MHC binding plane. "
        "These angles characterize the TCR binding geometry."
    )
    Input: ClassVar[type[BaseModel]] = AnalyzeAnglesInput

    # Tool metadata
    is_read_only: ClassVar[bool] = False
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 600.0

    async def call(self, args: AnalyzeAnglesInput, ctx: ToolContext) -> ToolResult:
        """Execute docking angle analysis."""
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
        required_chains = ['mhc_alpha', 'tcr_alpha', 'tcr_beta']
        missing = [c for c in required_chains if c not in args.chain_mapping]
        if missing:
            return ToolResult(
                is_error=True,
                error=f"Missing required chains in chain_mapping: {missing}"
            )

        await ctx.on_event({
            "type": "progress",
            "message": "Starting docking angle analysis..."
        })

        loop = asyncio.get_event_loop()
        try:
            result = await loop.run_in_executor(
                None,
                self._run_angle_analysis,
                args,
                ctx
            )
            return result
        except Exception as e:
            return ToolResult(
                is_error=True,
                content=f"Angle analysis failed: {str(e)}",
            )

    def _run_angle_analysis(
        self,
        args: AnalyzeAnglesInput,
        ctx: ToolContext
    ) -> ToolResult:
        """Run angle analysis (blocking operation)."""
        from immunoscope.analysis.angles.analyzer import DockingAngleAnalyzer
        from immunoscope.analysis.angles.angle_data_structures import DockingAngleInput
        import numpy as np

        # Prepare input
        angle_input = DockingAngleInput(
            trajectory_path=args.trajectory,
            topology_path=args.topology,
            structure_pdb=args.structure_pdb,
            mhc_alpha_chain=args.chain_mapping['mhc_alpha'],
            tcr_alpha_chain=args.chain_mapping['tcr_alpha'],
            tcr_beta_chain=args.chain_mapping['tcr_beta'],
            stride=args.stride
        )

        # Run analysis
        analyzer = DockingAngleAnalyzer()
        result = analyzer.calculate(angle_input)

        # Extract statistics
        crossing_angles = result.crossing_angles
        tilt_angles = result.tilt_angles

        mean_crossing = float(np.mean(crossing_angles))
        std_crossing = float(np.std(crossing_angles))
        mean_tilt = float(np.mean(tilt_angles))
        std_tilt = float(np.std(tilt_angles))

        # Save outputs
        output_files = {}
        if args.output_csv:
            result.to_dataframe().to_csv(args.output_csv, index=False)
            output_files['csv'] = args.output_csv

        if args.output_plot:
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            # Plot 1: Crossing angle distribution
            ax1 = axes[0]
            ax1.hist(crossing_angles, bins=30, alpha=0.7, edgecolor='black')
            ax1.axvline(mean_crossing, color='red', linestyle='--',
                       label=f'Mean: {mean_crossing:.1f}°')
            ax1.set_xlabel('Crossing Angle (°)')
            ax1.set_ylabel('Frequency')
            ax1.set_title('TCR Crossing Angle Distribution')
            ax1.legend()
            ax1.grid(alpha=0.3)

            # Plot 2: Tilt angle distribution
            ax2 = axes[1]
            ax2.hist(tilt_angles, bins=30, alpha=0.7, edgecolor='black', color='orange')
            ax2.axvline(mean_tilt, color='red', linestyle='--',
                       label=f'Mean: {mean_tilt:.1f}°')
            ax2.set_xlabel('Tilt Angle (°)')
            ax2.set_ylabel('Frequency')
            ax2.set_title('TCR Tilt Angle Distribution')
            ax2.legend()
            ax2.grid(alpha=0.3)

            plt.tight_layout()
            plt.savefig(args.output_plot, dpi=300, bbox_inches='tight')
            plt.close()
            output_files['plot'] = args.output_plot

        # Build a structured Markdown content payload (ToolResult only has
        # `content`, `is_error`, `truncated`, `saved_to` — no data/message fields).
        lines = [
            f"## Docking angle analysis — {len(crossing_angles)} frames",
            "",
            f"- **Crossing angle**: {mean_crossing:.2f}° ± {std_crossing:.2f}° "
            f"(range {float(np.min(crossing_angles)):.1f}°–{float(np.max(crossing_angles)):.1f}°)",
            f"- **Tilt angle**: {mean_tilt:.2f}° ± {std_tilt:.2f}° "
            f"(range {float(np.min(tilt_angles)):.1f}°–{float(np.max(tilt_angles)):.1f}°)",
        ]
        if output_files:
            lines.append("")
            lines.append("**Outputs:**")
            for k, v in output_files.items():
                lines.append(f"- {k}: `{v}`")

        return ToolResult(
            is_error=False,
            content="\n".join(lines),
        )

    def system_prompt_section(self) -> str:
        """Return system prompt section for this tool."""
        return """
## analyze_angles

Calculate TCR-pMHC docking angles (crossing angle and tilt angle).

**When to use:**
- To characterize TCR binding geometry
- To compare binding modes between different TCRs
- To understand TCR orientation relative to pMHC
- After structural analysis to quantify binding pose

**Input requirements:**
- Preprocessed trajectory
- Topology and reference PDB
- Chain mapping (mhc_alpha, tcr_alpha, tcr_beta required)

**Angle definitions:**

1. **Crossing Angle**: Rotation of TCR around MHC groove axis
   - Measures how TCR is rotated relative to peptide N→C direction
   - Range: 0-180°
   - Typical values: 20-70° (most TCRs bind diagonally)
   - 0°: TCR aligned with peptide N→C
   - 90°: TCR perpendicular to peptide
   - Canonical diagonal: ~45°

2. **Tilt Angle**: Inclination of TCR relative to MHC plane
   - Measures how much TCR tilts toward/away from MHC surface
   - Range: 0-90°
   - Typical values: 10-30°
   - 0°: TCR parallel to MHC surface
   - 90°: TCR perpendicular to MHC surface
   - Most TCRs: 15-25° (slight tilt)

**Output interpretation:**
- **Mean angles**: Average binding geometry
- **Std deviation**: Conformational flexibility
  - Low std (<5°): Rigid binding geometry
  - Moderate std (5-15°): Some flexibility
  - High std (>15°): Highly dynamic binding
- **Angle distributions**: Reveal conformational states
  - Single peak: One dominant geometry
  - Multiple peaks: Multiple binding modes

**Biological significance:**
- Crossing angle determines which CDRs contact peptide vs MHC
- Tilt angle affects binding interface size
- Conserved angles suggest conserved recognition mode
- Large angle changes indicate conformational transitions

**Best practices:**
- Use stride=10 for large trajectories
- Compare with crystal structure angles (if available)
- Correlate angle changes with RMSD/RMSF
- Large angle fluctuations may indicate unstable binding
- In pHLA-TCR: crossing angle ~45° is most common (diagonal binding)
"""
