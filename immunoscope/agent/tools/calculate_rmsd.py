"""Calculate RMSD tool for ImmunoScope Agent."""

from __future__ import annotations
import asyncio
from pathlib import Path
from typing import ClassVar, Optional

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class CalculateRMSDInput(BaseModel):
    """Input schema for calculate_rmsd tool."""

    trajectory: str = Field(
        ...,
        description="Path to trajectory file (.xtc, .trr)",
    )
    topology: str = Field(
        ...,
        description="Path to topology file (.tpr, .gro, .pdb)",
    )
    selection: str = Field(
        default="protein and name CA",
        description="Atom selection string (MDAnalysis syntax, default: 'protein and name CA')",
    )
    reference_frame: int = Field(
        default=0,
        ge=0,
        description="Reference frame index (default: 0 = first frame)",
    )
    output_file: Optional[str] = Field(
        default=None,
        description="Optional output file path for RMSD data (.csv or .xvg)",
    )
    method: str = Field(
        default="mdanalysis",
        description="Calculation method: 'mdanalysis' or 'gromacs' (default: mdanalysis)",
    )


@register_tool
class CalculateRMSDTool(Tool):
    """Calculate Root Mean Square Deviation (RMSD) of trajectory.

    RMSD measures structural deviation from a reference structure over time.
    It's a key metric for assessing simulation stability and convergence.

    Common selections:
    - 'protein and name CA': Backbone C-alpha atoms (default, most common)
    - 'protein': All protein atoms
    - 'backbone': Backbone atoms (N, CA, C, O)
    - 'resid 1-50': Specific residue range
    - 'segid A': Specific chain/segment

    The tool returns RMSD statistics and optionally saves time-series data.
    """

    name: ClassVar[str] = "calculate_rmsd"
    description: ClassVar[str] = (
        "Calculate Root Mean Square Deviation (RMSD) of trajectory relative to reference frame. "
        "RMSD measures structural stability over time. Returns mean, std, min, max RMSD values "
        "and optionally saves time-series data. Use after preprocessing to assess simulation convergence."
    )
    Input: ClassVar[type[BaseModel]] = CalculateRMSDInput

    is_read_only: ClassVar[bool] = False  # May write output file
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 300.0  # 5 minutes for large trajectories

    async def call(self, args: CalculateRMSDInput, ctx: ToolContext) -> ToolResult:
        """Execute RMSD calculation."""
        try:
            # Validate inputs
            traj_path = Path(args.trajectory).resolve()
            topo_path = Path(args.topology).resolve()

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

            # Prepare output file path
            output_path = None
            if args.output_file:
                output_path = Path(args.output_file).resolve()
                output_path.parent.mkdir(parents=True, exist_ok=True)

            # Report progress
            await ctx.on_event({
                "type": "progress",
                "message": f"Calculating RMSD for {traj_path.name}",
                "progress": 0.0,
            })

            # Run RMSD calculation in executor
            loop = asyncio.get_event_loop()
            result = await loop.run_in_executor(
                None,
                self._run_rmsd_calculation,
                str(traj_path),
                str(topo_path),
                args.selection,
                args.reference_frame,
                str(output_path) if output_path else None,
                args.method,
            )

            if result["success"]:
                rmsd_data = result["data"]

                return ToolResult(content={
                    "success": True,
                    "trajectory": str(traj_path),
                    "topology": str(topo_path),
                    "selection": args.selection,
                    "reference_frame": args.reference_frame,
                    "n_frames": rmsd_data["n_frames"],
                    "mean_rmsd_nm": rmsd_data["mean_rmsd"],
                    "std_rmsd_nm": rmsd_data["std_rmsd"],
                    "min_rmsd_nm": rmsd_data["min_rmsd"],
                    "max_rmsd_nm": rmsd_data["max_rmsd"],
                    "output_file": rmsd_data.get("output_file"),
                    "message": self._format_rmsd_message(rmsd_data),
                })
            else:
                return ToolResult(
                    content=f"RMSD calculation failed: {result['error']}",
                    is_error=True,
                )

        except Exception as e:
            return ToolResult(
                content=f"RMSD calculation failed: {str(e)}",
                is_error=True,
            )

    def _run_rmsd_calculation(
        self,
        trajectory: str,
        topology: str,
        selection: str,
        reference_frame: int,
        output_file: Optional[str],
        method: str,
    ) -> dict:
        """Run RMSD calculation (blocking operation)."""
        try:
            from immunoscope.analysis.trajectory.rmsd_refactored import (
                RMSDCalculator,
                RMSDInput,
            )

            # Create input
            rmsd_input = RMSDInput(
                topology=topology,
                trajectory=trajectory,
                selection=selection,
                reference_frame=reference_frame,
                output_file=output_file,
                method=method,
            )

            # Calculate RMSD
            calculator = RMSDCalculator()
            result = calculator.calculate(rmsd_input)

            if not result.success:
                return {
                    "success": False,
                    "error": result.error_message or "Unknown error",
                }

            # Save output if requested
            if output_file and result.times is not None and result.rmsd_values is not None:
                import pandas as pd
                df = pd.DataFrame({
                    'time_ps': result.times,
                    'rmsd_nm': result.rmsd_values,
                })

                if output_file.endswith('.csv'):
                    df.to_csv(output_file, index=False)
                elif output_file.endswith('.xvg'):
                    # GROMACS xvg format
                    with open(output_file, 'w') as f:
                        f.write("# RMSD data\n")
                        f.write("# Time (ps)  RMSD (nm)\n")
                        for t, r in zip(result.times, result.rmsd_values):
                            f.write(f"{t:.3f}  {r:.4f}\n")
                else:
                    df.to_csv(output_file, index=False)

            return {
                "success": True,
                "data": {
                    "n_frames": result.n_frames,
                    "mean_rmsd": round(result.mean_rmsd, 4),
                    "std_rmsd": round(result.std_rmsd, 4),
                    "min_rmsd": round(result.min_rmsd, 4),
                    "max_rmsd": round(result.max_rmsd, 4),
                    "output_file": output_file,
                },
            }

        except Exception as e:
            return {
                "success": False,
                "error": str(e),
            }

    def _format_rmsd_message(self, data: dict) -> str:
        """Format RMSD result message."""
        mean = data["mean_rmsd"]
        std = data["std_rmsd"]
        n_frames = data["n_frames"]

        # Assess stability
        if mean < 0.2:
            stability = "Excellent stability"
        elif mean < 0.3:
            stability = "Good stability"
        elif mean < 0.5:
            stability = "Moderate stability"
        else:
            stability = "High deviation"

        msg = f"{stability}. Mean RMSD: {mean:.3f} ± {std:.3f} nm "
        msg += f"(range: {data['min_rmsd']:.3f}-{data['max_rmsd']:.3f} nm). "
        msg += f"Analyzed {n_frames} frames."

        return msg

    def system_prompt_section(self) -> str:
        """Tool-specific guidance for the agent."""
        return """
## RMSD Calculation Guidelines

**When to calculate RMSD:**
- After preprocessing to assess simulation stability
- To check convergence before further analysis
- When user asks about structural stability or equilibration

**Interpreting RMSD values (for proteins):**
- < 0.2 nm: Excellent stability, well-equilibrated
- 0.2-0.3 nm: Good stability, typical for stable proteins
- 0.3-0.5 nm: Moderate stability, may need longer equilibration
- > 0.5 nm: High deviation, check for issues or large conformational changes

**Common selections:**
- 'protein and name CA': C-alpha atoms (default, most common)
- 'backbone': All backbone atoms (N, CA, C, O)
- 'protein': All protein atoms (more sensitive to side-chain motion)
- 'resid 1-50 and name CA': Specific domain or region

**Best practices:**
- Use C-alpha atoms for overall structural stability
- Use backbone for more detailed analysis
- For multi-domain proteins, analyze domains separately
- Save output file for plotting and further analysis
"""
