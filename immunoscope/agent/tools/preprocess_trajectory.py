"""Preprocess trajectory tool for ImmunoScope Agent."""

from __future__ import annotations
import asyncio
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class PreprocessTrajectoryInput(BaseModel):
    """Input schema for preprocess_trajectory tool."""

    trajectory: str = Field(
        ...,
        description="Path to input trajectory file (e.g., md.xtc, md.trr)",
    )
    topology: str = Field(
        ...,
        description="Path to topology file (e.g., md.tpr, md.gro)",
    )
    output_dir: str = Field(
        ...,
        description="Output directory for processed trajectory",
    )
    fit_group: str = Field(
        default="Backbone",
        description="Group for fitting (default: Backbone)",
    )
    output_group: str = Field(
        default="System",
        description="Group for output (default: System)",
    )


@register_tool
class PreprocessTrajectoryTool(Tool):
    """Preprocess MD trajectory with PBC correction.

    This tool performs periodic boundary condition (PBC) correction on MD trajectories.
    PBC correction is MANDATORY before any analysis to ensure molecules are not split
    across periodic boundaries. The process includes:

    1. Remove jumps across periodic boundaries (nojump)
    2. Fit rotation and translation to reference structure
    3. Center the system in the box

    Use this tool FIRST before any other analysis.
    """

    name: ClassVar[str] = "preprocess_trajectory"
    description: ClassVar[str] = (
        "Preprocess MD trajectory with PBC (Periodic Boundary Condition) correction. "
        "This is MANDATORY before any analysis. Removes jumps, fits rotation/translation, "
        "and centers the system. Input: trajectory file (.xtc/.trr), topology file (.tpr/.gro). "
        "Output: processed trajectory ready for analysis. Takes 1-5 minutes depending on size."
    )
    Input: ClassVar[type[BaseModel]] = PreprocessTrajectoryInput

    is_read_only: ClassVar[bool] = False  # Creates output files
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 600.0  # 10 minutes

    async def call(self, args: PreprocessTrajectoryInput, ctx: ToolContext) -> ToolResult:
        """Execute trajectory preprocessing."""
        try:
            # Validate inputs
            traj_path = Path(args.trajectory).resolve()
            topo_path = Path(args.topology).resolve()
            output_dir = Path(args.output_dir).resolve()

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

            # Create output directory
            output_dir.mkdir(parents=True, exist_ok=True)

            # Determine output filename
            output_filename = traj_path.stem + "_processed.xtc"
            output_path = output_dir / output_filename

            # Report progress
            await ctx.on_event({
                "type": "progress",
                "message": f"Starting PBC correction for {traj_path.name}",
                "progress": 0.0,
            })

            # Run preprocessing in executor (blocking operation)
            loop = asyncio.get_event_loop()
            result = await loop.run_in_executor(
                None,
                self._run_preprocessing,
                str(traj_path),
                str(topo_path),
                str(output_path),
                args.fit_group,
                args.output_group,
            )

            if result["success"]:
                return ToolResult(content={
                    "status": "success",
                    "output_trajectory": str(output_path),
                    "input_trajectory": str(traj_path),
                    "topology": str(topo_path),
                    "fit_group": args.fit_group,
                    "output_group": args.output_group,
                    "message": f"Trajectory preprocessed successfully: {output_filename}",
                })
            else:
                return ToolResult(
                    content=f"Preprocessing failed: {result['error']}",
                    is_error=True,
                )

        except Exception as e:
            return ToolResult(
                content=f"Preprocessing failed: {str(e)}",
                is_error=True,
            )

    def _run_preprocessing(
        self,
        trajectory: str,
        topology: str,
        output: str,
        fit_group: str,
        output_group: str,
    ) -> dict:
        """Run preprocessing (blocking operation)."""
        try:
            from immunoscope.analysis.trajectory.pbc import PBCProcessor

            processor = PBCProcessor(keep_temp_files=False)

            # Run PBC correction
            result_path = processor.remove_pbc_2step(
                trajectory=trajectory,
                topology=topology,
                output=output,
                fit_group=fit_group,
                output_group=output_group,
            )

            return {
                "success": True,
                "output_path": result_path,
            }

        except Exception as e:
            return {
                "success": False,
                "error": str(e),
            }

    def system_prompt_section(self) -> str:
        """Tool-specific guidance for the agent."""
        return """
## Preprocessing Guidelines

**ALWAYS preprocess trajectories first:**
- Raw MD trajectories have PBC artifacts (molecules split across boundaries)
- Preprocessing is MANDATORY before RMSD, RMSF, contacts, or any analysis
- Use default groups (Backbone for fit, System for output) unless user specifies

**When to preprocess:**
- User provides raw trajectory files (.xtc, .trr)
- Before any structural or interaction analysis
- When trajectory shows discontinuities or jumps

**Skip preprocessing if:**
- File is already named "*_processed.xtc"
- User explicitly says trajectory is already preprocessed
"""
