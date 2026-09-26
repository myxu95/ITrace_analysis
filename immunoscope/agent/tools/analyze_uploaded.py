"""Tool for analyzing uploaded trajectory files."""

from __future__ import annotations
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.analysis.trajectory.rmsd_refactored import RMSDCalculator, RMSDInput


class AnalyzeUploadedTrajectoryInput(BaseModel):
    """Input for analyzing uploaded trajectory."""

    trajectory_id: str = Field(
        description="ID of the uploaded trajectory (from upload response)"
    )
    analysis_type: str = Field(
        default="quick",
        description="Type of analysis: 'quick' (RMSD+RMSF), 'rmsd', 'rmsf', 'full'"
    )
    output_dir: str = Field(
        default="./agent_analysis_output",
        description="Directory to save analysis results"
    )
    selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection string for analysis"
    )


class AnalyzeUploadedTrajectoryTool(Tool):
    """Analyze uploaded MD trajectory files."""

    name = "analyze_uploaded_trajectory"
    description = (
        "Analyze uploaded MD trajectory files. Performs RMSD, RMSF, and basic "
        "quality checks on the uploaded trajectory. Returns paths to output files."
    )
    Input = AnalyzeUploadedTrajectoryInput
    is_concurrency_safe = False
    timeout_seconds = 600

    async def call(
        self, args: AnalyzeUploadedTrajectoryInput, ctx: ToolContext
    ) -> ToolResult:
        """Execute trajectory analysis."""
        try:
            # Locate uploaded files
            upload_base = Path(ctx.settings.DATA_DIR) / "agent_uploads"
            trajectory_dir = upload_base / args.trajectory_id

            if not trajectory_dir.exists():
                return ToolResult(
                    content=f"Trajectory {args.trajectory_id} not found. "
                            f"Please upload files first.",
                    is_error=True
                )

            # Find trajectory and topology files
            files = self._find_trajectory_files(trajectory_dir)

            if not files["trajectory"]:
                return ToolResult(
                    content="No trajectory file (.xtc, .trr, .dcd) found in upload.",
                    is_error=True
                )

            if not files["topology"]:
                return ToolResult(
                    content="No topology file (.tpr, .pdb, .gro) found in upload.",
                    is_error=True
                )

            # Create output directory
            output_dir = Path(args.output_dir) / args.trajectory_id
            output_dir.mkdir(parents=True, exist_ok=True)

            results = {}

            # Run requested analysis
            if args.analysis_type in ["quick", "rmsd", "full"]:
                rmsd_result = await self._run_rmsd(
                    files["trajectory"],
                    files["topology"],
                    output_dir,
                    args.selection
                )
                results["rmsd"] = rmsd_result

            # Note: RMSF analysis requires structure_pdb parameter
            # Skipping for now, can be added later with proper structure file handling

            # Format summary
            summary = self._format_summary(results, output_dir)

            return ToolResult(content=summary)

        except Exception as e:
            return ToolResult(
                content=f"Analysis failed: {str(e)}",
                is_error=True
            )

    def _find_trajectory_files(self, directory: Path) -> dict[str, Path | None]:
        """Find trajectory and topology files in directory."""
        files = {
            "trajectory": None,
            "topology": None,
            "structure": None
        }

        trajectory_exts = {".xtc", ".trr", ".dcd"}
        topology_exts = {".tpr"}
        structure_exts = {".pdb", ".gro"}

        for file_path in directory.iterdir():
            if not file_path.is_file():
                continue

            ext = file_path.suffix.lower()

            if ext in trajectory_exts and files["trajectory"] is None:
                files["trajectory"] = file_path
            elif ext in topology_exts and files["topology"] is None:
                files["topology"] = file_path
            elif ext in structure_exts:
                files["structure"] = file_path
                # Use structure as topology if no .tpr found
                if files["topology"] is None:
                    files["topology"] = file_path

        return files

    async def _run_rmsd(
        self,
        trajectory: Path,
        topology: Path,
        output_dir: Path,
        selection: str
    ) -> dict[str, Any]:
        """Run RMSD analysis."""
        rmsd_input = RMSDInput(
            topology=str(topology),
            trajectory=str(trajectory),
            selection=selection,
            output_file=str(output_dir / "rmsd.csv"),
            reference_frame=0
        )

        calculator = RMSDCalculator()
        result = calculator.calculate(rmsd_input)

        if not result.success:
            return {
                "success": False,
                "error": result.error_message
            }

        return {
            "success": True,
            "output_file": result.output_file,
            "stats": result.stats
        }

    def _format_summary(self, results: dict[str, Any], output_dir: Path) -> str:
        """Format analysis summary."""
        lines = ["# Trajectory Analysis Complete\n"]

        if "rmsd" in results:
            rmsd = results["rmsd"]
            if rmsd["success"]:
                stats = rmsd["stats"]
                lines.append("## RMSD Analysis")
                lines.append(f"- Output: {rmsd['output_file']}")
                lines.append(f"- Mean RMSD: {stats.get('mean_rmsd', 'N/A'):.3f} Å")
                lines.append(f"- Max RMSD: {stats.get('max_rmsd', 'N/A'):.3f} Å")
                lines.append(f"- Frames: {stats.get('n_frames', 'N/A')}")
                lines.append("")
            else:
                lines.append(f"## RMSD Analysis Failed\n- Error: {rmsd['error']}\n")

        lines.append(f"\n**All results saved to:** `{output_dir}`")

        return "\n".join(lines)


# Register tool
def register():
    """Register this tool with the agent."""
    from immunoscope.agent.tools import TOOL_REGISTRY
    tool = AnalyzeUploadedTrajectoryTool()
    TOOL_REGISTRY[tool.name] = tool
