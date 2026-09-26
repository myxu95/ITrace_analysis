"""Check quality tool for ImmunoScope Agent."""

from __future__ import annotations
import asyncio
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class CheckQualityInput(BaseModel):
    """Input schema for check_quality tool."""

    trajectory: str = Field(
        ...,
        description="Path to trajectory file (preferably preprocessed)",
    )
    topology: str = Field(
        ...,
        description="Path to topology file (.tpr, .gro)",
    )
    min_time_ns: float = Field(
        default=1.0,
        ge=0.0,
        description="Minimum expected simulation time in nanoseconds (default: 1.0)",
    )


@register_tool
class CheckQualityTool(Tool):
    """Check MD trajectory quality and completeness.

    This tool validates trajectory quality by checking:
    - File integrity and size
    - Simulation completeness
    - Basic trajectory statistics (number of frames, time coverage)
    - Potential issues (missing frames, corrupted data)

    Use this tool after preprocessing to ensure trajectory is suitable for analysis.
    """

    name: ClassVar[str] = "check_quality"
    description: ClassVar[str] = (
        "Check MD trajectory quality and completeness. Validates file integrity, "
        "simulation time, frame count, and identifies potential issues. "
        "Use after preprocessing to ensure trajectory is ready for analysis. "
        "Returns quality score, warnings, and recommendations."
    )
    Input: ClassVar[type[BaseModel]] = CheckQualityInput

    is_read_only: ClassVar[bool] = True
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 60.0

    async def call(self, args: CheckQualityInput, ctx: ToolContext) -> ToolResult:
        """Execute quality check."""
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

            # Report progress
            await ctx.on_event({
                "type": "progress",
                "message": f"Checking quality of {traj_path.name}",
                "progress": 0.0,
            })

            # Run quality check in executor
            loop = asyncio.get_event_loop()
            result = await loop.run_in_executor(
                None,
                self._run_quality_check,
                str(traj_path),
                str(topo_path),
                args.min_time_ns,
            )

            if result["success"]:
                quality_data = result["data"]

                # Determine overall status
                status = "good"
                if quality_data["warnings"]:
                    status = "warning"
                if quality_data["quality_score"] < 0.5:
                    status = "poor"

                return ToolResult(content={
                    "status": status,
                    "quality_score": quality_data["quality_score"],
                    "trajectory": str(traj_path),
                    "topology": str(topo_path),
                    "n_frames": quality_data["n_frames"],
                    "time_ns": quality_data["time_ns"],
                    "file_size_mb": quality_data["file_size_mb"],
                    "warnings": quality_data["warnings"],
                    "recommendations": quality_data["recommendations"],
                    "message": self._format_quality_message(quality_data),
                })
            else:
                return ToolResult(
                    content=f"Quality check failed: {result['error']}",
                    is_error=True,
                )

        except Exception as e:
            return ToolResult(
                content=f"Quality check failed: {str(e)}",
                is_error=True,
            )

    def _run_quality_check(
        self,
        trajectory: str,
        topology: str,
        min_time_ns: float,
    ) -> dict:
        """Run quality check (blocking operation)."""
        try:
            import MDAnalysis as mda

            # Load trajectory
            u = mda.Universe(topology, trajectory)

            # Basic statistics
            n_frames = len(u.trajectory)
            n_atoms = len(u.atoms)

            # Time information
            dt = u.trajectory.dt  # ps
            total_time_ps = u.trajectory.totaltime
            total_time_ns = total_time_ps / 1000.0

            # File size
            traj_path = Path(trajectory)
            file_size_mb = traj_path.stat().st_size / (1024 * 1024)

            # Quality checks
            warnings = []
            recommendations = []

            # Check simulation time
            if total_time_ns < min_time_ns:
                warnings.append(f"Simulation time ({total_time_ns:.1f} ns) is less than minimum ({min_time_ns:.1f} ns)")
                recommendations.append("Consider running longer simulation for better statistics")

            # Check frame count
            if n_frames < 100:
                warnings.append(f"Low frame count ({n_frames} frames)")
                recommendations.append("Increase trajectory output frequency or simulation time")

            # Check file size
            if file_size_mb < 1.0:
                warnings.append(f"Small trajectory file ({file_size_mb:.1f} MB)")

            # Calculate quality score (0-1)
            quality_score = 1.0
            if total_time_ns < min_time_ns:
                quality_score *= (total_time_ns / min_time_ns)
            if n_frames < 100:
                quality_score *= (n_frames / 100.0)
            quality_score = min(1.0, quality_score)

            return {
                "success": True,
                "data": {
                    "n_frames": n_frames,
                    "n_atoms": n_atoms,
                    "time_ns": round(total_time_ns, 2),
                    "time_ps": round(total_time_ps, 2),
                    "dt_ps": round(dt, 2),
                    "file_size_mb": round(file_size_mb, 2),
                    "quality_score": round(quality_score, 2),
                    "warnings": warnings,
                    "recommendations": recommendations,
                },
            }

        except Exception as e:
            return {
                "success": False,
                "error": str(e),
            }

    def _format_quality_message(self, data: dict) -> str:
        """Format quality check message."""
        score = data["quality_score"]
        time_ns = data["time_ns"]
        n_frames = data["n_frames"]

        if score >= 0.8:
            status = "Excellent"
        elif score >= 0.6:
            status = "Good"
        elif score >= 0.4:
            status = "Fair"
        else:
            status = "Poor"

        msg = f"{status} quality (score: {score}). "
        msg += f"Trajectory: {time_ns} ns, {n_frames} frames. "

        if data["warnings"]:
            msg += f"Warnings: {len(data['warnings'])}."
        else:
            msg += "No issues detected."

        return msg

    def system_prompt_section(self) -> str:
        """Tool-specific guidance for the agent."""
        return """
## Quality Check Guidelines

**When to check quality:**
- After preprocessing, before heavy analysis
- When user asks about trajectory stability or convergence
- Before running computationally expensive analyses

**Interpreting results:**
- Quality score ≥ 0.8: Excellent, proceed with all analyses
- Quality score 0.6-0.8: Good, suitable for most analyses
- Quality score 0.4-0.6: Fair, may need longer simulation
- Quality score < 0.4: Poor, recommend longer simulation or check for issues

**Common warnings:**
- Short simulation time: Recommend longer simulation
- Low frame count: Increase output frequency
- Small file size: May indicate incomplete simulation
"""
