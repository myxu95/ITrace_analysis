"""Sidechain chi-dihedral entropy pipeline node."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from ....analysis.geometry import ChiDihedralAnalyzer
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class ChiDihedralNode(PipelineNode):
    """Compute per-residue chi1/chi2 entropy and rotamer diversity."""

    def __init__(
        self,
        selection: str = "protein",
        stride: int = 1,
        name: Optional[str] = None,
    ):
        super().__init__(name=name or "ChiDihedralNode")
        self.selection = selection
        self.stride = stride

    def validate_inputs(self, context: PipelineContext):
        missing = []
        if not context.topology:
            missing.append("topology")
        if not (context.trajectory_processed or context.trajectory_raw):
            missing.append("trajectory_processed_or_raw")
        if missing:
            raise PipelineError(
                node_name=self.name,
                reason=f"Missing required inputs: {missing}",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        topology = context.topology
        trajectory = context.trajectory_processed or context.trajectory_raw

        try:
            analyzer = ChiDihedralAnalyzer(topology, trajectory)
            result = analyzer.calculate(
                selection=self.selection,
                stride=self.stride,
            )

            output_dir = Path(
                context.get_analysis_path("geometry/chi_dihedrals", "residue_chi_dihedrals.csv")
            ).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            csv_file = output_dir / "residue_chi_dihedrals.csv"
            json_file = output_dir / "residue_chi_dihedrals_summary.json"

            result.residue_frame.to_csv(csv_file, index=False)
            json_file.write_text(
                json.dumps(result.summary, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )

            context.results["chi_dihedrals"] = {
                "residue_csv": str(csv_file),
                "summary_json": str(json_file),
                "selection": self.selection,
                "stride": self.stride,
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"Chi dihedral analysis failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc
