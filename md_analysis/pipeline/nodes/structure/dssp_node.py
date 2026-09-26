"""Per-residue DSSP pipeline node."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from ....analysis.structure import DSSPAnalyzer
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class DSSPNode(PipelineNode):
    """Compute per-residue secondary structure (DSSP) over the trajectory."""

    def __init__(
        self,
        selection: str = "protein",
        stride: int = 1,
        name: Optional[str] = None,
    ):
        super().__init__(name=name or "DSSPNode")
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
            analyzer = DSSPAnalyzer(topology, trajectory)
            result = analyzer.calculate(
                selection=self.selection,
                stride=self.stride,
            )

            output_dir = Path(
                context.get_analysis_path("structure", "residue_secondary_structure.csv")
            ).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            csv_file = output_dir / "residue_secondary_structure.csv"
            json_file = output_dir / "residue_secondary_structure_summary.json"

            result.residue_frame.to_csv(csv_file, index=False)
            json_file.write_text(
                json.dumps(result.summary, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )

            context.results["dssp"] = {
                "residue_csv": str(csv_file),
                "summary_json": str(json_file),
                "selection": self.selection,
                "stride": self.stride,
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"DSSP analysis failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc
