"""cation-pi residue-pair calculation node."""

from __future__ import annotations

import json
from pathlib import Path

from ....analysis.interactions import CationPiPairAnalyzer
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class CationPiPairNode(PipelineNode):
    """Calculate residue-level cation-pi occupancy across the pHLA-TCR interface."""

    def __init__(self, stride: int = 1, distance_cutoff: float = 6.0, normal_angle_cutoff: float = 60.0, name: str | None = None):
        super().__init__(name=name or "CationPiPairNode")
        self.stride = stride
        self.distance_cutoff = distance_cutoff
        self.normal_angle_cutoff = normal_angle_cutoff

    def validate_inputs(self, context: PipelineContext):
        missing = []
        if not context.topology and not context.structure_pdb:
            missing.append("topology_or_structure_pdb")
        if not (context.trajectory_processed or context.trajectory_raw):
            missing.append("trajectory_processed_or_raw")
        if "chain_mapping" not in context.metadata:
            missing.append("chain_mapping")
        if missing:
            raise PipelineError(self.name, f"Missing required inputs: {missing}", {"system_id": context.system_id})

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)
        topology_file = context.topology or context.structure_pdb
        trajectory_file = context.trajectory_processed or context.trajectory_raw
        chain_mapping = context.metadata["chain_mapping"]
        try:
            analyzer = CationPiPairAnalyzer(
                topology_file,
                trajectory_file,
                reference_structure=context.structure_pdb,
            )
            pairs = analyzer.calculate_interface_pairs(
                chain_mapping=chain_mapping,
                stride=self.stride,
                distance_cutoff=self.distance_cutoff,
                normal_angle_cutoff=self.normal_angle_cutoff,
            )
            output_file = context.get_analysis_path("interactions/cation_pi_interactions", "residue_pair_cation_pi.csv")
            summary_file = context.get_analysis_path("interactions/cation_pi_interactions", "residue_pair_cation_pi_summary.json")
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            pairs.to_csv(output_file, index=False)
            with open(summary_file, "w", encoding="utf-8") as handle:
                json.dump(
                    {
                        "topology": topology_file,
                        "trajectory": trajectory_file,
                        "stride": self.stride,
                        "distance_cutoff_angstrom": self.distance_cutoff,
                        "normal_angle_cutoff_deg": self.normal_angle_cutoff,
                        "n_cation_pi_pairs": int(len(pairs)),
                    },
                    handle,
                    indent=2,
                )
            context.results["cation_pi_pairs"] = {
                "raw_pairs_file": output_file,
                "summary_file": summary_file,
                "n_cation_pi_pairs": int(len(pairs)),
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"Cation-pi pair calculation failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc
