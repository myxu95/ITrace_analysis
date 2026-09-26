"""Hydrogen-bond residue-pair calculation node."""

from __future__ import annotations

import json
from pathlib import Path

from ....analysis.interactions import HydrogenBondPairAnalyzer
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class HydrogenBondPairNode(PipelineNode):
    """Calculate residue-level hydrogen-bond occupancy across the pHLA-TCR interface."""

    def __init__(
        self,
        stride: int = 1,
        distance_cutoff: float = 3.5,
        angle_cutoff: float = 150.0,
        name: str | None = None,
    ):
        super().__init__(name=name or "HydrogenBondPairNode")
        self.stride = stride
        self.distance_cutoff = distance_cutoff
        self.angle_cutoff = angle_cutoff

    def validate_inputs(self, context: PipelineContext):
        missing = []
        if not context.topology and not context.structure_pdb:
            missing.append("topology_or_structure_pdb")
        if not (context.trajectory_processed or context.trajectory_raw):
            missing.append("trajectory_processed_or_raw")
        if "chain_mapping" not in context.metadata:
            missing.append("chain_mapping")
        if missing:
            raise PipelineError(
                node_name=self.name,
                reason=f"Missing required inputs: {missing}",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        # Prefer topology with bond/charge information for hydrogen-bond analysis.
        topology_file = context.topology or context.structure_pdb
        trajectory_file = context.trajectory_processed or context.trajectory_raw
        chain_mapping = context.metadata["chain_mapping"]

        try:
            analyzer = HydrogenBondPairAnalyzer(
                topology_file,
                trajectory_file,
                reference_structure=context.structure_pdb,
            )
            pairs = analyzer.calculate_interface_pairs(
                chain_mapping=chain_mapping,
                stride=self.stride,
                distance_cutoff=self.distance_cutoff,
                angle_cutoff=self.angle_cutoff,
            )

            output_file = context.get_analysis_path("interactions/hydrogen_bonds", "residue_pair_hbonds.csv")
            summary_file = context.get_analysis_path("interactions/hydrogen_bonds", "residue_pair_hbonds_summary.json")
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            pairs.to_csv(output_file, index=False)

            summary = {
                "topology": topology_file,
                "trajectory": trajectory_file,
                "stride": self.stride,
                "distance_cutoff_angstrom": self.distance_cutoff,
                "angle_cutoff_deg": self.angle_cutoff,
                "n_hbond_pairs": int(len(pairs)),
            }
            with open(summary_file, "w", encoding="utf-8") as handle:
                json.dump(summary, handle, indent=2)

            context.results["hbond_pairs"] = {
                "raw_pairs_file": output_file,
                "summary_file": summary_file,
                "n_hbond_pairs": int(len(pairs)),
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"Hydrogen bond pair calculation failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc
