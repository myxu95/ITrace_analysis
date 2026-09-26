"""Per-residue SASA pipeline node."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from ....analysis.interface import ResidueSasaAnalyzer
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class ResidueSasaNode(PipelineNode):
    """Calculate per-residue bound / unbound SASA and exposure metrics."""

    def __init__(
        self,
        selection_a: Optional[str] = None,
        selection_b: Optional[str] = None,
        probe_radius: float = 1.4,
        stride: int = 1,
        name: Optional[str] = None,
    ):
        super().__init__(name=name or "ResidueSasaNode")
        self.selection_a = selection_a
        self.selection_b = selection_b
        self.probe_radius = probe_radius
        self.stride = stride

    def validate_inputs(self, context: PipelineContext):
        missing = []
        if not context.topology:
            missing.append("topology")
        if not (context.trajectory_processed or context.trajectory_raw):
            missing.append("trajectory_processed_or_raw")
        if (self.selection_a is None or self.selection_b is None) and "chain_mapping" not in context.metadata:
            missing.append("chain_mapping_or_manual_selections")
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
        selection_a, selection_b = self._resolve_selections(context)

        try:
            analyzer = ResidueSasaAnalyzer(topology, trajectory)
            result = analyzer.calculate(
                selection_a=selection_a,
                selection_b=selection_b,
                probe_radius=self.probe_radius,
                stride=self.stride,
            )

            output_dir = Path(context.get_analysis_path("interface", "residue_sasa.csv")).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            csv_file = output_dir / "residue_sasa.csv"
            json_file = output_dir / "residue_sasa_summary.json"

            result.residue_frame.to_csv(csv_file, index=False)
            json_file.write_text(
                json.dumps(result.summary, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )

            context.results["residue_sasa"] = {
                "residue_csv": str(csv_file),
                "summary_json": str(json_file),
                "selection_a": selection_a,
                "selection_b": selection_b,
                "probe_radius": self.probe_radius,
                "stride": self.stride,
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"Residue SASA analysis failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc

    def _resolve_selections(self, context: PipelineContext) -> tuple[str, str]:
        if self.selection_a and self.selection_b:
            return self.selection_a, self.selection_b

        chain_mapping = context.metadata["chain_mapping"]
        phla_chains = [
            chain_mapping["mhc_alpha"],
            chain_mapping["b2m"],
            chain_mapping["peptide"],
        ]
        tcr_chains = [
            chain_mapping["tcr_alpha"],
            chain_mapping["tcr_beta"],
        ]
        selection_a = " or ".join(
            f"chainID {chain_id}" for chain_id in phla_chains if chain_id
        )
        selection_b = " or ".join(
            f"chainID {chain_id}" for chain_id in tcr_chains if chain_id
        )
        return selection_a, selection_b
