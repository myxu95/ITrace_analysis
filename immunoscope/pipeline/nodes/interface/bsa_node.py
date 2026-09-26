"""BSA analysis node."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt

from ....analysis.interface import BuriedSurfaceAreaAnalyzer
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class BSAAnalysisNode(PipelineNode):
    """Calculate pHLA-TCR buried surface area time series."""

    def __init__(
        self,
        selection_a: Optional[str] = None,
        selection_b: Optional[str] = None,
        probe_radius: float = 1.4,
        stride: int = 1,
        time_unit: str = "ps",
        name: Optional[str] = None,
    ):
        super().__init__(name=name or "BSAAnalysisNode")
        self.selection_a = selection_a
        self.selection_b = selection_b
        self.probe_radius = probe_radius
        self.stride = stride
        self.time_unit = time_unit

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

        topology_file = context.topology
        trajectory_file = context.trajectory_processed or context.trajectory_raw
        selection_a, selection_b = self._resolve_selections(context)

        try:
            analyzer = BuriedSurfaceAreaAnalyzer(topology_file, trajectory_file)
            result = analyzer.calculate(
                selection_a=selection_a,
                selection_b=selection_b,
                probe_radius=self.probe_radius,
                stride=self.stride,
                time_unit=self.time_unit,
            )

            output_dir = Path(context.get_analysis_path("interface", "bsa_timeseries.csv")).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            csv_file = output_dir / "bsa_timeseries.csv"
            json_file = output_dir / "interface_summary.json"
            plot_file = output_dir / "bsa_timeseries.png"

            result.to_timeseries_frame().to_csv(csv_file, index=False)
            json_file.write_text(
                json.dumps(result.to_summary(), indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            self._write_plot(result, plot_file)

            context.results["bsa"] = {
                "timeseries_file": str(csv_file),
                "summary_file": str(json_file),
                "plot_file": str(plot_file),
                "selection_a": selection_a,
                "selection_b": selection_b,
                "probe_radius": self.probe_radius,
                "stride": self.stride,
                "time_unit": self.time_unit,
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"BSA analysis failed: {exc}",
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
        selection_a = " or ".join(f"chainID {chain_id}" for chain_id in phla_chains)
        selection_b = " or ".join(f"chainID {chain_id}" for chain_id in tcr_chains)
        return selection_a, selection_b

    def _write_plot(self, result, output_file: Path) -> None:
        frame = result.to_timeseries_frame()
        time_column = f"time_{result.time_unit}"

        fig, axes = plt.subplots(2, 1, figsize=(9.6, 7.2), sharex=True)

        axes[0].plot(frame[time_column], frame["buried_surface_area"], color="#2f5d50", linewidth=1.8)
        axes[0].set_ylabel("BSA ($\\AA^2$)")
        axes[0].set_title("Buried surface area trajectory", fontsize=13, weight="bold")
        axes[0].grid(alpha=0.2, linestyle="--")

        axes[1].plot(frame[time_column], frame["interface_ratio"], color="#c37a3b", linewidth=1.8)
        axes[1].set_ylabel("Interface ratio")
        axes[1].set_xlabel(f"Time ({result.time_unit})")
        axes[1].set_title("Interface ratio trajectory", fontsize=13, weight="bold")
        axes[1].grid(alpha=0.2, linestyle="--")

        fig.tight_layout()
        fig.savefig(output_file, dpi=220, facecolor="white")
        plt.close(fig)
