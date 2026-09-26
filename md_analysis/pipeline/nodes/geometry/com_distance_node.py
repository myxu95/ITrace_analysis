"""COM distance analysis node."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt

from ....analysis.geometry import COMDistanceCalculator, COMDistanceInput
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class COMDistanceNode(PipelineNode):
    """Calculate the pHLA-TCR center-of-mass distance time series."""

    def __init__(
        self,
        tcr_selection: Optional[str] = None,
        mhc_selection: Optional[str] = None,
        stride: int = 1,
        time_unit: str = "ps",
        name: Optional[str] = None,
    ):
        super().__init__(name=name or "COMDistanceNode")
        self.tcr_selection = tcr_selection
        self.mhc_selection = mhc_selection
        self.stride = stride
        self.time_unit = time_unit

    def validate_inputs(self, context: PipelineContext):
        missing = []
        if not context.topology:
            missing.append("topology")
        if not (context.trajectory_processed or context.trajectory_raw):
            missing.append("trajectory_processed_or_raw")
        if (self.tcr_selection is None or self.mhc_selection is None) and "chain_mapping" not in context.metadata:
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
        tcr_selection, mhc_selection = self._resolve_selections(context)

        try:
            calculator = COMDistanceCalculator(topology_file, trajectory_file)
            result = calculator.calculate(
                COMDistanceInput(
                    tcr_selection=tcr_selection,
                    mhc_selection=mhc_selection,
                    stride=self.stride,
                    time_unit=self.time_unit,
                )
            )

            output_dir = Path(context.get_analysis_path("geometry", "com_distance_timeseries.csv")).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            csv_file = output_dir / "com_distance_timeseries.csv"
            summary_file = output_dir / "com_distance_summary.json"
            plot_file = output_dir / "com_distance_timeseries.png"

            result.to_timeseries_frame().to_csv(csv_file, index=False)
            summary_file.write_text(
                json.dumps(result.to_summary(), ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
            self._write_plot(result, plot_file)

            context.results["com_distance"] = {
                "timeseries_file": str(csv_file),
                "summary_file": str(summary_file),
                "plot_file": str(plot_file),
                "tcr_selection": tcr_selection,
                "mhc_selection": mhc_selection,
                "stride": self.stride,
                "time_unit": self.time_unit,
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"COM distance analysis failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc

    def _resolve_selections(self, context: PipelineContext) -> tuple[str, str]:
        if self.tcr_selection and self.mhc_selection:
            return self.tcr_selection, self.mhc_selection

        chain_mapping = context.metadata["chain_mapping"]
        mhc_selection = " or ".join(
            f"chainID {chain_id}"
            for chain_id in [
                chain_mapping.get("mhc_alpha"),
                chain_mapping.get("b2m"),
                chain_mapping.get("peptide"),
            ]
            if chain_id
        )
        tcr_selection = " or ".join(
            f"chainID {chain_id}"
            for chain_id in [
                chain_mapping.get("tcr_alpha"),
                chain_mapping.get("tcr_beta"),
            ]
            if chain_id
        )
        return tcr_selection, mhc_selection

    def _write_plot(self, result, output_file: Path) -> None:
        frame = result.to_timeseries_frame()
        time_column = f"time_{result.time_unit}"

        fig, ax = plt.subplots(figsize=(9.6, 4.8))
        ax.plot(frame[time_column], frame["com_distance_angstrom"], color="#315f93", linewidth=1.8)
        ax.set_title("TCR-MHC COM distance", fontsize=13, weight="bold")
        ax.set_xlabel(f"Time ({result.time_unit})")
        ax.set_ylabel("Distance (Å)")
        ax.grid(alpha=0.2, linestyle="--")
        fig.tight_layout()
        fig.savefig(output_file, dpi=220, facecolor="white")
        plt.close(fig)
