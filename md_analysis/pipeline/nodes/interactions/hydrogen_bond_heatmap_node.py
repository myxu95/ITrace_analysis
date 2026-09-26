"""Hydrogen-bond heatmap generation node."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ....analysis.topology import ContactHeatmapPlotter
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class HydrogenBondHeatmapNode(PipelineNode):
    """Generate region heatmaps from hydrogen-bond reports."""

    def __init__(self, name: str | None = None):
        super().__init__(name=name or "HydrogenBondHeatmapNode")

    def validate_inputs(self, context: PipelineContext):
        if "hbond_annotation" not in context.results:
            raise PipelineError(
                node_name=self.name,
                reason="hbond_annotation results not found in context",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        outputs = context.results["hbond_annotation"]
        report_file = outputs["hbond_report_file"]
        report_df = pd.read_csv(report_file)

        plotter = ContactHeatmapPlotter()
        heatmap_outputs = {
            "global": plotter.generate_bundle(report_df, Path(report_file).parent),
        }

        for region_name in ("cdr1", "cdr2", "cdr3", "non_cdr"):
            region_outputs = outputs.get(region_name)
            if not region_outputs:
                continue
            region_df = pd.read_csv(region_outputs["contacts_file"])
            heatmap_outputs[region_name] = plotter.generate_bundle(
                region_df,
                Path(region_outputs["heatmaps_dir"]),
            )

        context.results["hbond_heatmaps"] = heatmap_outputs
        return context
