"""Contact heatmap generation node."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ....analysis.topology import ContactHeatmapPlotter
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class ContactHeatmapNode(PipelineNode):
    """Generate contact heatmaps from contact reports."""

    def __init__(self, name: str | None = None):
        super().__init__(name=name or "ContactHeatmapNode")

    def validate_inputs(self, context: PipelineContext):
        if "contact_annotation" not in context.results:
            raise PipelineError(
                node_name=self.name,
                reason="contact_annotation results not found in context",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        outputs = context.results["contact_annotation"]
        report_file = outputs["contact_report_file"]
        report_df = pd.read_csv(report_file)

        plotter = ContactHeatmapPlotter()
        heatmap_outputs = {
            "global": plotter.generate_bundle(report_df, Path(report_file).parent),
        }

        for region_name in ("cdr1", "cdr2", "cdr3", "non_cdr"):
            region_outputs = outputs.get(region_name)
            if not region_outputs:
                continue
            contacts_file = region_outputs["contacts_file"]
            region_df = pd.read_csv(contacts_file)
            heatmap_outputs[region_name] = plotter.generate_bundle(
                region_df,
                Path(region_outputs["heatmaps_dir"]),
            )

        context.results["contact_heatmaps"] = heatmap_outputs
        return context
