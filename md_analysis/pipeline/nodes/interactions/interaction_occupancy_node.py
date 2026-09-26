"""
Generic interaction occupancy and stability analysis node.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from ....analysis.interactions import InteractionOccupancyAnalyzer
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class InteractionOccupancyNode(PipelineNode):
    """Analyze stability for an annotated report from one interaction family."""

    def __init__(
        self,
        source_result_key: str,
        report_file_key: str,
        output_subdir: str,
        family_name: str,
        result_key: str,
        name: str | None = None,
    ):
        super().__init__(name=name or f"{family_name.title()}OccupancyNode")
        self.source_result_key = source_result_key
        self.report_file_key = report_file_key
        self.output_subdir = output_subdir
        self.family_name = family_name
        self.result_key = result_key

    def validate_inputs(self, context: PipelineContext):
        if self.source_result_key not in context.results:
            raise PipelineError(
                node_name=self.name,
                reason=f"{self.source_result_key} results not found in context",
                context_state={"system_id": context.system_id},
            )
        source = context.results[self.source_result_key]
        if self.report_file_key not in source:
            raise PipelineError(
                node_name=self.name,
                reason=f"{self.report_file_key} not found under {self.source_result_key}",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)
        source = context.results[self.source_result_key]
        report_file = Path(source[self.report_file_key])
        report = pd.read_csv(report_file) if report_file.exists() else pd.DataFrame()

        analyzer = InteractionOccupancyAnalyzer()
        outputs = analyzer.build_all_outputs(report)
        pair_table = outputs["pair_table"]
        group_summary = outputs["interaction_class_contribution"]
        persistent_ranking = outputs["persistent_ranking"]
        occupancy_matrix = outputs["occupancy_matrix"]
        peptide_mhc_summary = outputs["peptide_vs_mhc_summary"]
        region_summary = outputs["region_summary"]

        output_dir = Path(context.get_analysis_path(self.output_subdir, "pair_stability.csv")).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        pair_file = output_dir / "pair_stability.csv"
        group_file = output_dir / "interaction_class_contribution.csv"
        persistent_file = output_dir / "persistent_interaction_ranking.csv"
        matrix_file = output_dir / "residue_pair_occupancy_matrix.csv"
        peptide_mhc_file = output_dir / "peptide_vs_mhc_summary.csv"
        region_file = output_dir / "region_occupancy_summary.csv"
        summary_file = output_dir / "occupancy_summary.json"

        pair_table.to_csv(pair_file, index=False)
        group_summary.to_csv(group_file, index=False)
        persistent_ranking.to_csv(persistent_file, index=False)
        occupancy_matrix.to_csv(matrix_file)
        peptide_mhc_summary.to_csv(peptide_mhc_file, index=False)
        region_summary.to_csv(region_file, index=False)
        analyzer.write_summary_json(pair_table, group_summary, summary_file)
        plots = analyzer.plot_outputs(pair_table, group_summary, output_dir, self.family_name)
        timeline_plot = analyzer.plot_timeline(pair_table, output_dir, self.family_name)

        context.results[self.result_key] = {
            "pair_stability_file": str(pair_file),
            "interaction_class_contribution_file": str(group_file),
            "persistent_interaction_ranking_file": str(persistent_file),
            "residue_pair_occupancy_matrix_file": str(matrix_file),
            "peptide_vs_mhc_summary_file": str(peptide_mhc_file),
            "region_occupancy_summary_file": str(region_file),
            "summary_file": str(summary_file),
            "timeline_plot": str(timeline_plot),
            **plots,
        }
        return context
