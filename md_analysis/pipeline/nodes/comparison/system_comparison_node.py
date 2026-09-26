"""Two-system comparison node."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import pandas as pd

from ....analysis.comparison import (
    SingleCaseLoader,
    SystemComparisonBuilder,
    write_flexibility_comparison_plot,
    write_interaction_family_comparison_plot,
    write_quality_interface_comparison_plot,
    write_rrcs_comparison_plot,
)
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class SystemComparisonNode(PipelineNode):
    """Load two single-system results and build comparison artifacts."""

    def __init__(
        self,
        case_a_root: str,
        case_b_root: str,
        label_a: str,
        label_b: str,
        comparison_mode: str = "generic",
        comparison_context: str = "",
        comparison_scope: str = "auto",
        alignment_selection: str = "phla_core_ca",
        residue_mapping: str = "auto",
        name: Optional[str] = None,
    ):
        super().__init__(name=name or "SystemComparisonNode")
        self.case_a_root = Path(case_a_root)
        self.case_b_root = Path(case_b_root)
        self.label_a = label_a
        self.label_b = label_b
        self.comparison_mode = comparison_mode
        self.comparison_context = comparison_context
        self.comparison_scope = comparison_scope
        self.alignment_selection = alignment_selection
        self.residue_mapping = residue_mapping

    def validate_inputs(self, context: PipelineContext) -> None:
        missing = []
        if not self.case_a_root.exists():
            missing.append(f"case_a_root={self.case_a_root}")
        if not self.case_b_root.exists():
            missing.append(f"case_b_root={self.case_b_root}")
        if missing:
            raise PipelineError(
                node_name=self.name,
                reason=f"Missing required inputs: {missing}",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        try:
            loader = SingleCaseLoader()
            case_a = loader.load(self.case_a_root, self.label_a)
            case_b = loader.load(self.case_b_root, self.label_b)

            builder = SystemComparisonBuilder()
            built = builder.build(
                case_a,
                case_b,
                comparison_mode=self.comparison_mode,
                comparison_context=self.comparison_context,
                comparison_scope=self.comparison_scope,
                alignment_selection=self.alignment_selection,
                residue_mapping=self.residue_mapping,
            )

            analysis_dir = Path(context.get_analysis_path("comparison", "comparison_table.csv")).parent
            analysis_dir.mkdir(parents=True, exist_ok=True)

            comparison_table_csv = analysis_dir / "comparison_table.csv"
            identity_csv = analysis_dir / "identity_comparison.csv"
            rmsf_region_csv = analysis_dir / "rmsf_region_comparison.csv"
            rrcs_region_csv = analysis_dir / "rrcs_region_comparison.csv"
            interaction_family_csv = analysis_dir / "interaction_family_comparison.csv"
            summary_json = analysis_dir / "comparison_summary.json"
            quality_plot = analysis_dir / "quality_interface_comparison.png"
            flexibility_plot = analysis_dir / "flexibility_comparison.png"
            rrcs_plot = analysis_dir / "rrcs_region_comparison.png"
            interaction_plot = analysis_dir / "interaction_family_comparison.png"

            pd.DataFrame(built.comparison_rows).to_csv(comparison_table_csv, index=False)
            pd.DataFrame(built.identity_rows).to_csv(identity_csv, index=False)
            pd.DataFrame(built.rmsf_region_rows).to_csv(rmsf_region_csv, index=False)
            pd.DataFrame(built.rrcs_region_rows).to_csv(rrcs_region_csv, index=False)
            pd.DataFrame(built.interaction_family_rows).to_csv(interaction_family_csv, index=False)
            summary_json.write_text(json.dumps(built.summary, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")

            comparison_table = pd.DataFrame(built.comparison_rows)
            write_quality_interface_comparison_plot(self.label_a, self.label_b, comparison_table, quality_plot)
            write_flexibility_comparison_plot(self.label_a, self.label_b, pd.DataFrame(built.rmsf_region_rows), flexibility_plot)
            write_rrcs_comparison_plot(self.label_a, self.label_b, pd.DataFrame(built.rrcs_region_rows), rrcs_plot)
            write_interaction_family_comparison_plot(
                self.label_a,
                self.label_b,
                pd.DataFrame(built.interaction_family_rows),
                interaction_plot,
            )

            context.results["comparison"] = {
                "summary": built.summary,
                "artifacts": {
                    "comparison_table_csv": str(comparison_table_csv),
                    "identity_comparison_csv": str(identity_csv),
                    "rmsf_region_comparison_csv": str(rmsf_region_csv),
                    "rrcs_region_comparison_csv": str(rrcs_region_csv),
                    "interaction_family_comparison_csv": str(interaction_family_csv),
                    "summary_json": str(summary_json),
                    "quality_interface_plot": str(quality_plot),
                    "flexibility_plot": str(flexibility_plot),
                    "rrcs_plot": str(rrcs_plot),
                    "interaction_plot": str(interaction_plot),
                },
                "tables": {
                    "comparison_rows": built.comparison_rows,
                    "identity_rows": built.identity_rows,
                    "rmsf_region_rows": built.rmsf_region_rows,
                    "rrcs_region_rows": built.rrcs_region_rows,
                    "interaction_family_rows": built.interaction_family_rows,
                },
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"System comparison failed: {exc}",
                context_state={
                    "system_id": context.system_id,
                    "case_a_root": str(self.case_a_root),
                    "case_b_root": str(self.case_b_root),
                },
            ) from exc
