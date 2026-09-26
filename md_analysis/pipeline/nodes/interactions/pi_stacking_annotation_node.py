"""Pair annotation and region-summary node for pi-pi results."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from ....analysis.topology import RegionInteractionSummaryBuilder, ResiduePairAnnotationAnnotator
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError
from .contact_annotation_node import ContactAnnotationNode


class PiStackingAnnotationNode(PipelineNode):
    """Add semantic annotation to residue-level pi-pi results and write region views."""

    def __init__(self, name: str | None = None):
        super().__init__(name=name or "PiStackingAnnotationNode")

    def validate_inputs(self, context: PipelineContext):
        if "pi_pi_pairs" not in context.results:
            raise PipelineError(self.name, "pi_pi_pairs results not found in context", {"system_id": context.system_id})
        if "chain_mapping" not in context.metadata:
            raise PipelineError(self.name, "chain_mapping metadata not found in context", {"system_id": context.system_id})

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)
        raw_file = context.results["pi_pi_pairs"]["raw_pairs_file"]
        pairs = pd.read_csv(raw_file)
        annotator = ResiduePairAnnotationAnnotator()
        annotated = annotator.annotate_pairs(
            pairs=pairs,
            chain_mapping=context.metadata["chain_mapping"],
            cdr_detection=context.metadata.get("cdr_detection"),
        )
        output_dir = Path(context.get_analysis_path("interactions/pi_interactions", "annotated_pi_pi_pairs.csv")).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        annotated_file = str(output_dir / "annotated_pi_pi_internal.csv")
        annotated.to_csv(annotated_file, index=False)
        report_table = ContactAnnotationNode._build_contact_report_table(annotated)
        report_file = str(output_dir / "pi_pi_report.csv")
        report_table.to_csv(report_file, index=False)
        groove_tcr = annotator.filter_pairs(annotated, interaction_class="hla_tcr", mhc_region="groove")
        groove_file = str(output_dir / "groove_tcr_pi_pi.csv")
        ContactAnnotationNode._build_contact_report_table(groove_tcr).to_csv(groove_file, index=False)
        region_builder = RegionInteractionSummaryBuilder(annotator=annotator)
        cdr1 = region_builder.write_region_bundle(output_dir, "cdr1", annotator.filter_pairs(annotated, tcr_region="CDR1"), ContactAnnotationNode._build_contact_report_table)
        cdr2 = region_builder.write_region_bundle(output_dir, "cdr2", annotator.filter_pairs(annotated, tcr_region="CDR2"), ContactAnnotationNode._build_contact_report_table)
        cdr3 = region_builder.write_region_bundle(output_dir, "cdr3", annotator.filter_pairs(annotated, tcr_region="CDR3"), ContactAnnotationNode._build_contact_report_table)
        non_cdr = region_builder.write_region_bundle(output_dir, "non_cdr", annotator.filter_pairs(annotated, tcr_region="non_cdr"), ContactAnnotationNode._build_contact_report_table)
        summary_file = str(output_dir / "pi_pi_annotation_summary.json")
        with open(summary_file, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "n_total_pairs": int(len(annotated)),
                    "n_cdr1_pairs": int(len(annotator.filter_pairs(annotated, tcr_region="CDR1"))),
                    "n_cdr2_pairs": int(len(annotator.filter_pairs(annotated, tcr_region="CDR2"))),
                    "n_cdr3_pairs": int(len(annotator.filter_pairs(annotated, tcr_region="CDR3"))),
                    "n_non_cdr_pairs": int(len(annotator.filter_pairs(annotated, tcr_region="non_cdr"))),
                },
                handle,
                indent=2,
            )
        context.results["pi_pi_annotation"] = {
            "annotated_pairs_file": annotated_file,
            "pi_pi_report_file": report_file,
            "groove_tcr_pi_pi_file": groove_file,
            "summary_file": summary_file,
            "cdr1": cdr1,
            "cdr2": cdr2,
            "cdr3": cdr3,
            "non_cdr": non_cdr,
        }
        return context
