"""Pair annotation and region-summary node for hydrogen-bond results."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from ....analysis.topology import RegionInteractionSummaryBuilder, ResiduePairAnnotationAnnotator
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError
from .contact_annotation_node import ContactAnnotationNode


class HydrogenBondAnnotationNode(PipelineNode):
    """Add semantic annotation to residue-level hydrogen-bond results and write region views."""

    def __init__(self, name: str | None = None):
        super().__init__(name=name or "HydrogenBondAnnotationNode")

    def validate_inputs(self, context: PipelineContext):
        if "hbond_pairs" not in context.results:
            raise PipelineError(
                node_name=self.name,
                reason="hbond_pairs results not found in context",
                context_state={"system_id": context.system_id},
            )
        if "chain_mapping" not in context.metadata:
            raise PipelineError(
                node_name=self.name,
                reason="chain_mapping metadata not found in context",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        raw_file = context.results["hbond_pairs"]["raw_pairs_file"]
        pairs = pd.read_csv(raw_file)

        annotator = ResiduePairAnnotationAnnotator()
        annotated = annotator.annotate_pairs(
            pairs=pairs,
            chain_mapping=context.metadata["chain_mapping"],
            cdr_detection=context.metadata.get("cdr_detection"),
        )

        output_dir = Path(context.get_analysis_path("interactions/hydrogen_bonds", "annotated_hbond_pairs.csv")).parent
        output_dir.mkdir(parents=True, exist_ok=True)

        annotated_file = str(output_dir / "annotated_hbond_internal.csv")
        annotated.to_csv(annotated_file, index=False)

        report_table = ContactAnnotationNode._build_contact_report_table(annotated)
        report_file = str(output_dir / "hbond_report.csv")
        report_table.to_csv(report_file, index=False)

        groove_tcr = annotator.filter_pairs(annotated, interaction_class="hla_tcr", mhc_region="groove")
        groove_file = str(output_dir / "groove_tcr_hbonds.csv")
        ContactAnnotationNode._build_contact_report_table(groove_tcr).to_csv(groove_file, index=False)

        cdr1_pairs = annotator.filter_pairs(annotated, tcr_region="CDR1")
        cdr2_pairs = annotator.filter_pairs(annotated, tcr_region="CDR2")
        cdr3_pairs = annotator.filter_pairs(annotated, tcr_region="CDR3")
        non_cdr_pairs = annotator.filter_pairs(annotated, tcr_region="non_cdr")

        region_builder = RegionInteractionSummaryBuilder(annotator=annotator)
        cdr1_outputs = region_builder.write_region_bundle(output_dir, "cdr1", cdr1_pairs, ContactAnnotationNode._build_contact_report_table)
        cdr2_outputs = region_builder.write_region_bundle(output_dir, "cdr2", cdr2_pairs, ContactAnnotationNode._build_contact_report_table)
        cdr3_outputs = region_builder.write_region_bundle(output_dir, "cdr3", cdr3_pairs, ContactAnnotationNode._build_contact_report_table)
        non_cdr_outputs = region_builder.write_region_bundle(output_dir, "non_cdr", non_cdr_pairs, ContactAnnotationNode._build_contact_report_table)

        summary = {
            "n_total_pairs": int(len(annotated)),
            "n_cdr1_pairs": int(len(cdr1_pairs)),
            "n_cdr2_pairs": int(len(cdr2_pairs)),
            "n_cdr3_pairs": int(len(cdr3_pairs)),
            "n_non_cdr_pairs": int(len(non_cdr_pairs)),
        }
        summary_file = str(output_dir / "hbond_annotation_summary.json")
        with open(summary_file, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2)

        context.results["hbond_annotation"] = {
            "annotated_pairs_file": annotated_file,
            "hbond_report_file": report_file,
            "groove_tcr_hbonds_file": groove_file,
            "summary_file": summary_file,
            "cdr1": cdr1_outputs,
            "cdr2": cdr2_outputs,
            "cdr3": cdr3_outputs,
            "non_cdr": non_cdr_outputs,
        }
        return context
