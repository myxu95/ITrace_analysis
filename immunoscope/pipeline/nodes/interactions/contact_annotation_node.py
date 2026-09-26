"""
Pipeline node for semantic annotation of residue contact tables.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from ....analysis.topology import RegionInteractionSummaryBuilder, ResiduePairAnnotationAnnotator
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class ContactAnnotationNode(PipelineNode):
    """Add pHLA/TCR and CDR annotations to contact frequency tables."""

    def __init__(self, name: str | None = None):
        super().__init__(name=name or "ContactAnnotationNode")

    def validate_inputs(self, context: PipelineContext):
        if "contact_frequency" not in context.results:
            raise PipelineError(
                node_name=self.name,
                reason="contact_frequency results not found in context",
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

        raw_file = context.results["contact_frequency"]["raw_contacts_file"]
        contacts = pd.read_csv(raw_file)

        annotator = ResiduePairAnnotationAnnotator()
        annotated = annotator.annotate_pairs(
            pairs=contacts,
            chain_mapping=context.metadata["chain_mapping"],
            cdr_detection=context.metadata.get("cdr_detection"),
        )

        output_dir = Path(context.get_analysis_path("contacts", "annotated_contact_frequencies.csv")).parent
        output_dir.mkdir(parents=True, exist_ok=True)

        annotated_file = str(output_dir / "annotated_contact_internal.csv")
        annotated.to_csv(annotated_file, index=False)

        report_table = self._build_contact_report_table(annotated)
        report_file = str(output_dir / "contact_report.csv")
        report_table.to_csv(report_file, index=False)

        annotator = ResiduePairAnnotationAnnotator()
        groove_tcr = annotator.filter_pairs(
            annotated,
            interaction_class="hla_tcr",
            mhc_region="groove",
        )
        groove_tcr_file = str(output_dir / "groove_tcr_contacts.csv")
        self._build_contact_report_table(groove_tcr).to_csv(groove_tcr_file, index=False)
        cdr1_contacts = annotator.filter_pairs(annotated, tcr_region="CDR1")
        cdr2_contacts = annotator.filter_pairs(annotated, tcr_region="CDR2")
        cdr3_contacts = annotator.filter_pairs(annotated, tcr_region="CDR3")
        non_cdr_contacts = annotator.filter_pairs(annotated, tcr_region="non_cdr")

        region_builder = RegionInteractionSummaryBuilder(annotator=annotator)
        cdr1_outputs = region_builder.write_region_bundle(output_dir, "cdr1", cdr1_contacts, self._build_contact_report_table)
        cdr2_outputs = region_builder.write_region_bundle(output_dir, "cdr2", cdr2_contacts, self._build_contact_report_table)
        cdr3_outputs = region_builder.write_region_bundle(output_dir, "cdr3", cdr3_contacts, self._build_contact_report_table)
        non_cdr_outputs = region_builder.write_region_bundle(output_dir, "non_cdr", non_cdr_contacts, self._build_contact_report_table)

        summary = {
            "n_total_contacts": int(len(annotated)),
            "n_cdr1_contacts": int(len(cdr1_contacts)),
            "n_cdr2_contacts": int(len(cdr2_contacts)),
            "n_cdr3_contacts": int(len(cdr3_contacts)),
            "n_non_cdr_contacts": int(len(non_cdr_contacts)),
        }
        summary_file = str(output_dir / "contact_annotation_summary.json")
        with open(summary_file, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2)

        context.results["contact_annotation"] = {
            "annotated_contacts_file": annotated_file,
            "contact_report_file": report_file,
            "groove_tcr_contacts_file": groove_tcr_file,
            "summary_file": summary_file,
            "cdr1": cdr1_outputs,
            "cdr2": cdr2_outputs,
            "cdr3": cdr3_outputs,
            "non_cdr": non_cdr_outputs,
        }
        return context

    @staticmethod
    def _build_contact_report_table(annotated: pd.DataFrame) -> pd.DataFrame:
        if annotated.empty:
            return pd.DataFrame(
                columns=[
                    "interaction_family",
                    "interaction_class",
                    "phla_region",
                    "mhc_region",
                    "mhc_subregion",
                    "phla_chain_id",
                    "phla_resid",
                    "phla_resname",
                    "phla_residue",
                    "tcr_chain_id",
                    "tcr_chain",
                    "tcr_resid",
                    "tcr_resname",
                    "tcr_residue",
                    "tcr_region",
                    "tcr_region_detailed",
                    "contact_frequency",
                    "contact_frames",
                    "total_frames",
                    "n_segments",
                    "max_consecutive_frames",
                    "max_consecutive_fraction",
                    "first_frame",
                    "last_frame",
                    "frame_segments",
                    "min_distance_angstrom",
                ]
            )

        rows = []
        for _, row in annotated.iterrows():
            if row["complex_side_1"] == "tcr":
                tcr_side = 1
                partner_side = 2
            elif row["complex_side_2"] == "tcr":
                tcr_side = 2
                partner_side = 1
            else:
                continue

            rows.append(
                {
                    "interaction_family": row.get("interaction_family", "contact"),
                    "interaction_class": row["interaction_class"],
                    "phla_region": row["phla_region"],
                    "mhc_region": row.get("mhc_region"),
                    "mhc_subregion": row.get("mhc_subregion"),
                    "phla_chain_id": row[f"chain_id_{partner_side}"],
                    "phla_resid": int(row[f"resid_{partner_side}"]),
                    "phla_resname": row[f"resname_{partner_side}"],
                    "phla_residue": row[f"residue_label_{partner_side}"],
                    "tcr_chain_id": row[f"chain_id_{tcr_side}"],
                    "tcr_chain": row["tcr_chain"],
                    "tcr_resid": int(row[f"resid_{tcr_side}"]),
                    "tcr_resname": row[f"resname_{tcr_side}"],
                    "tcr_residue": row[f"residue_label_{tcr_side}"],
                    "tcr_region": row["tcr_region"],
                    "tcr_region_detailed": row["tcr_region_detailed"],
                    "contact_frequency": row["contact_frequency"],
                    "contact_frames": int(row["contact_frames"]),
                    "total_frames": int(row["total_frames"]),
                    "n_segments": int(row["n_segments"]) if "n_segments" in row and pd.notna(row["n_segments"]) else None,
                    "max_consecutive_frames": int(row["max_consecutive_frames"]) if "max_consecutive_frames" in row and pd.notna(row["max_consecutive_frames"]) else None,
                    "max_consecutive_fraction": (
                        float(row["max_consecutive_frames"]) / float(row["total_frames"])
                        if "max_consecutive_frames" in row and pd.notna(row["max_consecutive_frames"]) and float(row["total_frames"]) > 0
                        else None
                    ),
                    "first_frame": int(row["first_frame"]) if "first_frame" in row and pd.notna(row["first_frame"]) else None,
                    "last_frame": int(row["last_frame"]) if "last_frame" in row and pd.notna(row["last_frame"]) else None,
                    "frame_segments": row.get("frame_segments"),
                    "min_distance_angstrom": row["min_distance_observed"],
                }
            )

        return pd.DataFrame(rows).sort_values(
            by=["contact_frequency", "phla_region", "tcr_chain", "tcr_region"],
            ascending=[False, True, True, True],
        ).reset_index(drop=True)
