"""
Region-level interaction filtering and summarization.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .residue_pair_annotation import ResiduePairAnnotationAnnotator


class RegionInteractionSummaryBuilder:
    """Build region-level outputs from annotated pair tables."""

    def __init__(self, annotator: ResiduePairAnnotationAnnotator | None = None):
        self.annotator = annotator or ResiduePairAnnotationAnnotator()

    def write_region_bundle(
        self,
        contacts_root: Path,
        region_name: str,
        region_pairs: pd.DataFrame,
        report_builder,
    ) -> dict:
        region_dir = contacts_root / region_name
        tables_dir = region_dir / "tables"
        summaries_dir = region_dir / "summaries"
        heatmaps_dir = region_dir / "heatmaps"

        tables_dir.mkdir(parents=True, exist_ok=True)
        summaries_dir.mkdir(parents=True, exist_ok=True)
        heatmaps_dir.mkdir(parents=True, exist_ok=True)

        region_report = report_builder(region_pairs)
        all_contacts_file = str(tables_dir / "contacts.csv")
        region_report.to_csv(all_contacts_file, index=False)

        peptide_tcr = self.annotator.filter_pairs(region_pairs, interaction_class="peptide_tcr")
        hla_tcr = self.annotator.filter_pairs(region_pairs, interaction_class="hla_tcr")
        groove_tcr = self.annotator.filter_pairs(
            region_pairs,
            interaction_class="hla_tcr",
            mhc_region="groove",
        )

        peptide_tcr_file = str(tables_dir / "peptide_tcr_contacts.csv")
        hla_tcr_file = str(tables_dir / "hla_tcr_contacts.csv")
        groove_tcr_file = str(tables_dir / "groove_tcr_contacts.csv")
        report_builder(peptide_tcr).to_csv(peptide_tcr_file, index=False)
        report_builder(hla_tcr).to_csv(hla_tcr_file, index=False)
        report_builder(groove_tcr).to_csv(groove_tcr_file, index=False)

        tcr_summary = self.annotator.summarize_pair_partners(region_pairs, "tcr_residue_label")
        partner_summary = self.annotator.summarize_pair_partners(region_pairs, "partner_residue_label")

        tcr_summary_file = str(summaries_dir / "tcr_residue_summary.csv")
        partner_summary_file = str(summaries_dir / "partner_residue_summary.csv")
        tcr_summary.to_csv(tcr_summary_file, index=False)
        partner_summary.to_csv(partner_summary_file, index=False)

        bundle_summary = {
            "region": region_name,
            "n_contacts": int(len(region_pairs)),
            "n_peptide_tcr_contacts": int(len(peptide_tcr)),
            "n_hla_tcr_contacts": int(len(hla_tcr)),
            "n_groove_tcr_contacts": int(len(groove_tcr)),
        }
        summary_file = str(summaries_dir / "summary.json")
        with open(summary_file, "w", encoding="utf-8") as handle:
            json.dump(bundle_summary, handle, indent=2)

        return {
            "region_dir": str(region_dir),
            "tables_dir": str(tables_dir),
            "summaries_dir": str(summaries_dir),
            "heatmaps_dir": str(heatmaps_dir),
            "contacts_file": all_contacts_file,
            "peptide_tcr_contacts_file": peptide_tcr_file,
            "hla_tcr_contacts_file": hla_tcr_file,
            "groove_tcr_contacts_file": groove_tcr_file,
            "tcr_residue_summary_file": tcr_summary_file,
            "partner_residue_summary_file": partner_summary_file,
            "summary_file": summary_file,
        }
