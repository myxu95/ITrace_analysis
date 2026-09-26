"""
Generic residue-pair semantic annotator.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import pandas as pd

from .complex_residue_semantics import ComplexResidueSemanticAnnotator


class ResiduePairAnnotationAnnotator:
    """
    Add pHLA/TCR semantic labels to any residue-pair result table.

    Currently used mainly for contact pairs; typed interactions should reuse the same fields.
    """

    EMPTY_COLUMNS = [
        "component_1",
        "complex_side_1",
        "phla_region_1",
        "mhc_region_1",
        "mhc_subregion_1",
        "tcr_chain_1",
        "tcr_region_1",
        "tcr_region_detailed_1",
        "component_2",
        "complex_side_2",
        "phla_region_2",
        "mhc_region_2",
        "mhc_subregion_2",
        "tcr_chain_2",
        "tcr_region_2",
        "tcr_region_detailed_2",
        "interaction_class",
        "tcr_chain",
        "tcr_region",
        "tcr_region_detailed",
        "phla_region",
        "mhc_region",
        "mhc_subregion",
        "tcr_residue_label",
        "partner_residue_label",
        "partner_component",
    ]

    def annotate_pairs(
        self,
        pairs: pd.DataFrame,
        chain_mapping: Dict[str, str],
        cdr_detection: Optional[Dict] = None,
    ) -> pd.DataFrame:
        annotated = pairs.copy()
        if "interaction_family" not in annotated.columns:
            annotated["interaction_family"] = "contact"
        if annotated.empty:
            for column in self.EMPTY_COLUMNS:
                annotated[column] = pd.Series(dtype="object")
            return annotated

        semantic_annotator = ComplexResidueSemanticAnnotator(
            chain_mapping=chain_mapping,
            cdr_detection=cdr_detection,
        )

        side1 = annotated.apply(
            lambda row: self._annotate_residue(
                chain_id=row["chain_id_1"],
                resid=int(row["resid_1"]),
                semantic_annotator=semantic_annotator,
            ),
            axis=1,
            result_type="expand",
        )
        side1.columns = [
            "component_1",
            "complex_side_1",
            "phla_region_1",
            "mhc_region_1",
            "mhc_subregion_1",
            "tcr_chain_1",
            "tcr_region_1",
            "tcr_region_detailed_1",
        ]

        side2 = annotated.apply(
            lambda row: self._annotate_residue(
                chain_id=row["chain_id_2"],
                resid=int(row["resid_2"]),
                semantic_annotator=semantic_annotator,
            ),
            axis=1,
            result_type="expand",
        )
        side2.columns = [
            "component_2",
            "complex_side_2",
            "phla_region_2",
            "mhc_region_2",
            "mhc_subregion_2",
            "tcr_chain_2",
            "tcr_region_2",
            "tcr_region_detailed_2",
        ]

        annotated = pd.concat([annotated, side1, side2], axis=1)

        interaction_info = annotated.apply(
            self._derive_interaction_semantics,
            axis=1,
            result_type="expand",
        )
        interaction_info.columns = [
            "interaction_class",
            "tcr_chain",
            "tcr_region",
            "tcr_region_detailed",
            "phla_region",
            "mhc_region",
            "mhc_subregion",
            "tcr_residue_label",
            "partner_residue_label",
            "partner_component",
        ]

        annotated = pd.concat([annotated, interaction_info], axis=1)
        return annotated

    @staticmethod
    def filter_pairs(
        pairs: pd.DataFrame,
        interaction_class: Optional[str] = None,
        tcr_region: Optional[str] = None,
        partner_component: Optional[str] = None,
        mhc_region: Optional[str] = None,
        mhc_subregion: Optional[str] = None,
    ) -> pd.DataFrame:
        filtered = pairs.copy()
        if interaction_class is not None:
            filtered = filtered[filtered["interaction_class"] == interaction_class]
        if tcr_region is not None:
            filtered = filtered[filtered["tcr_region"] == tcr_region]
        if partner_component is not None:
            filtered = filtered[filtered["partner_component"] == partner_component]
        if mhc_region is not None:
            filtered = filtered[filtered["mhc_region"] == mhc_region]
        if mhc_subregion is not None:
            filtered = filtered[filtered["mhc_subregion"] == mhc_subregion]
        return filtered.reset_index(drop=True)

    @staticmethod
    def summarize_pair_partners(
        pairs: pd.DataFrame,
        label_column: str,
    ) -> pd.DataFrame:
        if pairs.empty:
            return pd.DataFrame(
                columns=[label_column, "contact_frequency_sum", "contact_frequency_max", "n_pairs"]
            )

        return (
            pairs.groupby(label_column)["contact_frequency"]
            .agg(["sum", "max", "count"])
            .reset_index()
            .rename(
                columns={
                    "sum": "contact_frequency_sum",
                    "max": "contact_frequency_max",
                    "count": "n_pairs",
                }
            )
            .sort_values(
                by=["contact_frequency_sum", "contact_frequency_max"],
                ascending=[False, False],
            )
            .reset_index(drop=True)
        )

    @staticmethod
    def _annotate_residue(
        chain_id: str,
        resid: int,
        semantic_annotator: ComplexResidueSemanticAnnotator,
    ) -> Tuple[str, str, Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]]:
        semantics = semantic_annotator.annotate_residue(chain_id=chain_id, resid=resid)
        return (
            semantics.component,
            semantics.complex_side,
            semantics.phla_region,
            semantics.mhc_region,
            semantics.mhc_subregion,
            semantics.tcr_chain,
            semantics.tcr_region,
            semantics.tcr_region_detailed,
        )

    # Canonical naming order for non-TCR sub-interfaces: more-designable side first.
    # Drives both _component_key normalization and the f"{a}_{b}" class name spelling.
    _COMPONENT_KEY_ORDER: Tuple[str, ...] = ("peptide", "hla", "beta2m", "tcr", "other")

    @staticmethod
    def _component_key(component: Optional[str]) -> str:
        """Normalize a free-form component label to a taxonomy key."""
        if not component:
            return "other"
        c = component.lower()
        if "peptide" in c:
            return "peptide"
        if "hla" in c or "mhc" in c:
            return "hla"
        if "beta2m" in c or "b2m" in c:
            return "beta2m"
        if "tcr" in c:
            return "tcr"
        return "other"

    @staticmethod
    def _derive_interaction_semantics(row) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]]:
        """Classify a residue pair and emit TCR-centric labels when applicable.

        interaction_class taxonomy:
          - peptide_tcr / hla_tcr / beta2m_tcr / other_tcr
                pair contains exactly one TCR side; tcr_* and partner_* are populated.
          - peptide_hla / peptide_beta2m / hla_beta2m
                pHLA sub-interfaces (no TCR side); tcr_* are None, partner_* point
                at the non-subject side. Subject is the more-designable component
                per _COMPONENT_KEY_ORDER (peptide < hla < beta2m).
          - intra_peptide / intra_hla / intra_beta2m / intra_tcr
                both sides are the same component; partner_* are None.
          - other
                catch-all for unclassifiable pairs.

        Prior to 2026-05-27 non-TCR pairs returned all-None and were therefore
        invisible to every downstream consumer (data_formatter, contact_heatmap,
        comparison_builder, …). That blocked peptide-HLA / presentation-task
        analysis even though the underlying RRCS / contact tables already
        contained the data. Non-TCR pairs are now classified so the rows survive
        downstream filters, with tcr_* left None (we don't fabricate TCR labels
        for non-TCR data).
        """
        side1_is_tcr = row["complex_side_1"] == "tcr"
        side2_is_tcr = row["complex_side_2"] == "tcr"

        if side1_is_tcr and not side2_is_tcr:
            tcr_side = 1
            partner_side = 2
        elif side2_is_tcr and not side1_is_tcr:
            tcr_side = 2
            partner_side = 1
        else:
            return ResiduePairAnnotationAnnotator._derive_non_tcr_semantics(row)

        partner_component = row[f"component_{partner_side}"]
        partner_phla_region = row[f"phla_region_{partner_side}"]
        partner_mhc_region = row[f"mhc_region_{partner_side}"]
        partner_mhc_subregion = row[f"mhc_subregion_{partner_side}"]

        if partner_component == "peptide":
            interaction_class = "peptide_tcr"
            phla_region = "peptide"
        elif partner_component == "HLA_alpha":
            interaction_class = "hla_tcr"
            phla_region = "HLA_alpha"
        elif partner_component == "beta2m":
            interaction_class = "beta2m_tcr"
            phla_region = "beta2m"
        else:
            interaction_class = "other_tcr"
            phla_region = partner_component

        return (
            interaction_class,
            row[f"tcr_chain_{tcr_side}"],
            row[f"tcr_region_{tcr_side}"],
            row[f"tcr_region_detailed_{tcr_side}"],
            partner_phla_region or phla_region,
            partner_mhc_region,
            partner_mhc_subregion,
            row[f"residue_label_{tcr_side}"],
            row[f"residue_label_{partner_side}"],
            partner_component,
        )

    @staticmethod
    def _derive_non_tcr_semantics(row) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]]:
        """Classify pairs without a TCR residue (peptide-HLA, intra-X, etc.).

        These rows were previously dropped (all-None return). They are now
        emitted with an interaction_class so downstream filters and groupby
        aggregations can see them. tcr_* columns stay None (no TCR side exists);
        for asymmetric sub-interfaces (peptide_hla, peptide_beta2m, hla_beta2m)
        the partner_* columns point at the non-subject side.
        """
        order = ResiduePairAnnotationAnnotator._COMPONENT_KEY_ORDER
        k1 = ResiduePairAnnotationAnnotator._component_key(row.get("component_1"))
        k2 = ResiduePairAnnotationAnnotator._component_key(row.get("component_2"))

        if k1 == k2:
            interaction_class = f"intra_{k1}" if k1 != "other" else "other"
            return (interaction_class, None, None, None, None, None, None, None, None, None)

        i1, i2 = order.index(k1), order.index(k2)
        if i1 < i2:
            subject_side, partner_side = 1, 2
            first, second = k1, k2
        else:
            subject_side, partner_side = 2, 1
            first, second = k2, k1
        interaction_class = f"{first}_{second}"

        return (
            interaction_class,
            None,  # tcr_chain: no TCR side
            None,  # tcr_region
            None,  # tcr_region_detailed
            row.get(f"phla_region_{subject_side}"),
            row.get(f"mhc_region_{partner_side}"),
            row.get(f"mhc_subregion_{partner_side}"),
            None,  # tcr_residue_label: no TCR side
            row.get(f"residue_label_{partner_side}"),
            row.get(f"component_{partner_side}"),
        )


class ContactAnnotationAnnotator(ResiduePairAnnotationAnnotator):
    """Compatibility alias; typed interactions and contacts use the shared pair annotator."""

    def annotate_contacts(
        self,
        contacts: pd.DataFrame,
        chain_mapping: Dict[str, str],
        cdr_detection: Optional[Dict] = None,
    ) -> pd.DataFrame:
        return self.annotate_pairs(contacts, chain_mapping=chain_mapping, cdr_detection=cdr_detection)

    def filter_contacts(
        self,
        contacts: pd.DataFrame,
        interaction_class: Optional[str] = None,
        tcr_region: Optional[str] = None,
        partner_component: Optional[str] = None,
    ) -> pd.DataFrame:
        return self.filter_pairs(
            contacts,
            interaction_class=interaction_class,
            tcr_region=tcr_region,
            partner_component=partner_component,
        )

    def summarize_contact_partners(self, contacts: pd.DataFrame, label_column: str) -> pd.DataFrame:
        return self.summarize_pair_partners(contacts, label_column)
