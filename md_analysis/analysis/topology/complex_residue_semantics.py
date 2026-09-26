"""
Complex-wide residue semantic annotation utilities.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

from .tcr_residue_semantics import TCRResidueSemanticAnnotator, TCRResidueSemantics


@dataclass(frozen=True)
class ComplexResidueSemantics:
    """Semantic identity of a residue in the full pHLA-TCR complex."""

    chain_id: str
    resid: int
    component: str
    complex_side: str
    phla_region: Optional[str]
    mhc_region: Optional[str]
    mhc_subregion: Optional[str]
    tcr_chain: Optional[str]
    tcr_region: Optional[str]
    tcr_region_detailed: Optional[str]


class ComplexResidueSemanticAnnotator:
    """
    Unified residue semantic annotator for the full pHLA-TCR complex.
    """

    # Use a lightweight groove definition for class I MHC heavy chains to support region-level analysis.
    HLA_ALPHA1_RANGE = range(50, 87)
    HLA_ALPHA2_RANGE = range(140, 177)

    def __init__(
        self,
        chain_mapping: Dict[str, str],
        cdr_detection: Optional[Dict] = None,
    ):
        self.tcr_annotator = TCRResidueSemanticAnnotator(
            chain_mapping=chain_mapping,
            cdr_detection=cdr_detection,
        )

    def annotate_residue(self, chain_id: str, resid: int) -> ComplexResidueSemantics:
        tcr_semantics: TCRResidueSemantics = self.tcr_annotator.annotate_residue(chain_id, resid)
        component = tcr_semantics.component

        if component in {"TCR_alpha", "TCR_beta"}:
            complex_side = "tcr"
            phla_region = None
            mhc_region = None
            mhc_subregion = None
        elif component in {"HLA_alpha", "beta2m", "peptide"}:
            complex_side = "phla"
            phla_region = component
            mhc_region = self._infer_mhc_region(component=component, resid=int(resid))
            mhc_subregion = self._infer_mhc_subregion(component=component, resid=int(resid))
        else:
            complex_side = "unknown"
            phla_region = None
            mhc_region = None
            mhc_subregion = None

        return ComplexResidueSemantics(
            chain_id=chain_id,
            resid=int(resid),
            component=component,
            complex_side=complex_side,
            phla_region=phla_region,
            mhc_region=mhc_region,
            mhc_subregion=mhc_subregion,
            tcr_chain=tcr_semantics.tcr_chain,
            tcr_region=tcr_semantics.tcr_region,
            tcr_region_detailed=tcr_semantics.tcr_region_detailed,
        )

    @classmethod
    def _infer_mhc_region(cls, component: str, resid: int) -> Optional[str]:
        if component == "peptide":
            return "peptide"
        if component == "beta2m":
            return "beta2m"
        if component != "HLA_alpha":
            return None
        if resid in cls.HLA_ALPHA1_RANGE or resid in cls.HLA_ALPHA2_RANGE:
            return "groove"
        return "non_groove"

    @classmethod
    def _infer_mhc_subregion(cls, component: str, resid: int) -> Optional[str]:
        if component == "peptide":
            return "peptide"
        if component == "beta2m":
            return "beta2m"
        if component != "HLA_alpha":
            return None
        if resid in cls.HLA_ALPHA1_RANGE:
            return "alpha1_helix"
        if resid in cls.HLA_ALPHA2_RANGE:
            return "alpha2_helix"
        return "non_groove"
