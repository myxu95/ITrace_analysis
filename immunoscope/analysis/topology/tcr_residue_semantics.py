"""
TCR residue-level semantic annotation utilities.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple


@dataclass(frozen=True)
class TCRResidueSemantics:
    """Semantic identity of a residue in the pHLA-TCR complex."""

    chain_id: str
    resid: int
    component: str
    tcr_chain: Optional[str]
    tcr_region: Optional[str]
    tcr_region_detailed: Optional[str]


class TCRResidueSemanticAnnotator:
    """
    Annotate chain/residue pairs with TCR-specific semantics.

    Current scope:
    - component assignment: HLA_alpha, beta2m, peptide, TCR_alpha, TCR_beta
    - TCR chain assignment: alpha / beta
    - TCR region assignment: CDR1 / CDR2 / CDR3 / non_cdr
    """

    def __init__(
        self,
        chain_mapping: Dict[str, str],
        cdr_detection: Optional[Dict] = None,
    ):
        self.chain_mapping = chain_mapping
        self.cdr_detection = cdr_detection or {}
        self.component_chain_map = self._component_chain_map(chain_mapping)
        self.cdr_lookup = self._build_cdr_lookup(self.cdr_detection)

    def annotate_residue(self, chain_id: str, resid: int) -> TCRResidueSemantics:
        component_key = self.component_chain_map.get(chain_id)
        component = self._human_component_name(component_key)

        if component_key == "tcr_alpha":
            tcr_chain = "alpha"
        elif component_key == "tcr_beta":
            tcr_chain = "beta"
        else:
            tcr_chain = None

        tcr_region = None
        tcr_region_detailed = None
        if tcr_chain is not None:
            tcr_region, tcr_region_detailed = self.cdr_lookup.get(
                (chain_id, int(resid)),
                ("non_cdr", f"non_cdr_{tcr_chain}"),
            )

        return TCRResidueSemantics(
            chain_id=chain_id,
            resid=int(resid),
            component=component,
            tcr_chain=tcr_chain,
            tcr_region=tcr_region,
            tcr_region_detailed=tcr_region_detailed,
        )

    @staticmethod
    def _component_chain_map(chain_mapping: Dict[str, str]) -> Dict[str, str]:
        return {
            chain_id: component
            for component, chain_id in chain_mapping.items()
            if chain_id
        }

    @staticmethod
    def _build_cdr_lookup(cdr_detection: Dict) -> Dict[Tuple[str, int], Tuple[str, str]]:
        lookup: Dict[Tuple[str, int], Tuple[str, str]] = {}
        chain_name_map = {
            "TCR_alpha": "alpha",
            "TCR_beta": "beta",
        }

        for chain_name, chain_info in cdr_detection.get("chains", {}).items():
            suffix = chain_name_map.get(chain_name)
            if suffix is None:
                continue

            chain_id = chain_info.get("chain_id")
            if not chain_id:
                continue

            for cdr_num, cdr_data in chain_info.get("cdrs", {}).items():
                residue_range = cdr_data.get("residue_range")
                if not residue_range:
                    continue
                start, end = residue_range
                region = f"CDR{cdr_num}"
                detailed = f"CDR{cdr_num}_{suffix}"
                for resid in range(int(start), int(end) + 1):
                    lookup[(chain_id, resid)] = (region, detailed)

        return lookup

    @staticmethod
    def _human_component_name(component_key: Optional[str]) -> str:
        mapping = {
            "mhc_alpha": "HLA_alpha",
            "b2m": "beta2m",
            "peptide": "peptide",
            "tcr_alpha": "TCR_alpha",
            "tcr_beta": "TCR_beta",
        }
        return mapping.get(component_key, "unknown")
