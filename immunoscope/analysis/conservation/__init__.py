"""TCR germline evolutionary-conservation analysis (v1, 2026-05-28).

Computes a per-residue conservation score for the TCR alpha/beta chains
of a complex, grounded in the IMGT germline V-gene reference set bundled
with ANARCI. The IMGT-gapped germline sequences are pre-aligned by
position (the gaps *are* the alignment), so per-position conservation is
a column statistic — no MSA step is needed.

Scope (v1): germline V-REGION only. The reference panel is the whole
human germline V-gene set for the locus (all TRAV genes for alpha, all
TRBV for beta) — a phylogenetic conservation signal, NOT patient- or
allele-specific. Framework (FR) and CDR1/CDR2 are fully covered; the
germline V gene encodes the CDR3 *anchor* (IMGT ~105-107) but not the
somatic CDR3 *junction* (IMGT ~108-117), so junction residues are
reported `covered=False` rather than given a misleading score.

See `tcr_germline_conservation.py` for the algorithm.
"""

from __future__ import annotations

from .tcr_germline_conservation import (
    ConservationResidue,
    NoAnarciError,
    build_reference_profile,
    compute_chain_conservation,
    imgt_region,
)

__all__ = [
    "ConservationResidue",
    "NoAnarciError",
    "build_reference_profile",
    "compute_chain_conservation",
    "imgt_region",
]
