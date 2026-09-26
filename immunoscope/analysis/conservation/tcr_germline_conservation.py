"""Per-residue TCR germline conservation from the IMGT V-gene reference.

Algorithm
---------
1. Reference profile: take ANARCI's bundled human germline V genes for a
   locus (TRAV → chain class 'A', TRBV → 'B'). Each germline sequence is
   IMGT-gap-aligned to a fixed 128-column layout, so column *i* is IMGT
   position *i+1* for every gene. For each column we collect the non-gap
   amino acids across genes and score conservation as

       conservation = 1 - H / ln(20)

   where H is the Shannon entropy of the column's amino-acid frequency
   distribution. 1.0 = perfectly conserved, ~0 = uniformly variable.
   Columns with no germline coverage get score=None and are reported
   `covered=False`. The germline V gene encodes the CDR3 *anchor* (IMGT
   ~105-107, conserved) but not the somatic CDR3 junction (IMGT ~108-117),
   so junction residues are uncovered while the anchor is scored.

2. Per-residue mapping: ANARCI numbers the chain's own sequence, giving
   each real residue an IMGT position. We map structure residue →
   IMGT position → reference column score, and tag the IMGT region.

The entropy / profile maths are pure functions (no ANARCI needed) so
they are unit-testable offline; only the chain numbering and the
reference lookup touch ANARCI.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

log = logging.getLogger("immunoscope.analysis.conservation")

# Standard IMGT V-domain region boundaries (inclusive), unique numbering.
_FR1 = (1, 26)
_CDR1 = (27, 38)
_FR2 = (39, 55)
_CDR2 = (56, 65)
_FR3 = (66, 104)
_CDR3 = (105, 117)
_FR4 = (118, 128)

# TCR locus → ANARCI germline chain class.
_LOCUS_CLASS = {"alpha": "A", "beta": "B"}

# Conservation-class thresholds on the [0,1] score.
_CONSERVED_MIN = 0.80
_INTERMEDIATE_MIN = 0.50

# A reference column needs at least this fraction of non-gap residues to
# be considered "covered" (guards CDR3 / sparse tail columns).
_MIN_COLUMN_FILL = 0.50

_MAX_ENTROPY = math.log(20.0)


class NoAnarciError(RuntimeError):
    """Raised when ANARCI (or its germline data) is unavailable."""


@dataclass(frozen=True)
class ConservationResidue:
    """One structure residue's germline-conservation record."""

    chain_id: str
    resid: int
    resname: str
    imgt_position: Optional[int]
    region: str
    conservation_score: Optional[float]
    conservation_class: str  # conserved / intermediate / variable / uncovered
    covered: bool


# ---------------------------------------------------------------------------
# Pure helpers (no ANARCI) — unit-testable offline
# ---------------------------------------------------------------------------


def imgt_region(imgt_position: Optional[int]) -> str:
    """Return the IMGT V-domain region label for a position (or 'unknown')."""
    if imgt_position is None:
        return "unknown"
    p = imgt_position
    if _FR1[0] <= p <= _FR1[1]:
        return "FR1"
    if _CDR1[0] <= p <= _CDR1[1]:
        return "CDR1"
    if _FR2[0] <= p <= _FR2[1]:
        return "FR2"
    if _CDR2[0] <= p <= _CDR2[1]:
        return "CDR2"
    if _FR3[0] <= p <= _FR3[1]:
        return "FR3"
    if _CDR3[0] <= p <= _CDR3[1]:
        return "CDR3"
    if _FR4[0] <= p <= _FR4[1]:
        return "FR4"
    return "unknown"


def column_conservation(column_aas: Sequence[str]) -> Optional[float]:
    """Normalized-entropy conservation for one alignment column.

    `column_aas` is the list of amino acids observed at a column across
    the reference sequences (gaps already removed). Returns a score in
    [0, 1], or None when the column is empty.
    """
    n = len(column_aas)
    if n == 0:
        return None
    counts: Dict[str, int] = {}
    for aa in column_aas:
        counts[aa] = counts.get(aa, 0) + 1
    entropy = 0.0
    for c in counts.values():
        p = c / n
        entropy -= p * math.log(p)
    score = 1.0 - entropy / _MAX_ENTROPY
    # Clamp tiny negative/over-one drift from float error.
    return max(0.0, min(1.0, score))


def conservation_class(score: Optional[float]) -> str:
    """Bucket a score into conserved / intermediate / variable / uncovered."""
    if score is None:
        return "uncovered"
    if score >= _CONSERVED_MIN:
        return "conserved"
    if score >= _INTERMEDIATE_MIN:
        return "intermediate"
    return "variable"


def profile_from_matrix(
    sequences: Sequence[str],
) -> Tuple[Dict[int, Optional[float]], Dict[int, int]]:
    """Per-IMGT-position conservation from IMGT-gapped germline sequences.

    All `sequences` are expected to share one fixed length L (the
    IMGT-gapped germline layout), so column index `i` maps to IMGT position
    `i + 1`. Variable-length input is tolerated (a sequence shorter than a
    given column simply does not contribute to it); a column whose non-gap
    fill fraction falls below `_MIN_COLUMN_FILL` is scored None (uncovered).
    Returns (score_by_imgt_pos, n_observed_by_imgt_pos).
    """
    if not sequences:
        return {}, {}
    length = len(sequences[0])
    n_seqs = len(sequences)
    scores: Dict[int, Optional[float]] = {}
    observed: Dict[int, int] = {}
    for i in range(length):
        col = [s[i] for s in sequences if i < len(s) and s[i] not in ("-", ".")]
        imgt_pos = i + 1
        observed[imgt_pos] = len(col)
        if len(col) / n_seqs < _MIN_COLUMN_FILL:
            scores[imgt_pos] = None
        else:
            scores[imgt_pos] = column_conservation(col)
    return scores, observed


# ---------------------------------------------------------------------------
# ANARCI-backed: reference profile + chain numbering
# ---------------------------------------------------------------------------


def _germline_sequences(
    locus: str, species: str, dedupe_alleles: bool
) -> List[str]:
    """Return the IMGT-gapped germline V sequences for a locus/species.

    `dedupe_alleles` collapses `GENE*01/*02/...` to one entry per gene
    (keeping the first) so genes with many alleles don't dominate the
    column statistics.
    """
    try:
        from anarci import all_germlines
    except Exception as exc:  # noqa: BLE001
        raise NoAnarciError(f"ANARCI germline data unavailable: {exc}") from exc

    chain_class = _LOCUS_CLASS.get(locus)
    if chain_class is None:
        raise ValueError(f"Unknown locus {locus!r}; expected 'alpha'/'beta'.")
    try:
        gene_map: Dict[str, str] = all_germlines["V"][chain_class][species]
    except KeyError as exc:
        raise NoAnarciError(
            f"No germline V set for class={chain_class} species={species}"
        ) from exc

    if not dedupe_alleles:
        return list(gene_map.values())

    seen: Dict[str, str] = {}
    for gene_allele, seq in gene_map.items():
        base = gene_allele.split("*", 1)[0]
        if base not in seen:
            seen[base] = seq
    return list(seen.values())


def build_reference_profile(
    locus: str,
    *,
    species: str = "human",
    dedupe_alleles: bool = True,
) -> Tuple[Dict[int, Optional[float]], int]:
    """Build the per-IMGT-position conservation profile for a TCR locus.

    Returns (score_by_imgt_pos, n_reference_sequences).
    """
    seqs = _germline_sequences(locus, species, dedupe_alleles)
    scores, _observed = profile_from_matrix(seqs)
    return scores, len(seqs)


def _number_sequence(
    sequence: str,
) -> Tuple[Optional[str], List[Tuple[int, str, str]]]:
    """ANARCI-number a chain.

    Returns ``(chain_type, residues)`` where ``chain_type`` is ANARCI's
    detected class ('A'/'B'/'H'/...) or None on failure, and ``residues``
    is a list of ``(imgt_pos, insertion_code, aa)`` for real (non-gap)
    residues in order. The non-gap residues form a contiguous block — the
    V domain — within the (possibly longer) input.

    ANARCI's ``number()`` returns ``(False, False)`` when it cannot number
    the sequence (not an Ig/TCR V domain), so the falsy check below must
    catch ``False`` — ``result[0] is None`` would let it through and the
    iteration would raise ``TypeError: 'bool' object is not iterable``.
    """
    try:
        from anarci import number
    except Exception as exc:  # noqa: BLE001
        raise NoAnarciError(f"ANARCI unavailable: {exc}") from exc

    result = number(sequence)
    if not result or not result[0]:
        return None, []
    numbering, chain_type = result[0], result[1]
    out: List[Tuple[int, str, str]] = []
    for (imgt_pos, ins), aa in numbering:
        if aa == "-":
            continue
        out.append((imgt_pos, (ins or "").strip(), aa))
    return chain_type, out


def compute_chain_conservation(
    chain_id: str,
    locus: str,
    sequence: str,
    residue_ids: Sequence[int],
    resnames: Sequence[str],
    *,
    profile: Optional[Dict[int, Optional[float]]] = None,
    species: str = "human",
    dedupe_alleles: bool = True,
) -> List[ConservationResidue]:
    """Map one TCR chain's residues to germline-conservation records.

    `sequence`, `residue_ids`, and `resnames` are parallel sequences from
    the structure (e.g. from `PDBSequenceExtractor`). `profile` may be
    passed to reuse a pre-built reference profile across chains of the
    same locus; otherwise it is built here.

    Returns an empty list (with a logged warning) rather than raising when
    the chain cannot be numbered, or when ANARCI's detected chain type does
    not match the requested `locus` — scoring a mis-identified chain against
    the wrong germline reference would silently corrupt the output.
    """
    if profile is None:
        profile, _n = build_reference_profile(
            locus, species=species, dedupe_alleles=dedupe_alleles
        )

    expected_class = _LOCUS_CLASS.get(locus)
    chain_type, numbered = _number_sequence(sequence)
    if not numbered:
        log.warning(
            "conservation: ANARCI returned no numbering for chain %s "
            "(locus=%s) — skipping",
            chain_id, locus,
        )
        return []
    if expected_class and chain_type and chain_type != expected_class:
        log.warning(
            "conservation: chain %s numbered as ANARCI class %r but locus "
            "%s expects %r — skipping to avoid mis-scoring against the "
            "wrong germline reference",
            chain_id, chain_type, locus, expected_class,
        )
        return []

    # The non-gap numbered residues form a contiguous block (the V domain)
    # inside the possibly longer V+C input sequence. Locate that block by
    # exact substring match so trailing constant-region residues that
    # happen to repeat V-domain amino acids cannot desynchronize the
    # residue↔IMGT mapping (a per-amino-acid search can).
    numbered_aas = "".join(aa for _pos, _ins, aa in numbered)
    offset = sequence.find(numbered_aas)
    if offset < 0:
        log.warning(
            "conservation: could not align ANARCI numbering to the input "
            "sequence for chain %s — skipping",
            chain_id,
        )
        return []

    n = len(sequence)
    if len(residue_ids) != n or len(resnames) != n:
        log.warning(
            "conservation: residue_ids/resnames length (%d/%d) != sequence "
            "length (%d) for chain %s — some records may be dropped",
            len(residue_ids), len(resnames), n, chain_id,
        )

    records: List[ConservationResidue] = []
    for i, (imgt_pos, ins, _aa) in enumerate(numbered):
        seq_idx = offset + i
        if seq_idx >= len(residue_ids) or seq_idx >= len(resnames):
            # Cannot resolve this residue's structure identity — skip
            # rather than emit a (-1, 'UNK') sentinel record.
            continue
        resid = residue_ids[seq_idx]
        resname = resnames[seq_idx]

        if ins:
            # Insertion-coded residue (e.g. 111A): a somatic length variant
            # not represented in the IMGT-gapped germline grid, so it has no
            # germline column — report as uncovered rather than borrowing
            # the base position's conservation score.
            score = None
        else:
            score = profile.get(imgt_pos)
        records.append(
            ConservationResidue(
                chain_id=chain_id,
                resid=resid,
                resname=resname,
                imgt_position=imgt_pos,
                region=imgt_region(imgt_pos),
                conservation_score=score,
                conservation_class=conservation_class(score),
                covered=score is not None,
            )
        )
    return records


__all__ = [
    "ConservationResidue",
    "NoAnarciError",
    "build_reference_profile",
    "column_conservation",
    "compute_chain_conservation",
    "conservation_class",
    "imgt_region",
    "profile_from_matrix",
]
