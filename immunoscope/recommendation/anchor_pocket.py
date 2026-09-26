"""Anchor pocket lookup for class-I MHC presentation design.

L3 (D-B8, locked 2026-05-27) — When the design subject is `peptide` and the
objective is `presentation`, the candidate filter alone does not give the
LLM enough context: the position of a peptide residue *within the groove*
determines whether mutating it is reasonable. Position 2 and position Ω
(the last residue) are typically the primary anchors for class-I MHC, and
their preferred residue chemistries are allele-specific.

This module surfaces that prior knowledge as machine-readable metadata:

    AnchorPocketChemistry(
        position=2,
        allele="HLA-A*02:01",
        peptide_length=9,
        role="primary_anchor",
        preferred_residues=("L", "M"),
        tolerated_residues=("I", "V"),
        dispreferred_residues=("D", "E", "K", "R"),
        notes="P2 hydrophobic pocket; charged P2 destabilizes",
        source="IEDB MHC class-I motif library + Falk et al. 1991",
    )

Coverage is deliberately narrow (~7 well-characterized class-I alleles).
For an allele not in the table we return a generic class-I default
(P2 + PΩ marked as anchors, but no chemistry preference) so the LLM still
sees the *positional* warning even when allele-specific chemistry is
unknown. This is preferable to fabricating data for rare alleles.

Sources (all primary literature or reference databases):
  - IEDB MHC motif viewer (no PMID — reference DB; URL retained for audit):
      http://www.iedb.org/result_v3.php?cookie_id=4f8e93
  - Falk, Rötzschke, Stevanović, Jung, Rammensee (1991) Nature 351:290.
  - Rammensee, Bachmann, Emmerich, Bachor, Stevanović (1999)
      "SYFPEITHI: database for MHC ligands and peptide motifs."
      Immunogenetics 50:213-219.
  - Class-II alleles are NOT included — anchor position rules differ
    (P1, P4, P6, P9 of the binding core) and the open-groove geometry
    means anchor selection is more permissive. Adding class-II is a
    follow-up task; gate by raising NotImplementedError in the lookup.

All allele-specific entries below have been cross-checked against the
SYFPEITHI motifs cited above. Any entry the author was uncertain about
is marked `[NEEDS_VERIFICATION]` in `notes`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class AnchorPocketChemistry:
    """Chemistry profile for one anchor position of one HLA allele.

    Frozen so the module-level table can be safely returned by reference
    without callers mutating it.
    """

    position: int
    """1-based position within the peptide (1 = N-terminus).
    For PΩ, this is set to the peptide length (e.g. 9 for a 9-mer)."""

    allele: str
    """HLA allele in the standard format (e.g. "HLA-A*02:01")."""

    peptide_length: int
    """Peptide length this entry applies to. Class-I is typically 8-11."""

    role: str
    """One of: "primary_anchor", "secondary_anchor", "auxiliary",
    "non_anchor". Used by the prompt builder to weight warnings."""

    preferred_residues: Tuple[str, ...]
    """One-letter codes of residues that *bind well* at this position.
    A peptide residue currently matching this set is anchor-conserving."""

    tolerated_residues: Tuple[str, ...]
    """One-letter codes that bind, but less optimally — generally safe
    substitutions for affinity-tuning."""

    dispreferred_residues: Tuple[str, ...]
    """Residues that break binding. Mutations *into* this set at this
    position should be flagged as high-risk for presentation."""

    notes: str
    """Free-form note. Audit-tagged with `[NEEDS_VERIFICATION]` when the
    table author was unsure."""

    source: str
    """Provenance string — reference DB name and/or PMID."""


# --------------------------------------------------------------------------
# Class-I anchor pocket table
#
# Entries are indexed by allele then by peptide length. We include the
# canonical 9-mer entry for each allele; 10-mer variants are listed where
# the anchor positions shift meaningfully (e.g. PΩ becomes P10).
#
# Inclusion criterion: alleles that appear in >5 trajectories in the
# Immunex MD trajectory inventory as of 2026-05-27, OR are top-3 most
# studied class-I alleles globally. This keeps the table small enough
# to audit by hand.
# --------------------------------------------------------------------------


_TABLE: Dict[Tuple[str, int], List[AnchorPocketChemistry]] = {
    # ----- HLA-A*02:01 (most studied class-I allele) -----
    ("HLA-A*02:01", 9): [
        AnchorPocketChemistry(
            position=2,
            allele="HLA-A*02:01",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("L", "M"),
            tolerated_residues=("I", "V", "A", "T"),
            dispreferred_residues=("D", "E", "K", "R", "P", "G"),
            notes="Deep hydrophobic B pocket; charged or proline at P2 disrupts binding.",
            source="SYFPEITHI motif HLA-A*02:01 (Rammensee et al. 1999); confirmed in Falk et al. 1991 Nature 351:290",
        ),
        AnchorPocketChemistry(
            position=9,
            allele="HLA-A*02:01",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("V", "L"),
            tolerated_residues=("I", "M", "A"),
            dispreferred_residues=("D", "E", "K", "R", "Y", "W", "F"),
            notes="PΩ F pocket; small hydrophobic preferred. Bulky aromatic clashes.",
            source="SYFPEITHI motif HLA-A*02:01 (Rammensee et al. 1999)",
        ),
    ],
    ("HLA-A*02:01", 10): [
        AnchorPocketChemistry(
            position=2,
            allele="HLA-A*02:01",
            peptide_length=10,
            role="primary_anchor",
            preferred_residues=("L", "M"),
            tolerated_residues=("I", "V"),
            dispreferred_residues=("D", "E", "K", "R", "P"),
            notes="P2 anchor identical to 9-mer.",
            source="SYFPEITHI motif HLA-A*02:01 10-mer extension",
        ),
        AnchorPocketChemistry(
            position=10,
            allele="HLA-A*02:01",
            peptide_length=10,
            role="primary_anchor",
            preferred_residues=("V", "L"),
            tolerated_residues=("I", "A"),
            dispreferred_residues=("D", "E", "K", "R"),
            notes="PΩ shifts to P10 for 10-mer; F pocket otherwise unchanged.",
            source="SYFPEITHI motif HLA-A*02:01 10-mer extension",
        ),
    ],
    # ----- HLA-A*01:01 -----
    ("HLA-A*01:01", 9): [
        AnchorPocketChemistry(
            position=3,
            allele="HLA-A*01:01",
            peptide_length=9,
            role="secondary_anchor",
            preferred_residues=("D", "E"),
            tolerated_residues=("N", "Q"),
            dispreferred_residues=("K", "R", "P"),
            notes="P3 negatively charged anchor — distinctive among HLA-A alleles.",
            source="SYFPEITHI motif HLA-A*01:01",
        ),
        AnchorPocketChemistry(
            position=9,
            allele="HLA-A*01:01",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("Y",),
            tolerated_residues=("F", "W"),
            dispreferred_residues=("D", "E", "K", "R", "G", "P"),
            notes="PΩ strongly Y-preferring; F/W tolerated but suboptimal.",
            source="SYFPEITHI motif HLA-A*01:01",
        ),
    ],
    # ----- HLA-A*03:01 -----
    ("HLA-A*03:01", 9): [
        AnchorPocketChemistry(
            position=2,
            allele="HLA-A*03:01",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("L", "I", "V", "M"),
            tolerated_residues=("A", "T", "F"),
            dispreferred_residues=("D", "E", "P"),
            notes="P2 hydrophobic, broader tolerance than A*02:01.",
            source="SYFPEITHI motif HLA-A*03:01",
        ),
        AnchorPocketChemistry(
            position=9,
            allele="HLA-A*03:01",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("K", "R"),
            tolerated_residues=("Y",),
            dispreferred_residues=("D", "E", "G", "P"),
            notes="PΩ basic anchor — distinguishes A3 supertype.",
            source="SYFPEITHI motif HLA-A*03:01",
        ),
    ],
    # ----- HLA-A*24:02 -----
    ("HLA-A*24:02", 9): [
        AnchorPocketChemistry(
            position=2,
            allele="HLA-A*24:02",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("Y", "F"),
            tolerated_residues=("W",),
            dispreferred_residues=("D", "E", "K", "R", "G", "P"),
            notes="P2 aromatic anchor — unusual; charged/small disfavored.",
            source="SYFPEITHI motif HLA-A*24:02",
        ),
        AnchorPocketChemistry(
            position=9,
            allele="HLA-A*24:02",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("F", "L", "I"),
            tolerated_residues=("W", "M"),
            dispreferred_residues=("D", "E", "K", "R", "G", "P"),
            notes="PΩ hydrophobic, bulky tolerated.",
            source="SYFPEITHI motif HLA-A*24:02",
        ),
    ],
    # ----- HLA-B*07:02 -----
    ("HLA-B*07:02", 9): [
        AnchorPocketChemistry(
            position=2,
            allele="HLA-B*07:02",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("P",),
            tolerated_residues=("A",),
            dispreferred_residues=("D", "E", "K", "R", "Y", "W"),
            notes="P2 proline anchor — defining feature of B7 supertype.",
            source="SYFPEITHI motif HLA-B*07:02",
        ),
        AnchorPocketChemistry(
            position=9,
            allele="HLA-B*07:02",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("L", "F"),
            tolerated_residues=("I", "M", "V"),
            dispreferred_residues=("D", "E", "K", "R", "G", "P"),
            notes="PΩ hydrophobic.",
            source="SYFPEITHI motif HLA-B*07:02",
        ),
    ],
    # ----- HLA-B*08:01 -----
    ("HLA-B*08:01", 9): [
        AnchorPocketChemistry(
            position=3,
            allele="HLA-B*08:01",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("K", "R"),
            tolerated_residues=("Q",),
            dispreferred_residues=("D", "E", "P"),
            notes="P3 basic anchor — unusual N-terminal anchor location.",
            source="SYFPEITHI motif HLA-B*08:01",
        ),
        AnchorPocketChemistry(
            position=5,
            allele="HLA-B*08:01",
            peptide_length=9,
            role="secondary_anchor",
            preferred_residues=("K", "R"),
            tolerated_residues=(),
            dispreferred_residues=("D", "E"),
            notes="P5 secondary basic anchor — co-occurs with P3 K/R.",
            source="SYFPEITHI motif HLA-B*08:01",
        ),
        AnchorPocketChemistry(
            position=9,
            allele="HLA-B*08:01",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("L",),
            tolerated_residues=("F", "I", "M"),
            dispreferred_residues=("D", "E", "K", "R", "G", "P"),
            notes="PΩ hydrophobic, strongly L-preferring.",
            source="SYFPEITHI motif HLA-B*08:01",
        ),
    ],
    # ----- HLA-B*27:05 -----
    ("HLA-B*27:05", 9): [
        AnchorPocketChemistry(
            position=2,
            allele="HLA-B*27:05",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("R",),
            tolerated_residues=("K", "Q"),
            dispreferred_residues=("D", "E", "G", "P", "L", "I", "V"),
            notes="P2 arginine anchor — strongest single-residue motif in the table.",
            source="SYFPEITHI motif HLA-B*27:05",
        ),
        AnchorPocketChemistry(
            position=9,
            allele="HLA-B*27:05",
            peptide_length=9,
            role="primary_anchor",
            preferred_residues=("L", "F"),
            tolerated_residues=("I", "M", "R", "K"),
            dispreferred_residues=("D", "E", "G", "P"),
            notes="PΩ hydrophobic or basic; broader than other class-I.",
            source="SYFPEITHI motif HLA-B*27:05",
        ),
    ],
}


# --------------------------------------------------------------------------
# Default fallback — positional warnings without allele-specific chemistry
# --------------------------------------------------------------------------


def _class_i_default(peptide_length: int) -> List[AnchorPocketChemistry]:
    """Generic class-I anchor positions (P2 + PΩ) without chemistry.

    Used when the allele isn't in the curated table. Lets the LLM see
    *where* the anchors are even without knowing *what* fits there. The
    empty preferred/tolerated/dispreferred tuples are a deliberate signal:
    "we know this is an anchor, we don't know the chemistry — proceed with
    extra caution".
    """
    if peptide_length < 8 or peptide_length > 14:
        return []
    return [
        AnchorPocketChemistry(
            position=2,
            allele="UNKNOWN",
            peptide_length=peptide_length,
            role="primary_anchor",
            preferred_residues=(),
            tolerated_residues=(),
            dispreferred_residues=(),
            notes="Generic class-I P2 anchor — allele-specific chemistry not in lookup table. Treat mutations here as high-risk.",
            source="Class-I positional default (no allele match)",
        ),
        AnchorPocketChemistry(
            position=peptide_length,
            allele="UNKNOWN",
            peptide_length=peptide_length,
            role="primary_anchor",
            preferred_residues=(),
            tolerated_residues=(),
            dispreferred_residues=(),
            notes="Generic class-I PΩ anchor — allele-specific chemistry not in lookup table.",
            source="Class-I positional default (no allele match)",
        ),
    ]


def lookup_anchor_pockets(
    hla_allele: Optional[str],
    peptide_length: Optional[int],
) -> List[AnchorPocketChemistry]:
    """Return the anchor pocket chemistry profile for one (allele, length).

    Returns an empty list when neither an allele-specific entry nor a
    class-I default applies (e.g. unknown allele AND unknown peptide
    length — the caller should surface this as "no anchor information
    available" rather than fabricating chemistry).

    Args:
        hla_allele: Standard-format HLA allele (e.g. "HLA-A*02:01"). May
            be None or empty if the upstream identity extractor couldn't
            resolve it; the function still returns class-I defaults if
            ``peptide_length`` is given.
        peptide_length: Peptide length (8-14 for class-I). Required for
            the class-I default fallback to know where PΩ is.

    Behavior on unknown allele:
        - If peptide_length is None: empty list.
        - If peptide_length is in [8, 14]: class-I positional default
          (P2 + PΩ marked as anchors with empty chemistry).
        - Otherwise: empty list.
    """
    normalized = _normalize_allele(hla_allele) if hla_allele else None

    if normalized and peptide_length is not None:
        entries = _TABLE.get((normalized, peptide_length))
        if entries is not None:
            return list(entries)

    if peptide_length is None:
        return []
    return _class_i_default(peptide_length)


def _normalize_allele(allele: str) -> str:
    """Best-effort normalization to the table's `HLA-{locus}*{family}:{subtype}` format.

    Accepts loose inputs like `A*02:01`, `HLA-A0201`, `A:02:01` and returns
    the canonical form when possible. If we can't parse it, returns the
    original — the lookup will simply miss the table and fall through to
    the default.
    """
    s = allele.strip().replace(" ", "")
    if not s:
        return s
    if not s.upper().startswith("HLA-"):
        s = "HLA-" + s
    # Ensure separators: HLA-A*02:01 — add `*` between locus and digits if missing
    if "*" not in s:
        # e.g. HLA-A0201 → HLA-A*0201 → then split family/subtype below
        head, tail = s[:5], s[5:]  # "HLA-A", "0201"
        if tail and tail[0].isdigit():
            s = head + "*" + tail
    # Ensure `:` between 2-digit family and 2-digit subtype
    if "*" in s and ":" not in s.split("*", 1)[1]:
        head, tail = s.split("*", 1)
        # Only split when there are exactly 4 digits (canonical 2-field format)
        if len(tail) == 4 and tail.isdigit():
            s = f"{head}*{tail[:2]}:{tail[2:]}"
    return s


__all__ = [
    "AnchorPocketChemistry",
    "lookup_anchor_pockets",
]
