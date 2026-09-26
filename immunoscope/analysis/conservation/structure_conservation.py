"""Bridge structure → germline-conservation table (CSV).

Pulls the TCR alpha/beta chain sequences out of a structure PDB via
`PDBSequenceExtractor`, scores each residue against the IMGT germline
reference (`tcr_germline_conservation`), and writes the per-residue
`conservation_tcr.csv`. Kept separate from the pure scoring module so
that structure I/O (MDAnalysis) is not a hard import for the maths.

`build_conservation_table` / `write_conservation_csv` are runnable
standalone — given a PDB and the TCR chain IDs they produce the CSV
without needing the full MD preprocessing pipeline, which is handy for
back-filling cases or for tests.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from .tcr_germline_conservation import (
    ConservationResidue,
    build_reference_profile,
    compute_chain_conservation,
)

CSV_COLUMNS = [
    "chain_id",
    "resid",
    "resname",
    "locus",
    "imgt_position",
    "region",
    "conservation_score",
    "conservation_class",
    "covered",
    "reference_set",
    "n_ref_seqs",
]

_LOCUS_PREFIX = {"alpha": "TRAV", "beta": "TRBV"}


def build_conservation_table(
    structure_pdb: str,
    tcr_alpha_chain: Optional[str],
    tcr_beta_chain: Optional[str],
    *,
    species: str = "human",
    dedupe_alleles: bool = True,
) -> List[dict]:
    """Return per-residue conservation rows for the TCR chains of a PDB.

    Missing / unresolved chains are skipped. Raises NoAnarciError (from
    the scoring module) if ANARCI is unavailable — callers decide whether
    to degrade.
    """
    from immunoscope.analysis.structure import PDBSequenceExtractor

    extractor = PDBSequenceExtractor()
    chains = extractor.extract_sequences_from_pdb(structure_pdb)

    rows: List[dict] = []
    for locus, chain_id in (("alpha", tcr_alpha_chain), ("beta", tcr_beta_chain)):
        if not chain_id or chain_id not in chains:
            continue
        profile, n_ref = build_reference_profile(
            locus, species=species, dedupe_alleles=dedupe_alleles
        )
        chain = chains[chain_id]
        records = compute_chain_conservation(
            chain_id,
            locus,
            chain["sequence"],
            chain.get("residue_ids", []),
            chain.get("residue_names", []),
            profile=profile,
        )
        ref_label = f"IMGT {_LOCUS_PREFIX[locus]} germline ({species})"
        for r in records:
            rows.append(_record_to_row(r, locus, ref_label, n_ref))
    return rows


def _record_to_row(
    r: ConservationResidue, locus: str, ref_label: str, n_ref: int
) -> dict:
    return {
        "chain_id": r.chain_id,
        "resid": r.resid,
        "resname": r.resname,
        "locus": locus,
        "imgt_position": r.imgt_position if r.imgt_position is not None else "",
        "region": r.region,
        "conservation_score": (
            round(r.conservation_score, 4)
            if r.conservation_score is not None
            else ""
        ),
        "conservation_class": r.conservation_class,
        "covered": r.covered,
        "reference_set": ref_label,
        "n_ref_seqs": n_ref,
    }


def write_conservation_csv(
    structure_pdb: str,
    tcr_alpha_chain: Optional[str],
    tcr_beta_chain: Optional[str],
    out_csv: str,
    *,
    species: str = "human",
    dedupe_alleles: bool = True,
) -> tuple[str, List[dict]]:
    """Compute and write `conservation_tcr.csv`; return (path, rows).

    Always writes a CSV with the canonical header, even when there are
    zero rows (so downstream readers find a well-formed file).
    """
    import pandas as pd

    rows = build_conservation_table(
        structure_pdb,
        tcr_alpha_chain,
        tcr_beta_chain,
        species=species,
        dedupe_alleles=dedupe_alleles,
    )
    out_path = Path(out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=CSV_COLUMNS)
    df.to_csv(out_path, index=False)
    return str(out_path), rows


def write_empty_conservation_csv(out_csv: str) -> str:
    """Write a header-only CSV (used when ANARCI/inputs are unavailable)."""
    import pandas as pd

    out_path = Path(out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([], columns=CSV_COLUMNS).to_csv(out_path, index=False)
    return str(out_path)


__all__ = [
    "CSV_COLUMNS",
    "build_conservation_table",
    "write_conservation_csv",
    "write_empty_conservation_csv",
]
