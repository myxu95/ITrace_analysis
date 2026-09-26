"""
Example feature computers.

These are simple reference implementations that exercise the feature layer
end-to-end and verify that the core + locator + registry wiring is correct.
They are NOT the production Tier 1 features — those will be implemented in
PR 4 with full attention to physical correctness and edge cases.

Kept intentionally small so they serve as templates for future computers.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from .core.locator import CaseLocator
from .core.models import (
    FeatureSet,
    FeatureProvenance,
    ResidueFeatures,
    ResidueKey,
)
from .core.registry import FeatureComputer, register_feature


# ---------------------------------------------------------------------- #
# Helper: load the annotated RRCS pair table with robust column checking
# ---------------------------------------------------------------------- #

_RRCS_PAIR_FILENAME = "annotated_rrcs_pair_summary.csv"

# Columns the pair table is expected to provide. Missing columns are
# reported via provenance.missing_inputs so downstream consumers know the
# feature is partial.
_EXPECTED_COLUMNS = {
    "chain_id_1",
    "resid_1",
    "resname_1",
    "chain_id_2",
    "resid_2",
    "resname_2",
    "mean_rrcs",
    "rrcs_nonzero_fraction",
}


def _load_rrcs_pair_table(locator: CaseLocator) -> tuple[Optional[pd.DataFrame], Optional[Path]]:
    """
    Locate and load the annotated RRCS pair table.

    Returns (df, path). Either may be None if the table isn't available.
    """
    path = locator.find_file("rrcs", _RRCS_PAIR_FILENAME)
    if path is None:
        # Fall back to common report and pipeline layouts.
        path = locator.first_existing(
            Path("analysis/rrcs/analysis/interactions/rrcs") / _RRCS_PAIR_FILENAME,
            Path("analysis/interactions/rrcs") / _RRCS_PAIR_FILENAME,
            Path("overview/rrcs") / _RRCS_PAIR_FILENAME,
            Path("rrcs/analysis/interactions/rrcs") / _RRCS_PAIR_FILENAME,
        )
    if path is None or not path.exists():
        return None, None

    try:
        df = pd.read_csv(path)
    except Exception:
        return None, path

    return df, path


# ---------------------------------------------------------------------- #
# Example 1: contact_count — counts interface partners per residue
# ---------------------------------------------------------------------- #

@register_feature
class ContactCountComputer(FeatureComputer):
    """
    Counts distinct interface partners for each residue.

    This is a deliberately simple feature: for every residue that appears
    on either side of an RRCS pair, it counts how many distinct partner
    residues it has. It populates ResidueFeatures.contact_count.

    A partner is counted only once per residue regardless of how many
    frames it was observed in.
    """

    name = "contact_count"
    required_modules = ["rrcs"]
    produces = ["contact_count"]

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        df, path = _load_rrcs_pair_table(locator)

        if df is None:
            # Record missing input on the feature set itself.
            features.metadata.setdefault("missing_inputs", []).append(
                f"{_RRCS_PAIR_FILENAME} (rrcs module)"
            )
            return features

        missing_cols = _EXPECTED_COLUMNS - set(df.columns)
        source_rel = str(path.relative_to(locator.case_dir)) if path else ""

        # Build per-residue partner sets.
        partners: Dict[ResidueKey, set] = {}

        for _, row in df.iterrows():
            try:
                k1 = ResidueKey(
                    chain=str(row["chain_id_1"]),
                    resid=int(row["resid_1"]),
                    resname=str(row.get("resname_1", "")),
                )
                k2 = ResidueKey(
                    chain=str(row["chain_id_2"]),
                    resid=int(row["resid_2"]),
                    resname=str(row.get("resname_2", "")),
                )
            except (KeyError, ValueError, TypeError):
                continue

            partners.setdefault(k1, set()).add((k2.chain, k2.resid))
            partners.setdefault(k2, set()).add((k1.chain, k1.resid))

        # Write into FeatureSet.
        for key, partner_set in partners.items():
            prov = FeatureProvenance(
                sources=[source_rel] if source_rel else [],
                computer=self.name,
                missing_inputs=[c for c in missing_cols],
            )
            rf = ResidueFeatures(
                residue=key,
                contact_count=len(partner_set),
                provenance=prov,
            )
            features.upsert(rf)

        return features


# ---------------------------------------------------------------------- #
# Example 2: rrcs_contribution — sums mean_rrcs over partners per residue
# ---------------------------------------------------------------------- #

@register_feature
class RrcsContributionComputer(FeatureComputer):
    """
    Computes rrcs_contribution and rrcs_rank for each residue.

    For each residue, sums mean_rrcs across all pairs it participates in.
    Then ranks residues in descending order of contribution.

    Populates:
        - rrcs_contribution: sum of mean_rrcs across all pairs
        - rrcs_rank: 1-indexed rank, 1 = highest contribution
    """

    name = "rrcs_contribution"
    required_modules = ["rrcs"]
    produces = ["rrcs_contribution", "rrcs_rank"]

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        df, path = _load_rrcs_pair_table(locator)

        if df is None:
            features.metadata.setdefault("missing_inputs", []).append(
                f"{_RRCS_PAIR_FILENAME} (rrcs module)"
            )
            return features

        if "mean_rrcs" not in df.columns:
            features.metadata.setdefault("missing_inputs", []).append(
                "mean_rrcs column"
            )
            return features

        source_rel = str(path.relative_to(locator.case_dir)) if path else ""

        # Accumulate per-residue sum.
        contrib: Dict[ResidueKey, float] = {}

        for _, row in df.iterrows():
            try:
                mean_rrcs = float(row["mean_rrcs"])
            except (ValueError, TypeError):
                continue

            try:
                k1 = ResidueKey(
                    chain=str(row["chain_id_1"]),
                    resid=int(row["resid_1"]),
                    resname=str(row.get("resname_1", "")),
                )
                k2 = ResidueKey(
                    chain=str(row["chain_id_2"]),
                    resid=int(row["resid_2"]),
                    resname=str(row.get("resname_2", "")),
                )
            except (KeyError, ValueError, TypeError):
                continue

            contrib[k1] = contrib.get(k1, 0.0) + mean_rrcs
            contrib[k2] = contrib.get(k2, 0.0) + mean_rrcs

        # Rank (descending).
        ranked = sorted(contrib.items(), key=lambda kv: kv[1], reverse=True)

        for rank, (key, score) in enumerate(ranked, start=1):
            prov = FeatureProvenance(
                sources=[source_rel] if source_rel else [],
                computer=self.name,
            )
            rf = ResidueFeatures(
                residue=key,
                rrcs_contribution=round(score, 4),
                rrcs_rank=rank,
                provenance=prov,
            )
            features.upsert(rf)

        return features
