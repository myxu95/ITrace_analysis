"""Per-residue interaction-family count features.

Aggregates the 5 typed interaction families (hbond / saltbridge /
hydrophobic / pipi / cationpi) into per-residue counts and max occupancies.
Each `residue_pair_<family>.csv` row contributes +1 to BOTH residues it
connects. Missing families are recorded on FeatureProvenance.missing_inputs
so the feature degrades gracefully on cases that did not run every family.
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

import pandas as pd

from immunoscope.analysis.features.core.locator import CaseLocator
from immunoscope.analysis.features.core.models import (
    FeatureProvenance,
    FeatureSet,
    ResidueFeatures,
    ResidueKey,
)
from immunoscope.analysis.features.core.registry import (
    FeatureComputer,
    register_feature,
)


@register_feature
class InteractionFamilyCountsComputer(FeatureComputer):
    """Count partners per interaction family for every residue at the interface."""

    name = "interaction_family_counts"
    required_modules: List[str] = []  # graceful per-family fallback
    produces = [
        "hbond_count",
        "saltbridge_count",
        "hydrophobic_count",
        "pipi_count",
        "cationpi_count",
        "hbond_max_occupancy",
        "saltbridge_max_occupancy",
        "hydrophobic_max_occupancy",
        "pipi_max_occupancy",
        "cationpi_max_occupancy",
        "dominant_interaction_type",
        "interaction_diversity",
    ]

    # (module_key, csv_filename, count_field, occ_field)
    _FAMILIES: List[Tuple[str, str, str, str]] = [
        ("hbond",       "residue_pair_hbonds.csv",               "hbond_count",       "hbond_max_occupancy"),
        ("saltbridge",  "residue_pair_salt_bridges.csv",         "saltbridge_count",  "saltbridge_max_occupancy"),
        ("hydrophobic", "residue_pair_hydrophobic_contacts.csv", "hydrophobic_count", "hydrophobic_max_occupancy"),
        ("pipi",        "residue_pair_pi_pi.csv",                "pipi_count",        "pipi_max_occupancy"),
        ("cationpi",    "residue_pair_cation_pi.csv",            "cationpi_count",    "cationpi_max_occupancy"),
    ]

    _COUNT_FIELDS = (
        "hbond_count",
        "saltbridge_count",
        "hydrophobic_count",
        "pipi_count",
        "cationpi_count",
    )

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        per_residue: Dict[ResidueKey, Dict[str, Any]] = {}
        sources: List[str] = []
        missing: List[str] = []

        for module_key, csv_name, count_field, occ_field in self._FAMILIES:
            root = locator.get_module_root(module_key)
            if root is None:
                missing.append(module_key)
                continue
            csv_path = root / csv_name
            if not csv_path.exists():
                missing.append(f"{module_key}/{csv_name}")
                continue
            try:
                df = pd.read_csv(csv_path)
            except Exception as exc:
                missing.append(f"{module_key}/{csv_name}: {exc}")
                continue

            try:
                sources.append(str(csv_path.relative_to(locator.case_dir)))
            except ValueError:
                sources.append(str(csv_path))

            for row in df.itertuples(index=False):
                row_d = row._asdict()
                occ_raw = row_d.get("contact_frequency", 0.0)
                occ = float(occ_raw) if pd.notna(occ_raw) else 0.0
                for side in (1, 2):
                    chain = row_d.get(f"chain_id_{side}")
                    resid = row_d.get(f"resid_{side}")
                    resname = row_d.get(f"resname_{side}") or ""
                    if pd.isna(chain) or pd.isna(resid):
                        continue
                    key = ResidueKey(
                        chain=str(chain),
                        resid=int(resid),
                        resname=str(resname),
                    )
                    bucket = per_residue.setdefault(key, {})
                    bucket[count_field] = int(bucket.get(count_field, 0)) + 1
                    bucket[occ_field] = max(float(bucket.get(occ_field, 0.0)), occ)

        # Compute dominant_interaction_type + interaction_diversity per residue
        for bucket in per_residue.values():
            counts = {
                f.replace("_count", ""): bucket[f]
                for f in self._COUNT_FIELDS
                if f in bucket and bucket[f] > 0
            }
            bucket["interaction_diversity"] = len(counts)
            if not counts:
                bucket["dominant_interaction_type"] = "none"
                continue

            max_count = max(counts.values())
            tied = [k for k, c in counts.items() if c == max_count]
            if len(tied) == 1:
                bucket["dominant_interaction_type"] = tied[0]
            else:
                # Break tie on max occupancy
                occ_of = {k: bucket.get(f"{k}_max_occupancy", 0.0) for k in tied}
                top_occ = max(occ_of.values())
                still_tied = [k for k, v in occ_of.items() if v == top_occ]
                bucket["dominant_interaction_type"] = (
                    still_tied[0] if len(still_tied) == 1 else "mixed"
                )

        # Upsert into FeatureSet; merge with whatever was previously computed.
        provenance = FeatureProvenance(
            sources=sources,
            computer=self.name,
            missing_inputs=missing.copy(),
        )
        allowed_fields = set(ResidueFeatures.__dataclass_fields__.keys())

        for key, bucket in per_residue.items():
            payload: Dict[str, Any] = {"residue": key, "provenance": provenance}
            for k, v in bucket.items():
                if k in allowed_fields:
                    payload[k] = v
            features.upsert(ResidueFeatures(**payload))

        if missing:
            existing = features.metadata.setdefault("missing_inputs", [])
            existing.extend(f"interaction_family_counts:{m}" for m in missing)

        return features
