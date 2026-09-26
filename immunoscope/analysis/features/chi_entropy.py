"""Per-residue sidechain chi1/chi2 dihedral entropy features.

Reads `analysis/geometry/chi_dihedrals/residue_chi_dihedrals.csv` (produced
by `ims chi_dihedrals`) and populates rotamer-mobility fields on each
residue. For mutation design these complement RMSF: RMSF captures backbone
motion magnitude, chi entropy captures whether the sidechain locks into
one rotamer or hops between several — a key signal for choosing between
similar-volume substitutions.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List

import pandas as pd

from .core.locator import CaseLocator
from .core.models import (
    FeatureProvenance,
    FeatureSet,
    ResidueFeatures,
    ResidueKey,
)
from .core.registry import FeatureComputer, register_feature


_CSV_NAME = "residue_chi_dihedrals.csv"

_FLOAT_FIELDS = ("chi1_entropy", "chi2_entropy")


@register_feature
class ChiDihedralEntropyComputer(FeatureComputer):
    """Per-residue chi1/chi2 entropy aggregated from `ims chi_dihedrals`."""

    name = "chi_entropy"
    required_modules: List[str] = ["chi_dihedrals"]
    produces = [*_FLOAT_FIELDS, "rotamer_diversity"]

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        root = locator.get_module_root("chi_dihedrals")
        provenance = FeatureProvenance(computer=self.name)

        if root is None:
            provenance.missing_inputs.append("chi_dihedrals")
            self._record_missing(features, "chi_dihedrals module")
            return features

        csv_path = root / _CSV_NAME
        if not csv_path.exists():
            provenance.missing_inputs.append(f"chi_dihedrals/{_CSV_NAME}")
            self._record_missing(features, f"{_CSV_NAME} (run `ims chi_dihedrals`)")
            return features

        try:
            df = pd.read_csv(csv_path)
        except Exception as exc:
            provenance.missing_inputs.append(f"chi_dihedrals/{_CSV_NAME}: {exc}")
            self._record_missing(features, f"could not read {_CSV_NAME}: {exc}")
            return features

        try:
            provenance.sources.append(str(csv_path.relative_to(locator.case_dir)))
        except ValueError:
            provenance.sources.append(str(csv_path))

        allowed_fields = set(ResidueFeatures.__dataclass_fields__.keys())

        for row in df.itertuples(index=False):
            row_d = row._asdict()
            chain = row_d.get("chain_id")
            resid = row_d.get("resid")
            resname = row_d.get("resname") or ""
            if pd.isna(chain) or pd.isna(resid):
                continue

            key = ResidueKey(chain=str(chain), resid=int(resid), resname=str(resname))
            payload: Dict[str, Any] = {"residue": key, "provenance": provenance}

            for field in _FLOAT_FIELDS:
                if field not in allowed_fields:
                    continue
                value = row_d.get(field)
                if value is None or (isinstance(value, float) and math.isnan(value)):
                    continue
                payload[field] = float(value)

            rotamer_div = row_d.get("rotamer_diversity")
            if (
                "rotamer_diversity" in allowed_fields
                and rotamer_div is not None
                and not (isinstance(rotamer_div, float) and math.isnan(rotamer_div))
            ):
                payload["rotamer_diversity"] = int(rotamer_div)

            if len(payload) <= 2:
                continue

            features.upsert(ResidueFeatures(**payload))

        return features

    @staticmethod
    def _record_missing(features: FeatureSet, note: str) -> None:
        existing = features.metadata.setdefault("missing_inputs", [])
        existing.append(f"chi_entropy:{note}")
