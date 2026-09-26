"""Per-residue secondary structure features.

Reads `analysis/structure/residue_secondary_structure.csv` (produced by
`ims dssp`) and populates SS fields on each residue. Mutation design needs
to know whether a residue sits inside a rigid helix or sheet — where a
substitution can break the secondary element — or in a flexible loop where
exchanges are relatively safe.
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


_CSV_NAME = "residue_secondary_structure.csv"

_FLOAT_FIELDS = (
    "ss_propensity_helix",
    "ss_propensity_sheet",
    "ss_propensity_loop",
    "ss_stability",
)


@register_feature
class SecondaryStructureComputer(FeatureComputer):
    """Per-residue DSSP secondary structure aggregated from `ims dssp`."""

    name = "secondary_structure"
    required_modules: List[str] = ["structure"]
    produces = ["secondary_structure", *_FLOAT_FIELDS]

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        root = locator.get_module_root("structure")
        provenance = FeatureProvenance(computer=self.name)

        if root is None:
            provenance.missing_inputs.append("structure")
            self._record_missing(features, "structure module")
            return features

        csv_path = root / _CSV_NAME
        if not csv_path.exists():
            provenance.missing_inputs.append(f"structure/{_CSV_NAME}")
            self._record_missing(features, f"{_CSV_NAME} (run `ims dssp`)")
            return features

        try:
            df = pd.read_csv(csv_path)
        except Exception as exc:
            provenance.missing_inputs.append(f"structure/{_CSV_NAME}: {exc}")
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

            ss = row_d.get("secondary_structure")
            if (
                "secondary_structure" in allowed_fields
                and ss is not None
                and not (isinstance(ss, float) and math.isnan(ss))
                and str(ss) not in {"", "nan"}
            ):
                payload["secondary_structure"] = str(ss)

            for field in _FLOAT_FIELDS:
                if field not in allowed_fields:
                    continue
                value = row_d.get(field)
                if value is None or (isinstance(value, float) and math.isnan(value)):
                    continue
                payload[field] = float(value)

            if len(payload) <= 2:
                continue

            features.upsert(ResidueFeatures(**payload))

        return features

    @staticmethod
    def _record_missing(features: FeatureSet, note: str) -> None:
        existing = features.metadata.setdefault("missing_inputs", [])
        existing.append(f"secondary_structure:{note}")
