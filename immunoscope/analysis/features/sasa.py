"""Per-residue SASA features.

Reads `analysis/interface/residue_sasa.csv` (produced by `ims residue_sasa`)
and populates exposure fields on each residue. These signal how much room
exists at the interface for a larger / differently shaped sidechain — the
question RRCS and BSA at the chain level cannot answer.
"""

from __future__ import annotations

from typing import Any, Dict, List

import math
import pandas as pd

from .core.locator import CaseLocator
from .core.models import (
    FeatureProvenance,
    FeatureSet,
    ResidueFeatures,
    ResidueKey,
)
from .core.registry import FeatureComputer, register_feature


_CSV_NAME = "residue_sasa.csv"

_FIELDS = (
    "sasa_bound",
    "sasa_unbound",
    "delta_sasa",
    "relative_exposure",
    "sasa_sidechain_bound",
)


@register_feature
class ResidueSasaComputer(FeatureComputer):
    """Per-residue bound / unbound SASA and burial-state classification."""

    name = "residue_sasa"
    required_modules: List[str] = ["interface"]
    produces = [*_FIELDS, "burial_state"]

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        root = locator.get_module_root("interface")
        provenance = FeatureProvenance(computer=self.name)

        if root is None:
            provenance.missing_inputs.append("interface")
            self._record_missing(features, "interface module")
            return features

        csv_path = root / _CSV_NAME
        if not csv_path.exists():
            provenance.missing_inputs.append(f"interface/{_CSV_NAME}")
            self._record_missing(
                features,
                f"{_CSV_NAME} (run `ims residue_sasa`)",
            )
            return features

        try:
            df = pd.read_csv(csv_path)
        except Exception as exc:
            provenance.missing_inputs.append(f"interface/{_CSV_NAME}: {exc}")
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

            for field in _FIELDS:
                if field not in allowed_fields:
                    continue
                value = row_d.get(field)
                if value is None or (isinstance(value, float) and math.isnan(value)):
                    continue
                payload[field] = float(value)

            burial = row_d.get("burial_state")
            if (
                "burial_state" in allowed_fields
                and burial is not None
                and not (isinstance(burial, float) and math.isnan(burial))
                and str(burial) not in {"", "nan"}
            ):
                payload["burial_state"] = str(burial)

            if len(payload) <= 2:
                continue

            features.upsert(ResidueFeatures(**payload))

        return features

    @staticmethod
    def _record_missing(features: FeatureSet, note: str) -> None:
        existing = features.metadata.setdefault("missing_inputs", [])
        existing.append(f"residue_sasa:{note}")
