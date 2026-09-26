"""Backbone / sidechain RMSF split per residue.

Reads `analysis/rmsf/residue_rmsf_atom_sets.csv` (produced by
`ims rmsf --annotated --split-atom-sets`) and populates `rmsf_backbone`,
`rmsf_sidechain`, `rmsf_ratio` on each residue.

The ratio is the key design signal: high ratio means the backbone is rigid
while the sidechain explores rotamers — a strong indicator that the position
tolerates substitution. Low ratio means the sidechain is anchored on a moving
backbone — handle with care.
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


_CSV_NAME = "residue_rmsf_atom_sets.csv"


@register_feature
class RmsfBackboneSidechainComputer(FeatureComputer):
    """Per-residue backbone / sidechain RMSF and ratio."""

    name = "rmsf_backbone_sidechain"
    required_modules: List[str] = ["rmsf"]
    produces = ["rmsf_backbone", "rmsf_sidechain", "rmsf_ratio"]

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        root = locator.get_module_root("rmsf")
        provenance = FeatureProvenance(computer=self.name)

        if root is None:
            provenance.missing_inputs.append("rmsf")
            self._record_missing(features, "rmsf module")
            return features

        csv_path = root / _CSV_NAME
        if not csv_path.exists():
            provenance.missing_inputs.append(f"rmsf/{_CSV_NAME}")
            self._record_missing(
                features,
                f"{_CSV_NAME} (run `ims rmsf --annotated --split-atom-sets`)",
            )
            return features

        try:
            df = pd.read_csv(csv_path)
        except Exception as exc:
            provenance.missing_inputs.append(f"rmsf/{_CSV_NAME}: {exc}")
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

            for field in ("rmsf_backbone", "rmsf_sidechain", "rmsf_ratio"):
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
        existing.append(f"rmsf_backbone_sidechain:{note}")
