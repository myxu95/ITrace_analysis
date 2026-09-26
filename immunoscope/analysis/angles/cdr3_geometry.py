"""CDR3 loop geometry metrics for pHLA-TCR interfaces."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class CDR3GeometryResult:
    """Geometry summary for one CDR3 loop."""

    chain_id: str
    tip_residue_label: str
    loop_span_angstrom: float
    tip_to_peptide_distance_angstrom: float | None
    torsion_deg: float | None
    n_cdr3_residues: int

    def to_dict(self) -> dict:
        return {
            "chain_id": self.chain_id,
            "tip_residue_label": self.tip_residue_label,
            "loop_span_angstrom": self.loop_span_angstrom,
            "tip_to_peptide_distance_angstrom": self.tip_to_peptide_distance_angstrom,
            "torsion_deg": self.torsion_deg,
            "n_cdr3_residues": self.n_cdr3_residues,
        }


class CDR3GeometryAnalyzer:
    """
    Calculate CDR3 loop geometry from residue-level coordinates.

    Expected columns:
      - chain_id, resid, resname, region, x, y, z
      - component is optional but helps identify peptide residues
    """

    required_columns = {"chain_id", "resid", "resname", "region", "x", "y", "z"}

    def calculate_from_residue_coordinates(self, residue_coordinates: pd.DataFrame) -> pd.DataFrame:
        missing = self.required_columns - set(residue_coordinates.columns)
        if missing:
            raise ValueError(f"CDR3 geometry input missing columns: {sorted(missing)}")
        if residue_coordinates.empty:
            return self._empty_result()

        df = residue_coordinates.copy()
        for column in ("x", "y", "z"):
            df[column] = pd.to_numeric(df[column], errors="coerce")
        df = df.dropna(subset=["x", "y", "z"])
        if df.empty:
            return self._empty_result()

        peptide = self._peptide_coordinates(df)
        peptide_centroid = peptide.mean(axis=0) if len(peptide) else None

        results: list[CDR3GeometryResult] = []
        cdr3_df = df[df["region"].astype(str).str.contains("CDR3", case=False, na=False)].copy()
        for chain_id, chain_df in cdr3_df.groupby("chain_id", sort=True):
            chain_df = chain_df.sort_values("resid")
            coords = chain_df[["x", "y", "z"]].to_numpy(dtype=float)
            if len(coords) < 2:
                continue
            tip_index = len(coords) // 2
            tip = coords[tip_index]
            loop_span = float(np.linalg.norm(coords[-1] - coords[0]))
            tip_distance = (
                float(np.linalg.norm(tip - peptide_centroid))
                if peptide_centroid is not None
                else None
            )
            torsion = self._torsion(coords[0], tip, coords[-1], peptide_centroid)
            tip_row = chain_df.iloc[tip_index]
            results.append(
                CDR3GeometryResult(
                    chain_id=str(chain_id),
                    tip_residue_label=f"{tip_row['resname']}-{chain_id}{int(tip_row['resid'])}",
                    loop_span_angstrom=round(loop_span, 4),
                    tip_to_peptide_distance_angstrom=round(tip_distance, 4) if tip_distance is not None else None,
                    torsion_deg=round(torsion, 4) if torsion is not None else None,
                    n_cdr3_residues=int(len(chain_df)),
                )
            )

        if not results:
            return self._empty_result()
        return pd.DataFrame([item.to_dict() for item in results])

    def write_csv(self, residue_coordinates: pd.DataFrame, output_file: str | Path) -> pd.DataFrame:
        result = self.calculate_from_residue_coordinates(residue_coordinates)
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output_path, index=False)
        return result

    @staticmethod
    def _peptide_coordinates(df: pd.DataFrame) -> np.ndarray:
        if "component" in df.columns:
            peptide = df[df["component"].astype(str).str.contains("peptide", case=False, na=False)]
        else:
            peptide = df[df["region"].astype(str).str.contains("peptide", case=False, na=False)]
        return peptide[["x", "y", "z"]].to_numpy(dtype=float)

    @staticmethod
    def _torsion(
        start: np.ndarray,
        tip: np.ndarray,
        end: np.ndarray,
        peptide_centroid: np.ndarray | None,
    ) -> float | None:
        if peptide_centroid is None:
            return None
        v1 = start - tip
        v2 = end - tip
        v3 = peptide_centroid - tip
        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)
        if np.linalg.norm(n1) == 0 or np.linalg.norm(n2) == 0:
            return None
        n1 = n1 / np.linalg.norm(n1)
        n2 = n2 / np.linalg.norm(n2)
        x = float(np.clip(np.dot(n1, n2), -1.0, 1.0))
        return float(np.degrees(np.arccos(x)))

    @staticmethod
    def _empty_result() -> pd.DataFrame:
        return pd.DataFrame(
            columns=[
                "chain_id",
                "tip_residue_label",
                "loop_span_angstrom",
                "tip_to_peptide_distance_angstrom",
                "torsion_deg",
                "n_cdr3_residues",
            ]
        )
