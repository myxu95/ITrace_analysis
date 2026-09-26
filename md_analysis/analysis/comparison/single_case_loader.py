"""Single-system result loading and normalization."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .comparison_schema import SingleCaseArtifacts


class SingleCaseLoader:
    """Load the minimal result surface required by compare from a single-system result root."""

    def load(self, case_root: Path, label: str) -> SingleCaseArtifacts:
        case_root = Path(case_root).resolve()
        overview_root = self._resolve_overview_root(case_root)
        job_root = self._resolve_job_root(case_root, overview_root)

        identity_path = self._first_existing(
            overview_root / "identity" / "biological_identity.json",
            job_root / "analysis" / "identity" / "analysis" / "identity" / "biological_identity.json",
            job_root / "analysis" / "identity" / "biological_identity.json",
        )
        quality_path = self._first_existing(
            overview_root / "quality" / "preprocess_quality_report.json",
            job_root / "preparation" / "quality" / "preprocess_quality_report.json",
            job_root / "quality" / "preprocess_quality_report.json",
        )
        bsa_path = self._first_existing(
            overview_root / "bsa" / "interface_summary.json",
            job_root / "analysis" / "bsa" / "analysis" / "interface" / "interface_summary.json",
            job_root / "analysis" / "interface" / "interface_summary.json",
        )
        rmsf_summary_path = self._first_existing(
            overview_root / "rmsf" / "rmsf_summary.json",
            job_root / "analysis" / "rmsf" / "analysis" / "rmsf" / "rmsf_summary.json",
            job_root / "analysis" / "rmsf" / "rmsf_summary.json",
        )
        rmsf_region_path = self._first_existing(
            overview_root / "rmsf" / "region_rmsf_summary.csv",
            job_root / "analysis" / "rmsf" / "analysis" / "rmsf" / "region_rmsf_summary.csv",
            job_root / "analysis" / "rmsf" / "region_rmsf_summary.csv",
        )
        rrcs_summary_path = self._first_existing(
            overview_root / "rrcs" / "rrcs_summary.json",
            job_root / "analysis" / "rrcs" / "analysis" / "interactions" / "rrcs" / "rrcs_summary.json",
            job_root / "analysis" / "interactions" / "rrcs" / "rrcs_summary.json",
        )
        rrcs_region_path = self._first_existing(
            overview_root / "rrcs" / "rrcs_region_summary.csv",
            job_root / "analysis" / "rrcs" / "analysis" / "interactions" / "rrcs" / "rrcs_region_summary.csv",
            job_root / "analysis" / "interactions" / "rrcs" / "rrcs_region_summary.csv",
        )
        interaction_overview_path = self._first_existing(
            overview_root / "interaction_overview.csv",
        )
        run_summary_path = self._first_existing(job_root / "run_summary.json")
        run_manifest_path = self._first_existing(job_root / "run_manifest.json")

        return SingleCaseArtifacts(
            label=label,
            case_root=case_root,
            overview_root=overview_root,
            identity=self._load_json(identity_path),
            quality=self._load_json(quality_path),
            bsa=self._load_json(bsa_path),
            rmsf_summary=self._load_json(rmsf_summary_path),
            rmsf_regions=self._load_csv(rmsf_region_path),
            rrcs_summary=self._load_json(rrcs_summary_path),
            rrcs_regions=self._load_csv(rrcs_region_path),
            interaction_overview=self._load_csv(interaction_overview_path),
            run_summary=self._load_json(run_summary_path),
            run_manifest=self._load_json(run_manifest_path),
            source_paths={
                "identity": str(identity_path) if identity_path else "",
                "quality": str(quality_path) if quality_path else "",
                "bsa": str(bsa_path) if bsa_path else "",
                "rmsf_summary": str(rmsf_summary_path) if rmsf_summary_path else "",
                "rmsf_regions": str(rmsf_region_path) if rmsf_region_path else "",
                "rrcs_summary": str(rrcs_summary_path) if rrcs_summary_path else "",
                "rrcs_regions": str(rrcs_region_path) if rrcs_region_path else "",
                "interaction_overview": str(interaction_overview_path) if interaction_overview_path else "",
                "run_summary": str(run_summary_path) if run_summary_path else "",
                "run_manifest": str(run_manifest_path) if run_manifest_path else "",
            },
        )

    def _resolve_overview_root(self, case_root: Path) -> Path:
        overview_root = case_root / "overview"
        if overview_root.exists():
            return overview_root
        nested_overviews = sorted((case_root / "report").glob("*/overview"))
        for candidate in nested_overviews:
            if (candidate / "interaction_report_demo.html").exists() or (candidate / "interaction_overview.csv").exists():
                return candidate
        return case_root

    def _resolve_job_root(self, case_root: Path, overview_root: Path) -> Path:
        if (case_root / "run_summary.json").exists() or (case_root / "run_manifest.json").exists():
            return case_root
        if overview_root.name == "overview" and overview_root.parent.parent.name == "report":
            return overview_root.parent.parent.parent
        return case_root

    def _first_existing(self, *candidates: Path) -> Path | None:
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return None

    def _load_json(self, path: Path | None) -> dict:
        if not path:
            return {}
        return json.loads(path.read_text(encoding="utf-8"))

    def _load_csv(self, path: Path | None) -> pd.DataFrame:
        if not path:
            return pd.DataFrame()
        return pd.read_csv(path)

    def load_rmsd_timeseries(self, case_root: Path, selection: str = "backbone") -> Optional[np.ndarray]:
        """Load RMSD time series data.

        Args:
            case_root: Root directory of the case
            selection: Atom selection (backbone, ca, all)

        Returns:
            Array of RMSD values over time, or None if not found
        """
        case_root = Path(case_root).resolve()
        overview_root = self._resolve_overview_root(case_root)
        job_root = self._resolve_job_root(case_root, overview_root)

        rmsd_path = self._first_existing(
            overview_root / "rmsd" / f"rmsd_{selection}.csv",
            job_root / "analysis" / "rmsd" / f"rmsd_{selection}.csv",
        )

        if not rmsd_path:
            return None

        df = pd.read_csv(rmsd_path)
        if "rmsd" in df.columns:
            return df["rmsd"].values
        elif "RMSD" in df.columns:
            return df["RMSD"].values
        return None

    def load_rmsf_timeseries(self, case_root: Path) -> Optional[np.ndarray]:
        """Load per-residue RMSF values.

        Args:
            case_root: Root directory of the case

        Returns:
            Array of RMSF values per residue, or None if not found
        """
        case_root = Path(case_root).resolve()
        overview_root = self._resolve_overview_root(case_root)
        job_root = self._resolve_job_root(case_root, overview_root)

        rmsf_path = self._first_existing(
            overview_root / "rmsf" / "rmsf_per_residue.csv",
            job_root / "analysis" / "rmsf" / "rmsf_per_residue.csv",
        )

        if not rmsf_path:
            return None

        df = pd.read_csv(rmsf_path)
        if "rmsf" in df.columns:
            return df["rmsf"].values
        elif "RMSF" in df.columns:
            return df["RMSF"].values
        return None

    def load_distance_timeseries(self, case_root: Path, distance_name: str = "com_distance") -> Optional[np.ndarray]:
        """Load distance time series data.

        Args:
            case_root: Root directory of the case
            distance_name: Name of the distance metric

        Returns:
            Array of distance values over time, or None if not found
        """
        case_root = Path(case_root).resolve()
        overview_root = self._resolve_overview_root(case_root)
        job_root = self._resolve_job_root(case_root, overview_root)

        distance_path = self._first_existing(
            overview_root / "geometry" / f"{distance_name}.csv",
            job_root / "analysis" / "geometry" / f"{distance_name}.csv",
        )

        if not distance_path:
            return None

        df = pd.read_csv(distance_path)
        if "distance" in df.columns:
            return df["distance"].values
        elif "Distance" in df.columns:
            return df["Distance"].values
        return None

    def load_angle_timeseries(self, case_root: Path, angle_name: str) -> Optional[np.ndarray]:
        """Load angle time series data.

        Args:
            case_root: Root directory of the case
            angle_name: Name of the angle metric

        Returns:
            Array of angle values over time, or None if not found
        """
        case_root = Path(case_root).resolve()
        overview_root = self._resolve_overview_root(case_root)
        job_root = self._resolve_job_root(case_root, overview_root)

        angle_path = self._first_existing(
            overview_root / "angles" / f"{angle_name}.csv",
            job_root / "analysis" / "angles" / f"{angle_name}.csv",
        )

        if not angle_path:
            return None

        df = pd.read_csv(angle_path)
        if "angle" in df.columns:
            return df["angle"].values
        elif "Angle" in df.columns:
            return df["Angle"].values
        return None

    def load_contact_timeseries(self, case_root: Path) -> Optional[pd.DataFrame]:
        """Load contact formation time series.

        Args:
            case_root: Root directory of the case

        Returns:
            DataFrame with contact pairs and their occupancy over time, or None if not found
        """
        case_root = Path(case_root).resolve()
        overview_root = self._resolve_overview_root(case_root)
        job_root = self._resolve_job_root(case_root, overview_root)

        contact_path = self._first_existing(
            overview_root / "contacts" / "contact_timeseries.csv",
            job_root / "analysis" / "contacts" / "contact_timeseries.csv",
        )

        if not contact_path:
            return None

        return pd.read_csv(contact_path)
