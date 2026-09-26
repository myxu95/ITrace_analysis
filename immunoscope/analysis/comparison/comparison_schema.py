"""Data structures for the comparison module."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass(slots=True)
class SingleCaseArtifacts:
    """Single-system input surface for comparison."""

    label: str
    case_root: Path
    overview_root: Path
    identity: dict
    quality: dict
    bsa: dict
    rmsf_summary: dict
    rmsf_regions: pd.DataFrame = field(default_factory=pd.DataFrame)
    rrcs_summary: dict = field(default_factory=dict)
    rrcs_regions: pd.DataFrame = field(default_factory=pd.DataFrame)
    interaction_overview: pd.DataFrame = field(default_factory=pd.DataFrame)
    run_summary: dict = field(default_factory=dict)
    run_manifest: dict = field(default_factory=dict)
    source_paths: dict = field(default_factory=dict)

    # Optional time-series data for statistical analysis
    rmsd_timeseries: dict = field(default_factory=dict)  # {selection: array}
    rmsf_timeseries: dict = field(default_factory=dict)  # per-residue RMSF
    distance_timeseries: dict = field(default_factory=dict)  # {metric_name: array}
    angle_timeseries: dict = field(default_factory=dict)  # {angle_name: array}
    contact_timeseries: pd.DataFrame = field(default_factory=pd.DataFrame)


@dataclass(slots=True)
class ComparisonArtifacts:
    """Structured artifacts produced by the comparison module."""

    summary_json: Path
    comparison_table_csv: Path
    identity_comparison_csv: Path
    rmsf_region_comparison_csv: Path
    rrcs_region_comparison_csv: Path
    interaction_family_comparison_csv: Path
    quality_interface_plot: Path
    flexibility_plot: Path
    rrcs_plot: Path
    interaction_plot: Path
