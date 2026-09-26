"""Canonical D-B1 (locked 2026-05-26) feature → (layer, sub_flavor) registry.

This module is the **single source of truth** for the feature hierarchy
described in `paper_writting/methods_draft.md` Appendix D Table A1. It is
imported by:

  - `md_analysis/agent/presenters/` — to validate that the view → layer
    mapping (declared on each Presenter subclass) is consistent with the
    underlying metric → layer mapping.
  - `md_analysis/agent/tools/query_analysis_results.py` — to support the
    `spatial_layer=` / `sub_flavor=` discovery filters.
  - the §6.5.2 ablation harness — to enumerate which (metric, layer)
    cells to drop when running a hierarchy ablation.
  - documentation generators — Appendix D Table A1 in the paper is
    regenerated from this registry so the two never drift.

If you add a new metric to `md_analysis/analysis/*` you MUST add an
entry here in the same patch. The CI guard (planned, D-B1 follow-up)
will walk the analysis tree and fail the build when a metric has no
registry entry.

D-B1 vocabulary:
  spatial_layer ∈ {Complex, Interface, Region, Residue, Pair}
  sub_flavor    ∈ {static, dynamic}

Sub-flavor classification rules (auditable, from methods_draft.md §D):
  1. Categorical labels whose value depends on persistence or rotamer
     behavior are dynamic (e.g. chemistry_tags, risk_flags, persistence
     profile).
  2. Time-averaged scalars derived from a fluctuating quantity are
     static when the reported value is the scalar itself; dynamic when
     the reported value is a distribution shape, ensemble population,
     or temporal pattern.
  3. Internal utilities / on-demand artifacts are tagged with the
     `exposure` field (`INTERNAL_ONLY` / `RE_ANALYSIS_ONLY`) and do not
     enter the default prompt regardless of sub-flavor.

Counts: 82 metrics across Complex 23 / Interface 4 / Region 1 /
Residue 42 / Pair 12.  Sub-flavor totals: static 48 / dynamic 34.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Literal, Optional, Tuple

SpatialLayer = Literal["complex", "interface", "region", "residue", "pair"]
SubFlavor = Literal["static", "dynamic"]
Exposure = Literal["DEFAULT", "INTERNAL_ONLY", "RE_ANALYSIS_ONLY"]


@dataclass(frozen=True)
class FeatureEntry:
    """One row of Appendix D Table A1."""

    metric: str
    source_module: str
    spatial_layer: SpatialLayer
    sub_flavor: SubFlavor
    default_view: str
    exposure: Exposure = "DEFAULT"
    notes: str = ""


# ---------------------------------------------------------------------------
# Table A1.1 — Complex-layer metrics (23)
# ---------------------------------------------------------------------------
_COMPLEX_METRICS: Tuple[FeatureEntry, ...] = (
    FeatureEntry("backbone_rmsd_mean", "analysis/trajectory/rmsd_refactored.py", "complex", "static", "quality", notes="time-averaged Cα RMSD vs equilibrated reference"),
    FeatureEntry("backbone_rmsd_std", "analysis/trajectory/rmsd_refactored.py", "complex", "static", "quality", notes="trajectory-wide fluctuation magnitude"),
    FeatureEntry("backbone_rmsd_min", "analysis/trajectory/rmsd_refactored.py", "complex", "static", "quality"),
    FeatureEntry("backbone_rmsd_max", "analysis/trajectory/rmsd_refactored.py", "complex", "static", "quality"),
    FeatureEntry("convergence_grade", "analysis/quality/convergence_checker.py", "complex", "static", "quality", notes="categorical A / B / C from block-analysis"),
    FeatureEntry("convergence_time_estimate", "analysis/quality/convergence_checker.py", "complex", "static", "quality", notes="ps to equilibration"),
    FeatureEntry("energy_drift", "analysis/quality/energy_quality.py", "complex", "static", "quality", notes="total-energy drift % over trajectory"),
    FeatureEntry("temperature_stability", "analysis/quality/energy_quality.py", "complex", "static", "quality", notes="σ(T) (K)"),
    FeatureEntry("pressure_stability", "analysis/quality/energy_quality.py", "complex", "static", "quality", notes="σ(P) (bar), NPT only"),
    FeatureEntry("com_distance_mean", "analysis/geometry/com_distance.py", "complex", "static", "overview"),
    FeatureEntry("com_distance_std", "analysis/geometry/com_distance.py", "complex", "static", "overview"),
    FeatureEntry("com_distance_min", "analysis/geometry/com_distance.py", "complex", "static", "overview"),
    FeatureEntry("com_distance_max", "analysis/geometry/com_distance.py", "complex", "static", "overview"),
    FeatureEntry("docking_angles", "analysis/angles/analyzer.py", "complex", "static", "angles", notes="crossing / incident / tilt; mean & std each"),
    FeatureEntry("cdr3_loop_geometric_span", "analysis/angles/cdr3_geometry.py", "complex", "static", "angles", notes="end-to-end distance, Å"),
    FeatureEntry("tica_components", "analysis/landscape/tica.py", "complex", "dynamic", "clustering", notes="ensemble projection (default N=2)"),
    FeatureEntry("tica_eigenvalues", "analysis/landscape/tica.py", "complex", "dynamic", "clustering"),
    FeatureEntry("tica_implied_timescales", "analysis/landscape/tica.py", "complex", "dynamic", "clustering"),
    FeatureEntry("free_energy_landscape", "analysis/landscape/free_energy_landscape.py", "complex", "dynamic", "clustering", notes="Boltzmann-inverted density grid on chosen CVs"),
    FeatureEntry("fel_basin_populations", "analysis/landscape/free_energy_landscape.py", "complex", "dynamic", "clustering"),
    FeatureEntry("interface_cluster_assignments", "analysis/conformation/interface_clustering.py", "complex", "dynamic", "clustering"),
    FeatureEntry("interface_cluster_summary", "analysis/conformation/interface_clustering.py", "complex", "dynamic", "clustering"),
    FeatureEntry("feature_matrix", "analysis/landscape/feature_matrix.py", "complex", "dynamic", "INTERNAL_ONLY", exposure="INTERNAL_ONLY", notes="back-end TICA/FEL/clustering input"),
)

# ---------------------------------------------------------------------------
# Table A1.2 — Interface-layer metrics (4)
# ---------------------------------------------------------------------------
_INTERFACE_METRICS: Tuple[FeatureEntry, ...] = (
    FeatureEntry("total_bsa_mean", "analysis/interface/buried_surface_area.py", "interface", "static", "interface", notes="trajectory-mean Å²"),
    FeatureEntry("bsa_mean", "analysis/interface/buried_surface_area.py", "interface", "static", "interface", notes="redundant scalar for compat"),
    FeatureEntry("bsa_std", "analysis/interface/buried_surface_area.py", "interface", "static", "interface"),
    FeatureEntry("interface_ratio_mean", "analysis/interface/buried_surface_area.py", "interface", "static", "interface", notes="2·BSA / (SASA_a + SASA_b)"),
)

# ---------------------------------------------------------------------------
# Table A1.3 — Region-layer metrics (1)
# ---------------------------------------------------------------------------
_REGION_METRICS: Tuple[FeatureEntry, ...] = (
    FeatureEntry("per_region_rmsf", "analysis/trajectory/residue_rmsf.py", "region", "dynamic", "flexibility", notes="mean over residues per CDR / MHC helix / peptide / β2m / non-groove"),
)

# ---------------------------------------------------------------------------
# Table A1.4 — Residue-layer metrics (42)
# ---------------------------------------------------------------------------
_RESIDUE_METRICS: Tuple[FeatureEntry, ...] = (
    FeatureEntry("residue_backbone_rmsf", "analysis/features/rmsf_split.py", "residue", "dynamic", "flexibility"),
    FeatureEntry("residue_sidechain_rmsf", "analysis/features/rmsf_split.py", "residue", "dynamic", "flexibility"),
    FeatureEntry("sidechain_to_backbone_rmsf_ratio", "analysis/features/rmsf_split.py", "residue", "dynamic", "flexibility", notes="design-relevant tolerance signal"),
    FeatureEntry("residue_calpha_rmsf", "analysis/trajectory/residue_rmsf.py", "residue", "dynamic", "flexibility"),
    FeatureEntry("chi1_entropy", "analysis/features/chi_entropy.py", "residue", "dynamic", "dihedrals", notes="Shannon bits over rotamer bins"),
    FeatureEntry("chi2_entropy", "analysis/features/chi_entropy.py", "residue", "dynamic", "dihedrals"),
    FeatureEntry("rotamer_diversity", "analysis/features/chi_entropy.py", "residue", "dynamic", "dihedrals", notes="0–6 well count, wells ≥ 5%"),
    FeatureEntry("chi_time_series", "analysis/geometry/sidechain_dihedrals.py", "residue", "dynamic", "RE_ANALYSIS_ONLY", exposure="RE_ANALYSIS_ONLY", notes="per-frame χ1/χ2"),
    FeatureEntry("sasa_bound", "analysis/features/sasa.py", "residue", "static", "residue"),
    FeatureEntry("sasa_unbound", "analysis/features/sasa.py", "residue", "static", "residue"),
    FeatureEntry("delta_sasa", "analysis/features/sasa.py", "residue", "static", "residue"),
    FeatureEntry("relative_exposure", "analysis/features/sasa.py", "residue", "static", "residue"),
    FeatureEntry("sidechain_sasa_bound", "analysis/features/sasa.py", "residue", "static", "residue"),
    FeatureEntry("burial_state", "analysis/features/sasa.py", "residue", "static", "residue", notes="categorical: core / exposed / interface_core / interface_rim"),
    FeatureEntry("hbond_count", "analysis/features/interactions.py", "residue", "static", "hotspots", notes="unique partners"),
    FeatureEntry("hbond_max_occupancy", "analysis/features/interactions.py", "residue", "dynamic", "hotspots"),
    FeatureEntry("saltbridge_count", "analysis/features/interactions.py", "residue", "static", "hotspots"),
    FeatureEntry("saltbridge_max_occupancy", "analysis/features/interactions.py", "residue", "dynamic", "hotspots"),
    FeatureEntry("hydrophobic_count", "analysis/features/interactions.py", "residue", "static", "hotspots"),
    FeatureEntry("hydrophobic_max_occupancy", "analysis/features/interactions.py", "residue", "dynamic", "hotspots"),
    FeatureEntry("pipi_count", "analysis/features/interactions.py", "residue", "static", "hotspots"),
    FeatureEntry("pipi_max_occupancy", "analysis/features/interactions.py", "residue", "dynamic", "hotspots"),
    FeatureEntry("cationpi_count", "analysis/features/interactions.py", "residue", "static", "hotspots"),
    FeatureEntry("cationpi_max_occupancy", "analysis/features/interactions.py", "residue", "dynamic", "hotspots"),
    FeatureEntry("dominant_interaction_type", "analysis/features/interactions.py", "residue", "static", "hotspots", notes="categorical: hbond / saltbridge / hydrophobic / pipi / cationpi / mixed / none"),
    FeatureEntry("interaction_diversity", "analysis/features/interactions.py", "residue", "static", "hotspots", notes="0–5 distinct families"),
    FeatureEntry("dominant_secondary_structure", "analysis/features/secondary.py", "residue", "static", "residue", notes="DSSP H / E / L, frame-majority"),
    FeatureEntry("ss_helix_propensity", "analysis/features/secondary.py", "residue", "static", "residue"),
    FeatureEntry("ss_sheet_propensity", "analysis/features/secondary.py", "residue", "static", "residue"),
    FeatureEntry("ss_loop_propensity", "analysis/features/secondary.py", "residue", "static", "residue"),
    FeatureEntry("ss_stability", "analysis/features/secondary.py", "residue", "static", "residue"),
    FeatureEntry("chemistry_tags", "analysis/features/design.py", "residue", "dynamic", "residue", notes="categorical, set by occupancy/χ flexibility"),
    FeatureEntry("pocket_chemistry", "analysis/features/design.py", "residue", "static", "residue", notes="hydrophobic / charged / polar / aromatic / none"),
    FeatureEntry("partner_chemistry", "analysis/features/design.py", "residue", "static", "residue"),
    FeatureEntry("risk_flags", "analysis/features/design.py", "residue", "dynamic", "residue", notes="categorical, set by persistence behavior"),
    FeatureEntry("persistence_profile", "analysis/features/design.py", "residue", "dynamic", "residue", notes="categorical: stable_full / persistent / transient / early_lost / late_gained / oscillating"),
    FeatureEntry("persistence_mean", "analysis/features/design.py", "residue", "static", "residue"),
    FeatureEntry("persistence_max", "analysis/features/design.py", "residue", "static", "residue"),
    FeatureEntry("persistence_components", "analysis/features/design.py", "residue", "static", "residue", notes="dict: profile weights and segment counts"),
    FeatureEntry("cdr3_tip_distance_to_peptide", "analysis/features/design.py", "residue", "static", "angles", notes="Å (mirror of analysis/angles/cdr3_geometry.py)"),
    FeatureEntry("cdr3_loop_span_tip_residue", "analysis/features/design.py", "residue", "static", "angles"),
    FeatureEntry("cdr3_torsion_angle", "analysis/features/design.py", "residue", "static", "angles", notes="dihedral, degrees"),
)

# ---------------------------------------------------------------------------
# Table A1.5 — Pair-layer metrics (12)
# ---------------------------------------------------------------------------
_PAIR_METRICS: Tuple[FeatureEntry, ...] = (
    FeatureEntry("rrcs_mean", "analysis/interactions/rrcs.py", "pair", "static", "pair", notes="Wang et al. RRCS"),
    FeatureEntry("rrcs_median", "analysis/interactions/rrcs.py", "pair", "static", "pair"),
    FeatureEntry("rrcs_max", "analysis/interactions/rrcs.py", "pair", "static", "pair"),
    FeatureEntry("rrcs_occupancy_fraction", "analysis/interactions/rrcs.py", "pair", "dynamic", "pair", notes="fraction of frames with non-zero RRCS"),
    FeatureEntry("hbond_occupancy_pair", "analysis/interactions/hydrogen_bond_pairs.py", "pair", "dynamic", "pair", notes="distance, angle, occupancy per pair"),
    FeatureEntry("saltbridge_occupancy_pair", "analysis/interactions/salt_bridge_pairs.py", "pair", "dynamic", "pair"),
    FeatureEntry("hydrophobic_occupancy_pair", "analysis/interactions/hydrophobic_contact_pairs.py", "pair", "dynamic", "pair"),
    FeatureEntry("pi_interaction_occupancy_pair", "analysis/interactions/pi_interaction_pairs.py", "pair", "dynamic", "pair", notes="π–π and cation–π reported separately"),
    FeatureEntry("residue_contact_frequency_pair", "analysis/trajectory/residue_contacts.py", "pair", "dynamic", "pair"),
    FeatureEntry("occupancy_frame_segments", "analysis/interactions/occupancy_metrics.py", "pair", "dynamic", "RE_ANALYSIS_ONLY", exposure="RE_ANALYSIS_ONLY", notes="JSON list of contiguous frame ranges"),
    FeatureEntry("occupancy_n_segments", "analysis/interactions/occupancy_metrics.py", "pair", "dynamic", "RE_ANALYSIS_ONLY", exposure="RE_ANALYSIS_ONLY"),
    FeatureEntry("occupancy_max_consecutive", "analysis/interactions/occupancy_metrics.py", "pair", "dynamic", "RE_ANALYSIS_ONLY", exposure="RE_ANALYSIS_ONLY", notes="longest uninterrupted contact span"),
)


FEATURE_REGISTRY: Tuple[FeatureEntry, ...] = (
    _COMPLEX_METRICS
    + _INTERFACE_METRICS
    + _REGION_METRICS
    + _RESIDUE_METRICS
    + _PAIR_METRICS
)


def by_metric(name: str) -> Optional[FeatureEntry]:
    """Look up a registry row by metric name. Returns None when absent."""
    for entry in FEATURE_REGISTRY:
        if entry.metric == name:
            return entry
    return None


def by_layer(layer: SpatialLayer) -> List[FeatureEntry]:
    """Return all entries belonging to a spatial layer."""
    return [e for e in FEATURE_REGISTRY if e.spatial_layer == layer]


def by_sub_flavor(
    sub_flavor: SubFlavor,
    layer: Optional[SpatialLayer] = None,
) -> List[FeatureEntry]:
    """Return entries serving the given sub-flavor, optionally within a layer."""
    return [
        e for e in FEATURE_REGISTRY
        if e.sub_flavor == sub_flavor and (layer is None or e.spatial_layer == layer)
    ]


def by_default_view(view: str) -> List[FeatureEntry]:
    """Return entries whose primary surfacing view matches."""
    return [e for e in FEATURE_REGISTRY if e.default_view == view]


def hierarchy_counts() -> Dict[str, Dict[str, int]]:
    """Reproduce the 5 × 2 distribution table from methods_draft.md §D."""
    layers: Tuple[SpatialLayer, ...] = ("complex", "interface", "region", "residue", "pair")
    flavors: Tuple[SubFlavor, ...] = ("static", "dynamic")
    out: Dict[str, Dict[str, int]] = {}
    for layer in layers:
        row: Dict[str, int] = {f: 0 for f in flavors}
        for e in FEATURE_REGISTRY:
            if e.spatial_layer == layer:
                row[e.sub_flavor] += 1
        row["subtotal"] = sum(row[f] for f in flavors)
        out[layer] = row
    return out


def assert_complete() -> None:
    """Auditable check: the 5×2 grid covers every metric.

    Raises:
        AssertionError if any entry has an unknown spatial_layer or
        sub_flavor, or if the totals disagree with the locked counts.
    """
    valid_layers = {"complex", "interface", "region", "residue", "pair"}
    valid_flavors = {"static", "dynamic"}
    for entry in FEATURE_REGISTRY:
        assert entry.spatial_layer in valid_layers, f"{entry.metric}: bad layer {entry.spatial_layer}"
        assert entry.sub_flavor in valid_flavors, f"{entry.metric}: bad sub-flavor {entry.sub_flavor}"
    counts = hierarchy_counts()
    expected = {
        "complex":   {"static": 15, "dynamic": 8,  "subtotal": 23},
        "interface": {"static": 4,  "dynamic": 0,  "subtotal": 4},
        "region":    {"static": 0,  "dynamic": 1,  "subtotal": 1},
        "residue":   {"static": 26, "dynamic": 16, "subtotal": 42},
        "pair":      {"static": 3,  "dynamic": 9,  "subtotal": 12},
    }
    assert counts == expected, (
        "Feature registry counts have drifted from the locked D-B1 totals "
        f"(2026-05-26).\n  expected: {expected}\n  actual:   {counts}"
    )


__all__ = [
    "FeatureEntry",
    "FEATURE_REGISTRY",
    "by_metric",
    "by_layer",
    "by_sub_flavor",
    "by_default_view",
    "hierarchy_counts",
    "assert_complete",
]
