"""Design-oriented feature computers for mutation prioritization."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import json
from typing import Iterable

import pandas as pd

from .core.locator import CaseLocator
from .core.models import FeatureProvenance, FeatureSet, ResidueFeatures, ResidueKey
from .core.registry import FeatureComputer, register_feature
from .examples import _load_rrcs_pair_table


AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
POSITIVE = {"LYS", "ARG", "HIS"}
NEGATIVE = {"ASP", "GLU"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "HIS"}
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
SPECIAL_RISK = {"CYS", "PRO", "GLY"}


@dataclass
class _ResidueAggregate:
    """Internal accumulator built from pair-level analysis outputs."""

    key: ResidueKey
    rrcs_sum: float = 0.0
    max_rrcs: float = 0.0
    occupancy_sum: float = 0.0
    occupancy_count: int = 0
    occupancy_max: float = 0.0
    partners: set[tuple[str, int, str]] = field(default_factory=set)
    partner_components: set[str] = field(default_factory=set)
    regions: set[str] = field(default_factory=set)
    tags: set[str] = field(default_factory=set)
    risks: set[str] = field(default_factory=set)

    @property
    def mean_occupancy(self) -> float:
        """Mean occupancy over all pairs involving this residue.

        Note: the annotated RRCS pair table enumerates every candidate partner,
        so most pairs have occupancy ≈ 0. This mean is therefore heavily
        diluted and rarely reflects whether the residue has a strong contact.
        Prefer `max_occupancy` for design decisions.
        """
        if self.occupancy_count == 0:
            return 0.0
        return self.occupancy_sum / self.occupancy_count

    @property
    def max_occupancy(self) -> float:
        """Maximum occupancy across all partners — answers \"does this residue
        have at least one strong contact?\". This is the right signal for
        risk flags, persistence classification, and the mutability score
        occupancy component."""
        return self.occupancy_max


def _source_rel(path: Path | None, locator: CaseLocator) -> str:
    if path is None:
        return ""
    try:
        return str(path.relative_to(locator.case_dir))
    except ValueError:
        return str(path)


def _as_float(value: object, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _as_int(value: object) -> int | None:
    try:
        if pd.isna(value):
            return None
        return int(value)
    except (TypeError, ValueError):
        return None


def _clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _residue_key_from_row(row: pd.Series, side: int) -> ResidueKey | None:
    chain = _clean(row.get(f"chain_id_{side}"))
    resid = _as_int(row.get(f"resid_{side}"))
    resname = _clean(row.get(f"resname_{side}")).upper()
    if chain and resid is not None:
        return ResidueKey(chain=chain, resid=resid, resname=resname)

    prefix = "tcr" if side == 1 else "partner"
    chain = _clean(row.get(f"{prefix}_chain"))
    resid = _as_int(row.get(f"{prefix}_resid"))
    resname = _clean(row.get(f"{prefix}_resname")).upper()
    if chain and resid is not None:
        return ResidueKey(chain=chain, resid=resid, resname=resname)
    return None


def _row_region(row: pd.Series, side: int) -> str:
    if side == 1:
        for column in ("tcr_region_detailed", "tcr_region", "region_1", "region"):
            value = _clean(row.get(column))
            if value:
                return value
    for column in (f"region_{side}", "partner_region", "partner_component"):
        value = _clean(row.get(column))
        if value:
            return value
    return ""


def _partner_component(row: pd.Series, side: int) -> str:
    if side == 1:
        return _clean(row.get("partner_component")) or _clean(row.get("component_2"))
    return _clean(row.get("tcr_component")) or _clean(row.get("component_1")) or "TCR"


def _interaction_tags(row: pd.Series) -> set[str]:
    tags: set[str] = set()
    text = " ".join(
        _clean(row.get(column)).lower()
        for column in (
            "interaction_family",
            "interaction_type",
            "interaction_types",
            "interaction_category",
        )
    )
    if "hbond" in text or "hydrogen" in text:
        tags.add("hbond_network_hub")
    if "salt" in text or "ionic" in text:
        tags.add("salt_bridge_participant")
    if "hydrophobic" in text:
        tags.add("hydrophobic_contact")
    if "pi" in text or "π" in text:
        tags.add("pi_interaction")
    return tags


def _intrinsic_chemistry_tags(resname: str) -> set[str]:
    tags: set[str] = set()
    if resname in AROMATIC:
        tags.add("aromatic_residue")
    if resname in POSITIVE:
        tags.add("positive_residue")
    if resname in NEGATIVE:
        tags.add("negative_residue")
    if resname in POLAR:
        tags.add("polar_residue")
    if resname in HYDROPHOBIC:
        tags.add("hydrophobic_residue")
    return tags


def _is_peptide_anchor(key: ResidueKey, agg: _ResidueAggregate) -> bool:
    if not any("peptide" in component.lower() for component in agg.partner_components | agg.regions):
        return False
    return key.resid in {2, 9}


def _is_framework_region(region: str) -> bool:
    text = region.lower()
    return "framework" in text or "non_cdr" in text or text in {"fr", "fr1", "fr2", "fr3", "fr4"}


def _is_cdr_region(region: str) -> bool:
    return "cdr" in region.lower()


def _is_hla_facing(components: Iterable[str]) -> bool:
    joined = " ".join(component.lower() for component in components)
    return "hla" in joined or "mhc" in joined


def _is_peptide_facing(components: Iterable[str]) -> bool:
    return "peptide" in " ".join(component.lower() for component in components)


def _build_aggregates(df: pd.DataFrame) -> dict[ResidueKey, _ResidueAggregate]:
    aggregates: dict[ResidueKey, _ResidueAggregate] = {}

    for _, row in df.iterrows():
        k1 = _residue_key_from_row(row, 1)
        k2 = _residue_key_from_row(row, 2)
        if k1 is None or k2 is None:
            continue

        mean_rrcs = _as_float(row.get("mean_rrcs"))
        occupancy = _as_float(row.get("rrcs_nonzero_fraction", row.get("occupancy")))
        row_tags = _interaction_tags(row)

        for key, partner, side in ((k1, k2, 1), (k2, k1, 2)):
            agg = aggregates.setdefault(key, _ResidueAggregate(key=key))
            agg.rrcs_sum += mean_rrcs
            agg.max_rrcs = max(agg.max_rrcs, mean_rrcs)
            agg.occupancy_sum += occupancy
            agg.occupancy_count += 1
            agg.occupancy_max = max(agg.occupancy_max, occupancy)
            agg.partners.add((partner.chain, partner.resid, partner.resname))
            component = _partner_component(row, side)
            if component:
                agg.partner_components.add(component)
            region = _row_region(row, side)
            if region:
                agg.regions.add(region)
            agg.tags.update(row_tags)
            agg.tags.update(_intrinsic_chemistry_tags(key.resname))

    return aggregates


def _safe_residue_key_from_label(label: object) -> ResidueKey | None:
    text = _clean(label)
    if not text:
        return None
    try:
        return ResidueKey.parse(text)
    except Exception:
        return None


def _load_occupancy_pair_tables(locator: CaseLocator) -> list[tuple[pd.DataFrame, Path]]:
    tables: list[tuple[pd.DataFrame, Path]] = []
    seen: set[Path] = set()
    module_names = [
        "contact_occupancy",
        "hbond_occupancy",
        "saltbridge_occupancy",
        "hydrophobic_occupancy",
        "pipi_occupancy",
        "cationpi_occupancy",
    ]
    candidate_paths: list[Path] = []
    for module in module_names:
        root = locator.get_module_root(module)
        if root is not None:
            candidate_paths.extend([
                root / "pair_stability.csv",
                root / "persistent_interaction_ranking.csv",
            ])
    candidate_paths.extend(locator.case_dir.glob("analysis/**/pair_stability.csv"))
    candidate_paths.extend(locator.case_dir.glob("analysis/**/persistent_interaction_ranking.csv"))

    for path in candidate_paths:
        path = path.resolve()
        if path in seen or not path.exists():
            continue
        seen.add(path)
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if not df.empty:
            tables.append((df, path))
    return tables


def _as_segments(value: object) -> list[list[int]]:
    if isinstance(value, list):
        return value
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        return []
    if not isinstance(parsed, list):
        return []
    segments: list[list[int]] = []
    for segment in parsed:
        if not isinstance(segment, list) or len(segment) < 2:
            continue
        start = _as_int(segment[0])
        end = _as_int(segment[1])
        if start is None or end is None:
            continue
        segments.append([start, end])
    return segments


def _classify_persistence_from_row(row: pd.Series) -> tuple[str, dict[str, float]]:
    occupancy = _as_float(row.get("contact_frequency", row.get("occupancy")))
    n_segments = _as_float(row.get("n_segments"))
    max_consecutive = _as_float(row.get("max_consecutive_fraction"))
    first_frame = _as_float(row.get("first_frame"), default=-1.0)
    last_frame = _as_float(row.get("last_frame"), default=-1.0)
    total_frames = _as_float(row.get("total_frames"))

    segments = _as_segments(row.get("frame_segments"))
    if segments:
        first_frame = float(segments[0][0])
        last_frame = float(segments[-1][1])
        n_segments = float(len(segments))

    early = late = full = 0.0
    if total_frames > 0 and first_frame >= 0 and last_frame >= 0:
        early = 1.0 if last_frame <= total_frames * 0.45 else 0.0
        late = 1.0 if first_frame >= total_frames * 0.55 else 0.0
        full = 1.0 if first_frame <= total_frames * 0.20 and last_frame >= total_frames * 0.80 else 0.0

    if occupancy >= 0.80 and (full or max_consecutive >= 0.50 or n_segments <= 2):
        profile = "stable_full"
    elif early and occupancy >= 0.20:
        profile = "early_lost"
    elif late and occupancy >= 0.20:
        profile = "late_gained"
    elif n_segments >= 4 and occupancy >= 0.20:
        profile = "oscillating"
    elif occupancy >= 0.50:
        profile = "persistent"
    else:
        profile = "transient"

    return profile, {
        "occupancy": round(occupancy, 4),
        "n_segments": round(n_segments, 4),
        "max_consecutive_fraction": round(max_consecutive, 4),
        "covers_full_window": full,
    }


def _classify_persistence_from_aggregate(agg: _ResidueAggregate) -> str:
    occ = agg.max_occupancy
    if occ >= 0.80:
        return "stable_full"
    if occ >= 0.50:
        return "persistent"
    return "transient"


def _load_cdr3_geometry_table(locator: CaseLocator) -> tuple[pd.DataFrame | None, Path | None]:
    root = locator.get_module_root("cdr3_geometry")
    candidates: list[Path] = []
    if root is not None:
        candidates.extend([
            root / "cdr3_geometry.csv",
            root / "cdr3_geometry_summary.csv",
        ])
    candidates.extend([
        locator.case_dir / "analysis/angles/cdr3_geometry/cdr3_geometry.csv",
        locator.case_dir / "analysis/geometry/cdr3_geometry/cdr3_geometry.csv",
        locator.case_dir / "analysis/cdr3_geometry/cdr3_geometry.csv",
    ])
    for path in candidates:
        if not path.exists():
            continue
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        return df, path
    return None, None


def _normalize(value: float, max_value: float) -> float:
    if max_value <= 0:
        return 0.0
    return max(0.0, min(1.0, value / max_value))


def _risk_flags_for(agg: _ResidueAggregate) -> list[str]:
    risks: set[str] = set()
    regions = sorted(agg.regions)
    region_text = " ".join(regions)

    if any(_is_framework_region(region) for region in regions):
        risks.add("framework_or_non_cdr_region")
    if agg.key.resname in SPECIAL_RISK:
        risks.add(f"{agg.key.resname.lower()}_mutation_risk")
    if agg.max_rrcs >= 4.0 and agg.max_occupancy >= 0.8:
        risks.add("highly_persistent_core_contact")
    if _is_peptide_anchor(agg.key, agg):
        risks.add("possible_peptide_anchor_position")
    if _is_hla_facing(agg.partner_components) and not _is_peptide_facing(agg.partner_components):
        risks.add("hla_facing_specificity_risk")
    if "salt_bridge_participant" in agg.tags:
        risks.add("salt_bridge_disruption_risk")
    if not region_text and not agg.partner_components:
        risks.add("low_annotation_confidence")
    return sorted(risks)


@register_feature
class ChemistryTagsComputer(FeatureComputer):
    """Assign chemistry and interaction semantic tags to interface residues."""

    name = "chemistry_tags"
    required_modules = ["rrcs"]
    produces = ["chemistry_tags", "pocket_chemistry", "partner_chemistry"]

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        df, path = _load_rrcs_pair_table(locator)
        if df is None:
            features.metadata.setdefault("missing_inputs", []).append("rrcs pair table")
            return features

        source = _source_rel(path, locator)
        aggregates = _build_aggregates(df)
        for key, agg in aggregates.items():
            tags = sorted(agg.tags)
            pocket = self._dominant_pocket(tags)
            partner = ",".join(sorted(agg.partner_components)) or None
            features.upsert(
                ResidueFeatures(
                    residue=key,
                    chemistry_tags=tags,
                    pocket_chemistry=pocket,
                    partner_chemistry=partner,
                    is_interface=True,
                    provenance=FeatureProvenance(
                        sources=[source] if source else [],
                        computer=self.name,
                    ),
                )
            )
        return features

    @staticmethod
    def _dominant_pocket(tags: list[str]) -> str | None:
        if "hydrophobic_contact" in tags or "hydrophobic_residue" in tags:
            return "hydrophobic"
        if "salt_bridge_participant" in tags:
            return "charged"
        if "hbond_network_hub" in tags or "polar_residue" in tags:
            return "polar"
        if "pi_interaction" in tags or "aromatic_residue" in tags:
            return "aromatic"
        return None


@register_feature
class RiskFlagComputer(FeatureComputer):
    """Flag residues that should be handled conservatively in design."""

    name = "risk_flag"
    required_modules = ["rrcs"]
    produces = ["risk_flags"]

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        df, path = _load_rrcs_pair_table(locator)
        if df is None:
            features.metadata.setdefault("missing_inputs", []).append("rrcs pair table")
            return features

        source = _source_rel(path, locator)
        for key, agg in _build_aggregates(df).items():
            flags = _risk_flags_for(agg)
            if not flags:
                continue
            features.upsert(
                ResidueFeatures(
                    residue=key,
                    risk_flags=flags,
                    is_core="highly_persistent_core_contact" in flags,
                    provenance=FeatureProvenance(
                        sources=[source] if source else [],
                        computer=self.name,
                    ),
                )
            )
        return features


@register_feature
class MutabilityScoreComputer(FeatureComputer):
    """Compute a 0-1 design-priority score for interface residues.

    DEPRECATED (D-B5, 2026-05-26): This computer is no longer invoked from
    the LLM-facing pipeline. The composite formula (0.35·rrcs + 0.25·occ +
    0.15·redundancy + 0.25·region) × (1 − risk_penalty) was rejected at
    decision lock because feeding the agent a hand-tuned pre-ranked list
    collapses its reasoning into post-hoc rationalization of the heuristic
    order. Pure-RRCS-rank (Methods §6.2) is the replacement mono-signal
    baseline. The class is kept for now because legacy back-end consumers
    (`web/routers/design.py`, `recommendation/data_formatter.py`,
    `recommendation/prompts.py`) still import it; those callsites are
    scheduled for Phase-2 cleanup. Do NOT add new LLM-facing callers.

    Stored on ResidueFeatures as `design_priority_score`. The registry name
    `"mutability_score"` is kept as the stable string key for backwards
    compatibility with downstream callers of `compute_feature(...)`.
    """

    name = "mutability_score"
    required_modules = ["rrcs"]
    produces = ["design_priority_score", "priority_components"]

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        df, path = _load_rrcs_pair_table(locator)
        if df is None:
            features.metadata.setdefault("missing_inputs", []).append("rrcs pair table")
            return features

        source = _source_rel(path, locator)
        aggregates = _build_aggregates(df)
        max_rrcs = max((agg.rrcs_sum for agg in aggregates.values()), default=0.0)
        max_partner_count = max((len(agg.partners) for agg in aggregates.values()), default=1)

        for key, agg in aggregates.items():
            regions = sorted(agg.regions)
            rrcs_component = _normalize(agg.rrcs_sum, max_rrcs)
            occupancy_component = max(0.0, min(1.0, agg.max_occupancy))
            redundancy_component = _normalize(float(len(agg.partners)), float(max_partner_count))
            region_component = self._region_component(regions, agg.partner_components)
            risk_flags = set(_risk_flags_for(agg))
            risk_penalty = min(0.75, 0.15 * len(risk_flags))

            raw_score = (
                0.35 * rrcs_component
                + 0.25 * occupancy_component
                + 0.15 * redundancy_component
                + 0.25 * region_component
            )
            score = max(0.0, min(1.0, raw_score * (1.0 - risk_penalty)))

            features.upsert(
                ResidueFeatures(
                    residue=key,
                    design_priority_score=round(score, 4),
                    priority_components={
                        "rrcs": round(rrcs_component, 4),
                        "occupancy": round(occupancy_component, 4),
                        "redundancy": round(redundancy_component, 4),
                        "region": round(region_component, 4),
                        "risk_penalty": round(risk_penalty, 4),
                    },
                    contact_redundancy=round(redundancy_component, 4),
                    partner_diversity=len(agg.partner_components),
                    persistence_mean=round(agg.max_occupancy, 4),
                    region=regions[0] if regions else None,
                    is_interface=True,
                    risk_flags=sorted(risk_flags),
                    provenance=FeatureProvenance(
                        sources=[source] if source else [],
                        computer=self.name,
                        notes=[
                            "weights: rrcs=0.35, occupancy=0.25, redundancy=0.15, "
                            "region=0.25; risk flags apply multiplicative penalty"
                        ],
                    ),
                )
            )
        return features

    @staticmethod
    def _region_component(regions: list[str], partner_components: set[str]) -> float:
        if any(_is_peptide_facing([component]) for component in partner_components):
            peptide_bonus = 0.15
        else:
            peptide_bonus = 0.0

        if any("cdr3" in region.lower() for region in regions):
            return min(1.0, 0.95 + peptide_bonus)
        if any(_is_cdr_region(region) for region in regions):
            return min(1.0, 0.75 + peptide_bonus)
        if any("peptide" in region.lower() for region in regions):
            return 0.65
        if any(_is_framework_region(region) for region in regions):
            return 0.2
        return 0.45 + peptide_bonus


@register_feature
class PersistenceProfileComputer(FeatureComputer):
    """Classify whether each interface residue forms stable or time-local contacts."""

    name = "persistence_profile"
    required_modules = ["rrcs"]
    produces = [
        "persistence_profile",
        "persistence_components",
        "persistence_mean",
        "persistence_max",
    ]

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        occupancy_tables = _load_occupancy_pair_tables(locator)
        if occupancy_tables:
            return self._compute_from_occupancy_tables(locator, features, occupancy_tables)
        return self._compute_from_rrcs_fallback(locator, features)

    def _compute_from_occupancy_tables(
        self,
        locator: CaseLocator,
        features: FeatureSet,
        tables: list[tuple[pd.DataFrame, Path]],
    ) -> FeatureSet:
        residue_profiles: dict[ResidueKey, dict[str, float]] = {}
        residue_components: dict[ResidueKey, dict[str, float]] = {}
        sources = [_source_rel(path, locator) for _, path in tables]

        for df, _ in tables:
            for _, row in df.iterrows():
                profile, components = _classify_persistence_from_row(row)
                occupancy = components["occupancy"]
                for column in ("tcr_residue", "phla_residue"):
                    key = _safe_residue_key_from_label(row.get(column))
                    if key is None:
                        continue
                    weights = residue_profiles.setdefault(key, {})
                    weights[profile] = weights.get(profile, 0.0) + max(occupancy, 0.01)
                    comp = residue_components.setdefault(
                        key,
                        {
                            "occupancy_sum": 0.0,
                            "occupancy_count": 0.0,
                            "occupancy_max": 0.0,
                            "segment_sum": 0.0,
                        },
                    )
                    comp["occupancy_sum"] += occupancy
                    comp["occupancy_count"] += 1.0
                    comp["occupancy_max"] = max(comp["occupancy_max"], occupancy)
                    comp["segment_sum"] += components["n_segments"]

        for key, weights in residue_profiles.items():
            profile = max(weights.items(), key=lambda item: (item[1], item[0]))[0]
            comp = residue_components[key]
            count = max(comp["occupancy_count"], 1.0)
            features.upsert(
                ResidueFeatures(
                    residue=key,
                    persistence_profile=profile,
                    persistence_mean=round(comp["occupancy_sum"] / count, 4),
                    persistence_max=round(comp["occupancy_max"], 4),
                    persistence_components={
                        "mean_n_segments": round(comp["segment_sum"] / count, 4),
                        **{f"profile_weight_{name}": round(value, 4) for name, value in sorted(weights.items())},
                    },
                    is_interface=True,
                    provenance=FeatureProvenance(
                        sources=[source for source in sources if source],
                        computer=self.name,
                    ),
                )
            )
        return features

    def _compute_from_rrcs_fallback(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        df, path = _load_rrcs_pair_table(locator)
        if df is None:
            features.metadata.setdefault("missing_inputs", []).append("rrcs pair table")
            return features

        source = _source_rel(path, locator)
        for key, agg in _build_aggregates(df).items():
            profile = _classify_persistence_from_aggregate(agg)
            features.upsert(
                ResidueFeatures(
                    residue=key,
                    persistence_profile=profile,
                    # persistence_mean must be a real mean to match the
                    # occupancy-tables path (line ~641); previously this
                    # branch incorrectly aliased it to max_occupancy.
                    # Note: on full-size cases the annotated RRCS pair
                    # table enumerates every candidate partner, so this
                    # mean can be diluted by zero-occupancy rows — see
                    # _ResidueAggregate.mean_occupancy docstring. The
                    # `persistence_components` dict below preserves both
                    # signals so downstream consumers can pick.
                    persistence_mean=round(agg.mean_occupancy, 4),
                    persistence_max=round(agg.max_occupancy, 4),
                    persistence_components={
                        "max_pair_occupancy": round(agg.max_occupancy, 4),
                        "mean_pair_occupancy": round(agg.mean_occupancy, 4),
                        "n_partners": float(len(agg.partners)),
                    },
                    is_interface=True,
                    provenance=FeatureProvenance(
                        sources=[source] if source else [],
                        computer=self.name,
                        notes=[
                            "explicit occupancy pair-stability tables were not found; "
                            "used RRCS nonzero fraction as a persistence fallback"
                        ],
                    ),
                )
            )
        return features


@register_feature
class CDR3GeometryFeatureComputer(FeatureComputer):
    """Translate CDR3 loop geometry summaries into residue-level features."""

    name = "cdr3_geometry"
    required_modules = ["cdr3_geometry"]
    produces = ["cdr3_tip_distance", "cdr3_loop_span", "cdr3_torsion_deg"]

    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        df, path = _load_cdr3_geometry_table(locator)
        if df is None:
            features.metadata.setdefault("missing_inputs", []).append("cdr3_geometry table")
            return features

        source = _source_rel(path, locator)
        for _, row in df.iterrows():
            key = _safe_residue_key_from_label(row.get("tip_residue_label"))
            if key is None:
                continue
            features.upsert(
                ResidueFeatures(
                    residue=key,
                    cdr3_tip_distance=_as_float(row.get("tip_to_peptide_distance_angstrom"), default=None),
                    cdr3_loop_span=_as_float(row.get("loop_span_angstrom"), default=None),
                    cdr3_torsion_deg=_as_float(row.get("torsion_deg"), default=None),
                    region="CDR3",
                    is_interface=True,
                    provenance=FeatureProvenance(
                        sources=[source] if source else [],
                        computer=self.name,
                    ),
                )
            )
        return features
