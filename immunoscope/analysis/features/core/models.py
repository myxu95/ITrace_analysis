"""
Data models for the features layer.

The features layer turns raw analysis outputs into design-oriented features.
These features are the "adapter/translator" between:

  Analysis platform        →  Feature layer (this package)  →   Design platform
  (physical quantities:        (design language:                 (Agent prompts,
   RRCS, BSA, RMSF, ...)        mutability score,                 candidate ranking,
                                risk flags,                       design decisions)
                                chemistry tags, ...)

Everything in this module is pure data — no I/O, no computation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple


def _now_iso() -> str:
    """UTC timestamp in ISO format (timezone-aware)."""
    return datetime.now(timezone.utc).isoformat()


class FeatureError(Exception):
    """Raised when feature computation fails or required inputs are missing."""
    pass


# ---------------------------------------------------------------------- #
# Residue identity
# ---------------------------------------------------------------------- #

@dataclass(frozen=True)
class ResidueKey:
    """
    Canonical residue identifier used across all features.

    Kept hashable so it can be used as a dict key or in sets.

    Attributes:
        chain: Chain identifier (e.g. 'A', 'B', 'P' for peptide, ...).
        resid: Residue sequence number.
        resname: 3-letter residue name (e.g. 'TYR'). Optional — may be empty
            for a purely positional reference.
    """

    chain: str
    resid: int
    resname: str = ""

    def label(self) -> str:
        """Human-readable label: e.g. 'TYR-A95' or 'A95' if resname missing."""
        if self.resname:
            return f"{self.resname}-{self.chain}{self.resid}"
        return f"{self.chain}{self.resid}"

    @classmethod
    def parse(cls, text: str) -> "ResidueKey":
        """
        Parse a label of the form 'TYR-A95' or 'A95' back into a ResidueKey.

        This is intentionally permissive — extend as more formats appear.
        """
        text = text.strip()
        if "-" in text:
            resname, tail = text.split("-", 1)
        else:
            resname, tail = "", text

        # tail should be like 'A95' → chain 'A', resid 95
        i = 0
        while i < len(tail) and tail[i].isalpha():
            i += 1
        chain = tail[:i]
        try:
            resid = int(tail[i:])
        except ValueError as exc:
            raise FeatureError(f"Cannot parse residue label: {text!r}") from exc

        return cls(chain=chain, resid=resid, resname=resname)


# ---------------------------------------------------------------------- #
# Provenance
# ---------------------------------------------------------------------- #

@dataclass
class FeatureProvenance:
    """
    Tracks where a feature value came from.

    Provenance is critical because features are consumed by the Agent for
    design decisions — every derived number must be traceable back to the
    raw analysis file it came from.

    Attributes:
        sources: Relative paths (under case_dir) of files read to compute
            this feature. Always present; empty means the feature was
            synthesized without reading files.
        computer: Name of the FeatureComputer that produced the value.
        computed_at: ISO-format timestamp.
        missing_inputs: Names of expected inputs that were not found. A
            non-empty list means the feature is partial / best-effort.
        notes: Free-form notes (e.g. "used fallback formula because X missing").
    """

    sources: List[str] = field(default_factory=list)
    computer: str = ""
    computed_at: str = field(default_factory=_now_iso)
    missing_inputs: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def is_partial(self) -> bool:
        """True if any expected inputs were missing."""
        return bool(self.missing_inputs)


# ---------------------------------------------------------------------- #
# Per-residue features
# ---------------------------------------------------------------------- #

@dataclass
class ResidueFeatures:
    """
    All design-relevant features for a single residue.

    Fields are intentionally optional — a residue's feature set may be
    partially populated depending on which analysis modules ran and which
    feature computers were invoked.

    Tier 1 features (implemented in PR 4):
      - bsa_contribution, bsa_rank
      - rrcs_contribution, rrcs_rank
      - contact_redundancy
      - partner_diversity
      - persistence_mean
      - persistence_profile
      - design_priority_score (alias: mutability_score)
      - priority_components (alias: mutability_components)
      - risk_flags
      - chemistry_tags

    Tier 2/3 features (future):
      - pocket_chemistry, partner_chemistry
      - conformational_coupling
      - cdr3 geometry metrics
      - sasa, depth
      - evolutionary_conservation
    """

    residue: ResidueKey

    # --- Tier 1: contribution metrics --------------------------------- #
    bsa_contribution: Optional[float] = None   # Å² buried by this residue
    bsa_rank: Optional[int] = None             # 1 = largest contributor
    rrcs_contribution: Optional[float] = None  # summed RRCS score
    rrcs_rank: Optional[int] = None
    contact_count: Optional[int] = None        # # of interface partners

    # --- Tier 1: robustness / redundancy ------------------------------ #
    contact_redundancy: Optional[float] = None  # [0,1]: 0=unique, 1=many
    partner_diversity: Optional[int] = None     # # distinct partner regions

    # --- Tier 1: persistence ------------------------------------------ #
    persistence_mean: Optional[float] = None   # avg occupancy over frames
    persistence_max: Optional[float] = None
    persistence_profile: Optional[str] = None  # stable_full, early_lost, ...
    persistence_components: Dict[str, float] = field(default_factory=dict)

    # --- Tier 1: synthesized scores ----------------------------------- #
    # Design priority: higher = stronger candidate for mutation design.
    # Despite the legacy "mutability" name, this is a priority score
    # (contact strength + occupancy + redundancy + region weight, minus
    # risk penalty), NOT a measure of how easily the residue can be swapped.
    design_priority_score: Optional[float] = None   # [0,1]
    priority_components: Dict[str, float] = field(default_factory=dict)
    risk_flags: List[str] = field(default_factory=list)
    chemistry_tags: List[str] = field(default_factory=list)

    # --- Tier 2: chemistry (future) ----------------------------------- #
    pocket_chemistry: Optional[str] = None      # e.g. "hydrophobic", "polar"
    partner_chemistry: Optional[str] = None

    # --- Tier 3: dynamics (future) ------------------------------------ #
    rmsf_mean: Optional[float] = None
    rmsf_rank: Optional[int] = None
    cdr3_tip_distance: Optional[float] = None
    cdr3_loop_span: Optional[float] = None
    cdr3_torsion_deg: Optional[float] = None

    # --- Tier 1: backbone/sidechain RMSF split (Phase 2 / A2) --------- #
    # Per-residue mean RMSF over backbone (N/CA/C/O) and sidechain (heavy,
    # non-backbone) atom sets. rmsf_ratio = sidechain / backbone:
    #   > ~2.0  → rotamer-flexible on a rigid scaffold (ideal mutation target)
    #   < ~1.0  → anchored sidechain on a moving loop (handle with care)
    rmsf_backbone: Optional[float] = None
    rmsf_sidechain: Optional[float] = None
    rmsf_ratio: Optional[float] = None

    # --- Tier 1: per-residue SASA (Phase 2 / A3) ---------------------- #
    # All units Å². Bound = in full complex; unbound = own group alone.
    # delta_sasa = sasa_unbound − sasa_bound (≈ buried by partner).
    # relative_exposure = sasa_bound / sasa_unbound  (0 fully buried, 1 untouched).
    # burial_state: interface_core / interface_rim / exposed / core / unknown.
    sasa_bound: Optional[float] = None
    sasa_unbound: Optional[float] = None
    delta_sasa: Optional[float] = None
    relative_exposure: Optional[float] = None
    sasa_sidechain_bound: Optional[float] = None
    burial_state: Optional[str] = None

    # --- Tier 1: interaction-family per-residue counts ---------------- #
    # Aggregated from residue_pair_<family>.csv (one count per partner per family).
    # Missing values mean the family was not run for this case.
    hbond_count: Optional[int] = None
    saltbridge_count: Optional[int] = None
    hydrophobic_count: Optional[int] = None
    pipi_count: Optional[int] = None
    cationpi_count: Optional[int] = None

    # Per-family max occupancy (a single strong partner often matters more than many weak)
    hbond_max_occupancy: Optional[float] = None
    saltbridge_max_occupancy: Optional[float] = None
    hydrophobic_max_occupancy: Optional[float] = None
    pipi_max_occupancy: Optional[float] = None
    cationpi_max_occupancy: Optional[float] = None

    # Aggregate interaction signals
    dominant_interaction_type: Optional[str] = None   # hbond | saltbridge | hydrophobic | pipi | cationpi | mixed | none
    interaction_diversity: Optional[int] = None       # 0-5 distinct families

    # --- Tier 1: secondary structure (Phase 2 / A4) ------------------- #
    # Dominant DSSP code over the sampled trajectory (H = helix, E = sheet,
    # L = loop). ss_propensity_* are the per-frame fractions; ss_stability
    # is the share of frames matching the dominant code. Used in design to
    # weight against helix/sheet-breakers in stable secondary elements.
    secondary_structure: Optional[str] = None          # H | E | L
    ss_propensity_helix: Optional[float] = None
    ss_propensity_sheet: Optional[float] = None
    ss_propensity_loop: Optional[float] = None
    ss_stability: Optional[float] = None

    # --- Tier 1: chi1/chi2 dihedral entropy (Phase 2 / A5) ------------ #
    # Shannon entropy (bits) of the χ1 / χ2 angle distribution binned into
    # rotamer wells (60° wide; 6 bins). rotamer_diversity is the number of
    # rotamer wells visited ≥5% of the time — signals whether the sidechain
    # samples one rotamer or interchanges between several.
    chi1_entropy: Optional[float] = None
    chi2_entropy: Optional[float] = None
    rotamer_diversity: Optional[int] = None

    # --- Region / context --------------------------------------------- #
    region: Optional[str] = None                # e.g. "CDR3β", "peptide", ...
    is_interface: Optional[bool] = None
    is_core: Optional[bool] = None

    # --- Provenance --------------------------------------------------- #
    provenance: FeatureProvenance = field(default_factory=FeatureProvenance)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a plain dict (for JSON output, presenter input)."""
        out: Dict[str, Any] = {
            "residue": self.residue.label(),
            "chain": self.residue.chain,
            "resid": self.residue.resid,
            "resname": self.residue.resname,
        }
        # Copy over all set (non-None, non-empty) scalar fields.
        for name in (
            "bsa_contribution", "bsa_rank",
            "rrcs_contribution", "rrcs_rank", "contact_count",
            "contact_redundancy", "partner_diversity",
            "persistence_mean", "persistence_max",
            "persistence_profile", "design_priority_score",
            "pocket_chemistry", "partner_chemistry",
            "rmsf_mean", "rmsf_rank",
            "rmsf_backbone", "rmsf_sidechain", "rmsf_ratio",
            "sasa_bound", "sasa_unbound", "delta_sasa", "relative_exposure",
            "sasa_sidechain_bound", "burial_state",
            "secondary_structure", "ss_propensity_helix", "ss_propensity_sheet",
            "ss_propensity_loop", "ss_stability",
            "chi1_entropy", "chi2_entropy", "rotamer_diversity",
            "cdr3_tip_distance", "cdr3_loop_span", "cdr3_torsion_deg",
            "region", "is_interface", "is_core",
        ):
            val = getattr(self, name)
            if val is not None:
                out[name] = val

        # Legacy alias for downstream consumers that still read mutability_score
        if self.design_priority_score is not None:
            out["mutability_score"] = self.design_priority_score

        if self.risk_flags:
            out["risk_flags"] = list(self.risk_flags)
        if self.priority_components:
            out["priority_components"] = dict(self.priority_components)
            out["mutability_components"] = dict(self.priority_components)
        if self.persistence_components:
            out["persistence_components"] = dict(self.persistence_components)
        if self.chemistry_tags:
            out["chemistry_tags"] = list(self.chemistry_tags)

        out["provenance"] = {
            "computer": self.provenance.computer,
            "sources": list(self.provenance.sources),
            "missing_inputs": list(self.provenance.missing_inputs),
            "notes": list(self.provenance.notes),
            "is_partial": self.provenance.is_partial,
        }
        return out

    def merge(self, other: "ResidueFeatures") -> "ResidueFeatures":
        """
        Merge features from another ResidueFeatures into this one.

        Non-None values in `other` override None values in `self`. Lists
        (risk_flags, provenance.sources/notes/missing_inputs) are concatenated
        with de-duplication.

        Returns self (for chaining).
        """
        if other.residue != self.residue:
            raise FeatureError(
                f"Cannot merge features for different residues: "
                f"{self.residue.label()} vs {other.residue.label()}"
            )

        for name in (
            "bsa_contribution", "bsa_rank",
            "rrcs_contribution", "rrcs_rank", "contact_count",
            "contact_redundancy", "partner_diversity",
            "persistence_mean", "persistence_max",
            "persistence_profile", "design_priority_score",
            "pocket_chemistry", "partner_chemistry",
            "rmsf_mean", "rmsf_rank",
            "cdr3_tip_distance", "cdr3_loop_span", "cdr3_torsion_deg",
            "region", "is_interface", "is_core",
            # Interaction-family per-residue counts (Phase 1 / A1)
            "hbond_count", "saltbridge_count", "hydrophobic_count",
            "pipi_count", "cationpi_count",
            "hbond_max_occupancy", "saltbridge_max_occupancy",
            "hydrophobic_max_occupancy", "pipi_max_occupancy",
            "cationpi_max_occupancy",
            "dominant_interaction_type", "interaction_diversity",
            # Backbone/sidechain RMSF split (Phase 2 / A2)
            "rmsf_backbone", "rmsf_sidechain", "rmsf_ratio",
            # Per-residue SASA (Phase 2 / A3)
            "sasa_bound", "sasa_unbound", "delta_sasa", "relative_exposure",
            "sasa_sidechain_bound", "burial_state",
            # Secondary structure (Phase 2 / A4)
            "secondary_structure", "ss_propensity_helix", "ss_propensity_sheet",
            "ss_propensity_loop", "ss_stability",
            # Chi entropy (Phase 2 / A5)
            "chi1_entropy", "chi2_entropy", "rotamer_diversity",
        ):
            if getattr(self, name) is None:
                val = getattr(other, name)
                if val is not None:
                    setattr(self, name, val)

        # Merge risk_flags (preserve order, dedupe)
        for flag in other.risk_flags:
            if flag not in self.risk_flags:
                self.risk_flags.append(flag)

        for tag in other.chemistry_tags:
            if tag not in self.chemistry_tags:
                self.chemistry_tags.append(tag)

        for key, value in other.priority_components.items():
            self.priority_components.setdefault(key, value)

        for key, value in other.persistence_components.items():
            self.persistence_components.setdefault(key, value)

        # Merge provenance
        for src in other.provenance.sources:
            if src not in self.provenance.sources:
                self.provenance.sources.append(src)
        for note in other.provenance.notes:
            if note not in self.provenance.notes:
                self.provenance.notes.append(note)
        for missing in other.provenance.missing_inputs:
            if missing not in self.provenance.missing_inputs:
                self.provenance.missing_inputs.append(missing)

        # Prefer more specific computer name if current is empty
        if not self.provenance.computer:
            self.provenance.computer = other.provenance.computer

        return self

    # Legacy attribute aliases — let old consumers still read AND write
    # mutability_*. Read goes via __getattr__ (only triggered when the name
    # isn't a real attribute); write goes via __setattr__ which redirects
    # alias writes to the real field. Without the setter, `obj.mutability_score
    # = x` would silently create a *separate* attribute and the real
    # design_priority_score would stay None — a bug we hit in legacy
    # compute paths.
    _LEGACY_ATTRS = {
        "mutability_score": "design_priority_score",
        "mutability_components": "priority_components",
    }

    def __getattr__(self, name):
        target = type(self)._LEGACY_ATTRS.get(name) if hasattr(type(self), "_LEGACY_ATTRS") else None
        if target is not None:
            return self.__dict__[target] if target in self.__dict__ else None
        raise AttributeError(name)

    def __setattr__(self, name, value):
        # Redirect alias writes (e.g. mutability_score=...) to the canonical
        # field so old computers writing the legacy name still mutate the
        # real attribute. Use object.__setattr__ to avoid recursing through
        # this same method. dataclass-generated __init__ also goes through
        # this path — real field names just take the else-branch, so no
        # special handling is needed there.
        legacy_map = type(self).__dict__.get("_LEGACY_ATTRS")
        target = legacy_map.get(name) if legacy_map else None
        if target is not None:
            object.__setattr__(self, target, value)
        else:
            object.__setattr__(self, name, value)


# ---------------------------------------------------------------------- #
# Feature set (collection of per-residue features)
# ---------------------------------------------------------------------- #

@dataclass
class FeatureSet:
    """
    A collection of per-residue features for one case.

    Acts as a simple dict-like container plus metadata. Features can be
    built incrementally by multiple FeatureComputers and then queried
    by presenters / agent tools.
    """

    case_id: str
    residues: Dict[ResidueKey, ResidueFeatures] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    computed_at: str = field(default_factory=_now_iso)

    # ------------------------------------------------------------------ #
    # Basic accessors
    # ------------------------------------------------------------------ #

    def __len__(self) -> int:
        return len(self.residues)

    def __contains__(self, key: ResidueKey) -> bool:
        return key in self.residues

    def __iter__(self):
        return iter(self.residues.values())

    def get(self, key: ResidueKey) -> Optional[ResidueFeatures]:
        return self.residues.get(key)

    def upsert(self, features: ResidueFeatures) -> ResidueFeatures:
        """Insert a new ResidueFeatures or merge into existing one."""
        existing = self.residues.get(features.residue)
        if existing is None:
            self.residues[features.residue] = features
            return features
        existing.merge(features)
        return existing

    # ------------------------------------------------------------------ #
    # Queries
    # ------------------------------------------------------------------ #

    def top_by(
        self,
        field_name: str,
        n: int = 10,
        reverse: bool = True,
    ) -> List[ResidueFeatures]:
        """
        Return top-N residues sorted by a numeric feature.

        Residues where the field is None are excluded.
        """
        def key(f: ResidueFeatures) -> float:
            val = getattr(f, field_name, None)
            return float(val) if val is not None else float("-inf" if reverse else "inf")

        with_value = [f for f in self.residues.values() if getattr(f, field_name, None) is not None]
        return sorted(with_value, key=key, reverse=reverse)[:n]

    def filter(
        self,
        chain: Optional[str] = None,
        region: Optional[str] = None,
        interface_only: bool = False,
    ) -> List[ResidueFeatures]:
        """Simple filter over residues."""
        out: List[ResidueFeatures] = []
        for f in self.residues.values():
            if chain is not None and f.residue.chain != chain:
                continue
            if region is not None and f.region != region:
                continue
            if interface_only and f.is_interface is not True:
                continue
            out.append(f)
        return out

    # ------------------------------------------------------------------ #
    # Serialization
    # ------------------------------------------------------------------ #

    def to_dict(self) -> Dict[str, Any]:
        return {
            "case_id": self.case_id,
            "computed_at": self.computed_at,
            "metadata": dict(self.metadata),
            "residues": [rf.to_dict() for rf in self.residues.values()],
        }

    def summary(self) -> Dict[str, Any]:
        """
        Lightweight summary (residue count + which features are populated).

        Useful for debugging and for presenters that want to show what's
        available before querying detailed values.
        """
        populated: Dict[str, int] = {}
        for rf in self.residues.values():
            for name in (
                "bsa_contribution", "rrcs_contribution", "contact_redundancy",
                "partner_diversity", "persistence_mean", "design_priority_score",
                "persistence_profile", "pocket_chemistry", "rmsf_mean",
                "cdr3_tip_distance", "cdr3_loop_span", "cdr3_torsion_deg",
                "region",
            ):
                if getattr(rf, name, None) is not None:
                    populated[name] = populated.get(name, 0) + 1
            for name in (
                "risk_flags", "chemistry_tags", "priority_components",
                "persistence_components",
            ):
                if getattr(rf, name):
                    populated[name] = populated.get(name, 0) + 1

        return {
            "case_id": self.case_id,
            "residue_count": len(self.residues),
            "populated_features": populated,
            "computed_at": self.computed_at,
        }
