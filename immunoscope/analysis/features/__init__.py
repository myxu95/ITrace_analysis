"""
Feature engineering layer — the bridge between the analysis platform and
the design platform.

This layer sits between:
  - analysis/ : raw physical quantities (RRCS, BSA, RMSF, contacts, ...)
  - agent/    : LLM-driven design reasoning

A "feature" here is a design-relevant derived quantity computed from
analysis outputs. Examples:
  - bsa_contribution, bsa_rank
  - rrcs_contribution, rrcs_rank
  - contact_redundancy, partner_diversity
  - risk_flags, chemistry_tags, pocket_chemistry
  - design_priority_score (DEPRECATED per D-B5 2026-05-26; kept under the
    stable registry name "mutability_score" only for legacy back-end
    consumers — not exposed to the LLM)

Features are consumed by:
  - presenters/ : Markdown views for the Agent
  - agent/tools/: design-oriented tools (compute_design_features, ...)
  - reporter/   : structured JSON query answers

Public API:

    from immunoscope.analysis.features import (
        CaseLocator,           # shared with presenters/
        FeatureSet,            # collection of per-residue features
        ResidueKey,            # canonical residue identifier
        ResidueFeatures,       # per-residue feature record
        compute_feature,       # run one registered computer
        list_features,         # inspect registry
    )

Example:

    locator = CaseLocator(case_dir)
    fs = FeatureSet(case_id=locator.get_case_id())
    compute_feature("contact_count", locator, fs)
    compute_feature("rrcs_contribution", locator, fs)
    top = fs.top_by("rrcs_contribution", n=5)

This package exposes both lightweight reference computers and production
design features such as design priority scoring, risk flags, and chemistry tags.
"""

from .core.locator import CaseLocator, CaseLocatorError
from .core.models import (
    FeatureError,
    ResidueKey,
    ResidueFeatures,
    FeatureSet,
    FeatureProvenance,
)
from .core.registry import (
    FeatureComputer,
    FEATURE_REGISTRY,
    register_feature,
    compute_feature,
    list_features,
)

# Import examples to trigger their registration.
# Once real Tier 1 features land (PR 4), these may be moved or deprecated.
from . import examples as _examples  # noqa: F401
from . import design as _design  # noqa: F401
from . import interactions as _interactions  # noqa: F401
from . import rmsf_split as _rmsf_split  # noqa: F401
from . import sasa as _sasa  # noqa: F401
from . import secondary as _secondary  # noqa: F401
from . import chi_entropy as _chi_entropy  # noqa: F401

__all__ = [
    # Locator
    "CaseLocator",
    "CaseLocatorError",
    # Models
    "FeatureError",
    "ResidueKey",
    "ResidueFeatures",
    "FeatureSet",
    "FeatureProvenance",
    # Registry
    "FeatureComputer",
    "FEATURE_REGISTRY",
    "register_feature",
    "compute_feature",
    "list_features",
]
