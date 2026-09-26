"""Core types for the features layer."""

from .locator import CaseLocator, CaseLocatorError
from .models import (
    FeatureError,
    ResidueKey,
    ResidueFeatures,
    FeatureSet,
    FeatureProvenance,
)
from .registry import (
    FeatureComputer,
    FEATURE_REGISTRY,
    register_feature,
    compute_feature,
    list_features,
)

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
