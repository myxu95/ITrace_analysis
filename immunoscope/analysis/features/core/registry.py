"""
Feature computer registry.

A FeatureComputer reads raw analysis outputs (via CaseLocator) and populates
fields on ResidueFeatures objects within a FeatureSet.

Computers are registered by name and composed at call time. A typical flow:

    locator = CaseLocator(case_dir)
    fs = FeatureSet(case_id=locator.get_case_id())
    compute_feature("bsa_contribution", locator, fs)
    compute_feature("rrcs_contribution", locator, fs)
    compute_feature("mutability_score", locator, fs)  # depends on above

Dependencies between computers are expressed by ordering at the call site —
the registry does not currently auto-resolve a DAG. This keeps the layer
simple; if dependency graphs become painful to manage we can add resolution
later.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, List, Optional

from .locator import CaseLocator
from .models import FeatureError, FeatureSet


class FeatureComputer(ABC):
    """
    Base class for all feature computers.

    A computer reads analysis outputs and fills in one or more fields on
    the ResidueFeatures objects inside a FeatureSet. Computers should be
    idempotent: running twice on the same FeatureSet must not double-count.

    Subclasses declare:
      - name: unique identifier used in the registry
      - required_modules: list of module names whose output is needed
          (used for availability checks and clear error messages)
      - produces: list of ResidueFeatures field names this computer sets
    """

    name: str = ""
    required_modules: List[str] = []
    produces: List[str] = []

    @abstractmethod
    def compute(self, locator: CaseLocator, features: FeatureSet) -> FeatureSet:
        """
        Run the computation.

        Args:
            locator: Resolves analysis output paths for the case.
            features: FeatureSet to populate (modified in place and returned).

        Returns:
            The same FeatureSet (for chaining).

        Raises:
            FeatureError: Only for truly fatal errors. Missing-data cases
                should be handled by recording missing_inputs on the
                affected residues' provenance, not by raising.
        """

    # ------------------------------------------------------------------ #
    # Availability
    # ------------------------------------------------------------------ #

    def is_available(self, locator: CaseLocator) -> bool:
        """
        Whether this computer's required modules are present.

        Override if a computer has more nuanced availability rules
        (e.g. needs one of several modules, not all).
        """
        for module in self.required_modules:
            if locator.get_module_root(module) is None:
                return False
        return True

    def missing_modules(self, locator: CaseLocator) -> List[str]:
        """Return the subset of required_modules that are not available."""
        return [m for m in self.required_modules if locator.get_module_root(m) is None]


# ---------------------------------------------------------------------- #
# Registry
# ---------------------------------------------------------------------- #

FEATURE_REGISTRY: Dict[str, FeatureComputer] = {}


def register_feature(computer_cls: type[FeatureComputer]) -> type[FeatureComputer]:
    """
    Decorator to register a FeatureComputer in the global registry.

    Usage:

        @register_feature
        class BsaContributionComputer(FeatureComputer):
            name = "bsa_contribution"
            required_modules = ["bsa"]
            produces = ["bsa_contribution", "bsa_rank"]

            def compute(self, locator, features):
                ...

    Raises:
        FeatureError: If the class does not set `name` or if the name is
            already registered.
    """
    if not computer_cls.name:
        raise FeatureError(
            f"FeatureComputer {computer_cls.__name__} must define a non-empty 'name'"
        )

    if computer_cls.name in FEATURE_REGISTRY:
        raise FeatureError(
            f"Feature computer '{computer_cls.name}' is already registered"
        )

    FEATURE_REGISTRY[computer_cls.name] = computer_cls()
    return computer_cls


def compute_feature(
    name: str,
    locator: CaseLocator,
    features: FeatureSet,
    *,
    raise_on_missing_modules: bool = False,
) -> FeatureSet:
    """
    Run a registered feature computer by name.

    Args:
        name: Registered computer name.
        locator: Case locator.
        features: FeatureSet to populate.
        raise_on_missing_modules: If True, raise FeatureError when required
            modules are not available. Default False — the computer itself
            is responsible for recording missing-input provenance.

    Returns:
        The updated FeatureSet.

    Raises:
        FeatureError: If the computer is not registered (always), or if
            required modules are missing and raise_on_missing_modules=True.
    """
    computer = FEATURE_REGISTRY.get(name)
    if computer is None:
        available = ", ".join(sorted(FEATURE_REGISTRY.keys())) or "(none)"
        raise FeatureError(
            f"Unknown feature computer: {name!r}. Registered: {available}"
        )

    if raise_on_missing_modules and not computer.is_available(locator):
        missing = computer.missing_modules(locator)
        raise FeatureError(
            f"Feature '{name}' requires modules not found in case: {missing}"
        )

    return computer.compute(locator, features)


def list_features(locator: Optional[CaseLocator] = None) -> List[Dict]:
    """
    List all registered feature computers.

    Args:
        locator: If provided, annotates each entry with availability
            information for that specific case.

    Returns:
        List of {name, required_modules, produces, available?, missing?}.
    """
    out: List[Dict] = []
    for name, computer in sorted(FEATURE_REGISTRY.items()):
        entry: Dict = {
            "name": name,
            "required_modules": list(computer.required_modules),
            "produces": list(computer.produces),
        }
        if locator is not None:
            entry["available"] = computer.is_available(locator)
            entry["missing_modules"] = computer.missing_modules(locator)
        out.append(entry)
    return out
