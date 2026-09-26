"""
Case directory locator.

This module is a thin shim that re-exports the canonical implementation from
immunoscope.analysis.features.core.locator, so that the features layer and the
presenters layer share a single CaseLocator.

Existing imports of the form:

    from immunoscope.agent.presenters import CaseLocator
    from immunoscope.agent.presenters.locator import CaseLocator

continue to work.
"""

from immunoscope.analysis.features.core.locator import (
    CaseLocator,
    CaseLocatorError,
)

# Backward-compatible alias: old code imported PresenterError from .base,
# but if anyone catches CaseLocator errors directly they should use
# CaseLocatorError. We re-export both to be safe.
__all__ = ["CaseLocator", "CaseLocatorError"]
