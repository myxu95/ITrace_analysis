"""Compatibility exports: contact annotation now uses the generic pair annotator."""

from .residue_pair_annotation import ContactAnnotationAnnotator, ResiduePairAnnotationAnnotator

__all__ = ["ResiduePairAnnotationAnnotator", "ContactAnnotationAnnotator"]
