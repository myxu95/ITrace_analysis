"""Compatibility alias for the legacy flat node import path."""

from importlib import import_module
import sys

sys.modules[__name__] = import_module(".interactions.contact_annotation_node", package=__package__)
