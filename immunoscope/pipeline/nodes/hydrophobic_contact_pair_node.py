"""Compatibility alias for the legacy flat node import path."""

from importlib import import_module
import sys

sys.modules[__name__] = import_module(".interactions.hydrophobic_contact_pair_node", package=__package__)
