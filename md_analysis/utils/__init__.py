from .plotting import PlotManager
from .path_manager import PathManager
from .group_selector import GroupSelector
from .pdb_downloader import PDBDownloader
from .multimodel_concatenator import MultiModelConcatenator
from .selection_string_builder import SelectionStringBuilder


def __getattr__(name):
    """Lazily expose legacy structure/topology helpers without import cycles."""
    if name in {"PDBSequenceExtractor", "PDBChainStandardizer"}:
        from md_analysis.analysis.structure import PDBChainStandardizer, PDBSequenceExtractor
        return {
            "PDBSequenceExtractor": PDBSequenceExtractor,
            "PDBChainStandardizer": PDBChainStandardizer,
        }[name]
    if name == "IntelligentChainIdentifier":
        from md_analysis.analysis.topology import IntelligentChainIdentifier
        return IntelligentChainIdentifier
    raise AttributeError(name)

__all__ = [
    "PlotManager",
    "PathManager",
    "GroupSelector",
    "PDBDownloader",
    "MultiModelConcatenator",
    "SelectionStringBuilder",
    "PDBSequenceExtractor",
    "IntelligentChainIdentifier",
    "PDBChainStandardizer",
]
