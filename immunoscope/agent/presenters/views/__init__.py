"""View presenters package.

Per D-B1 (locked 2026-05-26) each presenter declares its (spatial_layer,
sub_flavors_served) on the class so that `query_analysis_results` can filter
or list views by spatial layer; see `feature_registry.py` in
`immunoscope/analysis/` for the canonical (metric → layer × sub_flavor) map.

The legacy `candidates` view (composite design_priority_score) was removed
per D-B5 (locked 2026-05-26) — Pure-RRCS-rank is the new mono-signal
heuristic baseline.
"""

# Import views to trigger registration
from .overview import OverviewPresenter
from .hotspots import HotspotsPresenter
from .interface import InterfacePresenter
from .flexibility import FlexibilityPresenter
from .quality import QualityPresenter
from .clustering import ClusteringPresenter
from .residue import ResiduePresenter
from .pair import PairPresenter
from .fingerprint import FingerprintPresenter
from .interface_comparison import InterfaceComparisonPresenter
from .angles import DockingAnglesPresenter
from .dihedrals import DihedralsPresenter
from .exposure import ExposurePresenter
from .conservation import ConservationPresenter

__all__ = [
    "OverviewPresenter",
    "HotspotsPresenter",
    "InterfacePresenter",
    "FlexibilityPresenter",
    "QualityPresenter",
    "ClusteringPresenter",
    "ResiduePresenter",
    "PairPresenter",
    "FingerprintPresenter",
    "InterfaceComparisonPresenter",
    "DockingAnglesPresenter",
    "DihedralsPresenter",
    "ExposurePresenter",
    "ConservationPresenter",
]
