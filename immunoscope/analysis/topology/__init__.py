"""Topology-aware analysis utilities and index generation contracts."""

from .index_generation import (
    IndexGenerator,
    IndexGenerationInput,
    IndexGenerationResult,
    IndexGenerationMethod,
    ComponentDefinition,
    ComponentIndexInfo,
    IndexGenerationError,
    InvalidInputError,
    TopologyFileError,
    GROMACSCommandError,
    ComponentNotFoundError,
    SideEffectTracker,
)
from .component_index_generator import ComponentIndexGenerator
from .intelligent_chain_identifier import IntelligentChainIdentifier, ChainIdentification
from .tpr_chain_extractor import TPRChainSequenceExtractor, UnifiedChainInfo
from .chain_identification_adapter import ChainIdentificationAdapter
from .topology_chain_identifier import TopologyChainIdentifier, TopologyChainInfo, TopologyChainMapping
from .shortest_chain_detector import ShortestChainDetector, create_shortest_chain_index
from .strict_shortest_chain_detector import StrictShortestChainDetector
from .index_manager import IndexManager, StandardComponent, GroupInfo
from .chain_based_index_generator import (
    ChainBasedIndexGenerator,
    generate_peptide_index_from_pdb,
    generate_peptide_index_from_tpr,
)
from .cdr_manager import CDRManager, ANARCIWrapper, CDRIndexGenerator, CDRMetadataManager
from .tcr_residue_semantics import TCRResidueSemantics, TCRResidueSemanticAnnotator
from .complex_residue_semantics import ComplexResidueSemantics, ComplexResidueSemanticAnnotator
from .residue_pair_annotation import ResiduePairAnnotationAnnotator, ContactAnnotationAnnotator
from .region_interaction_summary import RegionInteractionSummaryBuilder
from .contact_heatmap import ContactHeatmapArtifacts, ContactHeatmapPlotter
from .biological_identity import (
    HLAIdentityAnnotator,
    BiologicalIdentityAnnotator,
    derive_hla_reference_library,
    ensure_hla_reference_library,
)

__all__ = [
    "IndexGenerator",
    "IndexGenerationInput",
    "IndexGenerationResult",
    "IndexGenerationMethod",
    "ComponentDefinition",
    "ComponentIndexInfo",
    "IndexGenerationError",
    "InvalidInputError",
    "TopologyFileError",
    "GROMACSCommandError",
    "ComponentNotFoundError",
    "SideEffectTracker",
    "ComponentIndexGenerator",
    "IntelligentChainIdentifier",
    "ChainIdentification",
    "TPRChainSequenceExtractor",
    "UnifiedChainInfo",
    "ChainIdentificationAdapter",
    "TopologyChainIdentifier",
    "TopologyChainInfo",
    "TopologyChainMapping",
    "ShortestChainDetector",
    "create_shortest_chain_index",
    "StrictShortestChainDetector",
    "IndexManager",
    "StandardComponent",
    "GroupInfo",
    "ChainBasedIndexGenerator",
    "generate_peptide_index_from_pdb",
    "generate_peptide_index_from_tpr",
    "CDRManager",
    "ANARCIWrapper",
    "CDRIndexGenerator",
    "CDRMetadataManager",
    "TCRResidueSemantics",
    "TCRResidueSemanticAnnotator",
    "ComplexResidueSemantics",
    "ComplexResidueSemanticAnnotator",
    "ResiduePairAnnotationAnnotator",
    "ContactAnnotationAnnotator",
    "RegionInteractionSummaryBuilder",
    "ContactHeatmapArtifacts",
    "ContactHeatmapPlotter",
    "HLAIdentityAnnotator",
    "BiologicalIdentityAnnotator",
    "derive_hla_reference_library",
    "ensure_hla_reference_library",
]
