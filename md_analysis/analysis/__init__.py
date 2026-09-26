# Lazy import submodules to avoid MDAnalysis dependency
from . import quality
from . import interactions
from . import interface
from . import conformation
from . import comparison
from . import geometry
from . import landscape

# D-B1 (2026-05-26): canonical feature → (layer, sub_flavor) registry.
# Importing here so callers can `from md_analysis.analysis import
# feature_registry` without paying MDAnalysis import costs.
from . import feature_registry

# Import individual modules for convenience
from .trajectory import (
    PBCProcessor,
    RMSDCalculator,
    RMSDInput,
    RMSDResult,
    RMSDConvergenceAnalyzer,
    RMSFAnalyzer,
    ResidueRMSFAnalyzer,
    ResidueRMSFResult,
    ResidueContactFrequencyAnalyzer,
    extract_sequence_from_topology,
    find_subsequence_position,
    calculate_cdr3_rmsf,
    analyze_phla_tcr_rmsf,
    write_tcr_rmsf_profile,
    write_phla_rmsf_profile,
    write_region_rmsf_summary,
)

from .structure import (
    BFactorAnalyzer,
    ContactMapCalculator,
    GeometryAnalyzer,
    AtomInfoExtractor,
    PDBChainStandardizer,
    ChainInfo,
    StandardizationResult,
    PDBSequenceExtractor,
    extract_sequences_from_pdb,
    PDBStructureFixer,
    PDBDistanceTrimmer,
)

# Import quality analysis (no MDAnalysis dependency)
from .quality import (
    MDCompletenessChecker,
    StructureValidator,
    BatchTracker,
    QualityReporter
)

from .conformation import (
    InterfaceClusteringAnalyzer,
    InterfaceClusteringResult,
    write_cluster_id_vs_time,
    write_cluster_population_over_time,
)
from .comparison import (
    ComparisonArtifacts,
    SingleCaseArtifacts,
    SingleCaseLoader,
    ComparisonBuildResult,
    SystemComparisonBuilder,
    write_quality_interface_comparison_plot,
    write_flexibility_comparison_plot,
    write_interaction_family_comparison_plot,
)

# Import free energy analysis
from .free_energy import (
    FELCalculator,
    FELVisualizer
)

# Import topology-aware analysis contracts
from .topology import (
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
    IntelligentChainIdentifier,
    ChainIdentification,
    TPRChainSequenceExtractor,
    UnifiedChainInfo,
    ChainIdentificationAdapter,
    TopologyChainIdentifier,
    TopologyChainInfo,
    TopologyChainMapping,
    ShortestChainDetector,
    create_shortest_chain_index,
    StrictShortestChainDetector,
    IndexManager,
    StandardComponent,
    GroupInfo,
    ChainBasedIndexGenerator,
    generate_peptide_index_from_pdb,
    generate_peptide_index_from_tpr,
    CDRManager,
    ANARCIWrapper,
    CDRIndexGenerator,
    CDRMetadataManager,
    TCRResidueSemantics,
    TCRResidueSemanticAnnotator,
    ComplexResidueSemantics,
    ComplexResidueSemanticAnnotator,
    ContactAnnotationAnnotator,
    ContactHeatmapArtifacts,
    ContactHeatmapPlotter,
    HLAIdentityAnnotator,
    BiologicalIdentityAnnotator,
    derive_hla_reference_library,
    ensure_hla_reference_library,
)
from .interactions import (
    HydrogenBondPairAnalyzer,
    SaltBridgePairAnalyzer,
    HydrophobicContactPairAnalyzer,
    PiStackingPairAnalyzer,
    CationPiPairAnalyzer,
    RRCSAnalyzer,
    RRCSPairSpec,
    build_residue_heavy_atom_indices,
    parse_rrcs_pair_file,
    residue_chain_id,
)
from .interface import (
    BuriedSurfaceAreaAnalyzer,
    BuriedSurfaceAreaResult,
)

__all__ = [
    # Submodules
    "trajectory",
    "structure",
    "quality",
    "interactions",
    "interface",
    "conformation",
    "comparison",
    "free_energy",
    "geometry",
    "landscape",
    "topology",
    # Trajectory analysis modules (lazy loaded)
    "RMSDCalculator",
    "RMSDInput",
    "RMSDResult",
    "RMSDConvergenceAnalyzer",
    "PBCProcessor",
    "RMSFAnalyzer",
    "ResidueRMSFAnalyzer",
    "ResidueRMSFResult",
    "ResidueContactFrequencyAnalyzer",
    "extract_sequence_from_topology",
    "find_subsequence_position",
    "calculate_cdr3_rmsf",
    "analyze_phla_tcr_rmsf",
    "write_tcr_rmsf_profile",
    "write_phla_rmsf_profile",
    "write_region_rmsf_summary",
    # Structure analysis modules
    "BFactorAnalyzer",
    "ContactMapCalculator",
    "GeometryAnalyzer",
    "AtomInfoExtractor",
    "PDBChainStandardizer",
    "ChainInfo",
    "StandardizationResult",
    "PDBSequenceExtractor",
    "extract_sequences_from_pdb",
    "PDBStructureFixer",
    "PDBDistanceTrimmer",
    # Quality analysis modules
    "MDCompletenessChecker",
    "StructureValidator",
    "BatchTracker",
    "QualityReporter",
    "InterfaceClusteringAnalyzer",
    "InterfaceClusteringResult",
    "write_cluster_id_vs_time",
    "write_cluster_population_over_time",
    "ComparisonArtifacts",
    "SingleCaseArtifacts",
    "SingleCaseLoader",
    "ComparisonBuildResult",
    "SystemComparisonBuilder",
    "write_quality_interface_comparison_plot",
    "write_flexibility_comparison_plot",
    "write_interaction_family_comparison_plot",
    # Free energy analysis modules
    "FELCalculator",
    "FELVisualizer",
    # Topology/index manager
    "IndexManager",
    "StandardComponent",
    "GroupInfo",
    # Topology analysis modules
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
    "ContactAnnotationAnnotator",
    "ContactHeatmapArtifacts",
    "ContactHeatmapPlotter",
    "HLAIdentityAnnotator",
    "BiologicalIdentityAnnotator",
    "derive_hla_reference_library",
    "ensure_hla_reference_library",
    "HydrogenBondPairAnalyzer",
    "SaltBridgePairAnalyzer",
    "HydrophobicContactPairAnalyzer",
    "PiStackingPairAnalyzer",
    "CationPiPairAnalyzer",
    "RRCSAnalyzer",
    "RRCSPairSpec",
    "build_residue_heavy_atom_indices",
    "parse_rrcs_pair_file",
    "residue_chain_id",
    "BuriedSurfaceAreaAnalyzer",
    "BuriedSurfaceAreaResult",
    "COMDistanceCalculator",
    "COMDistanceInput",
    "COMDistanceResult",
    "FeatureMatrixBuilder",
    "FeatureMatrixInput",
    "FeatureMatrixResult",
    "LandscapeAnalyzer",
    "LandscapeInput",
    "LandscapeResult",
    "LandscapeVisualizer",
]

def __getattr__(name):
    """Lazy import for trajectory analysis classes that require MDAnalysis."""
    if name == "trajectory":
        from . import trajectory
        return trajectory
    elif name in [
        "RMSDCalculator",
        "RMSDInput",
        "RMSDResult",
        "RMSDConvergenceAnalyzer",
        "PBCProcessor",
        "RMSFAnalyzer",
        "ResidueRMSFAnalyzer",
        "ResidueRMSFResult",
        "ResidueContactFrequencyAnalyzer",
        "extract_sequence_from_topology",
        "find_subsequence_position",
        "calculate_cdr3_rmsf",
        "analyze_phla_tcr_rmsf",
        "write_tcr_rmsf_profile",
        "write_phla_rmsf_profile",
        "write_region_rmsf_summary",
    ]:
        from . import trajectory
        return getattr(trajectory, name)
    elif name == "structure":
        from . import structure
        return structure
    elif name == "conformation":
        from . import conformation
        return conformation
    elif name in [
        "InterfaceClusteringAnalyzer",
        "InterfaceClusteringResult",
        "write_cluster_id_vs_time",
        "write_cluster_population_over_time",
    ]:
        from . import conformation
        return getattr(conformation, name)
    elif name == "free_energy":
        from . import free_energy
        return free_energy
    elif name in ["FELCalculator", "FELVisualizer"]:
        from . import free_energy
        return getattr(free_energy, name)
    elif name == "geometry":
        from . import geometry
        return geometry
    elif name in ["COMDistanceCalculator", "COMDistanceInput", "COMDistanceResult"]:
        from . import geometry
        return getattr(geometry, name)
    elif name == "landscape":
        from . import landscape
        return landscape
    elif name in [
        "FeatureMatrixBuilder",
        "FeatureMatrixInput",
        "FeatureMatrixResult",
        "LandscapeAnalyzer",
        "LandscapeInput",
        "LandscapeResult",
        "LandscapeVisualizer",
    ]:
        from . import landscape
        return getattr(landscape, name)
    elif name == "topology":
        from . import topology
        return topology
    elif name == "interactions":
        from . import interactions
        return interactions
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
