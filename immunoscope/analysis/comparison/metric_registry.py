"""Comparison metric registry with unified hierarchical organization.

This module defines all comparison metrics with a two-tier priority system:
- Primary: Essential metrics shown by default (~12 metrics)
- Detailed: Complete set of metrics for in-depth analysis (~30 metrics)

The same metrics are used for all comparison types (same-system or cross-system).
Statistical tests are applied when time-series data is available.

Author: ImmunoScope Development Team
Date: 2026-05-06
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Literal


@dataclass
class MetricDefinition:
    """Definition of a single comparison metric."""

    name: str
    display_name: str
    category: str
    unit: str
    description: str
    priority: Literal["primary", "detailed"]
    show_in_overview: bool = False  # Show in hero section (top 4 cards)
    statistical_test: bool = True  # Apply statistical test if data available


class MetricRegistry:
    """Registry of all comparison metrics with hierarchical organization."""

    # Metric categories in display order
    CATEGORIES = [
        "identity",      # System identity (peptide, HLA, TCR)
        "structure",     # Structural metrics (BSA, interface, RMSD)
        "flexibility",   # Dynamics metrics (RMSF by region)
        "interaction",   # Interaction metrics (contacts, H-bonds, etc.)
        "hotspot",       # Hotspot analysis (RRCS)
        "quality",       # Simulation quality metrics
    ]

    CATEGORY_LABELS = {
        "identity": "System Identity",
        "structure": "Structure & Dynamics",
        "flexibility": "Regional Flexibility",
        "interaction": "Interactions",
        "hotspot": "Hotspot Analysis",
        "quality": "Simulation Quality",
    }

    # Define all metrics
    METRICS = [
        # ============================================================
        # IDENTITY METRICS - System identification
        # ============================================================
        MetricDefinition(
            name="peptide_sequence",
            display_name="Peptide Sequence",
            category="identity",
            unit="",
            description="Peptide amino acid sequence",
            priority="primary",
            statistical_test=False,
        ),
        MetricDefinition(
            name="peptide_length",
            display_name="Peptide Length",
            category="identity",
            unit="aa",
            description="Number of amino acids in peptide",
            priority="primary",
            statistical_test=False,
        ),
        MetricDefinition(
            name="hla_locus",
            display_name="HLA Locus",
            category="identity",
            unit="",
            description="HLA class I locus (A, B, or C)",
            priority="primary",
            statistical_test=False,
        ),
        MetricDefinition(
            name="tcr_alpha_genotype",
            display_name="TCR Alpha Genotype",
            category="identity",
            unit="",
            description="TCR alpha chain V-J gene usage",
            priority="detailed",
            statistical_test=False,
        ),
        MetricDefinition(
            name="tcr_beta_genotype",
            display_name="TCR Beta Genotype",
            category="identity",
            unit="",
            description="TCR beta chain V-J gene usage",
            priority="detailed",
            statistical_test=False,
        ),

        # ============================================================
        # STRUCTURE METRICS - Interface stability and geometry
        # ============================================================
        MetricDefinition(
            name="mean_bsa",
            display_name="Mean Buried Surface Area",
            category="structure",
            unit="Ų",
            description="Average buried surface area at TCR-pMHC interface",
            priority="primary",
            show_in_overview=True,
            statistical_test=True,
        ),
        MetricDefinition(
            name="mean_interface_ratio",
            display_name="Mean Interface Ratio",
            category="structure",
            unit="",
            description="Ratio of interface area to total surface area",
            priority="detailed",
            statistical_test=True,
        ),
        MetricDefinition(
            name="mean_rmsd",
            display_name="Mean RMSD",
            category="structure",
            unit="Å",
            description="Average RMSD from reference structure",
            priority="detailed",
            statistical_test=True,
        ),
        MetricDefinition(
            name="mean_com_distance",
            display_name="Mean COM Distance",
            category="structure",
            unit="Å",
            description="Average center-of-mass distance between TCR and pMHC",
            priority="detailed",
            statistical_test=True,
        ),

        # ============================================================
        # FLEXIBILITY METRICS - Regional dynamics
        # ============================================================
        MetricDefinition(
            name="cdr3_alpha_rmsf",
            display_name="CDR3α RMSF",
            category="flexibility",
            unit="Å",
            description="Average RMSF of CDR3 alpha loop",
            priority="primary",
            show_in_overview=True,
            statistical_test=True,
        ),
        MetricDefinition(
            name="cdr3_beta_rmsf",
            display_name="CDR3β RMSF",
            category="flexibility",
            unit="Å",
            description="Average RMSF of CDR3 beta loop",
            priority="primary",
            statistical_test=True,
        ),
        MetricDefinition(
            name="peptide_rmsf",
            display_name="Peptide RMSF",
            category="flexibility",
            unit="Å",
            description="Average RMSF of peptide",
            priority="primary",
            statistical_test=True,
        ),
        MetricDefinition(
            name="tcr_mean_rmsf",
            display_name="TCR Mean RMSF",
            category="flexibility",
            unit="Å",
            description="Average RMSF of entire TCR",
            priority="detailed",
            statistical_test=True,
        ),
        MetricDefinition(
            name="phla_mean_rmsf",
            display_name="pMHC Mean RMSF",
            category="flexibility",
            unit="Å",
            description="Average RMSF of pMHC (peptide-MHC)",
            priority="detailed",
            statistical_test=True,
        ),
        MetricDefinition(
            name="alpha1_helix_rmsf",
            display_name="MHC α1-helix RMSF",
            category="flexibility",
            unit="Å",
            description="Average RMSF of MHC alpha-1 helix",
            priority="detailed",
            statistical_test=True,
        ),
        MetricDefinition(
            name="alpha2_helix_rmsf",
            display_name="MHC α2-helix RMSF",
            category="flexibility",
            unit="Å",
            description="Average RMSF of MHC alpha-2 helix",
            priority="detailed",
            statistical_test=True,
        ),
        MetricDefinition(
            name="mean_rmsf",
            display_name="Overall Mean RMSF",
            category="flexibility",
            unit="Å",
            description="Average RMSF across entire complex",
            priority="detailed",
            statistical_test=True,
        ),

        # ============================================================
        # INTERACTION METRICS - Contact analysis
        # ============================================================
        MetricDefinition(
            name="total_contact_pairs",
            display_name="Total Contact Pairs",
            category="interaction",
            unit="pairs",
            description="Total number of residue pairs in contact",
            priority="primary",
            show_in_overview=True,
            statistical_test=True,
        ),
        MetricDefinition(
            name="hbond_pairs",
            display_name="Hydrogen Bond Pairs",
            category="interaction",
            unit="pairs",
            description="Number of hydrogen bond pairs",
            priority="primary",
            statistical_test=True,
        ),
        MetricDefinition(
            name="saltbridge_pairs",
            display_name="Salt Bridge Pairs",
            category="interaction",
            unit="pairs",
            description="Number of salt bridge pairs",
            priority="primary",
            statistical_test=True,
        ),
        MetricDefinition(
            name="hydrophobic_pairs",
            display_name="Hydrophobic Contact Pairs",
            category="interaction",
            unit="pairs",
            description="Number of hydrophobic contact pairs",
            priority="primary",
            statistical_test=True,
        ),
        MetricDefinition(
            name="pipi_pairs",
            display_name="π-π Stacking Pairs",
            category="interaction",
            unit="pairs",
            description="Number of π-π stacking pairs",
            priority="detailed",
            statistical_test=True,
        ),
        MetricDefinition(
            name="cationpi_pairs",
            display_name="Cation-π Pairs",
            category="interaction",
            unit="pairs",
            description="Number of cation-π interaction pairs",
            priority="detailed",
            statistical_test=True,
        ),

        # ============================================================
        # HOTSPOT METRICS - RRCS analysis
        # ============================================================
        MetricDefinition(
            name="peptide_tcr_rrcs",
            display_name="Peptide-TCR RRCS Sum",
            category="hotspot",
            unit="score",
            description="Sum of RRCS scores for peptide-TCR interactions",
            priority="primary",
            show_in_overview=True,
            statistical_test=True,
        ),
        MetricDefinition(
            name="hla_tcr_rrcs",
            display_name="HLA-TCR RRCS Sum",
            category="hotspot",
            unit="score",
            description="Sum of RRCS scores for HLA-TCR interactions",
            priority="primary",
            statistical_test=True,
        ),
        MetricDefinition(
            name="rrcs_nonzero_pairs",
            display_name="RRCS Active Pairs",
            category="hotspot",
            unit="pairs",
            description="Number of residue pairs with non-zero RRCS score",
            priority="detailed",
            statistical_test=True,
        ),
        MetricDefinition(
            name="rrcs_total_pairs",
            display_name="RRCS Total Pairs",
            category="hotspot",
            unit="pairs",
            description="Total number of residue pairs in RRCS analysis",
            priority="detailed",
            statistical_test=False,
        ),

        # ============================================================
        # QUALITY METRICS - Simulation quality and convergence
        # ============================================================
        MetricDefinition(
            name="time_span",
            display_name="Simulation Time",
            category="quality",
            unit="ns",
            description="Total simulation time span",
            priority="primary",
            statistical_test=False,
        ),
        MetricDefinition(
            name="n_frames",
            display_name="Number of Frames",
            category="quality",
            unit="frames",
            description="Total number of trajectory frames",
            priority="detailed",
            statistical_test=False,
        ),
        MetricDefinition(
            name="tail90_rmsd_variation",
            display_name="Tail 90% RMSD Variation",
            category="quality",
            unit="nm",
            description="RMSD variation in last 90% of trajectory (convergence indicator)",
            priority="detailed",
            statistical_test=True,
        ),
    ]

    @classmethod
    def get_metrics(
        cls,
        priority_filter: tuple[str, ...] = ("primary", "detailed"),
    ) -> list[MetricDefinition]:
        """Get metrics filtered by priority.

        Args:
            priority_filter: Which priority levels to include

        Returns:
            List of metrics in display order
        """
        return [
            metric for metric in cls.METRICS
            if metric.priority in priority_filter
        ]

    @classmethod
    def get_primary_metrics(cls) -> list[MetricDefinition]:
        """Get only primary metrics for standard report."""
        return cls.get_metrics(priority_filter=("primary",))

    @classmethod
    def get_overview_metrics(cls) -> list[MetricDefinition]:
        """Get metrics for hero section (top 4 cards)."""
        return [m for m in cls.METRICS if m.show_in_overview]

    @classmethod
    def get_metrics_by_category(
        cls,
        priority_filter: tuple[str, ...] = ("primary", "detailed"),
    ) -> dict[str, list[MetricDefinition]]:
        """Get metrics grouped by category.

        Args:
            priority_filter: Which priority levels to include

        Returns:
            Dictionary mapping category name to list of metrics
        """
        metrics = cls.get_metrics(priority_filter)
        result = {cat: [] for cat in cls.CATEGORIES}

        for metric in metrics:
            if metric.category in result:
                result[metric.category].append(metric)

        # Remove empty categories
        return {k: v for k, v in result.items() if v}

    @classmethod
    def get_metric(cls, name: str) -> MetricDefinition | None:
        """Get a specific metric by name."""
        for metric in cls.METRICS:
            if metric.name == name:
                return metric
        return None
