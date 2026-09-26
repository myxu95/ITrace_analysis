"""Question classification and result routing for reporter queries."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List


@dataclass
class QueryRoute:
    """Routing result for a report question."""

    query_type: str
    sources: List[Path]
    preferred_missing: List[str]


class QueryRouter:
    """Route a question to the minimum necessary result sources."""

    def __init__(self, case_dir: str | Path):
        self.case_dir = Path(case_dir).resolve()

    def classify(self, question: str) -> str:
        text = question.strip().lower()

        if any(token in text for token in ["what can you do", "capabilities", "available functions", "help"]):
            return "capability_overview"

        # Full-stage skills are intentionally routed only by explicit full-task wording.
        if any(token in text for token in ["full diagnostic", "system diagnostic", "complete analysis"]):
            return "system_diagnostic"

        if any(token in text for token in ["full design", "mutation design report", "design recommendation report"]):
            return "mutation_design"

        # Two-system comparison
        if any(token in text for token in ["compare", "comparison", "versus", "vs", "difference"]):
            return "system_compare"

        # Existing query types
        if any(token in text for token in ["cluster", "dominant", "state"]):
            return "cluster_status"
        if any(token in text for token in ["stability", "stable", "quality", "rmsd", "trajectory"]):
            return "quality_status"
        if any(token in text for token in ["fel", "free energy", "landscape", "energy basin", "pca landscape"]):
            return "landscape_status"
        if any(token in text for token in ["summarize interface", "interface summary", "interface"]):
            return "interface_summary"
        if any(token in text for token in ["mutation", "mutate"]):
            return "design_hint"
        if any(token in text for token in ["identify hotspots", "hotspots", "hotspot", "key site"]):
            return "hotspot_summary"
        if any(token in text for token in ["rrcs", "hotspot pair", "top pair", "strongest hotspot"]):
            return "top_pair"
        if any(token in text for token in ["persistent", "occupancy", "stable pair"]):
            return "persistent_pair"
        if any(token in text for token in ["region", "active region", "dominant region"]):
            return "top_region"
        if any(token in text for token in ["residue", "site", "key"]):
            return "top_residue"
        return "top_residue"

    def route(self, question: str) -> QueryRoute:
        query_type = self.classify(question)

        # New skill types
        if query_type == "system_diagnostic":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/quality_summary.json",
                    "overview/bsa/interface_summary.json",
                    "overview/interaction_overview.json",
                    "overview/rrcs/rrcs_summary.json",
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                    "overview/rrcs/rrcs_region_summary.csv",
                    "overview/rmsf/rmsf_summary.json",
                    "overview/cluster/interface_clustering_summary.json",
                    "overview/cluster/summary_table.csv",
                ]),
                preferred_missing=["overview directory with analysis summaries"],
            )

        if query_type == "mutation_design":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                    "overview/rrcs/rrcs_summary.json",
                    "overview/rrcs/rrcs_region_summary.csv",
                    "overview/interaction_overview.json",
                    "overview/bsa/interface_summary.json",
                ]),
                preferred_missing=["overview/rrcs/annotated_rrcs_pair_summary.csv"],
            )

        if query_type == "system_compare":
            # For comparison, we need two case directories
            # This will be handled specially in the handler
            return QueryRoute(
                query_type=query_type,
                sources=[],  # Will be populated by handler
                preferred_missing=["Two case directories for comparison"],
            )

        # Existing query types
        if query_type == "capability_overview":
            return QueryRoute(
                query_type=query_type,
                sources=[],
                preferred_missing=[],
            )

        if query_type == "cluster_status":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/cluster/interface_clustering_summary.json",
                    "overview/cluster/summary_table.csv",
                    "overview/cluster/cluster_feature_digest.csv",
                ]),
                preferred_missing=["overview/cluster/interface_clustering_summary.json", "overview/cluster/summary_table.csv"],
            )
        if query_type == "quality_status":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing(self._glob_candidates(["*quality*summary*.json", "*rmsd*summary*.json", "*tail*summary*.json"]))
                + self._existing([
                    "overview/rmsf/rmsf_summary.json",
                ]),
                preferred_missing=["quality summary or RMSD summary"],
            )
        if query_type == "landscape_status":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing(self._glob_candidates([
                    "*landscape_summary.json",
                    "*loadings.csv",
                    "*pca_coordinates.csv",
                ])),
                preferred_missing=["landscape_summary.json", "loadings.csv"],
            )
        if query_type == "interface_summary":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/bsa/interface_summary.json",
                    "overview/interaction_overview.json",
                    "overview/rrcs/rrcs_region_summary.csv",
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                ]),
                preferred_missing=["overview/bsa/interface_summary.json", "overview/interaction_overview.json"],
            )
        if query_type == "hotspot_summary":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                    "overview/rrcs/rrcs_summary.json",
                    "overview/rrcs/rrcs_region_summary.csv",
                    "overview/cluster/interface_clustering_summary.json",
                    "overview/cluster/summary_table.csv",
                ]),
                preferred_missing=["overview/rrcs/annotated_rrcs_pair_summary.csv", "overview/rrcs/rrcs_summary.json"],
            )
        if query_type == "design_hint":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                    "overview/rrcs/rrcs_summary.json",
                    "overview/rrcs/rrcs_region_summary.csv",
                    "overview/bsa/interface_summary.json",
                    "overview/interaction_overview.json",
                    "overview/cluster/interface_clustering_summary.json",
                    "overview/cluster/summary_table.csv",
                ]),
                preferred_missing=["overview/rrcs/annotated_rrcs_pair_summary.csv", "overview/rrcs/rrcs_summary.json"],
            )
        if query_type == "top_pair":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                    "overview/rrcs/rrcs_summary.json",
                    "overview/rrcs/rrcs_region_summary.csv",
                ]),
                preferred_missing=["overview/rrcs/annotated_rrcs_pair_summary.csv", "overview/rrcs/rrcs_summary.json"],
            )
        if query_type == "persistent_pair":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing(self._glob_candidates(["*occupancy*summary*.json", "*persistence*summary*.json", "*occupancy*region*.csv", "*persistent*pair*.csv"]))
                + self._existing([
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                    "overview/rrcs/rrcs_summary.json",
                ]),
                preferred_missing=["occupancy or persistence summary"],
            )
        if query_type == "top_region":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/rrcs/rrcs_region_summary.csv",
                    "overview/interaction_overview.json",
                    "overview/cluster/cluster_feature_digest.csv",
                ]),
                preferred_missing=["overview/rrcs/rrcs_region_summary.csv", "overview/interaction_overview.json"],
            )
        return QueryRoute(
            query_type="top_residue",
            sources=self._existing([
                "overview/rrcs/annotated_rrcs_pair_summary.csv",
                "overview/rrcs/rrcs_summary.json",
                "overview/rrcs/rrcs_region_summary.csv",
                "overview/interaction_overview.json",
            ]),
            preferred_missing=["RRCS annotated pair summary or RRCS summary"],
        )

    def _existing(self, candidates: list[str | Path]) -> list[Path]:
        paths: list[Path] = []
        seen: set[Path] = set()
        for candidate in candidates:
            path = candidate if isinstance(candidate, Path) else self.case_dir / candidate
            if path.exists() and path not in seen:
                paths.append(path)
                seen.add(path)
        return paths

    def _glob_candidates(self, patterns: list[str]) -> list[Path]:
        matched: list[Path] = []
        for pattern in patterns:
            matched.extend(sorted(self.case_dir.rglob(pattern)))
        return matched
