"""Answer builder for reporter query skills."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any

import pandas as pd

from .query_router import QueryRouter

try:
    from .llm_adapter import create_llm_adapter, LLMAdapter
    LLM_AVAILABLE = True
except ImportError:
    LLM_AVAILABLE = False
    LLMAdapter = None


@dataclass
class QueryAnswer:
    """Structured query answer."""

    query_type: str
    answer: str
    evidence: list[str]
    sources_used: list[str]
    top_items: list[dict[str, Any]]
    confidence: str
    skill: str = "reporter-query"
    response_panel: dict[str, Any] = field(default_factory=dict)
    evidence_panel: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class QueryAnswerBuilder:
    """Build lightweight answers from report results."""

    def __init__(
        self,
        case_dir: str | Path,
        use_llm: bool = False,
        llm_provider: str = "openai",
        llm_model: str | None = None,
        llm_api_key: str | None = None,
        llm_base_url: str | None = None,
        llm_temperature: float = 0.3,
        llm_max_tokens: int = 800,
    ):
        self.case_dir = Path(case_dir).resolve()
        self.router = QueryRouter(self.case_dir)
        self.use_llm = use_llm and LLM_AVAILABLE
        self.llm_adapter: LLMAdapter | None = None

        if self.use_llm:
            kwargs = {}
            if llm_model:
                kwargs["model"] = llm_model
            if llm_api_key:
                kwargs["api_key"] = llm_api_key
            if llm_base_url:
                kwargs["base_url"] = llm_base_url
            kwargs["temperature"] = llm_temperature
            kwargs["max_tokens"] = llm_max_tokens
            try:
                self.llm_adapter = create_llm_adapter(llm_provider, **kwargs)
            except Exception as e:
                print(f"Warning: Failed to initialize LLM adapter: {e}")
                self.use_llm = False

    def answer(self, question: str) -> QueryAnswer:
        route = self.router.route(question)
        handler = getattr(self, f"_handle_{route.query_type}", self._handle_top_residue)

        # Try LLM-enhanced answer first if enabled
        if self.use_llm and self.llm_adapter:
            try:
                return self._answer_with_llm(question, route.query_type, route.sources)
            except Exception as e:
                print(f"Warning: LLM answer failed, falling back to rule-based: {e}")

        # Fallback to rule-based answer
        return handler(route.sources, route.preferred_missing)

    def _answer_with_llm(
        self,
        question: str,
        query_type: str,
        sources: list[Path],
    ) -> QueryAnswer:
        """Generate answer using LLM with data context."""
        # Collect relevant data
        context = self._build_context_for_llm(query_type, sources)

        # Generate answer with LLM
        llm_response = self.llm_adapter.generate_answer(question, context, query_type)

        # Convert to QueryAnswer format
        return QueryAnswer(
            query_type=query_type,
            answer=llm_response.answer,
            evidence=llm_response.reasoning,
            sources_used=[str(s) for s in sources],
            top_items=[],
            confidence=llm_response.confidence,
            skill="reporter-query-llm",
        )

    def _build_context_for_llm(
        self,
        query_type: str,
        sources: list[Path],
    ) -> dict[str, Any]:
        """Build data context for LLM from available sources."""
        context = {"query_type": query_type, "data": {}}

        for source in sources:
            if not source.exists():
                continue

            try:
                if source.suffix == ".json":
                    data = json.loads(source.read_text(encoding="utf-8"))
                    context["data"][source.name] = data
                elif source.suffix == ".csv":
                    df = pd.read_csv(source)
                    # Convert to dict, limit to top 10 rows to avoid token overflow
                    context["data"][source.name] = df.head(10).to_dict(orient="records")
            except Exception as e:
                print(f"Warning: Failed to read {source}: {e}")

        return context

    def _handle_cluster_status(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        summary = self._read_json(self._find_name(sources, "interface_clustering_summary.json"))
        table = self._read_csv(self._find_name(sources, "summary_table.csv"))
        if not summary and table.empty:
            return self._missing("cluster_status", preferred_missing)

        dominant_cluster = None
        population = None
        descriptor = None
        top_items: list[dict[str, Any]] = []
        if not table.empty:
            ranked = table.sort_values("population_percent", ascending=False)
            row = ranked.iloc[0]
            dominant_cluster = int(row["cluster_id"])
            population = float(row["population_percent"])
            descriptor = row.get("main_structural_descriptor")
            top_items = ranked.head(3).to_dict(orient="records")
        elif summary:
            largest = summary.get("largest_cluster") or {}
            dominant_cluster = largest.get("cluster_id")
            population = largest.get("population_percent")
            descriptor = largest.get("main_structural_descriptor")

        answer = (
            f"The dominant cluster is cluster {dominant_cluster}." if population is None else
            f"The dominant cluster is cluster {dominant_cluster}, covering about {population:.1f}% of sampled frames."
        )
        evidence = [
            "Supported by interface clustering summary and cluster population table.",
        ]
        if summary:
            evidence.append(f"The clustering summary reports {summary.get('n_clusters')} interface states across {summary.get('n_frames')} sampled frames.")
        if descriptor:
            evidence.append(f"The dominant state is associated with: {descriptor}.")
        if population is not None:
            evidence.append(f"Cluster {dominant_cluster} ranks first by population share.")
        highlights = [
            {"label": "Dominant cluster", "value": f"cluster {dominant_cluster}"},
            {"label": "Population", "value": self._fmt_percent(population / 100) if population is not None else "n/a"},
            {"label": "State count", "value": str(summary.get("n_clusters", len(top_items)))},
        ]
        return self._finalize(
            "cluster_status",
            answer,
            evidence,
            sources,
            top_items,
            highlights=highlights,
            top_items_title="Cluster population snapshot",
            response_sections=[
                {
                    "title": "Direct readout",
                    "items": [
                        f"Dominant cluster: cluster {dominant_cluster}",
                        f"Population: {self._fmt_percent(population / 100) if population is not None else 'n/a'}",
                    ],
                },
                {
                    "title": "State interpretation",
                    "items": [descriptor] if descriptor else [],
                },
            ],
        )

    def _handle_capability_overview(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        capabilities = [
            "Query: answer focused questions about stability, interface summary, hotspots, RRCS pairs, residue ranking, regions, clusters, and landscape status.",
            "Diagnostic: combine available summaries into a single-system stability/interface diagnostic.",
            "Design: provide first-pass candidate mutation or validation sites from report evidence.",
            "Compare: compare two result folders when a second case is provided.",
        ]
        answer = "I can help read the current single-job ImmunoScope report, route questions to report skills, and return evidence-backed answers."
        return self._finalize(
            "capability_overview",
            answer,
            capabilities,
            sources,
            top_items=[
                {"skill": "reporter-query", "scope": "lightweight factual lookup"},
                {"skill": "reporter-diagnostic", "scope": "single-system diagnosis"},
                {"skill": "reporter-design", "scope": "candidate site suggestion"},
                {"skill": "reporter-compare", "scope": "multi-case comparison"},
            ],
            highlights=[
                {"label": "Current scope", "value": "Single trajectory report"},
                {"label": "Primary mode", "value": "Evidence-backed report Q&A"},
                {"label": "Data source", "value": "Generated ImmunoScope outputs"},
            ],
            top_items_title="Available skills",
            response_sections=[
                {
                    "title": "Available capabilities",
                    "items": capabilities,
                },
                {
                    "title": "Example questions",
                    "items": [
                        "Summarize interface",
                        "Identify hotspots",
                        "Diagnose stability",
                        "Suggest mutation sites",
                    ],
                },
            ],
            notes=["Capability questions are answered directly and are not routed to RRCS hotspot analysis."],
        )

    def _handle_top_pair(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        pair_df = self._read_csv(self._find_name(sources, "annotated_rrcs_pair_summary.csv"))
        summary = self._read_json(self._find_name(sources, "rrcs_summary.json"))
        if pair_df.empty and not summary:
            return self._missing("top_pair", preferred_missing)
        top_items: list[dict[str, Any]] = []
        if not pair_df.empty:
            ranked = self._rank_pairs(pair_df)
            row = ranked.iloc[0]
            pair_display = self._pair_display(row)
            answer = f"The strongest RRCS hotspot pair is {pair_display}."
            evidence = [
                "Supported by RRCS summary and annotated hotspot pair table.",
                "This pair ranks first by mean_rrcs.",
                f"It is observed in {self._fmt_percent(row.get('rrcs_nonzero_fraction'))} of sampled frames.",
            ]
            if row.get("interaction_class"):
                evidence.append(f"It belongs to the {row.get('interaction_class')} interaction class.")
            top_items = [self._pair_record(record) for record in ranked.head(5).to_dict(orient="records")]
            highlights = [
                {"label": "Top pair", "value": pair_display},
                {"label": "Mean RRCS", "value": self._fmt_float(row.get("mean_rrcs"), 3)},
                {"label": "Observed", "value": self._fmt_percent(row.get("rrcs_nonzero_fraction"))},
            ]
            return self._finalize(
                "top_pair",
                answer,
                evidence,
                sources,
                top_items,
                highlights=highlights,
                top_items_title="Top hotspot pairs",
                response_sections=[
                    {
                        "title": "Hotspot pair",
                        "items": [
                            pair_display,
                            f"Mean RRCS: {self._fmt_float(row.get('mean_rrcs'), 3)}",
                            f"Observed: {self._fmt_percent(row.get('rrcs_nonzero_fraction'))}",
                        ],
                    }
                ],
            )
        top_pair = (summary.get("top_pairs") or [{}])[0]
        answer = "The strongest RRCS hotspot pair is the top-ranked pair in the RRCS summary."
        evidence = [
            "Supported by RRCS summary.",
            f"The case reports {summary.get('n_nonzero_pairs')} identified nonzero RRCS pairs.",
        ]
        if top_pair:
            evidence.append(f"The top ranked pair entry is: {top_pair}.")
        return self._finalize("top_pair", answer, evidence, sources, top_items)

    def _handle_persistent_pair(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        pair_df = self._read_csv(self._find_name(sources, "annotated_rrcs_pair_summary.csv"))
        if pair_df.empty:
            return self._missing("persistent_pair", preferred_missing)
        ranked = pair_df.sort_values("rrcs_nonzero_fraction", ascending=False)
        row = ranked.iloc[0]
        pair_display = self._pair_display(row)
        answer = f"The most persistent hotspot pair is {pair_display}."
        evidence = [
            "Preferred persistence digest is missing; using RRCS pair persistence as the current fallback.",
            "This pair ranks first by rrcs_nonzero_fraction in the available annotated RRCS pair table.",
            f"Its mean RRCS is {self._fmt_float(row.get('mean_rrcs'), 3)}.",
            f"It is nonzero in {self._fmt_percent(row.get('rrcs_nonzero_fraction'))} of sampled frames.",
        ]
        top_items = [self._pair_record(record) for record in ranked.head(5).to_dict(orient="records")]
        highlights = [
            {"label": "Most persistent pair", "value": pair_display},
            {"label": "Observed", "value": self._fmt_percent(row.get("rrcs_nonzero_fraction"))},
            {"label": "Mean RRCS", "value": self._fmt_float(row.get("mean_rrcs"), 3)},
        ]
        return self._finalize(
            "persistent_pair",
            answer,
            evidence,
            sources,
            top_items,
            confidence="provisional",
            highlights=highlights,
            top_items_title="Persistent hotspot pairs",
            notes=["This answer currently uses RRCS persistence as a fallback because an explicit persistence digest is not available."],
            response_sections=[
                {
                    "title": "Fallback logic",
                    "items": [
                        "Using RRCS pair persistence as the current fallback.",
                        f"Leading pair: {pair_display}",
                    ],
                }
            ],
        )

    def _handle_top_region(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        region_df = self._read_csv(self._find_name(sources, "rrcs_region_summary.csv"))
        if region_df.empty:
            return self._missing("top_region", preferred_missing)
        ranked = region_df.sort_values("mean_rrcs_sum", ascending=False)
        row = ranked.iloc[0]
        region_name = self._format_region(row)
        answer = f"The most active interface region is {region_name}."
        evidence = [
            "Supported by the RRCS region summary.",
            "This region ranks first by mean_rrcs_sum.",
        ]
        if row.get("interaction_class"):
            evidence.append(f"It is classified as {row.get('interaction_class')}.")
        if row.get("rrcs_nonzero_fraction_mean") is not None:
            evidence.append(f"The mean nonzero RRCS fraction for this region is {self._fmt_percent(row.get('rrcs_nonzero_fraction_mean'), 2)}.")
        top_items = [self._region_record(record) for record in ranked.head(5).to_dict(orient="records")]
        highlights = [
            {"label": "Top region", "value": region_name},
            {"label": "Mean RRCS sum", "value": self._fmt_float(row.get("mean_rrcs_sum"), 3)},
            {"label": "Observed", "value": self._fmt_percent(row.get("rrcs_nonzero_fraction_mean"), 2)},
        ]
        return self._finalize(
            "top_region",
            answer,
            evidence,
            sources,
            top_items,
            highlights=highlights,
            top_items_title="Most active interface regions",
            response_sections=[
                {
                    "title": "Region focus",
                    "items": [
                        region_name,
                        f"Mean RRCS sum: {self._fmt_float(row.get('mean_rrcs_sum'), 3)}",
                    ],
                }
            ],
        )

    def _handle_quality_status(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        summary = self._read_json(self._find_name(sources, "rmsf_summary.json"))
        if not summary:
            summary = self._read_json(self._find_name_suffix(sources, "summary.json"))
        if not summary:
            return self._missing("quality_status", preferred_missing, confidence="provisional")
        label = "stable"
        mean_rmsf = summary.get("mean_rmsf_angstrom")
        max_rmsf = summary.get("max_rmsf_angstrom")
        if mean_rmsf is not None and float(mean_rmsf) > 3.0:
            label = "caution"
        answer = f"The current trajectory is classified as {label}."
        evidence = [
            "Supported by the available RMSF or quality-style summary output.",
        ]
        if mean_rmsf is not None:
            evidence.append(f"The mean RMSF is {self._fmt_float(mean_rmsf, 3)} Å.")
        if max_rmsf is not None:
            evidence.append(f"The maximum RMSF is {self._fmt_float(max_rmsf, 3)} Å.")
        if summary.get("n_frames") is not None:
            evidence.append(f"The summary covers {summary['n_frames']} sampled frames.")
        highlights = [
            {"label": "Verdict", "value": label},
            {"label": "Mean RMSF", "value": f"{self._fmt_float(mean_rmsf, 2)} Å" if mean_rmsf is not None else "n/a"},
            {"label": "Max RMSF", "value": f"{self._fmt_float(max_rmsf, 2)} Å" if max_rmsf is not None else "n/a"},
        ]
        return self._finalize(
            "quality_status",
            answer,
            evidence,
            sources,
            [],
            highlights=highlights,
            top_items_title="",
            response_sections=[
                {
                    "title": "Stability readout",
                    "items": [
                        f"Verdict: {label}",
                        f"Mean RMSF: {self._fmt_float(mean_rmsf, 2)} Å" if mean_rmsf is not None else "Mean RMSF unavailable",
                        f"Max RMSF: {self._fmt_float(max_rmsf, 2)} Å" if max_rmsf is not None else "Max RMSF unavailable",
                    ],
                }
            ],
        )

    def _handle_landscape_status(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        summary_path = self._find_name(sources, "landscape_summary.json")
        loadings_path = self._find_name(sources, "loadings.csv")
        pca_path = self._find_name(sources, "pca_coordinates.csv")
        used_sources = [path for path in [summary_path, loadings_path, pca_path] if path is not None]
        summary = self._read_json(summary_path)
        loadings = self._read_csv(loadings_path)
        if not summary and loadings.empty:
            return self._missing("landscape_status", preferred_missing, confidence="provisional")
        explained = summary.get("explained_variance") or []
        pc1 = explained[0] if len(explained) > 0 else None
        pc2 = explained[1] if len(explained) > 1 else None
        top_pc1 = summary.get("top_contributors_pc1") or []
        top_pc2 = summary.get("top_contributors_pc2") or []
        pc1_feature = top_pc1[0].get("feature") if top_pc1 else "n/a"
        pc2_feature = top_pc2[0].get("feature") if top_pc2 else "n/a"
        answer = f"FEL/PCA landscape is mainly driven by {pc1_feature} on PC1 and {pc2_feature} on PC2."
        evidence = [
            f"Frames: {summary.get('n_frames', 'n/a')}; features: {summary.get('n_features', 'n/a')}.",
            f"PC1 explained variance: {self._fmt_float((pc1 or 0) * 100, 1)}%.",
            f"PC2 explained variance: {self._fmt_float((pc2 or 0) * 100, 1)}%.",
            "Top contributors are read from landscape_summary.json.",
        ]
        top_items = [
            {"component": "PC1", "feature": item.get("feature"), "loading": item.get("loading")}
            for item in top_pc1[:3]
        ] + [
            {"component": "PC2", "feature": item.get("feature"), "loading": item.get("loading")}
            for item in top_pc2[:3]
        ]
        highlights = [
            {"label": "Frames", "value": summary.get("n_frames", "n/a")},
            {"label": "Features", "value": summary.get("n_features", "n/a")},
            {"label": "PC1 variance", "value": f"{self._fmt_float((pc1 or 0) * 100, 1)}%"},
            {"label": "PC2 variance", "value": f"{self._fmt_float((pc2 or 0) * 100, 1)}%"},
        ]
        return self._finalize(
            "landscape_status",
            answer,
            evidence,
            used_sources,
            top_items,
            confidence="high" if summary else "medium",
            highlights=highlights,
            top_items_title="Top FEL/PCA contributors",
            response_sections=[
                {
                    "title": "Landscape readout",
                    "items": [
                        f"PC1 driver: {pc1_feature}",
                        f"PC2 driver: {pc2_feature}",
                        f"PC1 + PC2 variance: {self._fmt_float(((pc1 or 0) + (pc2 or 0)) * 100, 1)}%",
                    ],
                }
            ],
        )

    def _handle_interface_summary(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        summary = self._read_json(self._find_name(sources, "interface_summary.json"))
        overview = self._read_json(self._find_name(sources, "interaction_overview.json"))
        region_df = self._read_csv(self._find_name(sources, "rrcs_region_summary.csv"))
        pair_df = self._read_csv(self._find_name(sources, "annotated_rrcs_pair_summary.csv"))
        if not summary and not overview and region_df.empty:
            return self._missing("interface_summary", preferred_missing)

        mean_bsa = ((summary.get("buried_surface_area") or {}).get("mean") if summary else None)
        mean_ratio = ((summary.get("interface_ratio") or {}).get("mean") if summary else None)
        families = [record for record in (overview.get("families") or []) if record.get("pair_rows", 0) > 0]
        family_items = [
            {
                "family": record.get("label") or record.get("family"),
                "pair_rows": record.get("pair_rows"),
                "heatmap_count": record.get("heatmap_count"),
            }
            for record in sorted(families, key=lambda item: item.get("pair_rows", 0), reverse=True)[:5]
        ]
        top_region = None
        if not region_df.empty:
            top_region = self._format_region(region_df.sort_values("mean_rrcs_sum", ascending=False).iloc[0])
        top_pair = None
        if not pair_df.empty:
            top_pair = self._pair_display(self._rank_pairs(pair_df).iloc[0])

        answer = f"The complex interface is characterized by a mean buried surface area of {self._fmt_float(mean_bsa, 1)} A^2 and {len(families) or 'n/a'} active interaction families."
        evidence = [
            "Supported by interface summary, interaction overview, and RRCS region digest.",
            f"Mean interface ratio: {self._fmt_float(mean_ratio, 3)}." if mean_ratio is not None else "Interface ratio is unavailable.",
            f"Most active RRCS region: {top_region}." if top_region else "RRCS region digest is unavailable.",
            f"Top hotspot pair: {top_pair}." if top_pair else "RRCS hotspot pair digest is unavailable.",
        ]
        highlights = [
            {"label": "Mean BSA", "value": f"{self._fmt_float(mean_bsa, 1)} Å²" if mean_bsa is not None else "n/a"},
            {"label": "Interface ratio", "value": self._fmt_float(mean_ratio, 3) if mean_ratio is not None else "n/a"},
            {"label": "Active families", "value": str(len(families))},
            {"label": "Top region", "value": top_region or "n/a"},
        ]
        return self._finalize(
            "interface_summary",
            answer,
            evidence,
            sources,
            family_items,
            highlights=highlights,
            top_items_title="Active interaction families",
            response_sections=[
                {
                    "title": "Interface summary",
                    "items": [
                        f"Mean BSA: {self._fmt_float(mean_bsa, 1)} Å²" if mean_bsa is not None else "Mean BSA unavailable",
                        f"Interface ratio: {self._fmt_float(mean_ratio, 3)}" if mean_ratio is not None else "Interface ratio unavailable",
                        f"Top region: {top_region}" if top_region else "Top region unavailable",
                    ],
                }
            ],
        )

    def _handle_hotspot_summary(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        summary = self._read_json(self._find_name(sources, "rrcs_summary.json"))
        pair_df = self._read_csv(self._find_name(sources, "annotated_rrcs_pair_summary.csv"))
        region_df = self._read_csv(self._find_name(sources, "rrcs_region_summary.csv"))
        cluster_summary = self._read_json(self._find_name(sources, "interface_clustering_summary.json"))
        if not summary and pair_df.empty and region_df.empty:
            return self._missing("hotspot_summary", preferred_missing)

        ranked_pairs = self._rank_pairs(pair_df) if not pair_df.empty else pd.DataFrame()
        top_pair = self._pair_display(ranked_pairs.iloc[0]) if not ranked_pairs.empty else None
        top_residue = self._rank_residues_from_pairs(pair_df)[0] if not pair_df.empty else None
        top_region_row = region_df.sort_values("mean_rrcs_sum", ascending=False).iloc[0] if not region_df.empty else None
        top_region = self._format_region(top_region_row) if top_region_row is not None else None

        answer_parts = []
        if top_pair:
            answer_parts.append(f"strongest hotspot pair: {top_pair}")
        if top_residue:
            answer_parts.append(f"top residue of interest: {top_residue.get('residue_label')}")
        if top_region:
            answer_parts.append(f"most active region: {top_region}")
        answer = "; ".join(answer_parts) + "." if answer_parts else "Hotspot interpretation requires an RRCS digest."
        evidence = [
            f"RRCS identifies {summary.get('n_nonzero_pairs')} nonzero hotspot pairs." if summary else "RRCS summary is unavailable.",
            f"The top RRCS pair is {top_pair}." if top_pair else "Top RRCS pair is unavailable.",
            f"The dominant hotspot region is {top_region}." if top_region else "Hotspot region digest is unavailable.",
            f"Dominant interface state: cluster {(cluster_summary.get('largest_cluster') or {}).get('cluster_id')} at {self._fmt_percent(((cluster_summary.get('largest_cluster') or {}).get('population_percent') or 0) / 100)}." if cluster_summary else "Cluster digest is unavailable.",
        ]
        top_items = [self._pair_record(record) for record in ranked_pairs.head(5).to_dict(orient="records")] if not ranked_pairs.empty else []
        highlights = [
            {"label": "Identified pairs", "value": str(summary.get("n_nonzero_pairs")) if summary else "n/a"},
            {"label": "Top pair", "value": top_pair or "n/a"},
            {"label": "Top residue", "value": top_residue.get("residue_label") if top_residue else "n/a"},
            {"label": "Top region", "value": top_region or "n/a"},
        ]
        return self._finalize(
            "hotspot_summary",
            answer,
            evidence,
            sources,
            top_items,
            highlights=highlights,
            top_items_title="Top hotspot pairs",
            response_sections=[
                {
                    "title": "Hotspot readout",
                    "items": [
                        f"Top pair: {top_pair}" if top_pair else "Top pair unavailable",
                        f"Top residue: {top_residue.get('residue_label')}" if top_residue else "Top residue unavailable",
                        f"Top region: {top_region}" if top_region else "Top region unavailable",
                    ],
                }
            ],
        )

    def _handle_design_hint(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        pair_df = self._read_csv(self._find_name(sources, "annotated_rrcs_pair_summary.csv"))
        region_df = self._read_csv(self._find_name(sources, "rrcs_region_summary.csv"))
        cluster_summary = self._read_json(self._find_name(sources, "interface_clustering_summary.json"))
        if pair_df.empty and region_df.empty:
            return self._missing("design_hint", preferred_missing, confidence="provisional")

        residue_rank = self._rank_residues_from_pairs(pair_df) if not pair_df.empty else []
        top_candidates = residue_rank[:5]
        candidate_labels = [item.get("residue_label") for item in top_candidates[:3] if item.get("residue_label")]
        top_region = None
        if not region_df.empty:
            top_region = self._format_region(region_df.sort_values("mean_rrcs_sum", ascending=False).iloc[0])
        top_pair = None
        if not pair_df.empty:
            top_pair = self._pair_display(self._rank_pairs(pair_df).iloc[0])
        answer = (
            f"First-pass candidate sites for manual review: {', '.join(candidate_labels)}. These sites come from RRCS hotspot and interface-state evidence and should be treated as screening hints only."
            if candidate_labels else
            "Only a preliminary screening hint is available; RRCS hotspot support should be completed first."
        )
        evidence = [
            f"Leading hotspot pair: {top_pair}." if top_pair else "Leading hotspot pair is unavailable.",
            f"Leading hotspot region: {top_region}." if top_region else "Hotspot region digest is unavailable.",
            (f"Dominant interface state: cluster {(cluster_summary.get('largest_cluster') or {}).get('cluster_id')} with descriptor {(cluster_summary.get('largest_cluster') or {}).get('main_structural_descriptor')}." if cluster_summary else "Cluster digest is unavailable."),
            "Treat these residues as screening candidates only; final mutation design still requires full pair tables, structure review, and manual validation.",
        ]
        highlights = [
            {"label": "Candidate count", "value": str(len(top_candidates))},
            {"label": "Top region", "value": top_region or "n/a"},
            {"label": "Top pair", "value": top_pair or "n/a"},
            {"label": "Confidence", "value": "provisional"},
        ]
        top_items = [
            {
                "candidate_residue": item.get("residue_label"),
                "region": item.get("region"),
                "best_mean_rrcs": item.get("best_mean_rrcs"),
                "observed": item.get("best_nonzero_fraction"),
            }
            for item in top_candidates
        ]
        return self._finalize(
            "design_hint",
            answer,
            evidence,
            sources,
            top_items,
            confidence="provisional",
            highlights=highlights,
            top_items_title="Screening candidate residues",
            notes=[
                "This output is intentionally conservative.",
                "Do not treat it as a validated mutation recommendation without structural review and full-table verification.",
            ],
            response_sections=[
                {
                    "title": "Screening scope",
                    "items": [
                        f"Prioritize manual review of: {', '.join(candidate_labels)}" if candidate_labels else "No residue shortlist available",
                        f"Primary hotspot region: {top_region}" if top_region else "Primary hotspot region unavailable",
                    ],
                },
                {
                    "title": "Boundary",
                    "items": [
                        "This is a screening hint, not a mutation prescription.",
                        "Manual structural review is still required before proposing sequence changes.",
                    ],
                },
            ],
        )

    def _handle_top_residue(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        pair_df = self._read_csv(self._find_name(sources, "annotated_rrcs_pair_summary.csv"))
        if pair_df.empty:
            return self._missing("top_residue", preferred_missing)
        ranked = self._rank_residues_from_pairs(pair_df)
        top = ranked[0]
        region_text = f" in the {top['region']} region" if top.get("region") else ""
        answer = f"The top hotspot residue of interest is {top['residue_label']}{region_text}."
        evidence = [
            "Primary residue ranking is unavailable, so the current answer is derived from RRCS hotspot support.",
            f"{top['residue_label']} appears in the strongest RRCS-supported hotspot pairs.",
            f"Its best supported mean RRCS is {self._fmt_float(top['best_mean_rrcs'], 3)}.",
            f"Its best nonzero RRCS fraction is {self._fmt_percent(top['best_nonzero_fraction'])}.",
        ]
        highlights = [
            {"label": "Top residue", "value": top['residue_label']},
            {"label": "Region", "value": top.get('region') or 'n/a'},
            {"label": "Best mean RRCS", "value": self._fmt_float(top['best_mean_rrcs'], 3)},
            {"label": "Observed", "value": self._fmt_percent(top['best_nonzero_fraction'])},
        ]
        return self._finalize(
            "top_residue",
            answer,
            evidence,
            sources,
            ranked[:5],
            confidence="provisional",
            highlights=highlights,
            top_items_title="Residue hotspot shortlist",
            notes=["This residue ranking is derived from RRCS hotspot support rather than a dedicated residue-scoring module."],
            response_sections=[
                {
                    "title": "Residue hotspot",
                    "items": [
                        f"Residue: {top['residue_label']}",
                        f"Region: {top.get('region') or 'n/a'}",
                    ],
                }
            ],
        )

    def _missing(self, query_type: str, preferred_missing: list[str], confidence: str = "low") -> QueryAnswer:
        missing_text = ", ".join(preferred_missing)
        return QueryAnswer(
            query_type=query_type,
            answer="The available summaries or digests are insufficient for a robust answer.",
            evidence=[
                f"Missing preferred source(s): {missing_text}.",
                "No stable fallback source was available for this query type.",
            ],
            sources_used=[],
            top_items=[],
            confidence=confidence,
            response_panel={
                "summary": "The preferred summary or digest is missing.",
                "sections": [
                    {
                        "title": "What is missing",
                        "items": [
                            f"Missing preferred source(s): {missing_text}.",
                            "No stable fallback source was available for this query type.",
                        ],
                    }
                ],
                "caution": "This answer is unavailable until the relevant summary or digest has been generated.",
            },
            evidence_panel={
                "title": "Evidence summary",
                "highlights": [],
                "bullets": [
                    f"Missing preferred source(s): {missing_text}.",
                    "No stable fallback source was available for this query type.",
                ],
                "top_items_title": "",
                "top_items": [],
                "sources": [],
                "notes": ["Check whether the corresponding summary or digest files were generated for this report."],
            },
        )

    def _finalize(
        self,
        query_type: str,
        answer: str,
        evidence: list[str],
        sources: list[Path],
        top_items: list[dict[str, Any]],
        confidence: str = "high",
        *,
        highlights: list[dict[str, Any]] | None = None,
        top_items_title: str = "Top items",
        notes: list[str] | None = None,
        response_sections: list[dict[str, Any]] | None = None,
    ) -> QueryAnswer:
        source_labels = [str(path.relative_to(self.case_dir)) for path in sources if path.exists()]
        evidence_list = [item for item in evidence if item][:5]
        panel = {
            "title": "Evidence summary",
            "highlights": highlights or [],
            "bullets": evidence_list,
            "top_items_title": top_items_title,
            "top_items": top_items[:5],
            "sources": source_labels,
            "notes": notes or [],
        }
        response_panel = {
            "summary": answer,
            "sections": [section for section in (response_sections or []) if section.get("items")],
            "caution": " ".join(notes or []) if notes else "",
        }
        return QueryAnswer(
            query_type=query_type,
            answer=answer,
            evidence=evidence_list,
            sources_used=source_labels,
            top_items=top_items[:5],
            confidence=confidence,
            response_panel=response_panel,
            evidence_panel=panel,
        )

    @staticmethod
    def _read_json(path: Path | None) -> dict[str, Any]:
        if path is None or not path.exists():
            return {}
        return json.loads(path.read_text(encoding="utf-8"))

    @staticmethod
    def _read_csv(path: Path | None) -> pd.DataFrame:
        if path is None or not path.exists():
            return pd.DataFrame()
        return pd.read_csv(path)

    @staticmethod
    def _find_name(paths: list[Path], name: str) -> Path | None:
        for path in paths:
            if path.name == name:
                return path
        return None

    @staticmethod
    def _find_name_suffix(paths: list[Path], suffix: str) -> Path | None:
        for path in paths:
            if path.name.endswith(suffix):
                return path
        return None

    @staticmethod
    def _fmt_float(value: Any, digits: int = 2) -> str:
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return "n/a"
        return f"{numeric:.{digits}f}"

    @staticmethod
    def _fmt_percent(value: Any, digits: int = 1) -> str:
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return "n/a"
        if numeric > 1.0:
            return f"{numeric:.{digits}f}%"
        return f"{numeric * 100:.{digits}f}%"

    @staticmethod
    def _pair_display(row: pd.Series | dict[str, Any]) -> str:
        if isinstance(row, pd.Series):
            row = row.to_dict()
        return row.get("pair_display") or f"{row.get('residue_label_1')} ↔ {row.get('residue_label_2')}"

    @staticmethod
    def _pair_record(record: dict[str, Any]) -> dict[str, Any]:
        return {
            "pair": record.get("pair_display") or f"{record.get('residue_label_1')} ↔ {record.get('residue_label_2')}",
            "interaction_class": record.get("interaction_class"),
            "mean_rrcs": record.get("mean_rrcs"),
            "max_rrcs": record.get("max_rrcs"),
            "observed": record.get("rrcs_nonzero_fraction"),
        }

    @staticmethod
    def _rank_pairs(pair_df: pd.DataFrame) -> pd.DataFrame:
        ranked = pair_df.copy()
        if "pair_display" not in ranked.columns:
            ranked["pair_display"] = ranked.apply(lambda row: f"{row.get('residue_label_1')} ↔ {row.get('residue_label_2')}", axis=1)
        return ranked.sort_values(["mean_rrcs", "rrcs_nonzero_fraction", "max_rrcs"], ascending=False)

    def _rank_residues_from_pairs(self, pair_df: pd.DataFrame) -> list[dict[str, Any]]:
        residue_scores: dict[str, dict[str, Any]] = {}
        for _, row in pair_df.iterrows():
            for side in (1, 2):
                label = row.get(f"residue_label_{side}") or f"{row.get(f'chain_id_{side}')}:{row.get(f'resid_{side}') }"
                rec = residue_scores.setdefault(label, {
                    "residue_label": label,
                    "best_mean_rrcs": 0.0,
                    "best_nonzero_fraction": 0.0,
                    "region": row.get(f"tcr_region_detailed_{side}") or row.get(f"mhc_subregion_{side}") or row.get(f"component_{side}"),
                })
                rec["best_mean_rrcs"] = max(rec["best_mean_rrcs"], float(row.get("mean_rrcs", 0.0) or 0.0))
                rec["best_nonzero_fraction"] = max(rec["best_nonzero_fraction"], float(row.get("rrcs_nonzero_fraction", 0.0) or 0.0))
        return sorted(residue_scores.values(), key=lambda item: (item["best_mean_rrcs"], item["best_nonzero_fraction"]), reverse=True)

    @staticmethod
    def _format_region(row: pd.Series | dict[str, Any]) -> str:
        if isinstance(row, pd.Series):
            row = row.to_dict()
        tcr_region = row.get("tcr_region")
        partner = row.get("partner_component")
        mhc_subregion = row.get("mhc_subregion")
        parts = []
        if tcr_region:
            parts.append(str(tcr_region))
        if partner:
            parts.append(str(partner))
        if mhc_subregion and str(mhc_subregion) != "nan":
            parts.append(str(mhc_subregion))
        return " ↔ ".join(parts) if parts else str(row.get("interaction_class", "top region"))

    @staticmethod
    def _region_record(record: dict[str, Any]) -> dict[str, Any]:
        return {
            "region": " ↔ ".join([str(item) for item in [record.get("tcr_region"), record.get("partner_component"), record.get("mhc_subregion")] if item and str(item) != "nan"]),
            "interaction_class": record.get("interaction_class"),
            "mean_rrcs_sum": record.get("mean_rrcs_sum"),
            "mean_rrcs_mean": record.get("mean_rrcs_mean"),
            "observed": record.get("rrcs_nonzero_fraction_mean"),
        }

    def _handle_system_diagnostic(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        """Handle full system diagnostic using SystemDiagnostic skill."""
        try:
            from skills.reporter_diagnostic import SystemDiagnostic

            diagnostic = SystemDiagnostic(self.case_dir)
            report = diagnostic.diagnose()

            # Format as QueryAnswer
            answer = report.overall_diagnosis
            evidence = []

            # Collect evidence from all sections
            for section in report.sections:
                evidence.append(f"{section.module}: {section.conclusion} (confidence: {section.confidence})")

            # Add hotspots if available
            if report.potential_hotspots:
                evidence.append(f"Identified {len(report.potential_hotspots)} potential hotspots")

            return QueryAnswer(
                query_type="system_diagnostic",
                answer=answer,
                evidence=evidence[:5],  # Top 5 evidence items
                sources_used=[str(s) for s in sources],
                top_items=report.potential_hotspots[:5],
                confidence=report.confidence,
                skill="reporter-diagnostic",
            )

        except Exception as e:
            return QueryAnswer(
                query_type="system_diagnostic",
                answer=f"System diagnostic failed: {str(e)}",
                evidence=["Error occurred during diagnostic analysis"],
                sources_used=[],
                top_items=[],
                confidence="low",
                skill="reporter-diagnostic",
            )

    def _handle_mutation_design(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        """Handle mutation design suggestions using MutationDesignAdvisor skill."""
        try:
            from skills.reporter_design import MutationDesignAdvisor

            advisor = MutationDesignAdvisor(self.case_dir)
            report = advisor.suggest_designs()

            # Format as QueryAnswer
            if report.candidates:
                top_candidate = report.candidates[0]
                answer = f"Top mutation candidate: {top_candidate.residue} ({top_candidate.chain}, {top_candidate.region}). {top_candidate.rationale}"
            else:
                answer = "No mutation candidates identified from available data."

            evidence = []
            for candidate in report.candidates[:3]:
                evidence.append(f"{candidate.residue}: {candidate.rationale} (confidence: {candidate.confidence})")

            # Convert candidates to dict format
            top_items = []
            for candidate in report.candidates[:5]:
                top_items.append({
                    "residue": candidate.residue,
                    "chain": candidate.chain,
                    "region": candidate.region,
                    "rationale": candidate.rationale,
                    "confidence": candidate.confidence,
                })

            return QueryAnswer(
                query_type="mutation_design",
                answer=answer,
                evidence=evidence,
                sources_used=[str(s) for s in sources],
                top_items=top_items,
                confidence=report.confidence,
                skill="reporter-design",
            )

        except Exception as e:
            return QueryAnswer(
                query_type="mutation_design",
                answer=f"Mutation design failed: {str(e)}",
                evidence=["Error occurred during design analysis"],
                sources_used=[],
                top_items=[],
                confidence="low",
                skill="reporter-design",
            )

    def _handle_system_compare(self, sources: list[Path], preferred_missing: list[str]) -> QueryAnswer:
        """Handle two-system comparison using SystemComparator skill."""
        # Note: This handler requires two case directories
        # For now, return a placeholder that explains the requirement
        return QueryAnswer(
            query_type="system_compare",
            answer="System comparison requires two case directories. Please use the SystemComparator class directly with two case_dir paths.",
            evidence=[
                "Comparison requires: case_dir_a and case_dir_b",
                "Use: SystemComparator(case_a, case_b, 'SystemA', 'SystemB').compare()",
            ],
            sources_used=[],
            top_items=[],
            confidence="medium",
            skill="reporter-compare",
        )
