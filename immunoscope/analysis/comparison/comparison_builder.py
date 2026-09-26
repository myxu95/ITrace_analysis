"""Build two-system comparison tables and summaries."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import pandas as pd

from .comparison_schema import SingleCaseArtifacts
from .comparison_preflight import ComparisonPreflightBuilder
from .statistical_tests import StatisticalComparator


@dataclass(slots=True)
class ComparisonBuildResult:
    """Comparison build result."""

    summary: dict
    comparison_rows: list[dict[str, Any]]
    identity_rows: list[dict[str, Any]]
    rmsf_region_rows: list[dict[str, Any]]
    rrcs_region_rows: list[dict[str, Any]]
    interaction_family_rows: list[dict[str, Any]]


class SystemComparisonBuilder:
    """Align two single-system results into unified comparison artifacts."""

    RMSF_REGIONS = (
        "CDR3_alpha",
        "CDR3_beta",
        "peptide",
        "alpha1_helix",
        "alpha2_helix",
    )

    FAMILY_ORDER = ("contact", "hbond", "saltbridge", "hydrophobic", "pipi", "cationpi")

    def __init__(self, enable_statistics: bool = True):
        """Initialize comparison builder.

        Args:
            enable_statistics: Whether to compute statistical tests for comparisons
        """
        self.enable_statistics = enable_statistics
        self.stat_comparator = StatisticalComparator() if enable_statistics else None

    def build(
        self,
        case_a: SingleCaseArtifacts,
        case_b: SingleCaseArtifacts,
        comparison_mode: str = "generic",
        comparison_context: str = "",
        comparison_scope: str = "auto",
        alignment_selection: str = "phla_core_ca",
        residue_mapping: str = "auto",
    ) -> ComparisonBuildResult:
        comparison_rows = []
        comparison_rows.extend(self._build_identity_metric_rows(case_a, case_b))
        comparison_rows.extend(self._build_quality_rows(case_a, case_b))
        comparison_rows.extend(self._build_interface_rows(case_a, case_b))
        comparison_rows.extend(self._build_flexibility_rows(case_a, case_b))
        comparison_rows.extend(self._build_rrcs_rows(case_a, case_b))
        comparison_rows.extend(self._build_interaction_rows(case_a, case_b))

        identity_rows = self._build_identity_rows(case_a, case_b)
        rmsf_region_rows = self._build_rmsf_region_rows(case_a, case_b)
        rrcs_region_rows = self._build_rrcs_region_rows(case_a, case_b)
        interaction_family_rows = self._build_interaction_family_rows(case_a, case_b)
        preflight = ComparisonPreflightBuilder().build(
            case_a=case_a,
            case_b=case_b,
            comparison_mode=comparison_mode,
            requested_scope=comparison_scope,
            alignment_selection=alignment_selection,
            residue_mapping=residue_mapping,
        )
        summary = self._build_summary(
            case_a,
            case_b,
            comparison_rows,
            rmsf_region_rows,
            rrcs_region_rows,
            interaction_family_rows,
            preflight=preflight.to_dict(),
            comparison_mode=comparison_mode,
            comparison_context=comparison_context,
        )

        return ComparisonBuildResult(
            summary=summary,
            comparison_rows=comparison_rows,
            identity_rows=identity_rows,
            rmsf_region_rows=rmsf_region_rows,
            rrcs_region_rows=rrcs_region_rows,
            interaction_family_rows=interaction_family_rows,
        )

    def _build_identity_metric_rows(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        peptide_a = case_a.identity.get("peptide_identity", {})
        peptide_b = case_b.identity.get("peptide_identity", {})
        hla_a = case_a.identity.get("hla_identity", {})
        hla_b = case_b.identity.get("hla_identity", {})
        tcr_a = case_a.identity.get("tcr_identity", {})
        tcr_b = case_b.identity.get("tcr_identity", {})

        return [
            self._row("identity", "peptide_sequence", "", peptide_a.get("sequence"), peptide_b.get("sequence")),
            self._row("identity", "peptide_length", "aa", peptide_a.get("length"), peptide_b.get("length"), numeric=True),
            self._row("identity", "hla_locus", "", hla_a.get("best_locus"), hla_b.get("best_locus")),
            self._row(
                "identity",
                "tcr_alpha_genotype",
                "",
                self._join_genotype(tcr_a.get("alpha_v_gene"), tcr_a.get("alpha_j_gene")),
                self._join_genotype(tcr_b.get("alpha_v_gene"), tcr_b.get("alpha_j_gene")),
            ),
            self._row(
                "identity",
                "tcr_beta_genotype",
                "",
                self._join_genotype(tcr_a.get("beta_v_gene"), tcr_a.get("beta_j_gene")),
                self._join_genotype(tcr_b.get("beta_v_gene"), tcr_b.get("beta_j_gene")),
            ),
        ]

    def _build_quality_rows(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        return [
            self._row("quality", "frames", "frames", case_a.quality.get("n_frames"), case_b.quality.get("n_frames"), numeric=True),
            self._row(
                "quality",
                "time_span_ps",
                "ps",
                self._time_span(case_a.quality),
                self._time_span(case_b.quality),
                numeric=True,
            ),
            self._row(
                "quality",
                "tail90_rmsd_variation_nm",
                "nm",
                case_a.quality.get("tail90_variation_nm"),
                case_b.quality.get("tail90_variation_nm"),
                numeric=True,
            ),
            self._row(
                "quality",
                "full_rmsd_variation_nm",
                "nm",
                case_a.quality.get("full_variation_nm"),
                case_b.quality.get("full_variation_nm"),
                numeric=True,
            ),
        ]

    def _build_interface_rows(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        return [
            self._row(
                "interface",
                "mean_bsa_angstrom2",
                "A^2",
                self._nested_mean(case_a.bsa, "buried_surface_area"),
                self._nested_mean(case_b.bsa, "buried_surface_area"),
                numeric=True,
            ),
            self._row(
                "interface",
                "mean_interface_ratio",
                "",
                self._nested_mean(case_a.bsa, "interface_ratio"),
                self._nested_mean(case_b.bsa, "interface_ratio"),
                numeric=True,
            ),
        ]

    def _build_flexibility_rows(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        rows = [
            self._row(
                "flexibility",
                "mean_rmsf_angstrom",
                "A",
                case_a.rmsf_summary.get("mean_rmsf_angstrom"),
                case_b.rmsf_summary.get("mean_rmsf_angstrom"),
                numeric=True,
            ),
            self._row(
                "flexibility",
                "tcr_mean_rmsf_angstrom",
                "A",
                case_a.rmsf_summary.get("tcr_mean_rmsf_angstrom"),
                case_b.rmsf_summary.get("tcr_mean_rmsf_angstrom"),
                numeric=True,
            ),
            self._row(
                "flexibility",
                "phla_mean_rmsf_angstrom",
                "A",
                case_a.rmsf_summary.get("phla_mean_rmsf_angstrom"),
                case_b.rmsf_summary.get("phla_mean_rmsf_angstrom"),
                numeric=True,
            ),
        ]
        for region in self.RMSF_REGIONS:
            rows.append(
                self._row(
                    "flexibility",
                    f"{region}_mean_rmsf_angstrom",
                    "A",
                    self._region_mean(case_a.rmsf_regions, region),
                    self._region_mean(case_b.rmsf_regions, region),
                    numeric=True,
                )
            )
        return rows

    def _build_interaction_rows(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        rows = [
            self._row(
                "interaction",
                "active_families",
                "count",
                self._active_family_count(case_a.interaction_overview),
                self._active_family_count(case_b.interaction_overview),
                numeric=True,
            ),
            self._row(
                "interaction",
                "total_residue_pairs",
                "count",
                self._total_pairs(case_a.interaction_overview),
                self._total_pairs(case_b.interaction_overview),
                numeric=True,
            ),
        ]
        for family in self.FAMILY_ORDER:
            rows.append(
                self._row(
                    "interaction",
                    f"{family}_pair_rows",
                    "count",
                    self._family_pair_rows(case_a.interaction_overview, family),
                    self._family_pair_rows(case_b.interaction_overview, family),
                    numeric=True,
                )
            )
        return rows

    def _build_rrcs_rows(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        return [
            self._row(
                "rrcs",
                "rrcs_total_pairs",
                "count",
                case_a.rrcs_summary.get("n_pairs"),
                case_b.rrcs_summary.get("n_pairs"),
                numeric=True,
            ),
            self._row(
                "rrcs",
                "rrcs_nonzero_pairs",
                "count",
                case_a.rrcs_summary.get("n_nonzero_pairs"),
                case_b.rrcs_summary.get("n_nonzero_pairs"),
                numeric=True,
            ),
            self._row(
                "rrcs",
                "peptide_tcr_rrcs_sum",
                "score",
                self._rrcs_interaction_sum(case_a.rrcs_regions, "peptide_tcr"),
                self._rrcs_interaction_sum(case_b.rrcs_regions, "peptide_tcr"),
                numeric=True,
            ),
            self._row(
                "rrcs",
                "hla_tcr_rrcs_sum",
                "score",
                self._rrcs_interaction_sum(case_a.rrcs_regions, "hla_tcr"),
                self._rrcs_interaction_sum(case_b.rrcs_regions, "hla_tcr"),
                numeric=True,
            ),
        ]

    def _build_identity_rows(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        tcr_a = case_a.identity.get("tcr_identity", {})
        tcr_b = case_b.identity.get("tcr_identity", {})
        peptide_a = case_a.identity.get("peptide_identity", {})
        peptide_b = case_b.identity.get("peptide_identity", {})
        hla_a = case_a.identity.get("hla_identity", {})
        hla_b = case_b.identity.get("hla_identity", {})
        return [
            {
                "field": "HLA locus",
                case_a.label: hla_a.get("best_locus", ""),
                case_b.label: hla_b.get("best_locus", ""),
            },
            {
                "field": "Top HLA candidate",
                case_a.label: hla_a.get("best_candidate_allele", ""),
                case_b.label: hla_b.get("best_candidate_allele", ""),
            },
            {
                "field": "Peptide",
                case_a.label: peptide_a.get("sequence", ""),
                case_b.label: peptide_b.get("sequence", ""),
            },
            {
                "field": "TCR alpha genotype",
                case_a.label: self._join_genotype(tcr_a.get("alpha_v_gene"), tcr_a.get("alpha_j_gene")),
                case_b.label: self._join_genotype(tcr_b.get("alpha_v_gene"), tcr_b.get("alpha_j_gene")),
            },
            {
                "field": "TCR beta genotype",
                case_a.label: self._join_genotype(tcr_a.get("beta_v_gene"), tcr_a.get("beta_j_gene")),
                case_b.label: self._join_genotype(tcr_b.get("beta_v_gene"), tcr_b.get("beta_j_gene")),
            },
            {
                "field": "CDR3 alpha",
                case_a.label: tcr_a.get("cdr3_alpha_sequence", ""),
                case_b.label: tcr_b.get("cdr3_alpha_sequence", ""),
            },
            {
                "field": "CDR3 beta",
                case_a.label: tcr_a.get("cdr3_beta_sequence", ""),
                case_b.label: tcr_b.get("cdr3_beta_sequence", ""),
            },
        ]

    def _build_rmsf_region_rows(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        rows = []
        for region in self.RMSF_REGIONS:
            value_a = self._region_mean(case_a.rmsf_regions, region)
            value_b = self._region_mean(case_b.rmsf_regions, region)
            rows.append(
                {
                    "region_group": region,
                    case_a.label: value_a,
                    case_b.label: value_b,
                    "delta": self._numeric_delta(value_a, value_b),
                }
            )
        return rows

    def _build_interaction_family_rows(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        merged = {}
        for family in self.FAMILY_ORDER:
            merged[family] = {
                "family": family,
                "label": self._family_label(case_a.interaction_overview, family) or self._family_label(case_b.interaction_overview, family) or family,
                case_a.label: self._family_pair_rows(case_a.interaction_overview, family),
                case_b.label: self._family_pair_rows(case_b.interaction_overview, family),
            }
            merged[family]["delta"] = self._numeric_delta(merged[family][case_a.label], merged[family][case_b.label])
        return list(merged.values())

    def _build_rrcs_region_rows(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        region_classes = ["peptide_tcr", "hla_tcr"]
        rows = []
        for interaction_class in region_classes:
            value_a = self._rrcs_interaction_sum(case_a.rrcs_regions, interaction_class)
            value_b = self._rrcs_interaction_sum(case_b.rrcs_regions, interaction_class)
            rows.append(
                {
                    "interaction_class": interaction_class,
                    case_a.label: value_a,
                    case_b.label: value_b,
                    "delta": self._numeric_delta(value_a, value_b),
                }
            )
        return rows

    def _build_summary(
        self,
        case_a: SingleCaseArtifacts,
        case_b: SingleCaseArtifacts,
        comparison_rows: list[dict[str, Any]],
        rmsf_region_rows: list[dict[str, Any]],
        rrcs_region_rows: list[dict[str, Any]],
        interaction_family_rows: list[dict[str, Any]],
        preflight: dict[str, Any],
        comparison_mode: str = "generic",
        comparison_context: str = "",
    ) -> dict:
        comparability = self._build_comparability(case_a, case_b)
        strongest_flex = max(
            rmsf_region_rows,
            key=lambda row: self._delta_rank(row.get("delta")),
            default=None,
        )
        strongest_family = max(
            interaction_family_rows,
            key=lambda row: self._delta_rank(row.get("delta")),
            default=None,
        )
        strongest_rrcs = max(
            rrcs_region_rows,
            key=lambda row: self._delta_rank(row.get("delta")),
            default=None,
        )
        takeaways = []
        flex_delta = self._to_float(strongest_flex.get("delta")) if strongest_flex else None
        if flex_delta is not None:
            delta = flex_delta
            if abs(delta) >= 1e-6:
                direction = "higher" if delta > 0 else "lower"
                takeaways.append(
                    f"Condition B shows {direction} flexibility in {strongest_flex['region_group']} relative to Condition A (Δ={delta:.3f} Å)."
                )
        family_delta = self._to_float(strongest_family.get("delta")) if strongest_family else None
        if family_delta is not None:
            delta = family_delta
            if abs(delta) >= 1e-6:
                direction = "more" if delta > 0 else "fewer"
                takeaways.append(
                    f"Condition B has {direction} {strongest_family['family']} residue pairs than Condition A (Δ={delta:.0f})."
                )
        rrcs_delta = self._to_float(strongest_rrcs.get("delta")) if strongest_rrcs else None
        if rrcs_delta is not None:
            delta = rrcs_delta
            if abs(delta) >= 1e-6:
                direction = "stronger" if delta > 0 else "weaker"
                takeaways.append(
                    f"Condition B shows {direction} RRCS coupling in {strongest_rrcs['interaction_class']} than Condition A (Δ={delta:.2f})."
                )
        quality_delta = self._numeric_delta(
            case_a.quality.get("tail90_variation_nm"),
            case_b.quality.get("tail90_variation_nm"),
        )
        if quality_delta is not None and abs(quality_delta) >= 1e-6:
            direction = "higher" if quality_delta > 0 else "lower"
            takeaways.append(
                f"Condition B has {direction} tail-90% RMSD variation than Condition A (Δ={quality_delta:.3f} nm)."
            )
        time_delta = self._numeric_delta(self._time_span(case_a.quality), self._time_span(case_b.quality))
        if time_delta is not None and abs(time_delta) >= 1e-6:
            direction = "longer" if time_delta > 0 else "shorter"
            takeaways.append(
                f"Condition B covers a {direction} analyzed time span than Condition A (Δ={time_delta:.1f} ps)."
            )
        bsa_delta = self._numeric_delta(
            self._nested_mean(case_a.bsa, "buried_surface_area"),
            self._nested_mean(case_b.bsa, "buried_surface_area"),
        )
        if bsa_delta is not None:
            if abs(float(bsa_delta)) >= 1e-6:
                direction = "higher" if bsa_delta > 0 else "lower"
                takeaways.append(
                    f"Condition B has a {direction} mean buried surface area than Condition A (Δ={bsa_delta:.1f} Å²)."
                )
        if not takeaways:
            takeaways.append("No material numeric differences were detected across the current comparison surface.")
        takeaways = self._prefix_takeaways_for_mode(takeaways, comparison_mode)
        return {
            "case_a": {
                "label": case_a.label,
                "case_root": str(case_a.case_root),
                "overview_root": str(case_a.overview_root),
                "source_paths": case_a.source_paths,
            },
            "case_b": {
                "label": case_b.label,
                "case_root": str(case_b.case_root),
                "overview_root": str(case_b.overview_root),
                "source_paths": case_b.source_paths,
            },
            "comparability": comparability,
            "preflight": preflight,
            "comparison_mode": comparison_mode,
            "comparison_context": comparison_context,
            "takeaways": takeaways,
            "metric_count": len(comparison_rows),
        }

    def _prefix_takeaways_for_mode(self, takeaways: list[str], comparison_mode: str) -> list[str]:
        prefix = {
            "generic": "",
            "mutation": "Mutation comparison: ",
            "sampling": "Sampling comparison: ",
            "replicate": "Replicate comparison: ",
        }.get(comparison_mode, "")
        if not prefix:
            return takeaways
        return [prefix + item for item in takeaways]

    def _build_comparability(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> dict:
        time_a = self._time_span(case_a.quality)
        time_b = self._time_span(case_b.quality)
        tail_a = case_a.quality.get("tail90_variation_nm")
        tail_b = case_b.quality.get("tail90_variation_nm")

        status = "Comparable with caution"
        reasons = []
        if pd.notna(time_a) and pd.notna(time_b):
            ratio = min(time_a, time_b) / max(time_a, time_b) if max(time_a, time_b) else 1.0
            if ratio >= 0.8:
                reasons.append("Time spans are within 20% of each other.")
            else:
                reasons.append("Time spans differ by more than 20%.")
        if pd.notna(tail_a) and pd.notna(tail_b):
            ratio = min(tail_a, tail_b) / max(tail_a, tail_b) if max(tail_a, tail_b) else 1.0
            if ratio >= 0.7:
                reasons.append("Tail-90% RMSD variation is in the same order of magnitude.")
            else:
                reasons.append("Tail-90% RMSD variation differs substantially.")
        if len(reasons) == 2 and "differ" not in " ".join(reasons):
            status = "Comparable"
        return {"status": status, "reasons": reasons}

    def _row(
        self,
        category: str,
        metric: str,
        unit: str,
        value_a: Any,
        value_b: Any,
        numeric: bool = False,
        timeseries_a: pd.Series | None = None,
        timeseries_b: pd.Series | None = None,
    ) -> dict[str, Any]:
        row = {
            "category": category,
            "metric": metric,
            "unit": unit,
            "case_a": value_a,
            "case_b": value_b,
            "delta": "",
            "relative_change_percent": "",
            "p_value": "",
            "effect_size": "",
            "significance": "",
        }
        if numeric:
            delta = self._numeric_delta(value_a, value_b)
            row["delta"] = delta
            numeric_a = self._to_float(value_a)
            if delta is not None and numeric_a not in (None, 0.0):
                row["relative_change_percent"] = (delta / numeric_a) * 100.0

            # Add statistical tests if time series data is available
            if (self.enable_statistics and self.stat_comparator is not None
                and timeseries_a is not None and timeseries_b is not None
                and len(timeseries_a) > 0 and len(timeseries_b) > 0):
                try:
                    stat_result = self.stat_comparator.compare_distributions(
                        timeseries_a, timeseries_b
                    )
                    row["p_value"] = stat_result.p_value
                    row["effect_size"] = stat_result.effect_size
                    row["significance"] = self.stat_comparator.significance_stars(stat_result.p_value)
                except Exception:
                    # Silently skip statistical tests if they fail
                    pass

        return row

    def _time_span(self, quality: dict) -> float | None:
        start = quality.get("tail90_start_time_ps", quality.get("trimmed_start_time_ps", quality.get("time_start")))
        end = quality.get("tail90_end_time_ps", quality.get("time_end"))
        numeric_start = self._to_float(start)
        numeric_end = self._to_float(end)
        if numeric_start is None or numeric_end is None:
            return None
        return numeric_end - numeric_start

    def _nested_mean(self, payload: dict, key: str) -> float | None:
        value = payload.get(key, {})
        if not isinstance(value, dict):
            return None
        return self._to_float(value.get("mean"))

    def _region_mean(self, frame: pd.DataFrame, region_group: str) -> float | None:
        if frame.empty or "region_group" not in frame.columns or "mean_rmsf_angstrom" not in frame.columns:
            return None
        matched = frame.loc[frame["region_group"] == region_group, "mean_rmsf_angstrom"]
        if matched.empty:
            return None
        numeric = pd.to_numeric(matched, errors="coerce").dropna()
        if numeric.empty:
            return None
        return float(numeric.iloc[0])

    def _active_family_count(self, frame: pd.DataFrame) -> int | None:
        if frame.empty or "pair_rows" not in frame.columns:
            return None
        pair_rows = pd.to_numeric(frame["pair_rows"], errors="coerce").fillna(0)
        return int((pair_rows > 0).sum())

    def _total_pairs(self, frame: pd.DataFrame) -> int | None:
        if frame.empty or "pair_rows" not in frame.columns:
            return None
        pair_rows = pd.to_numeric(frame["pair_rows"], errors="coerce").fillna(0)
        return int(pair_rows.sum())

    def _family_pair_rows(self, frame: pd.DataFrame, family: str) -> int | None:
        if frame.empty or "family" not in frame.columns or "pair_rows" not in frame.columns:
            return None
        matched = frame.loc[frame["family"] == family, "pair_rows"]
        if matched.empty:
            return 0
        numeric = pd.to_numeric(matched, errors="coerce").dropna()
        if numeric.empty:
            return 0
        return int(numeric.iloc[0])

    def _family_label(self, frame: pd.DataFrame, family: str) -> str | None:
        if frame.empty or "family" not in frame.columns or "label" not in frame.columns:
            return None
        matched = frame.loc[frame["family"] == family, "label"]
        if matched.empty:
            return None
        return str(matched.iloc[0])

    def _rrcs_interaction_sum(self, frame: pd.DataFrame, interaction_class: str) -> float | None:
        if frame.empty or "interaction_class" not in frame.columns or "mean_rrcs_sum" not in frame.columns:
            return None
        matched = frame.loc[frame["interaction_class"] == interaction_class, "mean_rrcs_sum"]
        if matched.empty:
            return None
        numeric = pd.to_numeric(matched, errors="coerce").dropna()
        if numeric.empty:
            return None
        return float(numeric.sum())

    def _numeric_delta(self, value_a: Any, value_b: Any) -> float | None:
        numeric_a = self._to_float(value_a)
        numeric_b = self._to_float(value_b)
        if numeric_a is None or numeric_b is None:
            return None
        return numeric_b - numeric_a

    def _delta_rank(self, value: Any) -> float:
        numeric = self._to_float(value)
        if numeric is None:
            return -1.0
        return abs(numeric)

    def _to_float(self, value: Any) -> float | None:
        if value in (None, ""):
            return None
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return None
        if math.isnan(numeric) or math.isinf(numeric):
            return None
        return numeric

    def _join_genotype(self, v_gene: str | None, j_gene: str | None) -> str:
        tokens = [token for token in (v_gene, j_gene) if token]
        return " / ".join(tokens)
