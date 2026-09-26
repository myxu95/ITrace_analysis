"""Comparison preflight checks."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .comparison_schema import SingleCaseArtifacts


@dataclass(slots=True)
class ComparisonPreflightResult:
    """Comparability assessment before compare execution."""

    requested_scope: str
    resolved_scope: str
    status: str
    alignment_policy: dict[str, Any]
    residue_mapping_policy: dict[str, Any]
    checks: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "requested_scope": self.requested_scope,
            "resolved_scope": self.resolved_scope,
            "status": self.status,
            "alignment_policy": self.alignment_policy,
            "residue_mapping_policy": self.residue_mapping_policy,
            "checks": self.checks,
            "warnings": self.warnings,
            "limitations": self.limitations,
        }


class ComparisonPreflightBuilder:
    """Determine comparison scope from identity annotations and job metadata."""

    SAME_SYSTEM_MODES = {"sampling", "replicate"}

    def build(
        self,
        case_a: SingleCaseArtifacts,
        case_b: SingleCaseArtifacts,
        comparison_mode: str,
        requested_scope: str = "auto",
        alignment_selection: str = "phla_core_ca",
        residue_mapping: str = "auto",
    ) -> ComparisonPreflightResult:
        resolved_scope = self._resolve_scope(requested_scope, comparison_mode)
        checks: list[dict[str, Any]] = []
        warnings: list[str] = []
        limitations: list[str] = []

        identity_checks = self._identity_checks(case_a, case_b)
        checks.extend(identity_checks)
        missing_identity = [check for check in identity_checks if check["status"] == "missing"]
        failed_identity = [check for check in identity_checks if check["status"] == "failed"]

        prepared_checks = self._prepared_input_checks(case_a, case_b)
        checks.extend(prepared_checks)

        if resolved_scope == "same-system":
            if failed_identity:
                limitations.append("Same-system comparison requested, but identity annotations differ.")
            if missing_identity:
                limitations.append("Same-system comparison cannot be fully verified because identity annotations are incomplete.")
            if not self._has_prepared_input(case_a) or not self._has_prepared_input(case_b):
                warnings.append("Prepared input metadata is unavailable for at least one case; coordinate alignment cannot be audited from job metadata.")
        else:
            limitations.append("Cross-system comparison should prioritize region-level metrics unless residue mapping is explicitly validated.")
            if residue_mapping == "auto":
                warnings.append("Residue-level deltas are descriptive only; explicit residue mapping was not provided.")

        alignment_policy = self._alignment_policy(
            case_a=case_a,
            case_b=case_b,
            resolved_scope=resolved_scope,
            alignment_selection=alignment_selection,
        )
        residue_mapping_policy = self._residue_mapping_policy(
            resolved_scope=resolved_scope,
            residue_mapping=residue_mapping,
            failed_identity=bool(failed_identity),
            missing_identity=bool(missing_identity),
        )

        status = "full"
        if limitations:
            status = "limited"
        if any(check["status"] == "failed" for check in prepared_checks if check.get("severity") == "blocking"):
            status = "blocked"

        return ComparisonPreflightResult(
            requested_scope=requested_scope,
            resolved_scope=resolved_scope,
            status=status,
            alignment_policy=alignment_policy,
            residue_mapping_policy=residue_mapping_policy,
            checks=checks,
            warnings=warnings,
            limitations=limitations,
        )

    def _resolve_scope(self, requested_scope: str, comparison_mode: str) -> str:
        if requested_scope in {"same-system", "cross-system"}:
            return requested_scope
        if comparison_mode in self.SAME_SYSTEM_MODES:
            return "same-system"
        return "cross-system"

    def _identity_checks(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        return [
            self._check_equal(
                name="peptide_sequence",
                label="Peptide sequence",
                value_a=case_a.identity.get("peptide_identity", {}).get("sequence"),
                value_b=case_b.identity.get("peptide_identity", {}).get("sequence"),
            ),
            self._check_equal(
                name="peptide_length",
                label="Peptide length",
                value_a=case_a.identity.get("peptide_identity", {}).get("length"),
                value_b=case_b.identity.get("peptide_identity", {}).get("length"),
            ),
            self._check_equal(
                name="hla_locus",
                label="HLA locus",
                value_a=case_a.identity.get("hla_identity", {}).get("best_locus"),
                value_b=case_b.identity.get("hla_identity", {}).get("best_locus"),
            ),
            self._check_equal(
                name="cdr3_alpha_sequence",
                label="CDR3 alpha sequence",
                value_a=case_a.identity.get("tcr_identity", {}).get("cdr3_alpha_sequence"),
                value_b=case_b.identity.get("tcr_identity", {}).get("cdr3_alpha_sequence"),
            ),
            self._check_equal(
                name="cdr3_beta_sequence",
                label="CDR3 beta sequence",
                value_a=case_a.identity.get("tcr_identity", {}).get("cdr3_beta_sequence"),
                value_b=case_b.identity.get("tcr_identity", {}).get("cdr3_beta_sequence"),
            ),
        ]

    def _prepared_input_checks(self, case_a: SingleCaseArtifacts, case_b: SingleCaseArtifacts) -> list[dict[str, Any]]:
        prepared_a = self._prepared_input(case_a)
        prepared_b = self._prepared_input(case_b)
        checks = []
        for key in ("structure", "topology", "trajectory"):
            checks.append(
                {
                    "name": f"prepared_{key}_available",
                    "label": f"Prepared {key} metadata",
                    "status": "passed" if prepared_a.get(key) and prepared_b.get(key) else "missing",
                    "case_a": prepared_a.get(key, ""),
                    "case_b": prepared_b.get(key, ""),
                    "severity": "informational",
                }
            )
        return checks

    def _alignment_policy(
        self,
        case_a: SingleCaseArtifacts,
        case_b: SingleCaseArtifacts,
        resolved_scope: str,
        alignment_selection: str,
    ) -> dict[str, Any]:
        prepared_a = self._prepared_input(case_a)
        prepared_b = self._prepared_input(case_b)
        reference = prepared_a.get("structure") or case_a.source_paths.get("identity") or ""
        return {
            "required": True,
            "executed_by_compare": False,
            "recommended_reference": reference,
            "recommended_fit_selection": alignment_selection,
            "recommended_analysis_selection": "full pHLA-TCR interface",
            "same_coordinate_system_required": resolved_scope == "same-system",
            "case_a_topology": prepared_a.get("topology", ""),
            "case_b_topology": prepared_b.get("topology", ""),
            "note": "Compare assumes upstream analyses were generated after a consistent alignment/preparation step.",
        }

    def _residue_mapping_policy(
        self,
        resolved_scope: str,
        residue_mapping: str,
        failed_identity: bool,
        missing_identity: bool,
    ) -> dict[str, Any]:
        if resolved_scope == "same-system" and not failed_identity and not missing_identity:
            level = "residue"
        else:
            level = "region"
        if residue_mapping == "region-only":
            level = "region"
        return {
            "requested": residue_mapping,
            "effective_level": level,
            "residue_level_allowed": level == "residue",
            "reason": "Residue-level comparison requires same-system identity checks to pass.",
        }

    def _check_equal(self, name: str, label: str, value_a: Any, value_b: Any) -> dict[str, Any]:
        if value_a in (None, "") or value_b in (None, ""):
            status = "missing"
        elif str(value_a) == str(value_b):
            status = "passed"
        else:
            status = "failed"
        return {
            "name": name,
            "label": label,
            "status": status,
            "case_a": value_a if value_a is not None else "",
            "case_b": value_b if value_b is not None else "",
            "severity": "warning",
        }

    def _has_prepared_input(self, case: SingleCaseArtifacts) -> bool:
        prepared = self._prepared_input(case)
        return bool(prepared.get("structure") and prepared.get("topology") and prepared.get("trajectory"))

    def _prepared_input(self, case: SingleCaseArtifacts) -> dict[str, Any]:
        summary_prepared = case.run_summary.get("prepared_input", {})
        if isinstance(summary_prepared, dict) and summary_prepared:
            return summary_prepared
        manifest_prepared = case.run_manifest.get("prepared_input", {})
        if isinstance(manifest_prepared, dict):
            return manifest_prepared
        return {}
