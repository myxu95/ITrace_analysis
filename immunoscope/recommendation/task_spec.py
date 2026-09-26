"""TaskSpec — first-class declaration of the mutation-design task.

D-B7 (2026-05-27): The recommendation engine and baseline harness were both
built assuming a single task — *TCR-side affinity optimization against a fixed
peptide-HLA*. The kernel (MD + spatial annotation + RRCS) is task-agnostic;
only the LLM-facing surface (candidate filter, prompt wording, baseline
selection) carries the TCR-affinity bias. To unblock peptide-presentation
design without forking the kernel, we represent the task as an explicit
field threaded from the caller all the way down to the formatter and prompt
builder.

Taxonomy (subject × objective):
    subject ∈ {tcr, peptide, hla}        — which chain we are mutating
    objective ∈ {affinity, stability, presentation}
                                        — what we are trying to optimize
    target_subinterface ∈ {tcr_pmhc, peptide_hla}
                                        — which interface drives candidate
                                          selection and prompt wording

Only two combinations are wired today (2026-05-27):
    1. (tcr, affinity, tcr_pmhc)        — legacy default, fully implemented
    2. (peptide, presentation, peptide_hla)
                                        — new presentation track (L3+ TODO)

Anything else raises ``NotImplementedError`` at construction time so the
caller hits the error early rather than mid-pipeline.

The dataclass is frozen so a TaskSpec can be safely cached and passed
through multiple layers without anyone mutating it under our feet.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import FrozenSet, Tuple


_VALID_SUBJECTS: FrozenSet[str] = frozenset({"tcr", "peptide", "hla"})
_VALID_OBJECTIVES: FrozenSet[str] = frozenset({"affinity", "stability", "presentation"})
_VALID_SUBINTERFACES: FrozenSet[str] = frozenset({"tcr_pmhc", "peptide_hla"})

# Wired combinations as of 2026-05-27. Adding a new combination requires
# updating:
#   - candidate filter in data_formatter._is_design_candidate_for_task
#   - prompt wording in prompts.build_user_prompt
#   - at least one baseline that can execute the task
#   - this allowlist
# Pinning the allowlist here (not at the call site) means we fail fast on
# unwired combinations even if a downstream layer happens to "work" on them
# by accident.
_WIRED_COMBINATIONS: FrozenSet[Tuple[str, str, str]] = frozenset({
    ("tcr", "affinity", "tcr_pmhc"),
    ("peptide", "presentation", "peptide_hla"),
})


@dataclass(frozen=True)
class TaskSpec:
    """The (subject, objective, target_subinterface) tuple for one design run.

    Construction validates both individual fields and the combination — an
    invalid subject raises ValueError, an unwired combination raises
    NotImplementedError. This lets call sites trust that `task` reaching the
    formatter / prompt / baseline is one of the supported combinations.
    """

    subject: str = "tcr"
    objective: str = "affinity"
    target_subinterface: str = "tcr_pmhc"

    def __post_init__(self) -> None:
        if self.subject not in _VALID_SUBJECTS:
            raise ValueError(
                f"TaskSpec.subject={self.subject!r} not in {sorted(_VALID_SUBJECTS)}"
            )
        if self.objective not in _VALID_OBJECTIVES:
            raise ValueError(
                f"TaskSpec.objective={self.objective!r} not in {sorted(_VALID_OBJECTIVES)}"
            )
        if self.target_subinterface not in _VALID_SUBINTERFACES:
            raise ValueError(
                f"TaskSpec.target_subinterface={self.target_subinterface!r} "
                f"not in {sorted(_VALID_SUBINTERFACES)}"
            )
        combo = (self.subject, self.objective, self.target_subinterface)
        if combo not in _WIRED_COMBINATIONS:
            raise NotImplementedError(
                f"TaskSpec{combo} is not wired yet. Wired combinations: "
                f"{sorted(_WIRED_COMBINATIONS)}. To add a new combination, "
                "update _WIRED_COMBINATIONS *and* the candidate filter / "
                "prompt builder / baseline that should execute it."
            )

    @property
    def slug(self) -> str:
        """Short identifier for filenames, log lines, and the results table."""
        return f"{self.subject}_{self.objective}_{self.target_subinterface}"

    @classmethod
    def tcr_affinity(cls) -> "TaskSpec":
        """Legacy default — TCR-side affinity optimization against pMHC."""
        return cls(subject="tcr", objective="affinity", target_subinterface="tcr_pmhc")

    @classmethod
    def peptide_presentation(cls) -> "TaskSpec":
        """Peptide-side presentation/stability optimization at peptide-HLA."""
        return cls(
            subject="peptide",
            objective="presentation",
            target_subinterface="peptide_hla",
        )


__all__ = ["TaskSpec"]
