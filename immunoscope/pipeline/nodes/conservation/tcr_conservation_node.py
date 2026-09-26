"""TCR germline-conservation pipeline node.

Sequence-only (no trajectory): extracts the TCR alpha/beta chain
sequences from the structure PDB, scores each residue against the IMGT
germline V-gene reference, and writes
`analysis/conservation/conservation_tcr.csv`.

Degrades gracefully: if ANARCI / its germline data is unavailable, or the
chain mapping does not name the TCR chains, it writes a header-only CSV
and records the gap rather than failing the pipeline.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from ....analysis.conservation.structure_conservation import (
    write_conservation_csv,
    write_empty_conservation_csv,
)
from ....analysis.conservation.tcr_germline_conservation import NoAnarciError
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class TcrConservationNode(PipelineNode):
    """Per-residue TCR germline conservation (FR / CDR1 / CDR2)."""

    def __init__(self, species: str = "human", name: Optional[str] = None):
        super().__init__(name=name or "TcrConservationNode")
        self.species = species

    def validate_inputs(self, context: PipelineContext):
        if not (context.structure_pdb or context.topology):
            raise PipelineError(
                node_name=self.name,
                reason="Missing structure_pdb / topology for sequence extraction",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        pdb = context.structure_pdb or context.topology
        out_csv = context.get_analysis_path("conservation", "conservation_tcr.csv")

        chain_mapping = context.metadata.get("chain_mapping", {}) or {}
        tcr_alpha = chain_mapping.get("tcr_alpha")
        tcr_beta = chain_mapping.get("tcr_beta")

        if not tcr_alpha and not tcr_beta:
            csv_file = write_empty_conservation_csv(out_csv)
            context.results["tcr_conservation"] = {
                "conservation_csv": csv_file,
                "status": "skipped",
                "reason": "chain_mapping has no tcr_alpha/tcr_beta",
            }
            self.logger.warning("conservation skipped: no TCR chains in mapping")
            return context

        try:
            csv_file, rows = write_conservation_csv(
                str(pdb), tcr_alpha, tcr_beta, out_csv, species=self.species
            )
        except NoAnarciError as exc:
            csv_file = write_empty_conservation_csv(out_csv)
            context.results["tcr_conservation"] = {
                "conservation_csv": csv_file,
                "status": "unavailable",
                "reason": f"ANARCI unavailable: {exc}",
            }
            self.logger.warning("conservation unavailable (ANARCI): %s", exc)
            return context
        except Exception as exc:  # noqa: BLE001
            raise PipelineError(
                node_name=self.name,
                reason=f"TCR conservation analysis failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc

        covered = sum(1 for r in rows if r.get("covered") is True)
        context.results["tcr_conservation"] = {
            "conservation_csv": csv_file,
            "status": "ok",
            "n_residues": len(rows),
            "n_covered": covered,
            "tcr_alpha_chain": tcr_alpha,
            "tcr_beta_chain": tcr_beta,
            "species": self.species,
        }
        return context


__all__ = ["TcrConservationNode"]
