"""
Pipeline node for residue-residue contact frequency calculation.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

from ....analysis import ResidueContactFrequencyAnalyzer
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class ContactFrequencyNode(PipelineNode):
    """Compute residue-level pHLA-TCR contact frequencies."""

    def __init__(
        self,
        cutoff: float = 4.5,
        stride: int = 1,
        min_frequency: float = 0.0,
        heavy_atoms_only: bool = True,
        name: str | None = None,
    ):
        super().__init__(name=name or "ContactFrequencyNode")
        self.cutoff = cutoff
        self.stride = stride
        self.min_frequency = min_frequency
        self.heavy_atoms_only = heavy_atoms_only

    def validate_inputs(self, context: PipelineContext):
        missing = []
        if not context.topology and not context.structure_pdb:
            missing.append("topology_or_structure_pdb")
        if not (context.trajectory_processed or context.trajectory_raw):
            missing.append("trajectory_processed_or_raw")
        if "chain_mapping" not in context.metadata:
            missing.append("chain_mapping")
        if missing:
            raise PipelineError(
                node_name=self.name,
                reason=f"Missing required inputs: {missing}",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        chain_mapping = context.metadata["chain_mapping"]
        topology_file = context.structure_pdb or context.topology
        trajectory_file = context.trajectory_processed or context.trajectory_raw

        required_chains = ["mhc_alpha", "b2m", "peptide", "tcr_alpha", "tcr_beta"]
        missing_chains = [chain for chain in required_chains if chain not in chain_mapping]
        if missing_chains:
            raise PipelineError(
                node_name=self.name,
                reason=f"chain_mapping missing required assignments: {missing_chains}",
                context_state={"system_id": context.system_id, "chain_mapping": chain_mapping},
            )

        try:
            analyzer = ResidueContactFrequencyAnalyzer(topology_file, trajectory_file)
            selection1 = self._build_selection(
                chain_mapping["mhc_alpha"],
                chain_mapping["b2m"],
                chain_mapping["peptide"],
            )
            selection2 = self._build_selection(
                chain_mapping["tcr_alpha"],
                chain_mapping["tcr_beta"],
            )

            contacts = analyzer.calculate_residue_contact_frequencies(
                selection1=selection1,
                selection2=selection2,
                cutoff=self.cutoff,
                stride=self.stride,
                heavy_atoms_only=self.heavy_atoms_only,
                min_frequency=self.min_frequency,
            )

            output_file = context.get_analysis_path("contacts", "residue_contact_frequencies.csv")
            summary_file = context.get_analysis_path("contacts", "residue_contact_summary.json")
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            contacts.to_csv(output_file, index=False)

            summary = {
                "topology": topology_file,
                "trajectory": trajectory_file,
                "selection_1": selection1,
                "selection_2": selection2,
                "cutoff_angstrom": self.cutoff,
                "stride": self.stride,
                "heavy_atoms_only": self.heavy_atoms_only,
                "min_frequency": self.min_frequency,
                "n_contact_pairs": int(len(contacts)),
            }
            with open(summary_file, "w", encoding="utf-8") as handle:
                json.dump(summary, handle, indent=2)

            context.results["contact_frequency"] = {
                "raw_contacts_file": output_file,
                "summary_file": summary_file,
                "selection_1": selection1,
                "selection_2": selection2,
                "n_contact_pairs": int(len(contacts)),
                "cutoff_angstrom": self.cutoff,
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"Contact frequency calculation failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc

    @staticmethod
    def _build_selection(*chain_ids: str) -> str:
        return " or ".join(f"chainID {chain_id}" for chain_id in chain_ids if chain_id)
