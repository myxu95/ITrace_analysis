"""Biological identity annotation node."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from ....analysis.topology import BiologicalIdentityAnnotator
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class BiologicalIdentityNode(PipelineNode):
    """Generate basic biological identity annotation for a single system."""

    def __init__(self, cdr_metadata_path: Optional[str] = None, name: Optional[str] = None):
        super().__init__(name=name or "BiologicalIdentityNode")
        self.cdr_metadata_path = cdr_metadata_path

    def validate_inputs(self, context: PipelineContext) -> None:
        if not context.structure_pdb:
            raise PipelineError(
                node_name=self.name,
                reason="structure_pdb not found in context",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        try:
            annotator = BiologicalIdentityAnnotator()
            cdr_metadata_path = self._resolve_cdr_metadata_path(context)
            payload = annotator.annotate(
                structure_pdb=Path(context.structure_pdb),
                cdr_metadata_path=Path(cdr_metadata_path) if cdr_metadata_path else None,
                chain_mapping=context.metadata.get("chain_mapping"),
                cdr_detection=context.metadata.get("cdr_detection"),
            )

            output_file = Path(context.get_analysis_path("identity", "biological_identity.json"))
            output_file.parent.mkdir(parents=True, exist_ok=True)
            output_file.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

            context.results["biological_identity"] = {
                "json_file": str(output_file),
                "complex_identity": payload.get("complex_identity", {}),
                "peptide_identity": payload.get("peptide_identity", {}),
                "tcr_identity": payload.get("tcr_identity", {}),
                "hla_identity": payload.get("hla_identity", {}),
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"Biological identity annotation failed: {exc}",
                context_state={"system_id": context.system_id, "structure_pdb": context.structure_pdb},
            ) from exc

    def _resolve_cdr_metadata_path(self, context: PipelineContext) -> Optional[str]:
        if self.cdr_metadata_path:
            return self.cdr_metadata_path
        return context.metadata.get("cdr_detection", {}).get("metadata_file")
