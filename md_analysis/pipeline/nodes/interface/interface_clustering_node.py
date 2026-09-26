"""Key interface-state clustering node."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from ....analysis.conformation import (
    InterfaceClusteringAnalyzer,
    write_cluster_id_vs_time,
    write_cluster_population_over_time,
)
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class InterfaceClusteringNode(PipelineNode):
    """Run key interface-state clustering."""

    def __init__(
        self,
        stride: int = 10,
        contact_cutoff_angstrom: float = 4.5,
        distance_cutoff: float = 0.35,
        linkage_method: str = "average",
        geometry_weight: float = 0.45,
        sidechain_weight: float = 0.25,
        interaction_weight: float = 0.30,
        name: Optional[str] = None,
    ):
        super().__init__(name=name or "InterfaceClusteringNode")
        self.stride = stride
        self.contact_cutoff_angstrom = contact_cutoff_angstrom
        self.distance_cutoff = distance_cutoff
        self.linkage_method = linkage_method
        self.geometry_weight = geometry_weight
        self.sidechain_weight = sidechain_weight
        self.interaction_weight = interaction_weight

    def validate_inputs(self, context: PipelineContext) -> None:
        missing = []
        if not context.topology:
            missing.append("topology")
        if not (context.trajectory_processed or context.trajectory_raw):
            missing.append("trajectory_processed/trajectory_raw")
        if not context.structure_pdb:
            missing.append("structure_pdb")
        if "chain_mapping" not in context.metadata:
            missing.append("chain_mapping")
        if "cdr_detection" not in context.metadata:
            missing.append("cdr_detection")
        if missing:
            raise PipelineError(
                node_name=self.name,
                reason=f"Missing required inputs: {missing}",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        try:
            analyzer = InterfaceClusteringAnalyzer(
                topology_file=context.topology,
                trajectory_file=context.trajectory_processed or context.trajectory_raw,
                structure_pdb=context.structure_pdb,
            )
            result = analyzer.calculate(
                chain_mapping=context.metadata["chain_mapping"],
                cdr_detection=context.metadata["cdr_detection"],
                stride=self.stride,
                contact_cutoff_angstrom=self.contact_cutoff_angstrom,
                distance_cutoff=self.distance_cutoff,
                linkage_method=self.linkage_method,
                geometry_weight=self.geometry_weight,
                sidechain_weight=self.sidechain_weight,
                interaction_weight=self.interaction_weight,
            )

            analysis_dir = Path(
                context.get_analysis_path("conformation/interface_clustering", "frame_cluster_assignments.csv")
            ).parent
            analysis_dir.mkdir(parents=True, exist_ok=True)
            self._cleanup_legacy_outputs(analysis_dir)

            frame_csv = analysis_dir / "frame_cluster_assignments.csv"
            summary_table_csv = analysis_dir / "summary_table.csv"
            difference_csv = analysis_dir / "cluster_difference_matrix.csv"
            feature_csv = analysis_dir / "cluster_feature_summary.csv"
            signature_csv = analysis_dir / "cluster_specific_signatures.csv"
            digest_csv = analysis_dir / "cluster_feature_digest.csv"
            summary_json = analysis_dir / "interface_clustering_summary.json"
            cluster_id_png = analysis_dir / "cluster_id_vs_time.png"
            population_png = analysis_dir / "state_population_over_time.png"
            structures_dir = analysis_dir / "structures"

            result.frame_assignments.to_csv(frame_csv, index=False)
            result.cluster_summary.to_csv(summary_table_csv, index=False)
            result.cluster_difference_matrix.to_csv(difference_csv, index_label="cluster")
            result.cluster_feature_summary.to_csv(feature_csv, index=False)
            result.cluster_signatures.to_csv(signature_csv, index=False)
            result.cluster_feature_digest.to_csv(digest_csv, index=False)
            summary_json.write_text(json.dumps(result.summary, ensure_ascii=False, indent=2), encoding="utf-8")

            structure_rows = analyzer.export_cluster_structures(
                chain_mapping=context.metadata["chain_mapping"],
                cdr_detection=context.metadata["cdr_detection"],
                frame_assignments=result.frame_assignments,
                cluster_summary=result.cluster_summary,
                stride=self.stride,
                output_dir=structures_dir,
            )
            write_cluster_id_vs_time(result.frame_assignments, cluster_id_png)
            write_cluster_population_over_time(result.frame_assignments, population_png)

            context.results["interface_clustering"] = {
                "frame_assignments_csv": str(frame_csv),
                "summary_table_csv": str(summary_table_csv),
                "cluster_difference_matrix_csv": str(difference_csv),
                "cluster_feature_summary_csv": str(feature_csv),
                "cluster_signatures_csv": str(signature_csv),
                "cluster_feature_digest_csv": str(digest_csv),
                "summary_json": str(summary_json),
                "cluster_id_vs_time_plot": str(cluster_id_png),
                "state_population_plot": str(population_png),
                "structures": structure_rows,
                "structures_dir": str(structures_dir),
                "distance_cutoff": self.distance_cutoff,
                "linkage_method": self.linkage_method,
                "contact_stride": int(result.summary.get("contact_stride", self.stride)),
                "weights": {
                    "geometry": self.geometry_weight,
                    "sidechain": self.sidechain_weight,
                    "interaction": self.interaction_weight,
                },
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"Interface clustering failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc

    @staticmethod
    def _cleanup_legacy_outputs(analysis_dir: Path) -> None:
        legacy_names = [
            "cluster_dendrogram.png",
            "cluster_population.png",
            "cluster_transition_timeline.png",
            "cluster_interface_signature.png",
            "cluster_summary.csv",
            "cluster_interface_signature.csv",
            "cluster_top_contacts.csv",
            "joint_distance_matrix.npz",
            "geometry_distance_matrix.npz",
            "sidechain_distance_matrix.npz",
            "interaction_distance_matrix.npz",
        ]
        for name in legacy_names:
            path = analysis_dir / name
            if path.exists():
                path.unlink()
