"""RRCS analysis node."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import MDAnalysis as mda
import pandas as pd

from ....analysis.interactions.rrcs import (
    RRCSAnalyzer,
    RRCSPairSpec,
    build_residue_heavy_atom_indices,
    parse_rrcs_pair_file,
    residue_chain_id,
)
from ....analysis.topology import ResiduePairAnnotationAnnotator
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class RRCSNode(PipelineNode):
    """Run RRCS calculation for key interface residue pairs."""

    def __init__(
        self,
        radius_min: float = 3.23,
        radius_max: float = 4.63,
        stride: int = 1,
        pair_scope: str = "interface",
        pair_file: str | None = None,
        name: Optional[str] = None,
    ):
        super().__init__(name=name or "RRCSNode")
        self.radius_min = radius_min
        self.radius_max = radius_max
        self.stride = stride
        self.pair_scope = pair_scope
        self.pair_file = pair_file

    def validate_inputs(self, context: PipelineContext) -> None:
        missing = []
        if not context.structure_pdb:
            missing.append("structure_pdb")
        if not (context.topology or context.structure_pdb):
            missing.append("topology_or_structure_pdb")
        if not (context.trajectory_processed or context.trajectory_raw):
            missing.append("trajectory_processed_or_raw")
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
            topology_file = context.topology or context.structure_pdb
            trajectory_file = context.trajectory_processed or context.trajectory_raw
            structure_pdb = context.structure_pdb
            chain_mapping = context.metadata["chain_mapping"]
            cdr_detection = context.metadata["cdr_detection"]

            structure_universe = mda.Universe(structure_pdb)
            pair_specs = self._build_pair_specs(
                structure_universe=structure_universe,
                chain_mapping=chain_mapping,
                cdr_detection=cdr_detection,
            )
            analyzer = RRCSAnalyzer(topology=topology_file, trajectory=trajectory_file)
            if analyzer.universe.atoms.n_atoms != structure_universe.atoms.n_atoms:
                raise ValueError(
                    "structure_pdb and topology/trajectory atom counts differ; "
                    "stable RRCS residue-pair mapping cannot be built"
                )
            timeseries = analyzer.calculate_pair_timeseries(
                pair_specs=pair_specs,
                radius_min=self.radius_min,
                radius_max=self.radius_max,
                stride=self.stride,
            )
            pair_summary = analyzer.summarize_pair_timeseries(timeseries)

            annotator = ResiduePairAnnotationAnnotator()
            annotated_summary = annotator.annotate_pairs(
                pair_summary,
                chain_mapping=chain_mapping,
                cdr_detection=cdr_detection,
            )
            region_summary = self._build_region_summary(annotated_summary)
            summary = self._build_summary_payload(
                context=context,
                pair_summary=pair_summary,
                region_summary=region_summary,
            )

            analysis_dir = Path(context.get_analysis_path("interactions/rrcs", "rrcs_timeseries.csv")).parent
            analysis_dir.mkdir(parents=True, exist_ok=True)

            timeseries_file = analysis_dir / "rrcs_timeseries.csv"
            pair_summary_file = analysis_dir / "rrcs_pair_summary.csv"
            annotated_file = analysis_dir / "annotated_rrcs_pair_summary.csv"
            region_summary_file = analysis_dir / "rrcs_region_summary.csv"
            summary_file = analysis_dir / "rrcs_summary.json"

            timeseries.to_csv(timeseries_file, index=False)
            pair_summary.to_csv(pair_summary_file, index=False)
            annotated_summary.to_csv(annotated_file, index=False)
            region_summary.to_csv(region_summary_file, index=False)
            summary_file.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

            context.results["rrcs"] = {
                "timeseries_file": str(timeseries_file),
                "pair_summary_file": str(pair_summary_file),
                "annotated_pair_summary_file": str(annotated_file),
                "region_summary_file": str(region_summary_file),
                "summary_file": str(summary_file),
                "radius_min": self.radius_min,
                "radius_max": self.radius_max,
                "stride": self.stride,
                "pair_scope": self.pair_scope,
                "pair_file": self.pair_file or "",
                "n_pairs": int(len(pair_summary)),
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"RRCS calculation failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc

    def _build_pair_specs(
        self,
        structure_universe: mda.Universe,
        chain_mapping: dict[str, str],
        cdr_detection: dict,
    ) -> list[RRCSPairSpec]:
        if self.pair_file:
            return self._build_pair_specs_from_file(
                structure_universe=structure_universe,
                pair_file=self.pair_file,
            )
        return self._build_default_pair_specs(
            structure_universe=structure_universe,
            chain_mapping=chain_mapping,
            cdr_detection=cdr_detection,
        )

    def _build_default_pair_specs(
        self,
        structure_universe: mda.Universe,
        chain_mapping: dict[str, str],
        cdr_detection: dict,
    ) -> list[RRCSPairSpec]:
        mhc_alpha_chain = chain_mapping.get("mhc_alpha")
        peptide_chain = chain_mapping.get("peptide")
        tcr_alpha_chain = chain_mapping.get("tcr_alpha")
        tcr_beta_chain = chain_mapping.get("tcr_beta")
        if not all([mhc_alpha_chain, peptide_chain, tcr_alpha_chain, tcr_beta_chain]):
            raise ValueError("chain_mapping is missing mhc_alpha/peptide/tcr_alpha/tcr_beta")

        residue_lookup = {
            (residue_chain_id(residue), int(residue.resid)): residue
            for residue in structure_universe.residues
        }
        groove_residue_keys = {
            (residue_chain_id(residue), int(residue.resid))
            for residue in structure_universe.residues
            if residue_chain_id(residue) == mhc_alpha_chain and (50 <= int(residue.resid) <= 86 or 140 <= int(residue.resid) <= 176)
        }
        peptide_residue_keys = {
            (residue_chain_id(residue), int(residue.resid))
            for residue in structure_universe.residues
            if residue_chain_id(residue) == peptide_chain
        }
        cdr_residue_keys = self._collect_tcr_cdr_keys(
            cdr_detection=cdr_detection,
            chain_mapping=chain_mapping,
        )
        all_tcr_keys = {
            (residue_chain_id(residue), int(residue.resid))
            for residue in structure_universe.residues
            if residue_chain_id(residue) in {tcr_alpha_chain, tcr_beta_chain}
        }

        source_keys, target_keys = self._resolve_scope_keys(
            pair_scope=self.pair_scope,
            cdr_residue_keys=cdr_residue_keys,
            peptide_residue_keys=peptide_residue_keys,
            groove_residue_keys=groove_residue_keys,
            all_tcr_keys=all_tcr_keys,
        )
        source_residues = [residue_lookup[key] for key in sorted(source_keys) if key in residue_lookup]
        target_residues = [residue_lookup[key] for key in sorted(target_keys) if key in residue_lookup]
        pair_specs: list[RRCSPairSpec] = []
        for source_residue in source_residues:
            atom_indices_1, atom_indices_1_no_backbone = build_residue_heavy_atom_indices(source_residue)
            if len(atom_indices_1) == 0:
                continue
            for target_residue in target_residues:
                atom_indices_2, atom_indices_2_no_backbone = build_residue_heavy_atom_indices(target_residue)
                if len(atom_indices_2) == 0:
                    continue
                pair_specs.append(
                    RRCSPairSpec(
                        chain_id_1=residue_chain_id(source_residue),
                        resid_1=int(source_residue.resid),
                        resname_1=str(source_residue.resname),
                        residue_label_1=f"{source_residue.resname}{int(source_residue.resid)}",
                        chain_id_2=residue_chain_id(target_residue),
                        resid_2=int(target_residue.resid),
                        resname_2=str(target_residue.resname),
                        residue_label_2=f"{target_residue.resname}{int(target_residue.resid)}",
                        atom_indices_1=atom_indices_1,
                        atom_indices_2=atom_indices_2,
                        atom_indices_1_no_backbone=atom_indices_1_no_backbone,
                        atom_indices_2_no_backbone=atom_indices_2_no_backbone,
                        selection_scope=self.pair_scope,
                    )
                )

        if not pair_specs:
            raise ValueError(
                f"Default RRCS pairs are empty; check pair_scope={self.pair_scope}, "
                "chain mapping, and CDR annotation"
            )
        return pair_specs

    def _build_pair_specs_from_file(
        self,
        structure_universe: mda.Universe,
        pair_file: str,
    ) -> list[RRCSPairSpec]:
        residue_lookup = {
            (residue_chain_id(residue), int(residue.resid)): residue
            for residue in structure_universe.residues
        }
        custom_pairs = parse_rrcs_pair_file(pair_file)
        pair_specs: list[RRCSPairSpec] = []
        missing_keys: list[str] = []

        for chain_1, resid_1, chain_2, resid_2 in custom_pairs:
            residue_1 = residue_lookup.get((chain_1, resid_1))
            residue_2 = residue_lookup.get((chain_2, resid_2))
            if residue_1 is None:
                missing_keys.append(f"{chain_1}:{resid_1}")
                continue
            if residue_2 is None:
                missing_keys.append(f"{chain_2}:{resid_2}")
                continue
            atom_indices_1, atom_indices_1_no_backbone = build_residue_heavy_atom_indices(residue_1)
            atom_indices_2, atom_indices_2_no_backbone = build_residue_heavy_atom_indices(residue_2)
            if len(atom_indices_1) == 0 or len(atom_indices_2) == 0:
                continue
            pair_specs.append(
                RRCSPairSpec(
                    chain_id_1=residue_chain_id(residue_1),
                    resid_1=int(residue_1.resid),
                    resname_1=str(residue_1.resname),
                    residue_label_1=f"{residue_1.resname}{int(residue_1.resid)}",
                    chain_id_2=residue_chain_id(residue_2),
                    resid_2=int(residue_2.resid),
                    resname_2=str(residue_2.resname),
                    residue_label_2=f"{residue_2.resname}{int(residue_2.resid)}",
                    atom_indices_1=atom_indices_1,
                    atom_indices_2=atom_indices_2,
                    atom_indices_1_no_backbone=atom_indices_1_no_backbone,
                    atom_indices_2_no_backbone=atom_indices_2_no_backbone,
                    selection_scope="custom_pair_file",
                )
            )

        if not pair_specs:
            detail = f"; missing residue examples: {', '.join(missing_keys[:10])}" if missing_keys else ""
            raise ValueError(f"RRCS custom pair file produced no valid pairs{detail}")
        return pair_specs

    @staticmethod
    def _collect_tcr_cdr_keys(
        cdr_detection: dict,
        chain_mapping: dict[str, str],
    ) -> dict[str, set[tuple[str, int]]]:
        region_keys = {"CDR1": set(), "CDR2": set(), "CDR3": set()}
        chain_id_map = {
            "TCR_alpha": chain_mapping.get("tcr_alpha"),
            "TCR_beta": chain_mapping.get("tcr_beta"),
        }
        for chain_name, chain_info in (cdr_detection.get("chains") or {}).items():
            chain_id = chain_info.get("chain_id") or chain_id_map.get(chain_name)
            if not chain_id:
                continue
            for cdr_num, cdr_meta in (chain_info.get("cdrs") or {}).items():
                residue_range = cdr_meta.get("residue_range")
                if not residue_range:
                    continue
                region_name = f"CDR{cdr_num}"
                if region_name not in region_keys:
                    continue
                start, end = int(residue_range[0]), int(residue_range[1])
                for resid in range(start, end + 1):
                    region_keys[region_name].add((chain_id, resid))
        return region_keys

    @staticmethod
    def _resolve_scope_keys(
        pair_scope: str,
        cdr_residue_keys: dict[str, set[tuple[str, int]]],
        peptide_residue_keys: set[tuple[str, int]],
        groove_residue_keys: set[tuple[str, int]],
        all_tcr_keys: set[tuple[str, int]],
    ) -> tuple[set[tuple[str, int]], set[tuple[str, int]]]:
        cdr_all = set().union(*cdr_residue_keys.values()) if cdr_residue_keys else set()
        cdr3 = set(cdr_residue_keys.get("CDR3", set()))
        if pair_scope == "interface":
            return cdr_all, peptide_residue_keys | groove_residue_keys
        if pair_scope == "cdr3_peptide":
            return cdr3, peptide_residue_keys
        if pair_scope == "cdr3_groove":
            return cdr3, groove_residue_keys
        if pair_scope == "tcr_peptide":
            return all_tcr_keys, peptide_residue_keys
        if pair_scope == "tcr_groove":
            return all_tcr_keys, groove_residue_keys
        if pair_scope == "tcr_interface":
            return all_tcr_keys, peptide_residue_keys | groove_residue_keys
        raise ValueError(f"Unsupported RRCS pair_scope: {pair_scope}")

    @staticmethod
    def _build_region_summary(annotated_summary: pd.DataFrame) -> pd.DataFrame:
        if annotated_summary.empty:
            return pd.DataFrame(
                columns=[
                    "interaction_class",
                    "tcr_chain",
                    "tcr_region",
                    "partner_component",
                    "mhc_subregion",
                    "n_pairs",
                    "mean_rrcs_sum",
                    "mean_rrcs_mean",
                    "median_rrcs_mean",
                    "max_rrcs_max",
                    "rrcs_nonzero_fraction_mean",
                ]
            )

        region_summary = (
            annotated_summary.groupby(
                ["interaction_class", "tcr_chain", "tcr_region", "partner_component", "mhc_subregion"],
                dropna=False,
                as_index=False,
            )
            .agg(
                n_pairs=("mean_rrcs", "size"),
                mean_rrcs_sum=("mean_rrcs", "sum"),
                mean_rrcs_mean=("mean_rrcs", "mean"),
                median_rrcs_mean=("median_rrcs", "mean"),
                max_rrcs_max=("max_rrcs", "max"),
                rrcs_nonzero_fraction_mean=("rrcs_nonzero_fraction", "mean"),
            )
        )
        return region_summary.sort_values(
            by=["mean_rrcs_sum", "mean_rrcs_mean", "n_pairs"],
            ascending=[False, False, False],
        ).reset_index(drop=True)

    def _build_summary_payload(
        self,
        context: PipelineContext,
        pair_summary: pd.DataFrame,
        region_summary: pd.DataFrame,
    ) -> dict[str, object]:
        top_pairs = []
        for _, row in pair_summary.head(10).iterrows():
            top_pairs.append(
                {
                    "pair": f"{row['chain_id_1']}:{int(row['resid_1'])}-{row['chain_id_2']}:{int(row['resid_2'])}",
                    "residues": f"{row['resname_1']}{int(row['resid_1'])} ↔ {row['resname_2']}{int(row['resid_2'])}",
                    "mean_rrcs": float(row["mean_rrcs"]),
                    "max_rrcs": float(row["max_rrcs"]),
                    "nonzero_fraction": float(row["rrcs_nonzero_fraction"]),
                }
            )

        top_regions = []
        for _, row in region_summary.head(10).iterrows():
            top_regions.append(
                {
                    "interaction_class": row["interaction_class"],
                    "tcr_region": row["tcr_region"],
                    "partner_component": row["partner_component"],
                    "mhc_subregion": row["mhc_subregion"],
                    "mean_rrcs_sum": float(row["mean_rrcs_sum"]),
                    "n_pairs": int(row["n_pairs"]),
                }
            )

        return {
            "system_id": context.system_id,
            "topology": context.topology or context.structure_pdb,
            "trajectory": context.trajectory_processed or context.trajectory_raw,
            "structure_pdb": context.structure_pdb,
            "radius_min_angstrom": self.radius_min,
            "radius_max_angstrom": self.radius_max,
            "stride": self.stride,
            "pair_scope": self.pair_scope,
            "pair_file": self.pair_file or "",
            "n_pairs": int(len(pair_summary)),
            "n_nonzero_pairs": int((pair_summary["rrcs_nonzero_frames"] > 0).sum()) if not pair_summary.empty else 0,
            "top_pairs": top_pairs,
            "top_regions": top_regions,
        }
