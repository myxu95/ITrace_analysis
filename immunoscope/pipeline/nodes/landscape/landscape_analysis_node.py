"""Energy-landscape analysis node."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ....analysis.angles.angle_data_structures import DockingAngleResult
from ....analysis.geometry import COMDistanceResult
from ....analysis.interface import BuriedSurfaceAreaResult
from ....analysis.landscape import (
    FeatureMatrixBuilder,
    FeatureMatrixInput,
    LandscapeAnalyzer,
    LandscapeInput,
    LandscapeVisualizer,
)
from ....analysis.trajectory import RMSDResult
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class LandscapeAnalysisNode(PipelineNode):
    """Integrate upstream results and generate feature matrices, embeddings, and free-energy landscapes."""

    def __init__(
        self,
        angle_result_key: str = "docking_angles",
        com_distance_result_key: str = "com_distance",
        bsa_result_key: str = "bsa",
        contact_result_keys: tuple[str, ...] = ("contact_annotation", "contact_frequency"),
        rmsd_result_keys: Optional[dict[str, str]] = None,
        hotspot_threshold: float = 0.5,
        reducer: str = "pca",
        n_components: int = 3,
        umap_n_neighbors: int = 15,
        umap_min_dist: float = 0.1,
        umap_metric: str = "euclidean",
        random_state: int | None = 42,
        tica_lag: int = 10,
        tica_regularization: float = 1e-6,
        tica_kinetic_map: bool = True,
        tica_stabilize_sign: bool = True,
        temperature: float = 300.0,
        bins: int | str | tuple[int, int] = "auto",
        include_interactive: bool = True,
        output_subdir: str = "landscape",
        result_key: str = "landscape",
        name: Optional[str] = None,
    ):
        super().__init__(name=name or "LandscapeAnalysisNode")
        self.angle_result_key = angle_result_key
        self.com_distance_result_key = com_distance_result_key
        self.bsa_result_key = bsa_result_key
        self.contact_result_keys = contact_result_keys
        self.rmsd_result_keys = rmsd_result_keys or {}
        self.hotspot_threshold = hotspot_threshold
        self.reducer = reducer
        self.n_components = n_components
        self.umap_n_neighbors = umap_n_neighbors
        self.umap_min_dist = umap_min_dist
        self.umap_metric = umap_metric
        self.random_state = random_state
        self.tica_lag = tica_lag
        self.tica_regularization = tica_regularization
        self.tica_kinetic_map = tica_kinetic_map
        self.tica_stabilize_sign = tica_stabilize_sign
        self.temperature = temperature
        self.bins = bins
        self.include_interactive = include_interactive
        self.output_subdir = output_subdir
        self.result_key = result_key

    def validate_inputs(self, context: PipelineContext):
        available = []
        if self.angle_result_key in context.results:
            available.append(self.angle_result_key)
        if self.com_distance_result_key in context.results:
            available.append(self.com_distance_result_key)
        if self.bsa_result_key in context.results:
            available.append(self.bsa_result_key)
        for key in self.contact_result_keys:
            if key in context.results:
                available.append(key)
                break
        for _, result_key in self.rmsd_result_keys.items():
            if result_key in context.results:
                available.append(result_key)

        if not available:
            raise PipelineError(
                node_name=self.name,
                reason="No upstream analysis results found for landscape construction",
                context_state={"system_id": context.system_id, "result_keys": list(context.results.keys())},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        try:
            contact_result = self._load_contact_result(context)
            angle_result = self._load_angle_result(context)
            com_distance_result = self._load_com_distance_result(context)
            bsa_result = self._load_bsa_result(context)
            rmsd_results = self._load_rmsd_results(context)

            feature_builder = FeatureMatrixBuilder()
            feature_result = feature_builder.build(
                FeatureMatrixInput(
                    contact_result=contact_result,
                    angle_result=angle_result,
                    com_distance_result=com_distance_result,
                    bsa_result=bsa_result,
                    rmsd_results=rmsd_results,
                    hotspot_threshold=self.hotspot_threshold,
                )
            )

            analyzer = LandscapeAnalyzer()
            landscape_result = analyzer.analyze(
                LandscapeInput(
                    feature_matrix=feature_result.feature_matrix,
                    feature_names=feature_result.feature_names,
                    times=feature_result.times,
                    reducer=self.reducer,
                    n_components=self.n_components,
                    umap_n_neighbors=self.umap_n_neighbors,
                    umap_min_dist=self.umap_min_dist,
                    umap_metric=self.umap_metric,
                    random_state=self.random_state,
                    tica_lag=self.tica_lag,
                    tica_regularization=self.tica_regularization,
                    tica_kinetic_map=self.tica_kinetic_map,
                    tica_stabilize_sign=self.tica_stabilize_sign,
                    temperature=self.temperature,
                    bins=self.bins,
                )
            )

            output_dir = Path(context.get_analysis_path(self.output_subdir, "feature_matrix.csv")).parent
            visualizer = LandscapeVisualizer(landscape_result)
            artifact_paths = visualizer.write_report_assets(
                output_dir=output_dir,
                feature_result=feature_result,
                include_interactive=self.include_interactive,
            )

            summary_path = artifact_paths.get("summary")
            summary = {}
            if summary_path and Path(summary_path).exists():
                summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))

            payload = {
                "output_dir": str(output_dir),
                **artifact_paths,
                "feature_names": feature_result.feature_names,
                "feature_groups": feature_result.get_feature_groups(),
                "n_frames": int(feature_result.feature_matrix.shape[0]),
                "n_features": int(feature_result.feature_matrix.shape[1]),
                "reducer": landscape_result.reducer,
                "coordinate_labels": landscape_result.coordinate_labels or [],
                "explained_variance": [float(value) for value in landscape_result.explained_variance],
                "reducer_metadata": landscape_result.reducer_metadata or {},
                "top_contributors_pc1": summary.get("top_contributors_pc1", []),
                "top_contributors_pc2": summary.get("top_contributors_pc2", []),
            }
            context.add_result(self.result_key, payload)
            context.metadata.setdefault("landscape", {}).update(
                {
                    "feature_names": feature_result.feature_names,
                    "output_dir": str(output_dir),
                    "reducer": landscape_result.reducer,
                    "explained_variance": payload["explained_variance"],
                    "reducer_metadata": payload["reducer_metadata"],
                }
            )
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"Landscape analysis failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc

    def _load_contact_result(self, context: PipelineContext) -> pd.DataFrame | None:
        for result_key in self.contact_result_keys:
            payload = context.get_result(result_key)
            if not payload:
                continue
            candidate_keys = ["contact_report_file", "annotated_contacts_file", "raw_contacts_file"]
            for candidate_key in candidate_keys:
                file_path = payload.get(candidate_key)
                if file_path and Path(file_path).exists():
                    return pd.read_csv(file_path)
        return None

    def _load_angle_result(self, context: PipelineContext) -> DockingAngleResult | None:
        payload = context.get_result(self.angle_result_key)
        if not payload:
            return None

        output_files = payload.get("output_files") or []
        csv_path = None
        for output_file in output_files:
            if str(output_file).endswith("docking_angles.csv") and Path(output_file).exists():
                csv_path = Path(output_file)
                break
        if csv_path is None:
            return None

        frame = pd.read_csv(csv_path)
        required_columns = ["Time(ps)", "Crossing(deg)", "Incident(deg)"]
        if not set(required_columns).issubset(frame.columns):
            return None

        return DockingAngleResult(
            success=True,
            times=frame["Time(ps)"].to_numpy(dtype=float),
            crossing_angles=frame["Crossing(deg)"].to_numpy(dtype=float),
            incident_angles=frame["Incident(deg)"].to_numpy(dtype=float),
            statistics=payload.get("statistics"),
            output_files=list(output_files),
            metadata=payload.get("metadata") or {},
        )

    def _load_com_distance_result(self, context: PipelineContext) -> COMDistanceResult | None:
        payload = context.get_result(self.com_distance_result_key)
        if not payload:
            return None

        timeseries_file = payload.get("timeseries_file")
        if not timeseries_file or not Path(timeseries_file).exists():
            return None

        frame = pd.read_csv(timeseries_file)
        time_column = next((column for column in frame.columns if column.startswith("time_")), None)
        if time_column is None or "com_distance_angstrom" not in frame.columns:
            return None

        distances = frame["com_distance_angstrom"].to_numpy(dtype=float)
        return COMDistanceResult(
            times=frame[time_column].to_numpy(dtype=float),
            distances=distances,
            mean=float(np.mean(distances)),
            std=float(np.std(distances)),
            min=float(np.min(distances)),
            max=float(np.max(distances)),
            tcr_selection=str(payload.get("tcr_selection", "")),
            mhc_selection=str(payload.get("mhc_selection", "")),
            stride=int(payload.get("stride", 1)),
            time_unit=str(payload.get("time_unit", time_column.replace("time_", ""))),
        )

    def _load_bsa_result(self, context: PipelineContext) -> BuriedSurfaceAreaResult | None:
        payload = context.get_result(self.bsa_result_key)
        if not payload:
            return None

        timeseries_file = payload.get("timeseries_file")
        if not timeseries_file or not Path(timeseries_file).exists():
            return None

        frame = pd.read_csv(timeseries_file)
        time_column = next((column for column in frame.columns if column.startswith("time_")), None)
        required_columns = [
            time_column,
            "buried_surface_area",
            "interface_ratio",
            "sasa_a",
            "sasa_b",
            "sasa_ab",
        ]
        if time_column is None or not all(column in frame.columns for column in required_columns if column):
            return None

        return BuriedSurfaceAreaResult(
            times=frame[time_column].to_numpy(dtype=float),
            sasa_a=frame["sasa_a"].to_numpy(dtype=float),
            sasa_b=frame["sasa_b"].to_numpy(dtype=float),
            sasa_ab=frame["sasa_ab"].to_numpy(dtype=float),
            buried_surface_area=frame["buried_surface_area"].to_numpy(dtype=float),
            interface_ratio=frame["interface_ratio"].to_numpy(dtype=float),
            selection_a=str(payload.get("selection_a", "")),
            selection_b=str(payload.get("selection_b", "")),
            probe_radius=float(payload.get("probe_radius", 1.4)),
            stride=int(payload.get("stride", 1)),
            time_unit=str(payload.get("time_unit", time_column.replace("time_", ""))),
        )

    def _load_rmsd_results(self, context: PipelineContext) -> dict[str, RMSDResult]:
        rmsd_results: dict[str, RMSDResult] = {}
        for feature_name, result_key in self.rmsd_result_keys.items():
            payload = context.get_result(result_key)
            if not payload:
                continue

            output_file = payload.get("output_file")
            if not output_file or not Path(output_file).exists():
                continue

            result = self._parse_rmsd_file(Path(output_file))
            if result is not None:
                rmsd_results[feature_name] = result

        return rmsd_results

    def _parse_rmsd_file(self, filepath: Path) -> RMSDResult | None:
        if filepath.suffix.lower() == ".csv":
            frame = pd.read_csv(filepath)
            time_column = next((column for column in frame.columns if "time" in column.lower()), None)
            value_column = next(
                (column for column in frame.columns if "rmsd" in column.lower() and column != time_column),
                None,
            )
            if time_column is None or value_column is None:
                return None
            times = frame[time_column].to_numpy(dtype=float)
            values = frame[value_column].to_numpy(dtype=float)
            return self._build_rmsd_result(filepath, times, values)

        times = []
        values = []
        with filepath.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line or line.startswith("#") or line.startswith("@"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    times.append(float(parts[0]))
                    values.append(float(parts[1]))
                except ValueError:
                    continue

        if not times:
            return None
        return self._build_rmsd_result(filepath, np.asarray(times, dtype=float), np.asarray(values, dtype=float))

    @staticmethod
    def _build_rmsd_result(filepath: Path, times: np.ndarray, values: np.ndarray) -> RMSDResult:
        return RMSDResult(
            success=True,
            output_file=str(filepath),
            times=times,
            rmsd_values=values,
            mean_rmsd=float(np.mean(values)),
            std_rmsd=float(np.std(values)),
            min_rmsd=float(np.min(values)),
            max_rmsd=float(np.max(values)),
            n_frames=int(values.size),
            metadata={"source_file": str(filepath)},
        )
