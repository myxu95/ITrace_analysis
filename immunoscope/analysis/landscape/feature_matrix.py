"""Energy-landscape feature matrix construction."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from ..angles.angle_data_structures import DockingAngleResult
from ..geometry import COMDistanceResult
from ..interface import BuriedSurfaceAreaResult
from ..trajectory import RMSDResult


@dataclass
class FeatureMatrixInput:
    """Feature matrix construction input."""

    contact_result: Any = None
    angle_result: DockingAngleResult | None = None
    com_distance_result: COMDistanceResult | None = None
    bsa_result: BuriedSurfaceAreaResult | None = None
    rmsd_result: RMSDResult | dict[str, RMSDResult] | None = None
    rmsd_results: dict[str, RMSDResult] | None = None
    extra_feature_series: dict[str, Any] | None = None
    hotspot_threshold: float = 0.5
    time_round_decimals: int = 6

    def validate(self) -> None:
        if not any(
            value is not None
            for value in [
                self.contact_result,
                self.angle_result,
                self.com_distance_result,
                self.bsa_result,
                self.rmsd_result,
                self.rmsd_results,
                self.extra_feature_series,
            ]
        ):
            raise ValueError("At least one result object is required to build a feature matrix")
        if not 0.0 <= self.hotspot_threshold <= 1.0:
            raise ValueError("hotspot_threshold must be within [0, 1]")


@dataclass(frozen=True)
class FeatureMatrixResult:
    """Unified feature matrix result."""

    feature_matrix: np.ndarray
    feature_names: list[str]
    times: np.ndarray
    feature_groups: dict[str, list[int]] = field(default_factory=dict)

    def to_dataframe(self) -> pd.DataFrame:
        frame = pd.DataFrame(self.feature_matrix, columns=self.feature_names)
        frame.insert(0, "time_ps", self.times)
        return frame

    def get_feature_groups(self) -> dict[str, list[int]]:
        return dict(self.feature_groups)


class FeatureMatrixBuilder:
    """Integrate multiple analysis outputs into one feature matrix."""

    def build(self, input_params: FeatureMatrixInput) -> FeatureMatrixResult:
        input_params.validate()

        dynamic_features: list[dict[str, Any]] = []
        static_features: list[dict[str, Any]] = []

        self._collect_contact_features(input_params, dynamic_features, static_features)
        self._collect_angle_features(input_params.angle_result, dynamic_features)
        self._collect_com_distance_features(input_params.com_distance_result, dynamic_features)
        self._collect_bsa_features(input_params.bsa_result, dynamic_features)
        self._collect_rmsd_features(input_params, dynamic_features)
        self._collect_extra_features(input_params.extra_feature_series, dynamic_features, static_features)

        common_times = self._resolve_common_times(dynamic_features, input_params.time_round_decimals)
        n_rows = int(common_times.size) if common_times.size else 1

        feature_columns: list[np.ndarray] = []
        feature_names: list[str] = []
        feature_groups: dict[str, list[int]] = {}

        for feature in dynamic_features:
            aligned_values = self._align_feature_to_times(
                feature["times"],
                feature["values"],
                common_times,
                input_params.time_round_decimals,
            )
            self._append_feature(
                feature_columns,
                feature_names,
                feature_groups,
                aligned_values,
                str(feature["name"]),
                str(feature["group"]),
            )

        for feature in static_features:
            broadcast_values = np.repeat(float(feature["value"]), n_rows)
            self._append_feature(
                feature_columns,
                feature_names,
                feature_groups,
                broadcast_values,
                str(feature["name"]),
                str(feature["group"]),
            )

        if not feature_columns:
            raise ValueError("No features could be parsed from the input results")

        feature_matrix = np.column_stack(feature_columns)
        if common_times.size == 0:
            common_times = np.asarray([0.0], dtype=float)

        return FeatureMatrixResult(
            feature_matrix=feature_matrix,
            feature_names=feature_names,
            times=common_times,
            feature_groups=feature_groups,
        )

    def _collect_contact_features(
        self,
        input_params: FeatureMatrixInput,
        dynamic_features: list[dict[str, Any]],
        static_features: list[dict[str, Any]],
    ) -> None:
        contact_result = input_params.contact_result
        if contact_result is None:
            return

        if isinstance(contact_result, pd.DataFrame):
            if self._is_timeseries_feature_frame(contact_result):
                time_column = self._infer_time_column(contact_result)
                for column in contact_result.columns:
                    if column == time_column:
                        continue
                    if pd.api.types.is_numeric_dtype(contact_result[column]):
                        dynamic_features.append(
                            {
                                "name": self._sanitize_feature_name(column),
                                "group": "contact",
                                "times": contact_result[time_column].to_numpy(dtype=float),
                                "values": contact_result[column].to_numpy(dtype=float),
                            }
                        )
                return

            self._collect_contact_summary_features(
                contact_result,
                static_features,
                hotspot_threshold=input_params.hotspot_threshold,
            )
            return

        if hasattr(contact_result, "to_dataframe") and callable(contact_result.to_dataframe):
            frame = contact_result.to_dataframe()
            self._collect_contact_features(
                FeatureMatrixInput(
                    contact_result=frame,
                    hotspot_threshold=input_params.hotspot_threshold,
                ),
                dynamic_features,
                static_features,
            )

    def _collect_contact_summary_features(
        self,
        contact_frame: pd.DataFrame,
        static_features: list[dict[str, Any]],
        hotspot_threshold: float,
    ) -> None:
        if "tcr_region_detailed" in contact_frame.columns and "contact_frequency" in contact_frame.columns:
            region_summary = (
                contact_frame.groupby("tcr_region_detailed", as_index=False)["contact_frequency"]
                .sum()
                .sort_values("tcr_region_detailed")
            )
            for _, row in region_summary.iterrows():
                static_features.append(
                    {
                        "name": f"contact_{self._sanitize_feature_name(str(row['tcr_region_detailed']))}_frequency_sum",
                        "group": "contact",
                        "value": float(row["contact_frequency"]),
                    }
                )

        if "phla_region" in contact_frame.columns and "contact_frequency" in contact_frame.columns:
            phla_summary = (
                contact_frame.groupby("phla_region", as_index=False)["contact_frequency"]
                .sum()
                .sort_values("phla_region")
            )
            for _, row in phla_summary.iterrows():
                static_features.append(
                    {
                        "name": f"contact_{self._sanitize_feature_name(str(row['phla_region']))}_frequency_sum",
                        "group": "contact",
                        "value": float(row["contact_frequency"]),
                    }
                )

        if "contact_frequency" not in contact_frame.columns:
            return

        distance_column = None
        for candidate in ["min_distance_observed", "min_distance_angstrom"]:
            if candidate in contact_frame.columns:
                distance_column = candidate
                break

        if distance_column is None:
            return

        hotspot_frame = contact_frame[contact_frame["contact_frequency"] >= hotspot_threshold].copy()
        if hotspot_frame.empty:
            return

        for _, row in hotspot_frame.iterrows():
            residue_1 = row.get("residue_label_1") or row.get("phla_residue") or row.get("resid_1") or "res1"
            residue_2 = row.get("residue_label_2") or row.get("tcr_residue") or row.get("resid_2") or "res2"
            static_features.append(
                {
                    "name": (
                        f"hotspot_{self._sanitize_feature_name(str(residue_1))}_"
                        f"{self._sanitize_feature_name(str(residue_2))}_min_distance"
                    ),
                    "group": "contact",
                    "value": float(row[distance_column]),
                }
            )

    def _collect_angle_features(
        self,
        angle_result: DockingAngleResult | None,
        dynamic_features: list[dict[str, Any]],
    ) -> None:
        if not angle_result or angle_result.times is None:
            return
        if angle_result.crossing_angles is not None:
            dynamic_features.append(
                {
                    "name": "crossing_angle",
                    "group": "geometry",
                    "times": np.asarray(angle_result.times, dtype=float),
                    "values": np.asarray(angle_result.crossing_angles, dtype=float),
                }
            )
        if angle_result.incident_angles is not None:
            dynamic_features.append(
                {
                    "name": "incident_angle",
                    "group": "geometry",
                    "times": np.asarray(angle_result.times, dtype=float),
                    "values": np.asarray(angle_result.incident_angles, dtype=float),
                }
            )

    def _collect_com_distance_features(
        self,
        com_distance_result: COMDistanceResult | None,
        dynamic_features: list[dict[str, Any]],
    ) -> None:
        if com_distance_result is None:
            return
        dynamic_features.append(
            {
                "name": "com_distance",
                "group": "geometry",
                "times": np.asarray(com_distance_result.times, dtype=float),
                "values": np.asarray(com_distance_result.distances, dtype=float),
            }
        )

    def _collect_bsa_features(
        self,
        bsa_result: BuriedSurfaceAreaResult | None,
        dynamic_features: list[dict[str, Any]],
    ) -> None:
        if bsa_result is None:
            return

        dynamic_features.append(
            {
                "name": "buried_surface_area",
                "group": "surface",
                "times": np.asarray(bsa_result.times, dtype=float),
                "values": np.asarray(bsa_result.buried_surface_area, dtype=float),
            }
        )
        dynamic_features.append(
            {
                "name": "interface_ratio",
                "group": "surface",
                "times": np.asarray(bsa_result.times, dtype=float),
                "values": np.asarray(bsa_result.interface_ratio, dtype=float),
            }
        )

        peptide_sasa = getattr(bsa_result, "peptide_exposed_sasa", None)
        if peptide_sasa is not None:
            dynamic_features.append(
                {
                    "name": "peptide_exposed_sasa",
                    "group": "surface",
                    "times": np.asarray(bsa_result.times, dtype=float),
                    "values": np.asarray(peptide_sasa, dtype=float),
                }
            )

    def _collect_rmsd_features(
        self,
        input_params: FeatureMatrixInput,
        dynamic_features: list[dict[str, Any]],
    ) -> None:
        rmsd_mapping: dict[str, RMSDResult] = {}
        if isinstance(input_params.rmsd_result, dict):
            rmsd_mapping.update(input_params.rmsd_result)
        elif input_params.rmsd_result is not None:
            rmsd_mapping["rmsd"] = input_params.rmsd_result
        if input_params.rmsd_results:
            rmsd_mapping.update(input_params.rmsd_results)

        for name, result in rmsd_mapping.items():
            if not result or result.times is None or result.rmsd_values is None:
                continue
            dynamic_features.append(
                {
                    "name": self._sanitize_feature_name(name),
                    "group": "conformation",
                    "times": np.asarray(result.times, dtype=float),
                    "values": np.asarray(result.rmsd_values, dtype=float),
                }
            )

    def _collect_extra_features(
        self,
        extra_feature_series: dict[str, Any] | None,
        dynamic_features: list[dict[str, Any]],
        static_features: list[dict[str, Any]],
    ) -> None:
        if not extra_feature_series:
            return

        for name, payload in extra_feature_series.items():
            feature_name = self._sanitize_feature_name(name)
            if isinstance(payload, dict):
                values = payload.get("values")
                times = payload.get("times")
                group = str(payload.get("group", "extra"))
                if values is None:
                    continue
                if times is None:
                    static_features.append(
                        {
                            "name": feature_name,
                            "group": group,
                            "value": float(np.asarray(values, dtype=float).reshape(-1)[0]),
                        }
                    )
                else:
                    dynamic_features.append(
                        {
                            "name": feature_name,
                            "group": group,
                            "times": np.asarray(times, dtype=float),
                            "values": np.asarray(values, dtype=float),
                        }
                    )
                continue

            values_array = np.asarray(payload, dtype=float)
            if values_array.ndim == 0 or values_array.size == 1:
                static_features.append(
                    {
                        "name": feature_name,
                        "group": "extra",
                        "value": float(values_array.reshape(-1)[0]),
                    }
                )

    @staticmethod
    def _append_feature(
        feature_columns: list[np.ndarray],
        feature_names: list[str],
        feature_groups: dict[str, list[int]],
        values: np.ndarray,
        feature_name: str,
        group_name: str,
    ) -> None:
        feature_columns.append(np.asarray(values, dtype=float))
        feature_names.append(feature_name)
        feature_groups.setdefault(group_name, []).append(len(feature_names) - 1)

    @staticmethod
    def _sanitize_feature_name(raw_name: str) -> str:
        sanitized = raw_name.strip().replace(" ", "_").replace("-", "_")
        while "__" in sanitized:
            sanitized = sanitized.replace("__", "_")
        return sanitized.lower()

    @staticmethod
    def _is_timeseries_feature_frame(frame: pd.DataFrame) -> bool:
        time_column = FeatureMatrixBuilder._infer_time_column(frame)
        if time_column is None:
            return False
        numeric_columns = [
            column
            for column in frame.columns
            if column != time_column and pd.api.types.is_numeric_dtype(frame[column])
        ]
        return bool(numeric_columns)

    @staticmethod
    def _infer_time_column(frame: pd.DataFrame) -> str | None:
        for candidate in ["time_ps", "Time(ps)", "time", "frame_time"]:
            if candidate in frame.columns:
                return candidate
        return None

    @staticmethod
    def _round_times(times: np.ndarray, decimals: int) -> np.ndarray:
        return np.round(np.asarray(times, dtype=float), decimals=decimals)

    def _resolve_common_times(
        self,
        dynamic_features: list[dict[str, Any]],
        decimals: int,
    ) -> np.ndarray:
        if not dynamic_features:
            return np.asarray([], dtype=float)

        rounded_sets = [self._round_times(feature["times"], decimals) for feature in dynamic_features]
        common_times = rounded_sets[0]
        for rounded_times in rounded_sets[1:]:
            common_times = np.intersect1d(common_times, rounded_times)
        return np.asarray(common_times, dtype=float)

    def _align_feature_to_times(
        self,
        times: np.ndarray,
        values: np.ndarray,
        common_times: np.ndarray,
        decimals: int,
    ) -> np.ndarray:
        time_array = np.asarray(times, dtype=float)
        value_array = np.asarray(values, dtype=float)
        if common_times.size == 0:
            return value_array

        rounded_time_array = self._round_times(time_array, decimals)
        index_map = {float(time_value): idx for idx, time_value in enumerate(rounded_time_array)}
        aligned_indices = [index_map[float(time_value)] for time_value in common_times if float(time_value) in index_map]
        if len(aligned_indices) != common_times.size:
            raise ValueError("Feature time series cannot be aligned to the common time axis")
        return value_array[np.asarray(aligned_indices, dtype=int)]
