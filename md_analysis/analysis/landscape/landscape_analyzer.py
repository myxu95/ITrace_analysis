"""Energy-landscape dimensionality reduction and free-energy calculation."""

from __future__ import annotations

from dataclasses import dataclass, field
import os
import tempfile
from typing import Any

import numpy as np
import pandas as pd

from ..free_energy import FELCalculator
from .tica import feature_embedding_correlations, run_internal_tica


@dataclass
class LandscapeInput:
    """Input parameters for energy-landscape analysis."""

    feature_matrix: np.ndarray
    feature_names: list[str]
    times: np.ndarray | None = None
    reducer: str = "pca"
    n_components: int = 2
    umap_n_neighbors: int = 15
    umap_min_dist: float = 0.1
    umap_metric: str = "euclidean"
    random_state: int | None = 42
    tica_lag: int = 10
    tica_regularization: float = 1e-6
    tica_kinetic_map: bool = True
    tica_stabilize_sign: bool = True
    temperature: float = 300.0
    bins: int | str | tuple[int, int] = "auto"
    smooth: bool = False
    sigma: float = 1.0
    vmax: float | None = None

    def validate(self) -> None:
        matrix = np.asarray(self.feature_matrix, dtype=float)
        if matrix.ndim != 2:
            raise ValueError("feature_matrix must be a 2D array")
        if matrix.shape[0] < 2:
            raise ValueError("At least 2 frames are required for energy-landscape analysis")
        if matrix.shape[1] < 2:
            raise ValueError("At least 2 features are required for dimensionality reduction")
        if len(self.feature_names) != matrix.shape[1]:
            raise ValueError("feature_names length must match the feature column count")
        reducer = self.reducer.lower()
        if reducer not in {"pca", "umap", "tica"}:
            raise ValueError("reducer must be pca, umap, or tica")
        if self.n_components < 2:
            raise ValueError("n_components must be >= 2")
        max_components = min(matrix.shape[0], matrix.shape[1])
        if reducer == "pca" and self.n_components > max_components:
            raise ValueError(
                f"n_components={self.n_components} exceeds the allowed maximum {max_components}"
            )
        if reducer == "umap" and self.n_components > matrix.shape[0]:
            raise ValueError(
                f"n_components={self.n_components} exceeds the UMAP maximum {matrix.shape[0]}"
            )
        if reducer == "tica":
            if self.tica_lag < 1:
                raise ValueError("tica_lag must be >= 1")
            if self.tica_lag >= matrix.shape[0]:
                raise ValueError(f"tica_lag={self.tica_lag} must be smaller than the frame count {matrix.shape[0]}")
            if self.n_components > max_components:
                raise ValueError(
                    f"n_components={self.n_components} exceeds the TICA maximum {max_components}"
                )
            if self.tica_regularization < 0:
                raise ValueError("tica_regularization must be >= 0")
        if self.umap_n_neighbors < 2:
            raise ValueError("umap_n_neighbors must be >= 2")
        if not 0.0 <= self.umap_min_dist <= 1.0:
            raise ValueError("umap_min_dist must be in the [0, 1] interval")
        if self.temperature <= 0:
            raise ValueError("temperature must be greater than 0")
        if isinstance(self.bins, int) and self.bins < 2:
            raise ValueError("integer bins must be >= 2")


@dataclass(frozen=True)
class LandscapeResult:
    """Energy-landscape analysis result."""

    pc_coordinates: np.ndarray
    explained_variance: np.ndarray
    loadings: np.ndarray
    free_energy: np.ndarray
    pc1_edges: np.ndarray
    pc2_edges: np.ndarray
    pc1_centers: np.ndarray
    pc2_centers: np.ndarray
    feature_names: list[str]
    n_components: int
    reducer: str = "pca"
    coordinate_labels: list[str] | None = None
    times: np.ndarray | None = None
    probability: np.ndarray | None = None
    histogram: np.ndarray | None = None
    reducer_metadata: dict[str, Any] = field(default_factory=dict)

    def to_pca_dataframe(self) -> pd.DataFrame:
        labels = self.coordinate_labels or [f"PC{i + 1}" for i in range(self.pc_coordinates.shape[1])]
        frame = pd.DataFrame(
            self.pc_coordinates,
            columns=labels,
        )
        if self.times is not None:
            frame.insert(0, "time_ps", np.asarray(self.times, dtype=float))
        return frame

    def to_loadings_dataframe(self) -> pd.DataFrame:
        labels = self.coordinate_labels or [f"PC{i + 1}" for i in range(self.loadings.shape[1])]
        return pd.DataFrame(
            self.loadings,
            index=self.feature_names,
            columns=labels[: self.loadings.shape[1]],
        ).reset_index(names="feature")

    def get_top_contributors(self, pc_index: int, top_n: int = 5) -> list[tuple[str, float]]:
        if pc_index < 0 or pc_index >= self.loadings.shape[1]:
            raise ValueError("pc_index is out of range")
        component = self.loadings[:, pc_index]
        ranked_indices = np.argsort(np.abs(component))[::-1][:top_n]
        return [(self.feature_names[idx], float(component[idx])) for idx in ranked_indices]


class LandscapeAnalyzer:
    """Run feature standardization, PCA/UMAP/TICA projection, and FEL calculation."""

    def __init__(self):
        self.pca = None
        self.scaler = None
        self.fel_calculator: FELCalculator | None = None

    def analyze(self, input_params: LandscapeInput) -> LandscapeResult:
        input_params.validate()

        try:
            from sklearn.preprocessing import StandardScaler
        except ImportError as exc:  # pragma: no cover - runtime dependency
            raise ImportError("LandscapeAnalyzer requires scikit-learn. Install it before running this module.") from exc

        feature_matrix = np.asarray(input_params.feature_matrix, dtype=float)
        self.scaler = StandardScaler()
        scaled_matrix = self.scaler.fit_transform(feature_matrix)

        reducer = input_params.reducer.lower()
        reducer_metadata: dict[str, Any] = {}
        if reducer == "pca":
            pc_coordinates, explained_variance, loadings, coordinate_labels = self._run_pca(
                scaled_matrix,
                input_params.n_components,
            )
        elif reducer == "umap":
            pc_coordinates, explained_variance, loadings, coordinate_labels = self._run_umap(
                scaled_matrix,
                input_params,
            )
        else:
            (
                pc_coordinates,
                explained_variance,
                loadings,
                coordinate_labels,
                reducer_metadata,
            ) = self._run_tica(scaled_matrix, input_params)

        self.fel_calculator = FELCalculator(temperature=input_params.temperature)
        fel_result: dict[str, Any] = self.fel_calculator.compute_2d_fel(
            pc_coordinates[:, 0],
            pc_coordinates[:, 1],
            bins=input_params.bins,
            smooth=input_params.smooth,
            sigma=input_params.sigma,
            vmax=input_params.vmax,
        )

        return LandscapeResult(
            pc_coordinates=pc_coordinates,
            explained_variance=explained_variance,
            loadings=loadings,
            free_energy=np.asarray(fel_result["free_energy"], dtype=float),
            pc1_edges=np.asarray(fel_result["cv1_edges"], dtype=float),
            pc2_edges=np.asarray(fel_result["cv2_edges"], dtype=float),
            pc1_centers=np.asarray(fel_result["cv1_centers"], dtype=float),
            pc2_centers=np.asarray(fel_result["cv2_centers"], dtype=float),
            feature_names=list(input_params.feature_names),
            n_components=int(input_params.n_components),
            reducer=reducer,
            coordinate_labels=coordinate_labels,
            times=None if input_params.times is None else np.asarray(input_params.times, dtype=float),
            probability=np.asarray(fel_result.get("probability"), dtype=float),
            histogram=np.asarray(fel_result.get("histogram"), dtype=float),
            reducer_metadata=reducer_metadata,
        )

    def _run_pca(
        self,
        scaled_matrix: np.ndarray,
        n_components: int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
        try:
            from sklearn.decomposition import PCA
        except ImportError as exc:  # pragma: no cover - runtime dependency
            raise ImportError("PCA reduction requires scikit-learn. Install it before running this module.") from exc

        self.pca = PCA(n_components=n_components)
        coordinates = self.pca.fit_transform(scaled_matrix)
        labels = [f"PC{i + 1}" for i in range(coordinates.shape[1])]
        return (
            coordinates,
            np.asarray(self.pca.explained_variance_ratio_, dtype=float),
            np.asarray(self.pca.components_.T, dtype=float),
            labels,
        )

    def _run_umap(
        self,
        scaled_matrix: np.ndarray,
        input_params: LandscapeInput,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
        os.environ.setdefault("NUMBA_CACHE_DIR", str(PathLikeTemp.numba_cache_dir()))
        os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
        try:
            import umap
        except ImportError as exc:  # pragma: no cover - runtime dependency
            raise ImportError("UMAP reduction requires umap-learn. Install with: conda install -c conda-forge umap-learn") from exc

        n_samples = int(scaled_matrix.shape[0])
        n_neighbors = min(int(input_params.umap_n_neighbors), max(n_samples - 1, 2))
        reducer = umap.UMAP(
            n_components=int(input_params.n_components),
            n_neighbors=n_neighbors,
            min_dist=float(input_params.umap_min_dist),
            metric=input_params.umap_metric,
            random_state=input_params.random_state,
        )
        coordinates = reducer.fit_transform(scaled_matrix)
        labels = [f"UMAP{i + 1}" for i in range(coordinates.shape[1])]
        correlations = self._feature_embedding_correlations(scaled_matrix, coordinates)
        return (
            np.asarray(coordinates, dtype=float),
            np.asarray([], dtype=float),
            correlations,
            labels,
        )

    def _run_tica(
        self,
        scaled_matrix: np.ndarray,
        input_params: LandscapeInput,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], dict[str, Any]]:
        """Implement lightweight TICA via a generalized covariance eigenproblem.

        TICA targets lagged slow variables without adding heavy dependencies.
        Coordinates are used for FEL; feature interpretation uses feature-TIC
        correlations, which are more suitable for report interpretation.
        """

        tica_result = run_internal_tica(
            scaled_matrix=scaled_matrix,
            lag=int(input_params.tica_lag),
            n_components=int(input_params.n_components),
            regularization=float(input_params.tica_regularization),
            kinetic_map=bool(input_params.tica_kinetic_map),
            stabilize_sign=bool(input_params.tica_stabilize_sign),
        )

        labels = [f"TIC{i + 1}" for i in range(tica_result.coordinates.shape[1])]
        return (
            tica_result.coordinates,
            np.asarray([], dtype=float),
            tica_result.correlations,
            labels,
            tica_result.metadata,
        )

    @staticmethod
    def _feature_embedding_correlations(
        scaled_matrix: np.ndarray,
        coordinates: np.ndarray,
    ) -> np.ndarray:
        return feature_embedding_correlations(scaled_matrix, coordinates)


class PathLikeTemp:
    """Small helper to keep optional numba cache under a writable temp path."""

    @staticmethod
    def numba_cache_dir() -> str:
        path = os.path.join(tempfile.gettempdir(), "itrace_numba_cache")
        os.makedirs(path, exist_ok=True)
        return path
