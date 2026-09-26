"""Lightweight internal TICA implementation used by landscape analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class InternalTICAResult:
    coordinates: np.ndarray
    eigenvalues: np.ndarray
    vectors: np.ndarray
    correlations: np.ndarray
    metadata: dict[str, Any]


def run_internal_tica(
    scaled_matrix: np.ndarray,
    lag: int,
    n_components: int,
    regularization: float = 1e-6,
    kinetic_map: bool = True,
    stabilize_sign: bool = True,
) -> InternalTICAResult:
    """Run a compact TICA projection on an already standardized feature matrix."""

    matrix = np.asarray(scaled_matrix, dtype=float)
    n_samples, n_features = matrix.shape
    if lag < 1:
        raise ValueError("tica_lag must be >= 1")
    if lag >= n_samples:
        raise ValueError(f"tica_lag={lag} must be smaller than n_frames={n_samples}")
    if regularization < 0:
        raise ValueError("tica_regularization must be >= 0")

    x0 = matrix[:-lag] - np.mean(matrix[:-lag], axis=0, keepdims=True)
    xt = matrix[lag:] - np.mean(matrix[lag:], axis=0, keepdims=True)

    n_pairs = x0.shape[0]
    denominator = max(n_pairs - 1, 1)
    c00 = (x0.T @ x0) / denominator
    ctt = (xt.T @ xt) / denominator
    c0t = (x0.T @ xt) / denominator
    covariance = 0.5 * (c00 + ctt)
    time_lagged_covariance = 0.5 * (c0t + c0t.T)

    if regularization > 0:
        covariance = covariance + np.eye(n_features, dtype=float) * float(regularization)

    eigenvalues, eigenvectors = solve_symmetric_generalized_eigenproblem(
        time_lagged_covariance,
        covariance,
    )
    order = rank_tica_eigenvalues(eigenvalues)
    component_count = min(int(n_components), eigenvectors.shape[1])
    selected = order[:component_count]

    selected_vectors = np.asarray(eigenvectors[:, selected], dtype=float)
    selected_values = np.asarray(eigenvalues[selected], dtype=float)
    if kinetic_map:
        selected_vectors = selected_vectors * selected_values.reshape(1, -1)

    coordinates = np.asarray(matrix @ selected_vectors, dtype=float)
    coordinates, selected_vectors = stabilize_tica_signs(
        matrix,
        coordinates,
        selected_vectors,
        enabled=stabilize_sign,
    )
    correlations = feature_embedding_correlations(matrix, coordinates)
    metadata = {
        "tica_lag": int(lag),
        "tica_regularization": float(regularization),
        "tica_kinetic_map": bool(kinetic_map),
        "tica_eigenvalues": [float(value) for value in selected_values],
        "tica_implied_timescales_frames": tica_implied_timescales(selected_values, int(lag)),
    }
    return InternalTICAResult(
        coordinates=coordinates,
        eigenvalues=selected_values,
        vectors=selected_vectors,
        correlations=correlations,
        metadata=metadata,
    )


def rank_tica_eigenvalues(eigenvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(eigenvalues, dtype=float)
    positive_indices = np.where(values > 0.0)[0]
    non_positive_indices = np.where(values <= 0.0)[0]
    positive_order = positive_indices[np.argsort(values[positive_indices])[::-1]]
    non_positive_order = non_positive_indices[np.argsort(np.abs(values[non_positive_indices]))[::-1]]
    return np.concatenate([positive_order, non_positive_order])


def stabilize_tica_signs(
    scaled_matrix: np.ndarray,
    coordinates: np.ndarray,
    vectors: np.ndarray,
    enabled: bool,
) -> tuple[np.ndarray, np.ndarray]:
    if not enabled:
        return coordinates, vectors

    adjusted_coordinates = np.asarray(coordinates, dtype=float).copy()
    adjusted_vectors = np.asarray(vectors, dtype=float).copy()
    correlations = feature_embedding_correlations(scaled_matrix, adjusted_coordinates)
    for coord_index in range(adjusted_coordinates.shape[1]):
        if correlations.shape[0] == 0:
            continue
        anchor_index = int(np.argmax(np.abs(correlations[:, coord_index])))
        if correlations[anchor_index, coord_index] < 0:
            adjusted_coordinates[:, coord_index] *= -1.0
            adjusted_vectors[:, coord_index] *= -1.0
    return adjusted_coordinates, adjusted_vectors


def tica_implied_timescales(eigenvalues: np.ndarray, lag: int) -> list[float | None]:
    timescales: list[float | None] = []
    for value in np.asarray(eigenvalues, dtype=float):
        if 0.0 < value < 1.0:
            timescales.append(float(-float(lag) / np.log(value)))
        else:
            timescales.append(None)
    return timescales


def solve_symmetric_generalized_eigenproblem(
    lhs: np.ndarray,
    rhs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    try:
        from scipy.linalg import eigh

        eigenvalues, eigenvectors = eigh(lhs, rhs)
        return np.asarray(eigenvalues, dtype=float), np.asarray(eigenvectors, dtype=float)
    except Exception:
        rhs_values, rhs_vectors = np.linalg.eigh(rhs)
        positive_mask = rhs_values > 1e-12
        if int(np.sum(positive_mask)) < 2:
            raise ValueError("TICA covariance matrix is rank-deficient")

        retained_vectors = rhs_vectors[:, positive_mask]
        retained_values = rhs_values[positive_mask]
        whitening = retained_vectors @ np.diag(1.0 / np.sqrt(retained_values))
        whitened_lhs = whitening.T @ lhs @ whitening
        eigenvalues, whitened_vectors = np.linalg.eigh(whitened_lhs)
        eigenvectors = whitening @ whitened_vectors
        return np.asarray(eigenvalues, dtype=float), np.asarray(eigenvectors, dtype=float)


def feature_embedding_correlations(
    scaled_matrix: np.ndarray,
    coordinates: np.ndarray,
) -> np.ndarray:
    correlations = np.zeros((scaled_matrix.shape[1], coordinates.shape[1]), dtype=float)
    for feature_index in range(scaled_matrix.shape[1]):
        feature_values = scaled_matrix[:, feature_index]
        if np.std(feature_values) == 0:
            continue
        for coord_index in range(coordinates.shape[1]):
            coord_values = coordinates[:, coord_index]
            if np.std(coord_values) == 0:
                continue
            correlations[feature_index, coord_index] = float(np.corrcoef(feature_values, coord_values)[0, 1])
    return correlations
