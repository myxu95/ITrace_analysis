"""Cluster representative conformation comparison utilities."""

from __future__ import annotations

import csv
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


@dataclass(slots=True)
class ClusterRepresentative:
    case_label: str
    cluster_id: int
    population_percent: float
    frames: int
    representative_frame: int
    descriptor: str
    source_pdb: Path
    report_pdb: Path | None = None


@dataclass(slots=True)
class ClusterCase:
    label: str
    root: Path
    summary_csv: Path
    summary_json: Path | None
    structures_dir: Path
    representatives: list[ClusterRepresentative]


def build_cluster_representative_comparison(
    case_a_root: Path,
    case_b_root: Path,
    case_a_label: str,
    case_b_label: str,
    output_dir: Path,
    top_n: int = 4,
) -> dict[str, Any]:
    """Build a report-ready comparison of cluster representative conformations."""
    output_dir = Path(output_dir)
    cluster_dir = output_dir / "analysis" / "comparison" / "cluster_representatives"
    cluster_dir.mkdir(parents=True, exist_ok=True)

    case_a = _load_cluster_case(case_a_root, case_a_label, top_n=top_n)
    case_b = _load_cluster_case(case_b_root, case_b_label, top_n=top_n)
    if not case_a or not case_b:
        return {
            "status": "missing",
            "reason": "Cluster representative artifacts were not found for both cases.",
            "case_a_found": bool(case_a),
            "case_b_found": bool(case_b),
        }

    _copy_representatives(case_a, cluster_dir / _safe_label(case_a.label))
    _copy_representatives(case_b, cluster_dir / _safe_label(case_b.label))

    rmsd_matrix = _compute_rmsd_matrix(case_a.representatives, case_b.representatives)
    matrix_csv = cluster_dir / "cluster_representative_rmsd_matrix.csv"
    matches_csv = cluster_dir / "cluster_representative_matches.csv"
    enriched_csv = cluster_dir / "sampling_enriched_states.csv"
    summary_json = cluster_dir / "cluster_representative_comparison.json"
    heatmap_png = cluster_dir / "cluster_representative_rmsd_heatmap.png"
    population_shift_png = cluster_dir / "cluster_population_shift.png"

    _write_matrix_csv(case_a, case_b, rmsd_matrix, matrix_csv)
    matches = _build_nearest_matches(case_a, case_b, rmsd_matrix)
    enriched_states = _build_sampling_enriched_states(case_a, case_b, matches)
    viewer_pairs = _write_viewer_pairs(case_a, case_b, matches, cluster_dir / "viewer_pairs")
    ensemble_models = _write_ensemble_models(case_a, case_b, cluster_dir / "ensemble_aligned")
    _write_matches_csv(matches, matches_csv)
    _write_enriched_states_csv(enriched_states, enriched_csv)
    _write_heatmap(case_a, case_b, rmsd_matrix, heatmap_png)
    _write_population_shift_plot(case_a, case_b, matches, population_shift_png)

    payload = {
        "status": "ready",
        "comparison_scope": "same_system_sampling" if _looks_like_sampling_pair(case_a.label, case_b.label) else "generic",
        "case_a": _case_payload(case_a),
        "case_b": _case_payload(case_b),
        "sampling_coverage": _sampling_coverage(case_a, case_b),
        "top_n": top_n,
        "rmsd_matrix_csv": str(matrix_csv.resolve()),
        "matches_csv": str(matches_csv.resolve()),
        "enriched_states_csv": str(enriched_csv.resolve()),
        "heatmap_png": str(heatmap_png.resolve()),
        "population_shift_png": str(population_shift_png.resolve()),
        "matches": matches,
        "enriched_states": enriched_states,
        "viewer_pairs": viewer_pairs,
        "ensemble_models": ensemble_models,
        "dominant_pair": matches[0] if matches else {},
    }
    summary_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    payload["summary_json"] = str(summary_json.resolve())
    return payload


def _load_cluster_case(case_root: Path, label: str, top_n: int) -> ClusterCase | None:
    root = _resolve_cluster_root(Path(case_root))
    if root is None:
        return None
    summary_csv = _first_existing(root / "summary_table.csv", root / "cluster_summary.csv")
    if summary_csv is None:
        return None
    structures_dir = _first_existing(root / "structures", root.parent / "structures")
    if structures_dir is None or not structures_dir.is_dir():
        return None

    frame = pd.read_csv(summary_csv)
    if frame.empty or "cluster_id" not in frame.columns:
        return None
    frame = frame.copy()
    frame["population_percent"] = pd.to_numeric(frame.get("population_percent", 0.0), errors="coerce").fillna(0.0)
    frame["frames"] = pd.to_numeric(frame.get("frames", 0), errors="coerce").fillna(0).astype(int)
    frame = frame.sort_values(["population_percent", "frames"], ascending=[False, False]).head(top_n)

    representatives: list[ClusterRepresentative] = []
    for _, row in frame.iterrows():
        cluster_id = int(row["cluster_id"])
        pdb = structures_dir / f"cluster_{cluster_id}_representative.pdb"
        if not pdb.exists():
            continue
        representatives.append(
            ClusterRepresentative(
                case_label=label,
                cluster_id=cluster_id,
                population_percent=float(row.get("population_percent", 0.0) or 0.0),
                frames=int(row.get("frames", 0) or 0),
                representative_frame=int(row.get("representative_frame", 0) or 0),
                descriptor=str(row.get("main_structural_descriptor", "") or ""),
                source_pdb=pdb,
            )
        )
    if not representatives:
        return None
    return ClusterCase(
        label=label,
        root=root,
        summary_csv=summary_csv,
        summary_json=_first_existing(root / "interface_clustering_summary.json"),
        structures_dir=structures_dir,
        representatives=representatives,
    )


def _resolve_cluster_root(case_root: Path) -> Path | None:
    candidates = [
        case_root / "analysis" / "conformation" / "interface_clustering",
        case_root / "analysis" / "inter_cluster" / "analysis" / "conformation" / "interface_clustering",
        case_root / "overview" / "cluster",
        case_root / "cluster",
    ]
    candidates.extend(sorted((case_root / "report").glob("*/overview/cluster")))
    for candidate in candidates:
        if (candidate / "summary_table.csv").exists() or (candidate / "cluster_summary.csv").exists():
            return candidate
    return None


def _first_existing(*candidates: Path) -> Path | None:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _copy_representatives(case: ClusterCase, target_dir: Path) -> None:
    target_dir.mkdir(parents=True, exist_ok=True)
    for representative in case.representatives:
        target = target_dir / f"cluster_{representative.cluster_id}_representative.pdb"
        shutil.copy2(representative.source_pdb, target)
        representative.report_pdb = target.resolve()


def _compute_rmsd_matrix(case_a: list[ClusterRepresentative], case_b: list[ClusterRepresentative]) -> np.ndarray:
    matrix = np.full((len(case_a), len(case_b)), np.nan, dtype=float)
    for i, left in enumerate(case_a):
        left_ca = _read_ca_coordinates(left.report_pdb or left.source_pdb)
        for j, right in enumerate(case_b):
            right_ca = _read_ca_coordinates(right.report_pdb or right.source_pdb)
            matrix[i, j] = _aligned_rmsd(left_ca, right_ca)
    return matrix


def _read_ca_coordinates(pdb_path: Path) -> np.ndarray:
    coords: list[tuple[float, float, float]] = []
    with Path(pdb_path).open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    if not coords:
        raise ValueError(f"No CA atoms found in {pdb_path}")
    return np.array(coords, dtype=float)


def _read_pdb_coordinates(pdb_path: Path) -> tuple[list[str], np.ndarray, list[int]]:
    lines = Path(pdb_path).read_text(encoding="utf-8", errors="ignore").splitlines()
    coords: list[tuple[float, float, float]] = []
    ca_indices: list[int] = []
    atom_index = 0
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError:
            coords.append((0.0, 0.0, 0.0))
        if line[12:16].strip() == "CA":
            ca_indices.append(atom_index)
        atom_index += 1
    return lines, np.array(coords, dtype=float), ca_indices


def _aligned_rmsd(left: np.ndarray, right: np.ndarray) -> float:
    n_atoms = min(len(left), len(right))
    if n_atoms < 3:
        return float("nan")
    left = left[:n_atoms]
    right = right[:n_atoms]
    rotation = _kabsch_rotation(right, left)
    right_aligned = (right - right.mean(axis=0)) @ rotation + left.mean(axis=0)
    diff = left - right_aligned
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _align_all_coordinates_to_reference(
    reference_ca: np.ndarray,
    mobile_ca: np.ndarray,
    mobile_all: np.ndarray,
) -> np.ndarray:
    n_atoms = min(len(reference_ca), len(mobile_ca))
    if n_atoms < 3:
        return mobile_all.copy()
    reference_ca = reference_ca[:n_atoms]
    mobile_ca = mobile_ca[:n_atoms]
    rotation = _kabsch_rotation(mobile_ca, reference_ca)
    return (mobile_all - mobile_ca.mean(axis=0)) @ rotation + reference_ca.mean(axis=0)


def _write_aligned_pdb(source_pdb: Path, output_pdb: Path, aligned_coordinates: np.ndarray) -> None:
    lines, _, _ = _read_pdb_coordinates(source_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    atom_index = 0
    output_lines: list[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and atom_index < len(aligned_coordinates):
            x, y, z = aligned_coordinates[atom_index]
            output_lines.append(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
            atom_index += 1
        else:
            output_lines.append(line)
    output_pdb.write_text("\n".join(output_lines) + "\n", encoding="utf-8")


def _kabsch_rotation(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    mobile_centered = mobile - mobile.mean(axis=0)
    reference_centered = reference - reference.mean(axis=0)
    covariance = mobile_centered.T @ reference_centered
    left, _, right_t = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right_t))
    return left @ correction @ right_t


def _write_matrix_csv(case_a: ClusterCase, case_b: ClusterCase, matrix: np.ndarray, output_file: Path) -> None:
    columns = [f"{case_b.label}_cluster_{item.cluster_id}" for item in case_b.representatives]
    rows = []
    for idx, representative in enumerate(case_a.representatives):
        row = {"cluster": f"{case_a.label}_cluster_{representative.cluster_id}"}
        row.update({column: matrix[idx, j] for j, column in enumerate(columns)})
        rows.append(row)
    pd.DataFrame(rows).to_csv(output_file, index=False)


def _build_nearest_matches(case_a: ClusterCase, case_b: ClusterCase, matrix: np.ndarray) -> list[dict[str, Any]]:
    matches: list[dict[str, Any]] = []
    for i, left in enumerate(case_a.representatives):
        if matrix.shape[1] == 0 or np.all(np.isnan(matrix[i])):
            continue
        j = int(np.nanargmin(matrix[i]))
        right = case_b.representatives[j]
        matches.append(
            {
                "case_a_cluster": int(left.cluster_id),
                "case_a_population_percent": float(left.population_percent),
                "case_a_representative_frame": int(left.representative_frame),
                "case_a_descriptor": left.descriptor,
                "case_b_cluster": int(right.cluster_id),
                "case_b_population_percent": float(right.population_percent),
                "case_b_representative_frame": int(right.representative_frame),
                "case_b_descriptor": right.descriptor,
                "representative_ca_rmsd_angstrom": float(matrix[i, j]),
                "population_delta_b_minus_a_percent": float(right.population_percent - left.population_percent),
            }
        )
    return sorted(matches, key=lambda item: item["case_a_population_percent"], reverse=True)


def _write_viewer_pairs(
    case_a: ClusterCase,
    case_b: ClusterCase,
    matches: list[dict[str, Any]],
    output_dir: Path,
) -> list[dict[str, Any]]:
    pairs: list[dict[str, Any]] = []
    left_by_cluster = {item.cluster_id: item for item in case_a.representatives}
    right_by_cluster = {item.cluster_id: item for item in case_b.representatives}
    for match in matches:
        left = left_by_cluster.get(int(match["case_a_cluster"]))
        right = right_by_cluster.get(int(match["case_b_cluster"]))
        if left is None or right is None:
            continue
        left_pdb = left.report_pdb or left.source_pdb
        right_pdb = right.report_pdb or right.source_pdb
        pair_dir = output_dir / f"{_safe_label(case_a.label)}_cluster_{left.cluster_id}__{_safe_label(case_b.label)}_cluster_{right.cluster_id}"
        pair_left = pair_dir / f"{_safe_label(case_a.label)}_cluster_{left.cluster_id}.pdb"
        pair_right = pair_dir / f"{_safe_label(case_b.label)}_cluster_{right.cluster_id}_aligned.pdb"
        pair_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(left_pdb, pair_left)

        left_ca = _read_ca_coordinates(left_pdb)
        right_lines, right_all, right_ca_indices = _read_pdb_coordinates(right_pdb)
        _ = right_lines
        right_ca = right_all[right_ca_indices]
        aligned_right = _align_all_coordinates_to_reference(left_ca, right_ca, right_all)
        _write_aligned_pdb(right_pdb, pair_right, aligned_right)
        pair = dict(match)
        pair.update(
            {
                "case_a_pdb": str(pair_left.resolve()),
                "case_b_aligned_pdb": str(pair_right.resolve()),
            }
        )
        pairs.append(pair)
    return pairs


def _write_ensemble_models(case_a: ClusterCase, case_b: ClusterCase, output_dir: Path) -> list[dict[str, Any]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    if not case_a.representatives:
        return []
    reference = case_a.representatives[0]
    reference_pdb = reference.report_pdb or reference.source_pdb
    reference_ca = _read_ca_coordinates(reference_pdb)
    models: list[dict[str, Any]] = []
    for side, case, palette in (
        ("case_a", case_a, ["#154f6f", "#2f6f89", "#6fa0b4", "#aac8d4"]),
        ("case_b", case_b, ["#9c581b", "#d08a3c", "#e2ad6e", "#f0d0a5"]),
    ):
        for rank, representative in enumerate(case.representatives):
            source_pdb = representative.report_pdb or representative.source_pdb
            target_pdb = output_dir / _safe_label(case.label) / f"cluster_{representative.cluster_id}_aligned_to_reference.pdb"
            if side == "case_a" and representative.cluster_id == reference.cluster_id:
                target_pdb.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source_pdb, target_pdb)
            else:
                _, all_coords, ca_indices = _read_pdb_coordinates(source_pdb)
                mobile_ca = all_coords[ca_indices]
                aligned = _align_all_coordinates_to_reference(reference_ca, mobile_ca, all_coords)
                _write_aligned_pdb(source_pdb, target_pdb, aligned)
            models.append(
                {
                    "side": side,
                    "label": case.label,
                    "cluster_id": representative.cluster_id,
                    "rank": int(rank + 1),
                    "population_percent": representative.population_percent,
                    "representative_frame": representative.representative_frame,
                    "descriptor": representative.descriptor,
                    "color": palette[min(rank, len(palette) - 1)],
                    "visible_by_default": bool(rank == 0),
                    "is_reference": bool(side == "case_a" and representative.cluster_id == reference.cluster_id),
                    "pdb": str(target_pdb.resolve()),
                }
            )
    return models


def _write_matches_csv(matches: list[dict[str, Any]], output_file: Path) -> None:
    fieldnames = [
        "case_a_cluster",
        "case_a_population_percent",
        "case_a_representative_frame",
        "case_a_descriptor",
        "case_b_cluster",
        "case_b_population_percent",
        "case_b_representative_frame",
        "case_b_descriptor",
        "representative_ca_rmsd_angstrom",
        "population_delta_b_minus_a_percent",
    ]
    with output_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(matches)


def _build_sampling_enriched_states(
    case_a: ClusterCase,
    case_b: ClusterCase,
    matches: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    matched_b = {int(match["case_b_cluster"]) for match in matches}
    rows = []
    for match in matches:
        delta = float(match["population_delta_b_minus_a_percent"])
        if delta <= 0:
            continue
        rows.append(
            {
                "state_type": "matched_rest2_enriched",
                "standard_cluster": match["case_a_cluster"],
                "rest2_cluster": match["case_b_cluster"],
                "standard_population_percent": match["case_a_population_percent"],
                "rest2_population_percent": match["case_b_population_percent"],
                "population_delta_percent": delta,
                "representative_ca_rmsd_angstrom": match["representative_ca_rmsd_angstrom"],
                "descriptor": match["case_b_descriptor"],
            }
        )
    for representative in case_b.representatives:
        if representative.cluster_id in matched_b:
            continue
        rows.append(
            {
                "state_type": "rest2_only_top_state",
                "standard_cluster": "",
                "rest2_cluster": representative.cluster_id,
                "standard_population_percent": 0.0,
                "rest2_population_percent": representative.population_percent,
                "population_delta_percent": representative.population_percent,
                "representative_ca_rmsd_angstrom": "",
                "descriptor": representative.descriptor,
            }
        )
    return sorted(rows, key=lambda item: float(item["population_delta_percent"]), reverse=True)


def _write_enriched_states_csv(rows: list[dict[str, Any]], output_file: Path) -> None:
    fieldnames = [
        "state_type",
        "standard_cluster",
        "rest2_cluster",
        "standard_population_percent",
        "rest2_population_percent",
        "population_delta_percent",
        "representative_ca_rmsd_angstrom",
        "descriptor",
    ]
    with output_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_heatmap(case_a: ClusterCase, case_b: ClusterCase, matrix: np.ndarray, output_file: Path) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.8, 5.8))
    image = ax.imshow(matrix, cmap="Blues_r")
    ax.set_xticks(range(len(case_b.representatives)))
    ax.set_xticklabels([f"{case_b.label} C{item.cluster_id}" for item in case_b.representatives], rotation=35, ha="right")
    ax.set_yticks(range(len(case_a.representatives)))
    ax.set_yticklabels([f"{case_a.label} C{item.cluster_id}" for item in case_a.representatives])
    ax.set_title("Cluster Representative Cα RMSD", fontsize=12, fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="#1e2d29", fontsize=9)
    cbar = fig.colorbar(image, ax=ax, shrink=0.86)
    cbar.set_label("Cα RMSD (Å)")
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def _write_population_shift_plot(
    case_a: ClusterCase,
    case_b: ClusterCase,
    matches: list[dict[str, Any]],
    output_file: Path,
) -> None:
    import matplotlib.pyplot as plt

    if not matches:
        return
    labels = [f"S C{match['case_a_cluster']} / R C{match['case_b_cluster']}" for match in matches]
    standard_values = [float(match["case_a_population_percent"]) for match in matches]
    rest2_values = [float(match["case_b_population_percent"]) for match in matches]
    deltas = [rest2 - standard for standard, rest2 in zip(standard_values, rest2_values)]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1.25, 0.95]})
    x = np.arange(len(labels))
    width = 0.38
    axes[0].bar(x - width / 2, standard_values, width=width, color="#2f6f89", label=case_a.label)
    axes[0].bar(x + width / 2, rest2_values, width=width, color="#d08a3c", label=case_b.label)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=28, ha="right")
    axes[0].set_ylabel("Cluster population (%)")
    axes[0].set_title("Matched Cluster Population", fontsize=11, fontweight="bold")
    axes[0].grid(alpha=0.18, axis="y")
    axes[0].legend(frameon=False)

    colors = ["#b55f3a" if value > 0 else "#5d7d74" for value in deltas]
    axes[1].barh(labels, deltas, color=colors)
    axes[1].axvline(0, color="#33433e", linewidth=1)
    axes[1].set_xlabel("REST2 - Standard MD population (%)")
    axes[1].set_title("Population Shift", fontsize=11, fontweight="bold")
    axes[1].grid(alpha=0.18, axis="x")
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def _case_payload(case: ClusterCase) -> dict[str, Any]:
    return {
        "label": case.label,
        "cluster_root": str(case.root.resolve()),
        "summary_csv": str(case.summary_csv.resolve()),
        "summary_json": str(case.summary_json.resolve()) if case.summary_json else "",
        "representatives": [
            {
                "cluster_id": item.cluster_id,
                "population_percent": item.population_percent,
                "frames": item.frames,
                "representative_frame": item.representative_frame,
                "descriptor": item.descriptor,
                "pdb": str((item.report_pdb or item.source_pdb).resolve()),
            }
            for item in case.representatives
        ],
    }


def _sampling_coverage(case_a: ClusterCase, case_b: ClusterCase) -> dict[str, Any]:
    population_a = np.array([item.population_percent for item in case_a.representatives], dtype=float) / 100.0
    population_b = np.array([item.population_percent for item in case_b.representatives], dtype=float) / 100.0
    return {
        "case_a_label": case_a.label,
        "case_b_label": case_b.label,
        "case_a_top_cluster_count": len(case_a.representatives),
        "case_b_top_cluster_count": len(case_b.representatives),
        "case_a_dominant_cluster": case_a.representatives[0].cluster_id if case_a.representatives else None,
        "case_b_dominant_cluster": case_b.representatives[0].cluster_id if case_b.representatives else None,
        "case_a_dominant_population_percent": case_a.representatives[0].population_percent if case_a.representatives else None,
        "case_b_dominant_population_percent": case_b.representatives[0].population_percent if case_b.representatives else None,
        "case_a_population_entropy": _entropy(population_a),
        "case_b_population_entropy": _entropy(population_b),
        "case_a_effective_state_count": _effective_state_count(population_a),
        "case_b_effective_state_count": _effective_state_count(population_b),
    }


def _entropy(population: np.ndarray) -> float:
    population = population[np.isfinite(population) & (population > 0)]
    if population.size == 0:
        return 0.0
    population = population / population.sum()
    return float(-np.sum(population * np.log(population)))


def _effective_state_count(population: np.ndarray) -> float:
    return float(np.exp(_entropy(population)))


def _looks_like_sampling_pair(label_a: str, label_b: str) -> bool:
    text = f"{label_a} {label_b}".lower()
    return "standard" in text and "rest2" in text


def _safe_label(label: str) -> str:
    safe = "".join(char.lower() if char.isalnum() else "_" for char in label.strip())
    return safe.strip("_") or "case"
