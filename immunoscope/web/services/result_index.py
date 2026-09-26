"""Build compact result indexes for web and assistant consumption."""

from __future__ import annotations

import json
import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _existing_path(value: Any) -> str | None:
    if not value:
        return None
    path = Path(str(value))
    return str(path) if path.exists() else None


def _read_csv_rows(path: Path, limit: int = 20) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = []
        for row in reader:
            rows.append({key: _coerce_scalar(value) for key, value in row.items()})
            if len(rows) >= limit:
                break
        return rows


def _coerce_scalar(value: Any) -> Any:
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return ""
    try:
        if "." in text or "e" in text.lower():
            return float(text)
        return int(text)
    except ValueError:
        return text


def _module_digest(name: str, root: Path | None, result: dict[str, Any]) -> dict[str, Any]:
    if not root or not root.exists():
        return {"summary": {}, "top_items": []}
    if name == "contact":
        return _contact_digest(root, result)
    if name == "rrcs":
        return _rrcs_digest(root)
    if name == "cluster":
        return _cluster_digest(root)
    if name == "landscape":
        return _landscape_digest(root)
    if name == "rmsf":
        return _json_summary_digest(root / "analysis" / "rmsf" / "rmsf_summary.json")
    if name == "bsa":
        return _json_summary_digest(root / "analysis" / "interface" / "interface_summary.json")
    return {"summary": {}, "top_items": []}


def _contact_digest(root: Path, result: dict[str, Any]) -> dict[str, Any]:
    report = Path(result["contact_report"]) if result.get("contact_report") else root / "analysis" / "contacts" / "contact_report.csv"
    rows = _read_csv_rows(report, limit=200)
    frequency_keys = ["frequency", "contact_frequency", "occupancy", "mean_occupancy"]
    def score(row: dict[str, Any]) -> float:
        for key in frequency_keys:
            if isinstance(row.get(key), (int, float)):
                return float(row[key])
        return 0.0
    top = sorted(rows, key=score, reverse=True)[:10]
    return {
        "summary": {
            "n_contact_pairs": result.get("n_contact_pairs") or len(rows),
            "cutoff_angstrom": result.get("cutoff_angstrom"),
            "min_frequency": result.get("min_frequency"),
        },
        "top_items": top,
    }


def _rrcs_digest(root: Path) -> dict[str, Any]:
    analysis = root / "analysis" / "interactions" / "rrcs"
    summary = _read_json(analysis / "rrcs_summary.json")
    pair_path = analysis / "annotated_rrcs_pair_summary.csv"
    if not pair_path.exists():
        pair_path = analysis / "rrcs_pair_summary.csv"
    pair_rows = _read_csv_rows(pair_path, limit=200)
    top = sorted(pair_rows, key=lambda row: float(row.get("mean_rrcs") or 0.0), reverse=True)[:10]
    return {
        "summary": {
            "pair_scope": summary.get("pair_scope"),
            "identified_pairs": summary.get("n_nonzero_pairs") or summary.get("n_pairs"),
            "n_frames": summary.get("n_frames"),
            "top_pairs": summary.get("top_pairs", [])[:5] if isinstance(summary.get("top_pairs"), list) else [],
        },
        "top_items": top,
    }


def _cluster_digest(root: Path) -> dict[str, Any]:
    analysis = root / "analysis" / "conformation" / "interface_clustering"
    summary_rows = _read_csv_rows(analysis / "summary_table.csv", limit=100)
    feature_rows = _read_csv_rows(analysis / "cluster_feature_digest.csv", limit=100)
    top = sorted(
        summary_rows,
        key=lambda row: float(row.get("population_percent") or row.get("population_fraction") or row.get("fraction") or 0.0),
        reverse=True,
    )[:8]
    return {
        "summary": {
            "n_clusters": len(summary_rows),
            "dominant_cluster": top[0].get("cluster_id") if top else None,
            "dominant_fraction": top[0].get("population_fraction") or top[0].get("fraction") if top else None,
            "dominant_percent": top[0].get("population_percent") if top else None,
            "feature_digest_rows": len(feature_rows),
        },
        "top_items": top,
    }


def _landscape_digest(root: Path) -> dict[str, Any]:
    summary_path = root / "landscape_summary.json"
    if not summary_path.exists():
        summary_path = root / "analysis" / "landscape" / "landscape_summary.json"
    summary = _read_json(summary_path)
    interactive = root / "analysis" / "landscape" / "landscape_interactive.html"
    static = root / "analysis" / "landscape" / "landscape_2d.png"
    return {
        "summary": {
            "reducer": summary.get("reducer"),
            "n_frames": summary.get("n_frames"),
            "coordinate_labels": summary.get("coordinate_labels", [])[:4] if isinstance(summary.get("coordinate_labels"), list) else [],
            "explained_variance": summary.get("explained_variance", [])[:4] if isinstance(summary.get("explained_variance"), list) else [],
            "landscape_summary": str(summary_path) if summary_path.exists() else None,
            "landscape_interactive": str(interactive) if interactive.exists() else None,
            "landscape_2d": str(static) if static.exists() else None,
        },
        "top_items": (summary.get("top_contributors_pc1", []) or [])[:5] if isinstance(summary.get("top_contributors_pc1"), list) else [],
    }


def _json_summary_digest(path: Path) -> dict[str, Any]:
    summary = _read_json(path)
    return {"summary": summary, "top_items": []}


def build_job_result_index(job_dir: Path) -> dict[str, Any]:
    """Create a stable, compact index over an ImmunoScope run directory."""
    run_dir = job_dir / "run"
    summary_path = run_dir / "run_summary.json"
    manifest_path = run_dir / "run_manifest.json"
    summary = _read_json(summary_path)
    manifest = _read_json(manifest_path)
    module_results = summary.get("module_results", {}) if isinstance(summary.get("module_results"), dict) else {}

    modules = []
    for name in summary.get("modules", manifest.get("modules", [])) or []:
        result = module_results.get(name, {}) if isinstance(module_results.get(name), dict) else {}
        root = Path(result["root"]) if result.get("root") else None
        digest = _module_digest(name, root, result)
        artifacts = {}
        for key in ["root", "html", "contact_report"]:
            path = _existing_path(result.get(key))
            if path:
                artifacts[key] = path
        modules.append(
            {
                "name": name,
                "status": result.get("status", "unavailable"),
                "artifacts": artifacts,
                "metrics": {
                    key: value
                    for key, value in result.items()
                    if key not in {"status", "exit_code", "argv", "root", "html", "contact_report", "error"}
                },
                "summary": digest.get("summary", {}),
                "top_items": digest.get("top_items", []),
                "error": result.get("error"),
            }
        )

    report_html = _existing_path(summary.get("report_html"))
    index = {
        "schema_version": "immunoscope.web.result_index.v1",
        "job_id": summary.get("job_id") or job_dir.name,
        "status": summary.get("status", manifest.get("status", "unknown")),
        "run_dir": str(run_dir),
        "manifest": str(manifest_path) if manifest_path.exists() else None,
        "summary": str(summary_path) if summary_path.exists() else None,
        "report_html": report_html,
        "prepared_input": summary.get("prepared_input", manifest.get("prepared_input", {})),
        "modules": modules,
        "updated_at": datetime.now(timezone.utc).isoformat(),
    }
    output = job_dir / "job_result_index.json"
    output.write_text(json.dumps(index, ensure_ascii=False, indent=2), encoding="utf-8")
    return index
