"""Structured result/evidence schema for ImmunoScope assistant skills."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass
class EvidenceCard:
    """Compact evidence item that can be cited by an assistant response."""

    module: str
    title: str
    summary: str
    metrics: dict[str, Any] = field(default_factory=dict)
    top_items: list[dict[str, Any]] = field(default_factory=list)
    artifacts: dict[str, str] = field(default_factory=dict)
    confidence: str = "medium"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def build_agent_result_schema(result_index: dict[str, Any]) -> dict[str, Any]:
    """Convert a web result index into an agent-facing schema."""
    modules = result_index.get("modules", []) if isinstance(result_index.get("modules"), list) else []
    evidence_cards = [_module_to_evidence_card(module).to_dict() for module in modules if isinstance(module, dict)]
    completed_modules = [module.get("name") for module in modules if module.get("status") == "completed"]
    failed_modules = [module.get("name") for module in modules if module.get("status") == "failed"]
    return {
        "schema_version": "immunoscope.agent.result_schema.v1",
        "job_id": result_index.get("job_id"),
        "status": result_index.get("status"),
        "available_modules": completed_modules,
        "failed_modules": failed_modules,
        "report_html": result_index.get("report_html"),
        "prepared_input": result_index.get("prepared_input", {}),
        "evidence_cards": evidence_cards,
        "suggested_skills": suggest_skills(completed_modules),
    }


def answer_from_agent_schema(question: str, schema: dict[str, Any], *, skill: str | None = None) -> dict[str, Any] | None:
    """Answer simple single-job skill requests directly from agent result schema."""
    intent = _normalize_intent(skill or question)
    cards = schema.get("evidence_cards", []) if isinstance(schema.get("evidence_cards"), list) else []
    by_module = {card.get("module"): card for card in cards if isinstance(card, dict)}
    if intent == "explain_fel":
        return _answer_fel(schema, by_module.get("landscape"))
    if intent == "summarize_existing_job":
        return _answer_summary(schema, cards)
    if intent == "identify_hotspots":
        return _answer_hotspots(schema, by_module)
    return None


def suggest_skills(completed_modules: list[str]) -> list[dict[str, str]]:
    """Return assistant skills that are meaningful for the available outputs."""
    available = set(completed_modules)
    skills: list[dict[str, str]] = [
        {
            "name": "summarize_existing_job",
            "kind": "query",
            "description": "Summarize completed modules and available evidence.",
        }
    ]
    if {"contact", "rrcs"} & available:
        skills.append(
            {
                "name": "identify_hotspots",
                "kind": "query",
                "description": "Rank interface residues or pairs from contact and RRCS evidence.",
            }
        )
    if {"contact", "bsa", "rmsf", "landscape"} & available:
        skills.append(
            {
                "name": "diagnose_stability",
                "kind": "diagnostic",
                "description": "Combine interface, flexibility, BSA, and landscape outputs into a stability diagnosis.",
            }
        )
    if "landscape" in available:
        skills.append(
            {
                "name": "explain_fel",
                "kind": "query",
                "description": "Explain the free-energy landscape result and link FEL artifacts.",
            }
        )
    if {"contact", "rrcs", "rmsf"} & available:
        skills.append(
            {
                "name": "suggest_mutation_sites",
                "kind": "design",
                "description": "Suggest conservative candidate residues for manual mutation review.",
            }
        )
    return skills


def _answer_fel(schema: dict[str, Any], card: dict[str, Any] | None) -> dict[str, Any] | None:
    if not card:
        return None
    metrics = card.get("metrics", {}) if isinstance(card.get("metrics"), dict) else {}
    artifacts = {
        key: value
        for key, value in {
            "landscape_summary": metrics.get("landscape_summary"),
            "landscape_interactive": metrics.get("landscape_interactive"),
            "landscape_2d": metrics.get("landscape_2d"),
        }.items()
        if value
    }
    reducer = metrics.get("reducer") or "unknown"
    n_frames = metrics.get("n_frames") or "n/a"
    labels = metrics.get("coordinate_labels") or []
    variance = metrics.get("explained_variance") or []
    evidence = [
        f"FEL reducer: {reducer}.",
        f"FEL sampled frames in the reduced landscape: {n_frames}.",
    ]
    if labels:
        evidence.append(f"Coordinate labels: {', '.join(map(str, labels[:4]))}.")
    if variance:
        evidence.append(f"Explained variance snapshot: {', '.join(_fmt_number(v) for v in variance[:4])}.")
    if artifacts:
        evidence.append("FEL plot and summary artifacts are available from the job-level result index.")
    return {
        "query_type": "landscape_status",
        "answer": f"The FEL module completed with {reducer} reduction over {n_frames} sampled frames.",
        "evidence": evidence,
        "sources_used": list(artifacts.values()),
        "sources": list(artifacts.values()),
        "top_items": card.get("top_items", [])[:5],
        "confidence": card.get("confidence", "medium"),
        "skill": "explain_fel",
        "response_panel": {
            "summary": f"FEL completed with {reducer}; use the linked 2D and interactive artifacts for basin inspection.",
            "sections": [{"title": "FEL readout", "items": evidence}],
            "caution": "Basin interpretation should be reviewed with the plotted landscape and trajectory context.",
        },
        "evidence_panel": {
            "title": "FEL evidence",
            "highlights": [
                {"label": "Reducer", "value": reducer},
                {"label": "Frames", "value": n_frames},
            ],
            "bullets": evidence,
            "top_items_title": "Top contributors",
            "top_items": card.get("top_items", [])[:5],
            "sources": list(artifacts.values()),
            "notes": [],
        },
    }


def _answer_summary(schema: dict[str, Any], cards: list[dict[str, Any]]) -> dict[str, Any]:
    available = schema.get("available_modules", [])
    failed = schema.get("failed_modules", [])
    evidence = [f"{card.get('module')}: {card.get('summary')}" for card in cards[:6]]
    answer = f"This job has {len(available)} completed modules: {', '.join(map(str, available))}."
    if failed:
        answer += f" Failed modules: {', '.join(map(str, failed))}."
    return {
        "query_type": "job_summary",
        "answer": answer,
        "evidence": evidence,
        "sources_used": [],
        "sources": [],
        "top_items": cards[:6],
        "confidence": "high" if available else "low",
        "skill": "summarize_existing_job",
        "response_panel": {
            "summary": answer,
            "sections": [{"title": "Available evidence", "items": evidence}],
            "caution": "This summary is based on generated result indexes, not raw trajectory inspection.",
        },
        "evidence_panel": {
            "title": "Job evidence",
            "highlights": [
                {"label": "Completed modules", "value": len(available)},
                {"label": "Failed modules", "value": len(failed)},
            ],
            "bullets": evidence,
            "top_items_title": "Evidence cards",
            "top_items": cards[:6],
            "sources": [],
            "notes": [],
        },
    }


def _answer_hotspots(schema: dict[str, Any], by_module: dict[str, dict[str, Any]]) -> dict[str, Any] | None:
    rrcs = by_module.get("rrcs")
    contact = by_module.get("contact")
    if not rrcs and not contact:
        return None
    top_items = []
    if rrcs:
        top_items.extend(rrcs.get("top_items", [])[:5])
    if contact:
        top_items.extend(contact.get("top_items", [])[:5])
    evidence = []
    if rrcs:
        evidence.append(rrcs.get("summary", "RRCS evidence is available."))
    if contact:
        evidence.append(contact.get("summary", "Contact evidence is available."))
    top = top_items[0] if top_items else {}
    label = (
        top.get("pair")
        or _pair_label(top)
        or top.get("residue")
        or top.get("residue_i")
        or top.get("residue_j")
        or "available ranked evidence"
    )
    return {
        "query_type": "hotspot_summary",
        "answer": f"Hotspot screening should start from {label}.",
        "evidence": evidence,
        "sources_used": [],
        "sources": [],
        "top_items": top_items[:8],
        "confidence": "high" if rrcs and top_items else "medium",
        "skill": "identify_hotspots",
        "response_panel": {
            "summary": f"Top hotspot evidence starts from {label}.",
            "sections": [{"title": "Hotspot evidence", "items": evidence}],
            "caution": "Treat hotspot ranking as a first-pass screen; inspect structures before mutation decisions.",
        },
        "evidence_panel": {
            "title": "Hotspot evidence",
            "highlights": [
                {"label": "RRCS items", "value": len(rrcs.get("top_items", [])) if rrcs else 0},
                {"label": "Contact items", "value": len(contact.get("top_items", [])) if contact else 0},
            ],
            "bullets": evidence,
            "top_items_title": "Ranked hotspot evidence",
            "top_items": top_items[:8],
            "sources": [],
            "notes": [],
        },
    }


def _module_to_evidence_card(module: dict[str, Any]) -> EvidenceCard:
    name = str(module.get("name") or "module")
    status = module.get("status") or "unavailable"
    summary = module.get("summary", {}) if isinstance(module.get("summary"), dict) else {}
    metrics = module.get("metrics", {}) if isinstance(module.get("metrics"), dict) else {}
    top_items = module.get("top_items", []) if isinstance(module.get("top_items"), list) else []
    artifacts = module.get("artifacts", {}) if isinstance(module.get("artifacts"), dict) else {}

    if name == "contact":
        title = "Interface contact evidence"
        text = f"{summary.get('n_contact_pairs', 'n/a')} contact pairs detected"
        confidence = "high" if status == "completed" and summary.get("n_contact_pairs") else "medium"
    elif name == "rrcs":
        title = "RRCS hotspot evidence"
        text = f"{summary.get('identified_pairs', 'n/a')} RRCS pairs available"
        confidence = "high" if status == "completed" and top_items else "medium"
    elif name == "cluster":
        title = "Interface conformation evidence"
        dominant = summary.get("dominant_cluster")
        percent = summary.get("dominant_percent") or summary.get("dominant_fraction")
        text = f"{summary.get('n_clusters', 'n/a')} clusters; dominant cluster={dominant}, population={percent}"
        confidence = "high" if status == "completed" and summary.get("n_clusters") else "medium"
    elif name == "landscape":
        title = "FEL evidence"
        text = f"Reducer={summary.get('reducer', 'n/a')}; frames={summary.get('n_frames', 'n/a')}"
        confidence = "high" if status == "completed" and summary.get("landscape_summary") else "medium"
    elif name == "bsa":
        title = "Buried surface area evidence"
        text = "BSA summary is available" if summary else "BSA summary is not available"
        confidence = "medium"
    elif name == "rmsf":
        title = "Flexibility evidence"
        text = "RMSF summary is available" if summary else "RMSF summary is not available"
        confidence = "medium"
    else:
        title = f"{name} evidence"
        text = f"Status={status}"
        confidence = "medium" if status == "completed" else "low"

    if status != "completed":
        confidence = "low"
        text = module.get("error") or text

    return EvidenceCard(
        module=name,
        title=title,
        summary=text,
        metrics={**summary, **metrics},
        top_items=top_items[:8],
        artifacts=artifacts,
        confidence=confidence,
    )


def _normalize_intent(value: str) -> str:
    text = str(value or "").strip().lower().replace("_", "-")
    if any(token in text for token in ["fel", "free energy", "landscape", "explain-fel"]):
        return "explain_fel"
    if any(token in text for token in ["identify-hotspots", "hotspot", "hotspots"]):
        return "identify_hotspots"
    if any(token in text for token in ["summarize-existing-job", "summarize", "summary"]):
        return "summarize_existing_job"
    return text.replace("-", "_")


def _fmt_number(value: Any) -> str:
    try:
        return f"{float(value):.3f}"
    except (TypeError, ValueError):
        return str(value)


def _pair_label(item: dict[str, Any]) -> str | None:
    left = item.get("residue_label_1") or item.get("tcr_residue_label")
    right = item.get("residue_label_2") or item.get("partner_residue_label")
    if left and right:
        return f"{left} - {right}"
    return None
