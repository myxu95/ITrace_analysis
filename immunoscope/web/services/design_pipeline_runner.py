"""Background runner: the multi-agent mutation-design pipeline as a design "run".

Bridges a web Design task to ``immunoscope.agent.multi_agent``. The pipeline is
IO-bound on LLM calls (like the chat agent), so it runs in the FastAPI event
loop via ``asyncio.create_task`` — no thread offload needed (cf. job_runner,
which offloads CPU-bound MD analysis with ``asyncio.to_thread``).

Decision A (pipeline-led): the multi-agent pipeline is the recommender. This
runner:
  - streams stage progress to ``{task_id}_pipeline_status.json`` (frontend polls)
  - writes the full trace to ``{task_id}_pipeline.json`` (4 briefs + candidate
    table + recommendation + critic verdict/findings + deterministic findings +
    numeric audit)
  - folds the recommended sites into the task's drafts (``source="multi_agent"``)

The chat agent later reads the artifact and only *explains* it (see
``design_context``); it does not generate competing recommendations.

Design helpers from ``immunoscope.web.routers.design`` are imported lazily
inside functions to avoid a circular import (design.py imports this module's
``start_design_pipeline``).
"""
from __future__ import annotations

import asyncio
import json
import logging
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from immunoscope.web.config import settings

logger = logging.getLogger(__name__)

# Same directory design.py persists tasks/drafts into. Defined here (rather than
# imported) so path resolution never depends on the design router module.
DESIGN_DIR = settings.JOBS_DIR.parent / "design"

_RUNNING: dict[str, asyncio.Task] = {}

_DEFAULT_INTENT = (
    "Improve the binding affinity of this TCR for its peptide-HLA target "
    "while preserving specificity."
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _status_path(task_id: str) -> Path:
    return DESIGN_DIR / f"{task_id}_pipeline_status.json"


def _artifact_path(task_id: str) -> Path:
    return DESIGN_DIR / f"{task_id}_pipeline.json"


def read_status(task_id: str) -> dict[str, Any]:
    """Current pipeline-run status for a task (``{"status": "idle"}`` if none)."""
    p = _status_path(task_id)
    if not p.exists():
        return {"status": "idle"}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {"status": "unknown"}


def read_artifact(task_id: str) -> Optional[dict[str, Any]]:
    p = _artifact_path(task_id)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None


def _write_status(task_id: str, state: dict[str, Any]) -> None:
    state["updated_at"] = _now()
    DESIGN_DIR.mkdir(parents=True, exist_ok=True)
    _status_path(task_id).write_text(
        json.dumps(state, ensure_ascii=False, indent=2), encoding="utf-8"
    )


# Rough stage → percentage map, so the frontend bar advances monotonically.
_STAGE_PROGRESS = {
    "start": 5,
    "round1_done": 60,
    "round2": 78,
    "critic": 88,
    "revision": 92,
    "done": 100,
}


def _progress_for(event: dict, state: dict) -> int:
    stage = event.get("stage")
    cur = int(state.get("progress", 0))
    if stage == "round1":
        # four readers; nudge ~12% each as they land, capped below round1_done.
        cur = min(58, cur + 13)
        return cur
    return max(cur, _STAGE_PROGRESS.get(stage, cur))


# ---------------------------------------------------------------------------
# Artifact serialization
# ---------------------------------------------------------------------------

def _ref_to_dict(r: Any) -> dict[str, Any]:
    return {
        "ref_id": getattr(r, "ref_id", ""),
        "kind": getattr(r, "kind", ""),
        "excerpt": getattr(r, "excerpt", "") or getattr(r, "content", ""),
        "msg_id": getattr(r, "msg_id", ""),
        "pmid": getattr(r, "pmid", ""),
        "path": getattr(r, "path", ""),
    }


def _message_to_dict(run_result: Any) -> Optional[dict[str, Any]]:
    if run_result is None:
        return None
    msg = getattr(run_result, "message", None)
    if msg is None:
        return None
    audit = getattr(run_result, "numeric_audit", None)
    return {
        "msg_id": getattr(msg, "msg_id", ""),
        "narrative": getattr(msg, "narrative", "") or "",
        "structured_payload": getattr(msg, "structured_payload", None),
        "evidence_refs": [_ref_to_dict(r) for r in getattr(msg, "evidence_refs", []) or []],
        "turns_used": getattr(run_result, "turns_used", None),
        "numeric_audit": audit.summary() if audit and hasattr(audit, "summary") else None,
    }


def _serialize_result(result: Any, ctx: Any) -> dict[str, Any]:
    """Flatten a PipelineResult into a JSON-safe artifact for the frontend."""
    briefs = {}
    for key, attr in (
        ("bio", "bio_result"),
        ("interaction", "interaction_result"),
        ("conformation", "conformation_result"),
        ("interface_exposure", "interface_exposure_result"),
    ):
        briefs[key] = _message_to_dict(getattr(result, attr, None))

    crit = getattr(result, "critique_result", None)
    critique = None
    if crit is not None:
        cm = getattr(crit, "message", None)
        payload = getattr(cm, "structured_payload", None) or {}
        critique = {
            "verdict": payload.get("verdict"),
            "narrative": getattr(cm, "narrative", "") or "",
            "findings": payload.get("findings", []),
            "turns_used": getattr(crit, "turns_used", None),
        }

    det = []
    for f in getattr(result, "deterministic_findings", []) or []:
        det.append({
            "issue_type": getattr(f, "issue_type", ""),
            "severity": getattr(f, "severity", ""),
            "target": getattr(f, "target", ""),
            "detail": getattr(f, "detail", ""),
        })

    return {
        "succeeded": bool(getattr(result, "succeeded", False)),
        "errors": {k: str(v) for k, v in (getattr(result, "errors", {}) or {}).items()},
        "briefs": briefs,
        "candidate_table": getattr(ctx, "recommendation_preamble", None),
        "recommendation": _message_to_dict(getattr(result, "recommendation_result", None)),
        "critique": critique,
        "deterministic_findings": det,
        "numeric_audit_summary": getattr(ctx, "numeric_audit_summary", None),
        "generated_at": _now(),
    }


# ---------------------------------------------------------------------------
# Recommendation → drafts
# ---------------------------------------------------------------------------

# Three-letter -> one-letter residue codes. The recommendation agent labels
# sites either one-letter ("β-A99") or three-letter ("α-ASP93"); both encode
# the wild-type residue, so we normalize to a one-letter `current_aa`.
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M",
}

# chain prefix (Greek) + residue token (1 or 3 letters) + position. Handles
# "α-ASP93", "β-A99", "α-V92", "GLY95", "CYS23 (D)".
_RES_RE = re.compile(r"^\s*([αβγδ])?\s*[-\s]?\s*([A-Za-z]{1,3})(\d+)")


def _split_residue(label: str) -> tuple[str, str]:
    """Best-effort (chain, current_aa) from a label like 'α-ASP93' / 'β-A99' /
    'CYS23 (D)'.

    Recognizes both one-letter and three-letter wild-type codes and normalizes
    `current_aa` to one letter. Returns ("", "") only when nothing parses — the
    residue string itself is always kept verbatim, so a failed split never loses
    information."""
    if not label:
        return "", ""
    m = _RES_RE.match(label)
    aa = ""
    chain = ""
    if m:
        chain_tok, res_tok, _ = m.groups()
        if chain_tok in {"α", "β", "γ", "δ"}:
            chain = chain_tok
        res_tok = res_tok.upper()
        if len(res_tok) == 3:
            aa = _THREE_TO_ONE.get(res_tok, "")
        elif len(res_tok) == 1:
            aa = res_tok
    # Trailing "(D)" chain hint, if present and no Greek chain was found.
    if not chain:
        tail = re.search(r"\(([A-Za-z0-9])\)", label)
        if tail:
            chain = tail.group(1)
    return chain, aa


# Human-facing chain labels: the recommender prefixes residues with a Greek
# TCR-domain letter (α/β); show "TCRβ" rather than a bare "β" on the cards.
_FRIENDLY_CHAIN = {"α": "TCRα", "β": "TCRβ", "γ": "TCRγ", "δ": "TCRδ"}


def _friendly_chain(residue: str, chain: str) -> str:
    """'β' -> 'TCRβ'. Falls back to the raw chain token for non-TCR chains."""
    for tok in (chain, residue[:1] if residue else ""):
        if tok in _FRIENDLY_CHAIN:
            return _FRIENDLY_CHAIN[tok]
    return chain or ""


def _clean_region(region: str) -> str:
    """'CDR3β' -> 'CDR3' (the chain is already shown separately)."""
    if not region:
        return ""
    return re.sub(r"[αβγδ]$", "", region.strip())


def _res_key(label: str) -> str:
    """Normalized join key for a residue across label formats:
    'β-A99' / 'β-ARG99' / 'CDR3β-A99' all collapse to 'β99'. The recommender
    and the reader briefs spell residues differently (1- vs 3-letter, with/out
    a region prefix), so a chain+number key is the only reliable join."""
    if not label:
        return ""
    s = str(label)
    gm = re.search(r"[αβγδ]", s)
    # The RESIDUE number is the last digit run — a leading "CDR3" would
    # otherwise hijack the key ("CDR3β-A99" must be β99, not β3).
    nums = re.findall(r"\d+", s)
    if gm and nums:
        return gm.group(0) + nums[-1]
    if nums:
        return nums[-1]
    return s.upper()


def _build_evidence_index(interaction_payload: Optional[dict],
                          bio_payload: Optional[dict]) -> dict:
    """res_key -> {'region', 'hotspot', 'prior'} joined from the reader briefs.

    The recommendation payload carries only labels + a strategy class; the
    quantitative evidence (RRCS / occupancy / dominant partner) lives in the
    InteractionReader hotspots and the literature precedents in the BioAgent's
    `prior_engineering`. We index both by `_res_key` so each card can show real
    'Contact environment' content instead of opaque internal refs."""
    index: dict = {}
    for h in ((interaction_payload or {}).get("hotspots") or []):
        k = _res_key(h.get("residue") or "")
        if not k:
            continue
        index.setdefault(k, {})["hotspot"] = h
        reg = str(h.get("region") or "").strip()
        if reg:
            index[k]["region"] = reg
    for p in ((bio_payload or {}).get("prior_engineering") or []):
        k = _res_key(p.get("residue") or "")
        if not k:
            continue
        index.setdefault(k, {}).setdefault("prior", []).append(p)
    return index


def _hotspot_md_line(hs: dict) -> str:
    """Human-readable contact-environment summary from a hotspot dict."""
    bits = []
    rr = hs.get("rrcs_mean")
    if isinstance(rr, (int, float)):
        bits.append(f"RRCS {rr:.1f}")
    occ = hs.get("occupancy_max")
    if isinstance(occ, (int, float)):
        bits.append(f"{round(occ * 100)}% occupancy")
    partner = hs.get("dominant_partner_residue")
    if partner:
        pc = hs.get("dominant_partner_chain")
        bits.append(f"partner {partner}" + (f" ({pc})" if pc else ""))
    role = hs.get("role")
    if role:
        bits.append(str(role))
    return " · ".join(bits)


def _build_supporting_evidence(rec: dict, ev: dict) -> list[dict]:
    """md_data (MD contact environment) + literature (prior engineering PMIDs)
    in the shape the card renders. Falls back to the brief refs only when no
    grounded evidence joined."""
    out: list[dict] = []
    hs = ev.get("hotspot")
    if hs:
        line = _hotspot_md_line(hs)
        if line:
            out.append({"type": "md_data", "data": line})
    seen_pmids: set = set()
    for p in ev.get("prior", []) or []:
        pmid = str(p.get("pmid") or "").strip()
        sub = str(p.get("substitution") or "").strip()
        outcome = str(p.get("outcome") or "").strip()
        data = (f"{sub}: {outcome}" if sub or outcome else "prior engineering")
        if pmid and pmid not in seen_pmids:
            seen_pmids.add(pmid)
            out.append({"type": "literature", "pmid": pmid, "data": data})
        elif not pmid:
            out.append({"type": "md_data", "data": data})
    if not out:
        refs = ", ".join(rec.get("supporting_brief_refs", []) or [])
        if refs:
            out.append({"type": "multi_agent", "data": refs})
    return out


def _recommendation_to_drafts(rec_payload: dict,
                              evidence_index: Optional[dict] = None) -> list[dict]:
    evidence_index = evidence_index or {}
    drafts = []
    for rec in (rec_payload or {}).get("recommendations", []) or []:
        residue = rec.get("residue", "") or ""
        chain, current_aa = _split_residue(residue)
        ev = evidence_index.get(_res_key(residue), {})
        handoff = rec.get("handoff_to_physchem", {}) or {}
        direction = handoff.get("design_hint_class", "")
        cand_aas = [str(a).strip().upper() for a in (handoff.get("candidate_residues") or [])
                    if str(a).strip()]
        # You cannot "mutate" a residue to itself: drop the wild-type from the
        # proposed substitutions (it is frequently listed as a "keep" option),
        # de-duplicating while preserving order. A site left with nothing is a
        # preserve / do-not-mutate call, flagged so the UI never shows "R→R".
        wt = (current_aa or "").upper()
        seen: set = set()
        muts: list = []
        for c in cand_aas:
            if c == wt or c in seen:
                continue
            seen.add(c)
            muts.append(c)
        # Preserve == nothing real left to substitute to. A "preserve"
        # direction that still offers conservative swaps (e.g. R→K/H) is NOT a
        # pure preserve — show those rather than hiding them.
        is_preserve = not muts
        # Region: authoritative from the interaction hotspot (joined via
        # _res_key); else a "(CDR3β)" already embedded in the label; else blank
        # — never the strategy role.
        region_full = ev.get("region") or ""
        if not region_full:
            m = re.search(r"\(([^)]+)\)", residue)
            if m and re.match(r"(CDR|FR|FW|HV)", m.group(1), re.I):
                region_full = m.group(1)
        drafts.append({
            "residue": residue,
            "chain": _friendly_chain(residue, chain),
            "region": _clean_region(region_full),
            "current_aa": current_aa,
            "suggested_mutations": muts,
            "candidate_residues": cand_aas,
            "direction": direction,
            "is_preserve": is_preserve,
            "priority": rec.get("priority", "medium") or "medium",
            "confidence": "medium",
            "rationale": rec.get("rationale_brief", "") or "",
            "risks": rec.get("conflicts", []) or [],
            "constraints": handoff.get("constraints", []) or [],
            "supporting_evidence": _build_supporting_evidence(rec, ev),
            "source": "multi_agent",
            "rec_id": rec.get("rec_id", ""),
        })
    return drafts


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

async def run_design_pipeline(task_id: str) -> None:
    """Execute the multi-agent pipeline for a design task and persist results."""
    from immunoscope.web.routers.design import _load_task, append_draft

    task = _load_task(task_id)
    if not task:
        _write_status(task_id, {
            "task_id": task_id, "status": "failed", "started_at": _now(),
            "error": f"design task not found: {task_id}", "progress": 100,
        })
        return

    run_dir = task.get("source_run_dir")
    if not run_dir or not Path(run_dir).exists():
        _write_status(task_id, {
            "task_id": task_id, "status": "failed", "started_at": _now(),
            "error": f"source run dir unavailable: {run_dir}", "progress": 100,
        })
        return

    design_intent = " ".join(task.get("design_goals") or []).strip() or _DEFAULT_INTENT
    case_id = task.get("source_job_name") or task.get("source_job_id") or task_id

    state: dict[str, Any] = {
        "task_id": task_id,
        "status": "running",
        "started_at": _now(),
        "finished_at": None,
        "stage": "start",
        "progress": 5,
        "events": [],
        "error": None,
        "succeeded": None,
        "n_drafts_added": 0,
    }
    _write_status(task_id, state)

    # progress_cb runs synchronously inside the pipeline; keep it cheap.
    def cb(event: dict) -> None:
        state["stage"] = event.get("stage", state["stage"])
        state["progress"] = _progress_for(event, state)
        state["events"].append(event)
        if event.get("stage") == "done":
            state["succeeded"] = bool(event.get("succeeded"))
        _write_status(task_id, state)

    try:
        from immunoscope.agent.config import get_settings
        from immunoscope.agent.llm import make_llm_client
        from immunoscope.agent.multi_agent import (
            MutationDesignContext, build_default_pipeline,
        )

        llm = make_llm_client(get_settings())
        pipeline = build_default_pipeline(llm)
        ctx = MutationDesignContext(
            case_dir=str(run_dir), case_id=str(case_id), design_intent=design_intent,
        )
        result = await pipeline.run(ctx, progress_cb=cb)

        # Persist the full trace artifact.
        artifact = _serialize_result(result, ctx)
        _artifact_path(task_id).write_text(
            json.dumps(artifact, ensure_ascii=False, indent=2), encoding="utf-8"
        )

        # Fold recommended sites into the task drafts (decision A). The CDR
        # loop + quantitative contact evidence per residue live in the reader
        # briefs (InteractionReader hotspots, BioAgent prior_engineering), not
        # the recommendation, so build an evidence index to annotate the cards.
        def _payload_of(attr):
            r = getattr(result, attr, None)
            if r is not None and getattr(r, "message", None) is not None:
                return getattr(r.message, "structured_payload", None)
            return None
        evidence_index = _build_evidence_index(
            _payload_of("interaction_result"), _payload_of("bio_result"),
        )
        rec = getattr(result, "recommendation_result", None)
        n_added = 0
        if rec is not None and getattr(rec, "message", None) is not None:
            payload = getattr(rec.message, "structured_payload", None) or {}
            for draft in _recommendation_to_drafts(payload, evidence_index):
                try:
                    append_draft(task_id, draft)
                    n_added += 1
                except Exception as exc:  # noqa: BLE001
                    logger.warning("append_draft failed for %s: %s", task_id, exc)

        state.update({
            "status": "completed" if result.succeeded else "completed_with_errors",
            "finished_at": _now(),
            "progress": 100,
            "stage": "done",
            "succeeded": bool(result.succeeded),
            "n_drafts_added": n_added,
            "errors": {k: str(v) for k, v in (result.errors or {}).items()},
        })
        _write_status(task_id, state)
        logger.info("design pipeline %s finished: succeeded=%s drafts+=%d",
                    task_id, result.succeeded, n_added)

    except Exception as exc:  # noqa: BLE001
        logger.exception("design pipeline %s failed", task_id)
        state.update({
            "status": "failed", "finished_at": _now(), "progress": 100,
            "error": str(exc),
        })
        _write_status(task_id, state)
    finally:
        _RUNNING.pop(task_id, None)


def start_design_pipeline(task_id: str) -> dict[str, str]:
    """Start (or report already-running) the pipeline for a design task."""
    task = _RUNNING.get(task_id)
    if task and not task.done():
        return {"status": "already_running"}
    new_task = asyncio.create_task(run_design_pipeline(task_id))
    _RUNNING[task_id] = new_task
    return {"status": "started"}


def is_running(task_id: str) -> bool:
    t = _RUNNING.get(task_id)
    return bool(t and not t.done())
