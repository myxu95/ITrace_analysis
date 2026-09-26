"""Design Copilot router - conversational mutation design via Agent.

Replaces the old batch recommendation engine. Each design task is now a
persistent conversation session: the user describes intent, the Agent
investigates with tools (presenters, RAG), commits recommendations to a
draft, and the user iterates.
"""

from __future__ import annotations

import json
import logging
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field, field_validator

from immunoscope.recommendation.task_spec import TaskSpec
from immunoscope.web.config import settings
from immunoscope.web.services.job_store import job_store

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/design", tags=["design"])

# Design tasks are stored as JSON files (lightweight, no need for new DB table)
DESIGN_DIR = settings.JOBS_DIR.parent / "design"
DESIGN_DIR.mkdir(parents=True, exist_ok=True)

# Task key → TaskSpec factory. Keys are the *only* values accepted on the
# wire; the wider (subject, objective, target_subinterface) triple is an
# internal detail held by TaskSpec. Keep this map in sync with the wired
# combinations in `immunoscope.recommendation.task_spec._WIRED_COMBINATIONS`.
# D-B7 (2026-05-27): adding a new task here also requires:
#   - a candidate-extraction branch in _extract_design_context
#   - a system-prompt branch in agent.design_context.system_prompt_for_task
#   - a baseline routed for the task in evaluation.benchmark.task_routing
SUPPORTED_TASK_KEYS: dict[str, "type[TaskSpec] | None"] = {
    "tcr_affinity": None,        # → TaskSpec.tcr_affinity()
    "peptide_presentation": None,  # → TaskSpec.peptide_presentation()
}


def _task_spec_from_key(key: str) -> TaskSpec:
    """Resolve a wire-level task key into a TaskSpec.

    Raises HTTPException(400) on unknown keys so the API surface gives a
    structured 4xx rather than a Python ValueError.
    """
    if key == "tcr_affinity":
        return TaskSpec.tcr_affinity()
    if key == "peptide_presentation":
        return TaskSpec.peptide_presentation()
    raise HTTPException(
        status_code=400,
        detail=(
            f"Unknown task key {key!r}. Supported keys: "
            f"{sorted(SUPPORTED_TASK_KEYS.keys())}"
        ),
    )


def _task_key_for_spec(task: TaskSpec) -> str:
    """Reverse of _task_spec_from_key — wire-level key for a TaskSpec.

    Keeps the wire-level key as the canonical externalized form even when
    the rest of the pipeline holds a TaskSpec. Falls back to ``task.slug``
    on an unrecognized spec rather than raising, because reaching here with
    an unwired spec would mean TaskSpec.__post_init__ allowed it — that's
    a kernel bug, not an API-surface bug, so we surface the slug instead
    of a 500.
    """
    if task.subject == "tcr" and task.objective == "affinity":
        return "tcr_affinity"
    if task.subject == "peptide" and task.objective == "presentation":
        return "peptide_presentation"
    return task.slug


# ----------------------------------------------------------------------------
# Models
# ----------------------------------------------------------------------------

class DesignTaskRequest(BaseModel):
    """Request to create a new Design Copilot task."""
    source_job_id: str = Field(..., description="ID of the source analysis job")
    task: str = Field(
        default="tcr_affinity",
        description=(
            "Design task key. One of: "
            "'tcr_affinity' (legacy default — TCR mutations against fixed "
            "pMHC), 'peptide_presentation' (peptide mutations to tune HLA "
            "presentation). Mirrors `TaskSpec` wired combinations."
        ),
    )
    design_goals: list[str] = Field(
        default_factory=list,
        description="Optional design goals (free text or known keywords)",
    )
    name: str = Field(default="", description="Optional human-readable name")
    notes: str = Field(default="")
    constraints: dict[str, Any] = Field(
        default_factory=dict,
        description="User-defined constraints: must_include, must_exclude, focus_region, notes",
    )

    @field_validator("task")
    @classmethod
    def _validate_task(cls, v: str) -> str:
        if v not in SUPPORTED_TASK_KEYS:
            raise ValueError(
                f"task={v!r} not supported. Choose one of: "
                f"{sorted(SUPPORTED_TASK_KEYS.keys())}"
            )
        return v


class DesignTaskSummary(BaseModel):
    id: str
    name: str
    status: str  # active | archived
    source_job_id: str
    source_job_name: str
    # Task is "tcr_affinity" by default for both new tasks and legacy task
    # files that predate the field — see _load_task for the back-compat fill.
    task: str = "tcr_affinity"
    design_goals: list[str] = Field(default_factory=list)
    created_at: str
    updated_at: str
    n_drafts: int = 0
    n_messages: int = 0
    notes: str = ""


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

def _task_path(task_id: str) -> Path:
    return DESIGN_DIR / f"{task_id}.json"


def _drafts_path(task_id: str) -> Path:
    return DESIGN_DIR / f"{task_id}_drafts.json"


def _save_task(task: dict[str, Any]) -> None:
    path = _task_path(task["id"])
    path.write_text(json.dumps(task, indent=2, default=str), encoding="utf-8")


def _load_task(task_id: str) -> dict[str, Any] | None:
    path = _task_path(task_id)
    if not path.exists():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    # Back-compat: tasks created before 2026-05-27 have no `task` key.
    # All of those were TCR-affinity designs — that was the only path that
    # existed — so fill the default rather than raising on missing field.
    data.setdefault("task", "tcr_affinity")
    return data


def _load_drafts(task_id: str) -> list[dict[str, Any]]:
    path = _drafts_path(task_id)
    if not path.exists():
        return []
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return []


def _save_drafts(task_id: str, drafts: list[dict[str, Any]]) -> None:
    path = _drafts_path(task_id)
    path.write_text(json.dumps(drafts, indent=2, ensure_ascii=False), encoding="utf-8")


def _list_tasks() -> list[dict[str, Any]]:
    tasks = []
    for path in sorted(DESIGN_DIR.glob("*.json"), key=lambda p: p.stat().st_mtime, reverse=True):
        # Skip drafts files and legacy report files
        if path.stem.endswith("_drafts") or path.stem.endswith("_report"):
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            if "id" not in data or "status" not in data:
                continue
            # Annotate with current draft count
            data["n_drafts"] = len(_load_drafts(data["id"]))
            # Mirror _load_task back-compat: legacy task files are all
            # TCR-affinity. _list_tasks reads the JSON itself so we apply
            # the fill here too.
            data.setdefault("task", "tcr_affinity")
            tasks.append(data)
        except (json.JSONDecodeError, OSError) as e:
            logger.warning(f"Failed to load task {path}: {e}")
    return tasks


def _extract_design_context(
    source_job_run_dir: Path,
    design_goals: list[str],
    constraints: dict,
    task_id: str,
    task: TaskSpec | None = None,
) -> dict[str, Any]:
    """Build the design context dict from the source job's analysis data.

    This will be injected into the Agent session's design_context so the
    system prompt knows about the source system.

    The ``task`` argument controls which side of the interface gets surfaced
    as design candidates and whether anchor-pocket metadata is rendered:
        - TaskSpec.tcr_affinity()        → TCR-side candidate filter (default,
                                           matches the legacy behavior of this
                                           function before D-B7).
        - TaskSpec.peptide_presentation() → peptide-side candidate filter +
                                           anchor pocket chemistry lookup
                                           (so the LLM sees P2/PΩ before it
                                           picks a mutation site).
    """
    if task is None:
        task = TaskSpec.tcr_affinity()

    context: dict[str, Any] = {
        "task_id": task_id,
        # Expose the task key on the context so downstream code (system
        # prompt, frontend) doesn't have to re-derive it from the
        # (subject, objective, target_subinterface) triple.
        "task": _task_key_for_spec(task),
        "task_subject": task.subject,
        "task_objective": task.objective,
        "design_goals": design_goals,
        "constraints": constraints or {},
        # IMPORTANT: case_dir for the Agent's query_analysis_results tool.
        # This tells the Agent exactly which directory to query — no need to
        # guess or list files.
        "case_dir": str(source_job_run_dir),
    }

    # Source job info
    run_summary_path = source_job_run_dir / "run_summary.json"
    if run_summary_path.exists():
        try:
            run_summary = json.loads(run_summary_path.read_text(encoding="utf-8"))
            context["system_id"] = run_summary.get("job_id", "unknown")
        except json.JSONDecodeError:
            pass

    # Identity (peptide / HLA / TCR)
    identity_file = source_job_run_dir / "analysis" / "identity" / "analysis" / "identity" / "biological_identity.json"
    if identity_file.exists():
        try:
            identity = json.loads(identity_file.read_text(encoding="utf-8"))
            context["peptide"] = identity.get("peptide", {}).get("sequence")
            hla = identity.get("hla", {})
            if hla:
                context["hla"] = hla.get("allele") or hla.get("name")
            tcr = identity.get("tcr", {})
            if tcr:
                a = tcr.get("alpha_v") or tcr.get("trav")
                b = tcr.get("beta_v") or tcr.get("trbv")
                if a or b:
                    context["tcr_genes"] = f"{a or '?'} / {b or '?'}"
        except (json.JSONDecodeError, OSError):
            pass

    # Quality
    quality_file = source_job_run_dir / "preparation" / "quality" / "preprocess_quality_report.json"
    if quality_file.exists():
        try:
            q = json.loads(quality_file.read_text(encoding="utf-8"))
            context["n_frames"] = q.get("n_frames")
            context["duration_ns"] = round(q.get("tail90_end_time_ps", 0) / 1000, 1)
            context["rmsd"] = round(q.get("tail90_mean_rmsd_nm", 0), 3)
            std = q.get("tail90_std_rmsd_nm", 0)
            context["convergence"] = (
                "excellent" if std < 0.05 else
                "good" if std < 0.1 else
                "moderate" if std < 0.15 else "unstable"
            )
        except (json.JSONDecodeError, OSError):
            pass

    # Top hotspots (from RRCS annotated pair summary) +
    # build a chain mapping from component label → PDB chain ID
    # so the 3D viewer can look up residues in the correct chain.
    rrcs_csv = source_job_run_dir / "analysis" / "rrcs" / "analysis" / "interactions" / "rrcs" / "annotated_rrcs_pair_summary.csv"
    chain_mapping: dict[str, str] = {}
    if rrcs_csv.exists():
        try:
            import csv
            rows = []
            with rrcs_csv.open(encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    rows.append(row)
            # Build chain map from any column pair (chain_id_1/component_1 and *_2)
            for row in rows:
                for side in ("1", "2"):
                    ch = row.get(f"chain_id_{side}", "").strip()
                    comp = row.get(f"component_{side}", "").strip()
                    if ch and comp and comp not in chain_mapping:
                        chain_mapping[comp] = ch
                    # Also map by tcr_chain ("alpha"/"beta") when component is TCR_*
                    if comp.startswith("TCR_") and ch:
                        sub = comp.replace("TCR_", "")  # "alpha" or "beta"
                        if sub and sub not in chain_mapping:
                            chain_mapping[sub] = ch

            # Sort by mean_rrcs and take top 10 unique TCR residues
            rows.sort(key=lambda r: float(r.get("mean_rrcs", 0) or 0), reverse=True)
            seen = set()
            hotspots = []
            for row in rows:
                tcr_res = row.get("tcr_residue_label", "")
                if not tcr_res or tcr_res in seen:
                    continue
                seen.add(tcr_res)
                hotspots.append({
                    "residue": tcr_res,
                    "region": row.get("tcr_region_detailed") or row.get("tcr_region", "-"),
                    "rrcs": round(float(row.get("mean_rrcs", 0) or 0), 2),
                    "occupancy": round(float(row.get("rrcs_nonzero_fraction", 0) or 0), 2),
                    "partner": row.get("partner_component", "-"),
                    "partner_residue": row.get("partner_residue_label", "-"),
                    # Add PDB chain info for 3D viewer
                    "pdb_chain": row.get("chain_id_1", ""),
                    "pdb_resid": row.get("resid_1", ""),
                    "partner_pdb_chain": row.get("chain_id_2", ""),
                    "partner_pdb_resid": row.get("resid_2", ""),
                })
                if len(hotspots) >= 10:
                    break
            context["top_hotspots"] = hotspots
        except Exception as e:
            logger.warning(f"Failed to load hotspots: {e}")
    context["chain_mapping"] = chain_mapping

    # Build a direct residue label → PDB chain/resid lookup table.
    # This is the most reliable way for the 3D viewer to resolve any
    # residue label that appears in MD analysis output (it lists every
    # residue that participates in any interface contact).
    residue_lookup: dict[str, dict[str, Any]] = {}
    if rrcs_csv.exists():
        try:
            with rrcs_csv.open(encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    for side in ("1", "2"):
                        label = row.get(f"residue_label_{side}", "").strip()
                        chain = row.get(f"chain_id_{side}", "").strip()
                        resid = row.get(f"resid_{side}", "").strip()
                        comp = row.get(f"component_{side}", "").strip()
                        if label and chain and resid and label not in residue_lookup:
                            residue_lookup[label] = {
                                "pdb_chain": chain,
                                "pdb_resid": int(resid) if resid.isdigit() else resid,
                                "component": comp,
                            }
        except Exception as e:
            logger.warning(f"Failed to build residue lookup: {e}")
    context["residue_lookup"] = residue_lookup

    # Already-saved drafts (so Agent knows what's been committed)
    context["current_drafts"] = _load_drafts(task_id)

    # Multi-agent pipeline result (decision A: the pipeline is the recommender).
    # When a run has produced its artifact, surface the authoritative picks +
    # critic verdict so the chat agent EXPLAINS them rather than re-deriving its
    # own recommendations. A trimmed view keeps the prompt bounded.
    pipeline_path = DESIGN_DIR / f"{task_id}_pipeline.json"
    if pipeline_path.exists():
        try:
            pj = json.loads(pipeline_path.read_text(encoding="utf-8"))
            rec_msg = pj.get("recommendation") or {}
            rec_payload = rec_msg.get("structured_payload") or {}
            crit = pj.get("critique") or {}
            context["pipeline_result"] = {
                "succeeded": pj.get("succeeded"),
                "design_strategy_in_use": rec_payload.get("design_strategy_in_use"),
                "recommendations": rec_payload.get("recommendations", []),
                "skipped_candidates": rec_payload.get("skipped_candidates", []),
                "global_caveats": rec_payload.get("global_caveats", []),
                "critic_verdict": crit.get("verdict"),
                "critic_findings": crit.get("findings", []),
                "deterministic_findings": pj.get("deterministic_findings", []),
                "recommendation_narrative": rec_msg.get("narrative", ""),
            }
        except (json.JSONDecodeError, OSError) as e:
            logger.warning(f"Failed to load pipeline artifact for {task_id}: {e}")

    # Docking angle summary (if angles module was run)
    angles_summary_path = source_job_run_dir / "analysis" / "angles" / "docking_angles_summary.json"
    if angles_summary_path.exists():
        try:
            ang = json.loads(angles_summary_path.read_text(encoding="utf-8"))
            cross = ang.get("crossing_angle") or {}
            inc = ang.get("incident_angle") or {}
            context["docking_angles"] = {
                "crossing_mean": round(cross.get("mean", 0), 2) if cross else None,
                "crossing_std": round(cross.get("std", 0), 2) if cross else None,
                "incident_mean": round(inc.get("mean", 0), 2) if inc else None,
                "incident_std": round(inc.get("std", 0), 2) if inc else None,
            }
        except (json.JSONDecodeError, OSError):
            pass

    # Design candidates ranked by RRCS contribution (single deterministic
    # signal, matching the Pure-RRCS baseline). Per D-B5 (2026-05-26) the
    # composite design_priority_score is no longer surfaced to the LLM;
    # chemistry tags and risk flags are categorical evidence rather than
    # weights in a scalar.
    #
    # The candidate filter branches on task.subject:
    #   - tcr     → keep CDR/framework/non_cdr residues + cross-interface
    #               TCR↔(hla|peptide) contacts. Mirrors the inline filter
    #               that lived here before D-B7.
    #   - peptide → keep peptide-side residues whose partner is HLA/β2m.
    #               Excludes peptide-TCR contacts (those are TCR-affinity
    #               territory, not presentation).
    # Centralising the filter on AnalysisDataFormatter would be cleaner, but
    # the formatter's helper consumes `ResidueFeatures` items and applies
    # the same `region` / `partner_chemistry` checks — duplicating the
    # 6-line predicate here is acceptable while we still have two callers
    # (the agent path and the legacy batch engine kept around as fallback).
    try:
        from immunoscope.analysis.features import (
            CaseLocator as _FeatureLocator,
            FeatureSet as _FeatureSet,
            compute_feature as _compute_feature,
        )

        _fl = _FeatureLocator(source_job_run_dir)
        _fs = _FeatureSet(case_id=_fl.get_case_id())
        for _name in (
            "contact_count",
            "rrcs_contribution",
            "chemistry_tags",
            "risk_flag",
            "persistence_profile",
        ):
            _compute_feature(_name, _fl, _fs)

        def _accept(_item) -> bool:
            _region = (_item.region or "").lower()
            _partner = (_item.partner_chemistry or "").lower()
            if task.subject == "tcr":
                if "cdr" in _region or "framework" in _region or "non_cdr" in _region:
                    return True
                # Reject anything whose region clearly marks a non-TCR component.
                if any(tag in _region for tag in ("peptide", "hla", "mhc", "beta2m", "b2m")):
                    return False
                return ("hla" in _partner or "peptide" in _partner) and "tcr" not in _partner
            if task.subject == "peptide":
                if "peptide" not in _region:
                    return False
                # Peptide↔TCR is a presentation-irrelevant contact.
                if "tcr" in _partner:
                    return False
                return True
            return False

        _candidates_out: list[dict[str, Any]] = []
        for _item in _fs.top_by("rrcs_contribution", n=30):
            if not _accept(_item):
                continue
            _candidates_out.append({
                "residue": f"{_item.residue.resname}{_item.residue.resid}",
                "chain": _item.residue.chain,
                "region": _item.region or "",
                "rrcs": round(_item.rrcs_contribution or 0.0, 2),
                "occupancy": round(_item.persistence_mean or 0.0, 2),
                "partner": _item.partner_chemistry or "",
                "chemistry": list(_item.chemistry_tags[:3]),
                "risk_flags": list(_item.risk_flags[:2]),
            })
            if len(_candidates_out) >= 10:
                break
        if _candidates_out:
            context["candidates"] = _candidates_out
    except Exception as e:
        logger.warning(f"Failed to load design candidates: {e}")

    # Anchor pocket chemistry — presentation track only.
    # The LLM cannot reason about peptide mutations safely without knowing
    # which positions are P2/PΩ anchors and what residues are tolerated
    # there. We piggyback on the same lookup table used by the (now-retired)
    # batch formatter so both paths surface identical chemistry to the LLM.
    if task.subject == "peptide" and task.objective == "presentation":
        try:
            from immunoscope.recommendation.anchor_pocket import (
                lookup_anchor_pockets,
            )

            hla_allele = context.get("hla")
            peptide_seq = (context.get("peptide") or "").strip() or None
            peptide_length = len(peptide_seq) if peptide_seq else None
            entries = lookup_anchor_pockets(hla_allele, peptide_length)
            context["anchor_pockets"] = [
                {
                    "position": e.position,
                    "allele": e.allele,
                    "peptide_length": e.peptide_length,
                    "role": e.role,
                    "preferred_residues": list(e.preferred_residues),
                    "tolerated_residues": list(e.tolerated_residues),
                    "dispreferred_residues": list(e.dispreferred_residues),
                    "notes": e.notes,
                    "source": e.source,
                    "current_residue": (
                        peptide_seq[e.position - 1]
                        if peptide_seq and 1 <= e.position <= len(peptide_seq)
                        else None
                    ),
                }
                for e in entries
            ]
        except Exception as e:
            logger.warning(f"Failed to look up anchor pockets: {e}")
            context["anchor_pockets"] = []

    # Backbone Ramachandran summary (if dihedrals module was run)
    dih_summary_path = source_job_run_dir / "analysis" / "dihedrals" / "dihedrals_summary.json"
    if dih_summary_path.exists():
        try:
            dih = json.loads(dih_summary_path.read_text(encoding="utf-8"))
            context["dihedrals"] = {
                "mean_angular_spread_deg": dih.get("mean_angular_spread_deg"),
                "n_residues": dih.get("n_residues_analyzed"),
                # Map residue label → angular spread + dominant region for
                # quick lookup during mutation site reasoning.
                "residue_index": {},
            }
            # Build per-residue index for the LLM (most-flexible 20)
            residue_csv = source_job_run_dir / "analysis" / "dihedrals" / "residue_dihedrals.csv"
            if residue_csv.exists():
                import csv as _csv
                with residue_csv.open(encoding="utf-8") as f:
                    rows = list(_csv.DictReader(f))
                for row in rows:
                    label = f"{row.get('resname','')}{row.get('resid','')}".strip()
                    if not label:
                        continue
                    try:
                        context["dihedrals"]["residue_index"][label] = {
                            "angular_spread_deg": round(float(row.get("angular_spread_deg", 0)), 1),
                            "dominant_region": row.get("dominant_region", ""),
                            "phi_mean": round(float(row.get("phi_mean_deg", 0)), 1),
                            "psi_mean": round(float(row.get("psi_mean_deg", 0)), 1),
                        }
                    except (ValueError, TypeError):
                        pass
        except (json.JSONDecodeError, OSError):
            pass

    return context


# ----------------------------------------------------------------------------
# Routes
# ----------------------------------------------------------------------------

@router.get("/jobs")
async def list_design_tasks() -> list[DesignTaskSummary]:
    """List all design tasks."""
    tasks = _list_tasks()
    return [DesignTaskSummary(**t) for t in tasks]


@router.post("/jobs")
async def create_design_task(request: DesignTaskRequest) -> DesignTaskSummary:
    """Create a new Design Copilot task (conversational session)."""
    source_job = await job_store.get(request.source_job_id)
    if not source_job:
        raise HTTPException(status_code=404, detail=f"Source job not found: {request.source_job_id}")

    status_val = source_job.status.value if hasattr(source_job.status, "value") else str(source_job.status)
    if status_val != "completed":
        raise HTTPException(
            status_code=400,
            detail=f"Source job must be completed (current: {status_val})",
        )

    # Resolve task key now so a malformed request fails fast at create
    # time, not deep in the agent path. The resolved TaskSpec isn't
    # serialized here — we persist the wire-level key and let downstream
    # consumers re-resolve via _task_spec_from_key.
    _task_spec_from_key(request.task)

    task_id = uuid.uuid4().hex[:12]
    now = datetime.now(timezone.utc).isoformat()
    task = {
        "id": task_id,
        "name": request.name or f"design_{task_id}",
        "status": "active",
        "source_job_id": request.source_job_id,
        "source_job_name": source_job.name,
        "source_run_dir": source_job.run_dir,
        "task": request.task,
        "design_goals": request.design_goals,
        "constraints": request.constraints,
        "notes": request.notes,
        "created_at": now,
        "updated_at": now,
        "n_messages": 0,
    }
    _save_task(task)
    _save_drafts(task_id, [])  # initialize empty drafts file

    return DesignTaskSummary(
        id=task_id,
        name=task["name"],
        status=task["status"],
        source_job_id=task["source_job_id"],
        source_job_name=task["source_job_name"],
        task=task["task"],
        design_goals=task["design_goals"],
        created_at=now,
        updated_at=now,
        n_drafts=0,
        n_messages=0,
        notes=task["notes"],
    )


@router.get("/jobs/{task_id}")
async def get_design_task(task_id: str) -> dict[str, Any]:
    """Get design task details, including its context and saved drafts."""
    task = _load_task(task_id)
    if not task:
        raise HTTPException(status_code=404, detail=f"Design task not found: {task_id}")

    drafts = _load_drafts(task_id)
    task["n_drafts"] = len(drafts)
    return {
        "task": task,
        "drafts": drafts,
    }


@router.get("/jobs/{task_id}/context")
async def get_design_context(task_id: str) -> dict[str, Any]:
    """Return the design context that will be injected into the Agent session.

    This is what the Agent sees as pre-loaded MD evidence about the source
    system. Frontend can use this to render the left-side context panel.
    """
    task = _load_task(task_id)
    if not task:
        raise HTTPException(status_code=404, detail=f"Design task not found: {task_id}")

    run_dir = task.get("source_run_dir")
    if not run_dir:
        source_job = await job_store.get(task["source_job_id"])
        if source_job:
            run_dir = source_job.run_dir
    if not run_dir:
        raise HTTPException(status_code=500, detail="Source run directory not available")

    # Resolve the persisted wire-level key back to a TaskSpec. _load_task
    # already filled the "tcr_affinity" default for legacy task files, so
    # this never sees a missing key.
    task_spec = _task_spec_from_key(task["task"])
    context = _extract_design_context(
        Path(run_dir),
        task.get("design_goals", []),
        task.get("constraints", {}),
        task_id,
        task=task_spec,
    )
    return context


@router.get("/jobs/{task_id}/drafts")
async def get_drafts(task_id: str) -> list[dict[str, Any]]:
    """Get current draft recommendations for this task."""
    task = _load_task(task_id)
    if not task:
        raise HTTPException(status_code=404, detail=f"Design task not found: {task_id}")
    return _load_drafts(task_id)


@router.delete("/jobs/{task_id}/drafts/{draft_idx}")
async def delete_draft(task_id: str, draft_idx: int) -> dict[str, Any]:
    """Remove a single draft recommendation."""
    drafts = _load_drafts(task_id)
    if draft_idx < 0 or draft_idx >= len(drafts):
        raise HTTPException(status_code=404, detail=f"Draft index out of range")
    drafts.pop(draft_idx)
    _save_drafts(task_id, drafts)
    return {"status": "ok", "n_drafts": len(drafts)}


def _resolve_structure_pdb(run_dir: Path) -> Path | None:
    """Locate the analysis structure PDB for 3D viewing across run layouts.

    Web-prepared runs put it at ``preparation/analysis_input/analysis_structure.pdb``;
    ``--skip-preparation`` runs (e.g. imported benchmark analyses) put it at
    ``input/prepared_structure.pdb``. The authoritative path is whatever the run
    actually used, recorded in ``run_summary.json``'s ``prepared_input``. Try the
    layout-relative candidates first, then fall back to the recorded absolute
    paths, so any run with a structure on disk renders rather than 404-ing."""
    candidates = [
        run_dir / "preparation" / "analysis_input" / "analysis_structure.pdb",
        run_dir / "input" / "prepared_structure.pdb",
        run_dir / "input" / "raw_topology.pdb",
    ]
    summ = run_dir / "run_summary.json"
    if summ.exists():
        try:
            prepared = json.loads(summ.read_text(encoding="utf-8")).get("prepared_input", {}) or {}
            for key in ("structure", "prepared_structure", "topology"):
                val = prepared.get(key)
                if val:
                    candidates.append(Path(val))
        except (json.JSONDecodeError, OSError):
            pass
    return next((c for c in candidates if c.exists()), None)


@router.get("/jobs/{task_id}/structure")
async def get_design_structure(task_id: str):
    """Return the source job's prepared analysis structure (PDB) for 3D viewing."""
    from fastapi.responses import FileResponse
    task = _load_task(task_id)
    if not task:
        raise HTTPException(status_code=404, detail=f"Design task not found: {task_id}")
    run_dir = task.get("source_run_dir")
    if not run_dir:
        source_job = await job_store.get(task["source_job_id"])
        if source_job:
            run_dir = source_job.run_dir
    if not run_dir:
        raise HTTPException(status_code=500, detail="Source run directory not available")
    pdb_path = _resolve_structure_pdb(Path(run_dir))
    if pdb_path is None:
        raise HTTPException(
            status_code=404,
            detail=f"No analysis structure PDB found under {run_dir}",
        )
    return FileResponse(str(pdb_path), media_type="chemical/x-pdb", filename="structure.pdb")


@router.delete("/jobs/{task_id}")
async def delete_design_task(task_id: str) -> dict[str, str]:
    """Delete a design task and its drafts."""
    task = _load_task(task_id)
    if not task:
        raise HTTPException(status_code=404, detail=f"Design task not found: {task_id}")

    for p in [_task_path(task_id), _drafts_path(task_id),
              DESIGN_DIR / f"{task_id}_report.json",
              DESIGN_DIR / f"{task_id}_pipeline.json",
              DESIGN_DIR / f"{task_id}_pipeline_status.json"]:
        if p.exists():
            p.unlink()
    return {"status": "deleted", "id": task_id}


# ----------------------------------------------------------------------------
# Multi-agent pipeline run (decision A: the pipeline is the recommender; the
# chat agent only explains its output). These three endpoints drive a
# background run of immunoscope.agent.multi_agent over the task's source MD
# evidence, expose poll-able progress, and return the full trace artifact.
# The runner is imported lazily to avoid a circular import.
# ----------------------------------------------------------------------------

@router.post("/jobs/{task_id}/run_pipeline")
async def run_design_pipeline_endpoint(task_id: str) -> dict[str, Any]:
    """Kick off (or report already-running) the multi-agent pipeline."""
    task = _load_task(task_id)
    if not task:
        raise HTTPException(status_code=404, detail=f"Design task not found: {task_id}")

    run_dir = task.get("source_run_dir")
    if not run_dir:
        source_job = await job_store.get(task["source_job_id"])
        if source_job:
            run_dir = source_job.run_dir
    if not run_dir or not Path(run_dir).exists():
        raise HTTPException(
            status_code=400,
            detail="Source analysis run directory is unavailable; cannot run pipeline.",
        )

    from immunoscope.web.services.design_pipeline_runner import (
        is_running, start_design_pipeline,
    )
    if is_running(task_id):
        return {"status": "already_running", "task_id": task_id}
    result = start_design_pipeline(task_id)
    return {"status": result["status"], "task_id": task_id}


@router.get("/jobs/{task_id}/run_pipeline/status")
async def get_design_pipeline_status(task_id: str) -> dict[str, Any]:
    """Poll the multi-agent pipeline run status (stage / progress / events)."""
    from immunoscope.web.services.design_pipeline_runner import read_status
    return read_status(task_id)


@router.get("/jobs/{task_id}/pipeline")
async def get_design_pipeline_artifact(task_id: str) -> dict[str, Any]:
    """Return the full multi-agent trace artifact (briefs / candidate table /
    recommendation / critic / deterministic findings / numeric audit)."""
    from immunoscope.web.services.design_pipeline_runner import (
        read_artifact, read_status,
    )
    artifact = read_artifact(task_id)
    if artifact is None:
        return {"available": False, "status": read_status(task_id).get("status", "idle")}
    return {"available": True, "artifact": artifact}


# ----------------------------------------------------------------------------
# Helper for the save_recommendation tool to append drafts
# ----------------------------------------------------------------------------

def append_draft(task_id: str, recommendation: dict[str, Any]) -> int:
    """Append a recommendation to the task's draft list.

    Called by the `save_recommendation` Agent tool. Returns the new total
    count of drafts.
    """
    drafts = _load_drafts(task_id)
    # Add timestamp
    recommendation["saved_at"] = datetime.now(timezone.utc).isoformat()
    drafts.append(recommendation)
    _save_drafts(task_id, drafts)

    # Update task updated_at
    task = _load_task(task_id)
    if task:
        task["updated_at"] = recommendation["saved_at"]
        task["n_drafts"] = len(drafts)
        _save_task(task)

    return len(drafts)
