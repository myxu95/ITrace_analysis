"""Run the two-round mutation-design pipeline end-to-end (v1, 2026-05-28).

Round 1 is a parallel fan-out of four agents:

    BioAgent              — literature-grounded biological brief.
    InteractionReader     — hotspots / pair / fingerprint.
    ConformationReader    — flexibility / dihedrals.
    InterfaceExposureReader — interface (BSA) + per-residue SASA / cavities.

Before the fan-out the orchestrator pre-renders the `overview` and
`quality` views once and parks them in `ctx.system_preamble`, so every
agent shares a baseline picture of the system + trajectory quality
without each calling those views.

Round 2 is a single integrator, the RecommendationAgent. After the
round-1 fan-out completes, the orchestrator deterministically
pre-aggregates a candidate-site table from the four `structured_payload`s
(`_aggregate_candidates` → `ctx.recommendation_preamble`) and runs the
RecommendationAgent, which folds the four briefs + the candidate table
into one recommendation `AgentMessage` (site + direction class, no
specific amino acid).

Failure handling:
  - Any round-1 agent that raises is captured in
    `PipelineResult.errors[<name>]` and the rest of the run continues —
    one bad agent should not blank out the other three briefs.
  - Round 2 runs even if some round-1 briefs are missing; the candidate
    aggregation simply uses whatever payloads are present.
  - `succeeded` is true only if all four round-1 agents AND the
    round-2 agent produced a message with no recorded errors.
  - Validation failures (`MessageValidationError`) inside an agent
    bubble up as that agent's error — the agent prompt is what needs
    fixing.
"""

from __future__ import annotations

import asyncio
import logging
from dataclasses import dataclass, field
from typing import Callable, Dict, Optional

from immunoscope.agent.llm import LLMClient
from immunoscope.agent.presenters import CaseLocator, PresenterError, render

from . import critic_checks
from .base import AgentRunResult, BaseAgent
from .context import MutationDesignContext
from .residue_id import ChainAliasMap, parse_residue_label

log = logging.getLogger("immunoscope.agent.multi_agent.orchestrator")


_PREAMBLE_VIEWS = ("overview", "quality")


class EmptyRecommendationError(RuntimeError):
    """The pipeline finished but the recommendation carries zero sites.

    Recorded in ``PipelineResult.errors`` so an empty plan reads as a failed
    run (``succeeded == False``) instead of masquerading as a clean success.
    """


@dataclass
class PipelineResult:
    """Outcome of one two-round pipeline run."""

    context: MutationDesignContext
    bio_result: Optional[AgentRunResult] = None
    interaction_result: Optional[AgentRunResult] = None
    conformation_result: Optional[AgentRunResult] = None
    interface_exposure_result: Optional[AgentRunResult] = None
    recommendation_result: Optional[AgentRunResult] = None
    critique_result: Optional[AgentRunResult] = None
    # Deterministic, model-free findings on the final recommendation.
    deterministic_findings: list = field(default_factory=list)
    errors: Dict[str, BaseException] = field(default_factory=dict)
    # True when the LLM critic ran but never produced a usable verdict (after
    # retries) — the run was NOT cleanly approved, so callers shouldn't treat a
    # silent critic as an endorsement.
    critic_unverified: bool = False

    @property
    def succeeded(self) -> bool:
        return (
            not self.errors
            and self.bio_result is not None
            and self.interaction_result is not None
            and self.conformation_result is not None
            and self.interface_exposure_result is not None
            and self.recommendation_result is not None
        )

    def recommendation_count(self) -> int:
        """Number of sites in the final recommendation payload (0 if none)."""
        rec = self.recommendation_result
        msg = getattr(rec, "message", None) if rec is not None else None
        payload = (getattr(msg, "structured_payload", None) or {}) if msg else {}
        recs = payload.get("recommendations") or []
        return len(recs) if isinstance(recs, list) else 0


class MutationDesignPipeline:
    """Holds one instance of each round-1 agent and runs them in parallel.

    Construction is decoupled from execution: the same pipeline object
    can run against multiple `MutationDesignContext` instances, which
    matters when a calling script wants to sweep design intents over
    the same cached LLM client.
    """

    def __init__(
        self,
        *,
        bio_agent: BaseAgent,
        interaction_reader: BaseAgent,
        conformation_reader: BaseAgent,
        interface_exposure_reader: BaseAgent,
        recommendation_agent: BaseAgent,
        design_critic: Optional[BaseAgent] = None,
    ) -> None:
        self._bio = bio_agent
        self._interaction = interaction_reader
        self._conformation = conformation_reader
        self._interface_exposure = interface_exposure_reader
        self._recommendation = recommendation_agent
        # Optional round-2.5 reviewer. When present, it critiques the
        # recommendation and may trigger one bounded revision.
        self._critic = design_critic

    async def run(
        self,
        ctx: MutationDesignContext,
        progress_cb: Optional[Callable[[dict], None]] = None,
    ) -> PipelineResult:
        """Run the full pipeline.

        ``progress_cb`` (optional) is called with a small dict at each stage
        boundary so a caller (e.g. the web runner) can stream live progress.
        Every callback is wrapped in a guard — a broken callback never
        interferes with the pipeline. Events:
            {"stage":"start","case_id":...}
            {"stage":"round1","agent":<slot>,"status":"done"|"failed"}
            {"stage":"round1_done"}
            {"stage":"round2","status":"done"|"failed"}
            {"stage":"critic","pass":N,"llm_verdict":...,"det_blocking":bool,
             "decision":"revise"|"approve"}
            {"stage":"revision","status":"done"|"failed"}
            {"stage":"done","succeeded":bool}
        """
        result = PipelineResult(context=ctx)
        log.info("multi_agent round-1 start: case_id=%s", ctx.case_id)
        self._emit(progress_cb, {"stage": "start", "case_id": ctx.case_id})

        # Best-effort preamble. If a view can't render (missing module,
        # etc.) we skip it rather than aborting — agents have their own
        # views to fall back on.
        ctx.system_preamble = self._build_preamble(ctx)

        agents = (
            ("bio", self._bio),
            ("interaction", self._interaction),
            ("conformation", self._conformation),
            ("interface_exposure", self._interface_exposure),
        )

        outcomes = await asyncio.gather(
            *(agent.run(ctx) for _, agent in agents),
            return_exceptions=True,
        )

        for (slot, agent), outcome in zip(agents, outcomes):
            if isinstance(outcome, BaseException):
                log.error("agent %s failed: %s", agent.name, outcome)
                result.errors[agent.name] = outcome
                self._emit(progress_cb,
                           {"stage": "round1", "agent": slot, "status": "failed"})
                continue
            if slot == "bio":
                result.bio_result = outcome
                ctx.bio_msg_id = outcome.msg_id
            elif slot == "interaction":
                result.interaction_result = outcome
                ctx.interaction_msg_id = outcome.msg_id
            elif slot == "conformation":
                result.conformation_result = outcome
                ctx.conformation_msg_id = outcome.msg_id
            elif slot == "interface_exposure":
                result.interface_exposure_result = outcome
                ctx.interface_exposure_msg_id = outcome.msg_id
            log.info("agent %s emitted msg_id=%s", agent.name, outcome.msg_id)
            self._emit(progress_cb,
                       {"stage": "round1", "agent": slot, "status": "done"})

        self._emit(progress_cb, {"stage": "round1_done"})

        # Deterministic numeric-fidelity flags from round-1, handed to the
        # critic so it can escalate a fabricated number to a blocking finding.
        ctx.numeric_audit_summary = self._build_numeric_audit_summary(result)

        # ---- round 2: candidate aggregation + single integrator --------
        # Runs even if some round-1 briefs are missing; the aggregation
        # uses whatever payloads landed in ctx.message_history.
        log.info("multi_agent round-2 start: case_id=%s", ctx.case_id)
        alias_map = ChainAliasMap.from_case_dir(ctx.case_dir)
        ctx.recommendation_preamble = self._aggregate_candidates(ctx, alias_map)
        # The integrator occasionally burns its whole turn budget without
        # emitting (a transient structured-output miss). Retry once before
        # giving up so a single stochastic stall doesn't blank the run.
        for attempt in range(1, self._ROUND2_ATTEMPTS + 1):
            try:
                rec = await self._recommendation.run(ctx)
                result.recommendation_result = rec
                ctx.recommendation_msg_id = rec.msg_id
                result.errors.pop(self._recommendation.name, None)  # clear prior fail
                log.info(
                    "agent %s emitted msg_id=%s (attempt %d/%d)",
                    self._recommendation.name, rec.msg_id, attempt,
                    self._ROUND2_ATTEMPTS,
                )
                self._emit(progress_cb, {"stage": "round2", "status": "done"})
                break
            except BaseException as exc:  # noqa: BLE001
                log.error("agent %s failed (attempt %d/%d): %s",
                          self._recommendation.name, attempt,
                          self._ROUND2_ATTEMPTS, exc)
                result.errors[self._recommendation.name] = exc
                self._emit(progress_cb, {
                    "stage": "round2", "status": "failed", "attempt": attempt,
                })

        # ---- round 2.5: deterministic checks + adversarial critic ---------
        # Each pass: run model-free constraint/coverage checks on the current
        # recommendation, then the LLM critic (which is told about both the
        # deterministic findings and the numeric-audit flags). A revision is
        # triggered when EITHER the LLM critic says `revise` OR a deterministic
        # `blocking` finding exists — so a hard-constraint violation forces a
        # fix even if the same-model critic missed it. Bounded to
        # `_MAX_CRITIC_PASSES` so a revision is re-checked but cannot oscillate.
        if self._critic is not None and result.recommendation_result is not None:
            log.info("multi_agent round-2.5 critic start: case_id=%s", ctx.case_id)
            for cpass in range(self._MAX_CRITIC_PASSES):
                det = critic_checks.run_deterministic_checks(ctx, alias_map)
                result.deterministic_findings = det
                ctx.deterministic_findings_summary = critic_checks.summarize(det)
                det_blocking = critic_checks.has_blocking(det)
                if det_blocking:
                    log.info("deterministic blocking finding(s): %s",
                             critic_checks.summarize(det))
                crit, verdict, crit_exc = await self._run_critic_for_verdict(ctx)
                if crit is not None:
                    result.critique_result = crit
                    ctx.critique_msg_id = crit.msg_id
                if crit_exc is not None:
                    # Critic threw — advisory, must not blank a recommendation
                    # that already succeeded. Defer to deterministic checks.
                    result.errors[self._critic.name] = crit_exc
                # A critic that EMITTED but produced no usable verdict (after
                # retries) is NOT an approval — flag it and force a review pass
                # instead of letting silence rubber-stamp the plan.
                emitted_no_verdict = crit is not None and verdict is None
                if emitted_no_verdict:
                    result.critic_unverified = True
                should_revise = (verdict == "revise") or det_blocking or emitted_no_verdict
                log.info("critic pass %d: llm_verdict=%s det_blocking=%s -> %s",
                         cpass + 1, verdict, det_blocking,
                         "revise" if should_revise else "approve")
                self._emit(progress_cb, {
                    "stage": "critic", "pass": cpass + 1,
                    "llm_verdict": verdict, "det_blocking": det_blocking,
                    "decision": "revise" if should_revise else "approve",
                })
                if not should_revise:
                    break
                if cpass == self._MAX_CRITIC_PASSES - 1:
                    log.info("still revise after %d passes — accepting latest",
                             self._MAX_CRITIC_PASSES)
                    break
                if crit is None and not det_blocking:
                    break  # nothing actionable to revise against
                try:
                    rev = await self._recommendation.run(ctx)  # revision mode
                    result.recommendation_result = rev
                    ctx.recommendation_msg_id = rev.msg_id
                    log.info("recommendation revised: msg_id=%s", rev.msg_id)
                    self._emit(progress_cb,
                               {"stage": "revision", "status": "done"})
                except BaseException as exc:  # noqa: BLE001
                    log.error("revision failed: %s", exc)
                    result.errors[self._recommendation.name] = exc
                    self._emit(progress_cb,
                               {"stage": "revision", "status": "failed"})
                    break

        # ---- empty-plan guard ---------------------------------------------
        # If the run produced a recommendation object but it carries zero sites
        # (even after the revision pass), record it as a failure so it can't
        # masquerade as a clean success downstream (`succeeded` -> False).
        if (result.recommendation_result is not None
                and result.recommendation_count() == 0):
            log.error("recommendation is empty (0 sites) for case_id=%s — "
                      "marking run failed", ctx.case_id)
            result.errors.setdefault(
                self._recommendation.name,
                EmptyRecommendationError(
                    "recommendation produced zero sites after "
                    f"{self._MAX_CRITIC_PASSES} critic pass(es)"
                ),
            )
            self._emit(progress_cb, {"stage": "round2", "status": "empty"})

        self._emit(progress_cb, {"stage": "done", "succeeded": result.succeeded})
        return result

    @staticmethod
    def _emit(progress_cb: Optional[Callable[[dict], None]], event: dict) -> None:
        """Best-effort progress callback. A broken callback must never break
        the pipeline, so every call is guarded."""
        if progress_cb is None:
            return
        try:
            progress_cb(event)
        except Exception as exc:  # noqa: BLE001
            log.debug("progress_cb raised (ignored): %s", exc)

    # Max critic passes per run (a revision is run between consecutive passes).
    _MAX_CRITIC_PASSES = 2
    # Attempts for the round-2 integrator (1 retry) and for getting a usable
    # verdict out of the critic within a single pass (1 retry).
    _ROUND2_ATTEMPTS = 2
    _CRITIC_ATTEMPTS = 2

    _VALID_VERDICTS = ("approve", "revise")

    async def _run_critic_for_verdict(self, ctx):
        """Run the LLM critic, retrying once if it emits no usable verdict.

        Returns ``(crit, verdict, exc)`` where ``verdict`` is a normalized
        member of ``_VALID_VERDICTS`` or ``None`` (critic emitted but gave no
        parseable verdict), and ``exc`` is set only if the critic *threw*. The
        retry targets the transient structured-output miss that previously let
        a silent critic default to approval.
        """
        crit = None
        for attempt in range(1, self._CRITIC_ATTEMPTS + 1):
            try:
                crit = await self._critic.run(ctx)
            except BaseException as exc:  # noqa: BLE001
                log.error("agent %s failed (attempt %d/%d): %s",
                          self._critic.name, attempt, self._CRITIC_ATTEMPTS, exc)
                return None, None, exc
            raw = (crit.message.structured_payload or {}).get("verdict")
            verdict = str(raw).strip().lower() if raw is not None else None
            if verdict in self._VALID_VERDICTS:
                return crit, verdict, None
            log.info("critic emitted no usable verdict (raw=%r, attempt %d/%d)",
                     raw, attempt, self._CRITIC_ATTEMPTS)
        return crit, None, None  # emitted, but never a usable verdict

    @staticmethod
    def _build_numeric_audit_summary(result: PipelineResult) -> Optional[str]:
        """Collect each round-1 reader's ungrounded-decimal flags into a short
        block for the critic. Returns None when every reader was clean."""
        lines: list[str] = []
        for name in (
            "bio_result", "interaction_result",
            "conformation_result", "interface_exposure_result",
        ):
            r = getattr(result, name)
            audit = getattr(r, "numeric_audit", None) if r else None
            if audit is None or not getattr(audit, "ungrounded_decimals", None):
                continue
            flagged = ", ".join(
                f"{c.raw} ({c.where})" for c in audit.ungrounded_decimals[:8]
            )
            lines.append(f"- {audit.agent}: {flagged}")
        if not lines:
            return None
        return "\n".join(lines)

    @staticmethod
    def _build_preamble(ctx: MutationDesignContext) -> Optional[str]:
        """Render `overview` + `quality` once and concatenate them.

        Returns None if both views fail — callers treat None the same as
        an empty preamble and the agents just don't get that context.
        """
        try:
            locator = CaseLocator(ctx.case_dir)
        except Exception as exc:  # noqa: BLE001
            log.warning("preamble: CaseLocator(%s) failed: %s", ctx.case_dir, exc)
            return None

        chunks: list[str] = []
        for view in _PREAMBLE_VIEWS:
            try:
                markdown = render(view, locator, filters={})
            except PresenterError as exc:
                log.info("preamble: view %s unavailable: %s", view, exc)
                continue
            except Exception as exc:  # noqa: BLE001
                log.warning("preamble: view %s raised: %s", view, exc)
                continue
            chunks.append(f"### {view}\n\n{markdown.strip()}")

        if not chunks:
            return None
        return "\n\n".join(chunks)

    @staticmethod
    def _aggregate_candidates(
        ctx: MutationDesignContext,
        alias_map: Optional[ChainAliasMap] = None,
    ) -> str:
        """Build the deterministic candidate-site table for round 2.

        Walks `ctx.message_history`, reads each round-1 agent's
        `structured_payload` by sender name, and applies fixed rules to
        assemble a per-residue candidate record (merged roles + tags +
        the source msg_id(s) to cite). Bio anchor positions are excluded
        from candidates and surfaced under a separate "avoid" list.
        Missing senders / fields are skipped silently — round 2 degrades
        gracefully when a round-1 brief is absent.

        Residue labels are **canonicalized** (`residue_id.parse_residue_label`)
        before keying, so the same physical residue cited under different
        reader dialects (`α-ASP92` vs `ASP92 (D)`) merges into one row
        instead of splitting. The greek↔chain-letter map is loaded from the
        case (`alias_map`); when absent, greek labels still collapse among
        themselves.

        Returns a markdown blob: a candidate table, an avoid list, and a
        conformation-tag note. Returns a "no candidates" marker when no
        usable payload is present.
        """
        alias = alias_map or ChainAliasMap()

        # sender name -> structured_payload (last one wins if duplicated)
        payloads: Dict[str, dict] = {}
        msg_ids: Dict[str, str] = {}
        for msg in ctx.message_history:
            if isinstance(msg.structured_payload, dict):
                payloads[msg.sender] = msg.structured_payload
                msg_ids[msg.sender] = msg.msg_id

        # canonical key (chain, resid) | raw-fallback string -> record.
        # Each record keeps a display label, the raw labels seen (for audit),
        # and the merged roles / tags / sources.
        candidates: Dict[object, Dict[str, object]] = {}

        def _key_for(raw: str) -> tuple[object, str]:
            """Return (dedupe_key, display_label) for a raw reader label.

            Parseable labels key on the canonical (chain, resid) tuple and
            display via the canonical form. Unparseable labels fall back to
            keying on the raw string so they are never silently dropped."""
            res = parse_residue_label(raw, alias)
            if res is None:
                return (f"raw::{raw.strip()}", str(raw).strip())
            return (res.key, res.display())

        def _touch(raw: str, role: str, source: Optional[str]) -> Optional[dict]:
            if not raw:
                return None
            key, display = _key_for(raw)
            rec = candidates.get(key)
            if rec is None:
                rec = {
                    "display": display,
                    "raw_labels": {str(raw).strip()},
                    "roles": set(),
                    "tags": [],
                    "sources": set(),
                }
                candidates[key] = rec
            else:
                rec["raw_labels"].add(str(raw).strip())  # type: ignore[union-attr]
                # Prefer a display label that carries a resname + chain.
                if "(" in display and "(" not in str(rec["display"]):
                    rec["display"] = display
            if role:
                rec["roles"].add(role)  # type: ignore[union-attr]
            if source:
                rec["sources"].add(source)  # type: ignore[union-attr]
            return rec

        def _find(raw: str) -> Optional[dict]:
            """Look up an existing candidate by canonical key (for tagging)."""
            if not raw:
                return None
            key, _ = _key_for(raw)
            return candidates.get(key)

        # --- Interaction reader: hotspots -> hotspot_strengthen ---------
        ia = payloads.get("interaction_reader", {})
        ia_src = msg_ids.get("interaction_reader")
        for h in ia.get("hotspots", []) or []:
            res = h.get("residue")
            rec = _touch(res, "hotspot_strengthen", ia_src)
            if rec is None:
                continue
            role = h.get("role")
            if role:
                rec["tags"].append(f"interaction.role={role}")  # type: ignore[union-attr]

        # --- Interface & exposure: cavity / rim / burial ----------------
        ie = payloads.get("interface_exposure_reader", {})
        ie_src = msg_ids.get("interface_exposure_reader")
        for c in ie.get("cavity_candidates", []) or []:
            res = c.get("residue")
            rec = _touch(res, "cavity_fill", ie_src)
            if rec is None:
                continue
            sev = c.get("severity")
            hint = c.get("design_hint")
            if sev:
                rec["tags"].append(f"cavity.severity={sev}")  # type: ignore[union-attr]
            if hint:
                rec["tags"].append(f"cavity.hint={hint}")  # type: ignore[union-attr]
        for r in ie.get("rim_specificity_knobs", []) or []:
            _touch(r.get("residue"), "rim_tune", ie_src)
        for b in ie.get("notable_burial_events", []) or []:
            _touch(b.get("residue"), "caution", ie_src)

        # --- Bio: anchor primaries are excluded (avoid list) ------------
        bio = payloads.get("bio_agent", {})
        bio_src = msg_ids.get("bio_agent")
        anchors = bio.get("anchor_positions", {}) or {}
        avoid = [a for a in (anchors.get("primary") or []) if a]
        for a in avoid:
            key, _ = _key_for(a)
            candidates.pop(key, None)

        # --- Bio: germline conservation tags (canonical-key match) ------
        conservation = bio.get("conservation", {}) or {}
        for c in conservation.get("highly_conserved", []) or []:
            rec = _find(c.get("residue", ""))
            if rec is not None:
                rec["tags"].append("conservation=conserved")  # type: ignore[union-attr]
                if bio_src:
                    rec["sources"].add(bio_src)  # type: ignore[union-attr]
        for c in conservation.get("variable_tolerant", []) or []:
            rec = _find(c.get("residue", ""))
            if rec is not None:
                rec["tags"].append("conservation=variable")  # type: ignore[union-attr]
                if bio_src:
                    rec["sources"].add(bio_src)  # type: ignore[union-attr]

        # --- Conformation: tag existing candidates ----------------------
        cf = payloads.get("conformation_reader", {})
        cf_src = msg_ids.get("conformation_reader")
        for pr in cf.get("per_residue_conformation", []) or []:
            rec = _find(pr.get("residue", ""))
            if rec is None:
                continue
            cls = pr.get("class")
            sig = pr.get("design_signal")
            if cls:
                rec["tags"].append(f"conformation.class={cls}")  # type: ignore[union-attr]
            if sig:
                rec["tags"].append(f"conformation.signal={sig}")  # type: ignore[union-attr]
            if cf_src:
                rec["sources"].add(cf_src)  # type: ignore[union-attr]

        return MutationDesignPipeline._render_candidate_table(candidates, avoid)

    @staticmethod
    def _render_candidate_table(
        candidates: Dict[object, Dict[str, object]],
        avoid: list,
    ) -> str:
        lines: list[str] = []
        if not candidates:
            lines.append("(no candidate sites extracted from round-1 payloads)")
        else:
            lines.append("| residue | roles | tags | source msg_ids |")
            lines.append("|---|---|---|---|")
            # Sort by display label for stable, human-readable ordering.
            ordered = sorted(
                candidates.values(), key=lambda r: str(r["display"])
            )
            for rec in ordered:
                display = str(rec["display"])
                # De-dupe tags while preserving order (merged rows can repeat).
                seen: set = set()
                tag_list = [
                    t for t in rec["tags"]  # type: ignore[union-attr]
                    if not (t in seen or seen.add(t))
                ]
                roles = ", ".join(sorted(rec["roles"])) or "—"  # type: ignore[arg-type]
                tags = "; ".join(tag_list) or "—"
                sources = ", ".join(sorted(rec["sources"])) or "—"  # type: ignore[arg-type]
                lines.append(f"| `{display}` | {roles} | {tags} | {sources} |")
        if avoid:
            lines.append("")
            lines.append(
                "**Avoid (Bio anchor primaries — do not recommend without "
                "flagged justification):** " + ", ".join(f"`{a}`" for a in avoid)
            )
        return "\n".join(lines)


def make_pipeline(
    llm: LLMClient,
    *,
    bio_agent_cls: type[BaseAgent],
    interaction_reader_cls: type[BaseAgent],
    conformation_reader_cls: type[BaseAgent],
    interface_exposure_reader_cls: type[BaseAgent],
    recommendation_agent_cls: type[BaseAgent],
    design_critic_cls: Optional[type[BaseAgent]] = None,
) -> MutationDesignPipeline:
    """Convenience constructor: instantiate each agent against one shared LLM."""
    return MutationDesignPipeline(
        bio_agent=bio_agent_cls(llm),
        interaction_reader=interaction_reader_cls(llm),
        conformation_reader=conformation_reader_cls(llm),
        interface_exposure_reader=interface_exposure_reader_cls(llm),
        recommendation_agent=recommendation_agent_cls(llm),
        design_critic=design_critic_cls(llm) if design_critic_cls else None,
    )


__all__ = [
    "MutationDesignPipeline",
    "PipelineResult",
    "make_pipeline",
]
