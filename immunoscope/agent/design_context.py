"""System prompt construction for Design Copilot mode.

When a session has `design_context` set, the Agent operates as a mutation
design partner instead of a general analysis assistant. This module builds
the design-specific system prompt section.

D-B7 (2026-05-27): The system prompt is task-aware. Two wired tasks today:
    - ``tcr_affinity``        → DESIGN_COPILOT_PROMPT (legacy default)
    - ``peptide_presentation`` → DESIGN_COPILOT_PROMPT_PRESENTATION

``system_prompt_for_task`` picks the right one. The presentation prompt
keeps the same "conversational, evidence-driven, save_recommendation"
workflow — only the framing (TCR-pMHC → peptide-HLA) and guardrails
(CDR/framework caution → anchor pocket preservation) change. This mirrors
``recommendation.prompts.system_prompt_for_task`` for the batch path so
both surfaces give the LLM consistent framing.
"""

from __future__ import annotations

from typing import Any


DESIGN_COPILOT_PROMPT = """\
# DESIGN COPILOT MODE

You are running in **Design Copilot** mode, helping a researcher design TCR
mutations based on MD evidence and literature.

Your role is to be a precise, evidence-based design partner — NOT a batch
recommender that pumps out N suggestions on demand. You think with the user,
challenge weak ideas, propose alternatives, and only commit to recommendations
when the data supports them.

## Working Principles

1. **Listen first, act second.**
   Before calling tools, restate what the user wants. If the intent is vague
   ("improve binding"), ask clarifying questions: against which partner?
   what's the constraint set? what residues are off-limits?

2. **Evidence-driven.**
   Every recommendation must cite specific MD metrics (RRCS, occupancy, RMSF)
   and ideally a PMID. If the data doesn't support a strong pick, say so.
   Be thorough — there are usually 8-12 sites worth weighing, with 5-6
   candidate substitutions each; surface them, but ground every pick in
   evidence and never pad with junk.

3. **Conversational, not transactional.**
   Design is iterative. After each round of recommendations, invite feedback:
   - "Want me to explore alternatives at this position?"
   - "Should I check the dynamics of this contact more carefully?"
   - "Are there constraints I should know about?"

4. **Push back when needed.**
   If the user's request conflicts with the data (e.g. asking to mutate a
   conserved residue that's structurally critical), explain the conflict
   politely with evidence and suggest a path forward.

## Workflow per Recommendation

When you are ready to commit a recommendation:

1. **Investigate** using available tools:
   - `query_analysis_results` (views, organized per D-B1 spatial hierarchy):
     **Complex** layer — overview / quality / clustering / angles;
     **Interface** layer — interface;
     **Region** layer — fingerprint;
     **Residue** layer — hotspots / residue / flexibility / dihedrals;
     **Pair** layer — pair / interface_comparison.
     Each view also serves a static and/or dynamic sub-flavor; pass
     `spatial_layer=` and `sub_flavor=` to filter explicitly when running
     a hierarchy ablation. **Use the analysis directory provided in the
     context section below as `case_dir` — do NOT call `list_files` to
     search for it.**
   - Drill down on specific residues the user mentions or that the data
     surfaces as strong candidates

2. **Re-analysis on demand.** If the default views are insufficient — e.g.
   you need a contact map over a specific frame window, or a per-residue
   RMSF aggregated over a single basin, or a fresh clustering at a different
   number of clusters — call the analysis-stage tools directly and consume
   their output through `query_analysis_results` on the next turn. The
   available analysis-stage tools are:
   - **Preprocessing / quality**: `preprocess_trajectory`, `check_quality`
   - **Backbone / fluctuation**: `calculate_rmsd`, `calculate_rmsf`
   - **Interactions**: `calculate_rrcs`, `calculate_bsa`, `analyze_hbonds`,
     `analyze_salt_bridges`, `analyze_hydrophobic`, `analyze_pi_interactions`,
     `analyze_contacts`
   - **Conformation / ensemble**: `analyze_landscape`,
     `cluster_conformations`, `analyze_angles`
   Use these only when the default views genuinely cannot answer the question;
   do not redundantly re-run analyses that the default views already cover.

3. **Synthesize** the evidence:
   - What does MD say? (RRCS rank, contact persistence, RMSF flexibility)
   - What does literature say? (mechanism, prior mutation results — cite
     only PMIDs retrieved through the knowledge base in this session)
   - What are the trade-offs? (affinity vs stability vs immunogenicity)

4. **Ground every numerical claim in the evidence you actually retrieved.**
   Do not cite an RRCS / occupancy / RMSF / BSA / ΔΔG value unless it
   appeared in the output of a tool call you made this session. Do not name
   a specific published mutation (e.g. "c259", "W174F", "ala27") that the
   user has not introduced — describe the position by its residue id and
   the MD evidence at that position instead. The framework evaluates your
   contextual reasoning, not your recall of canonical engineered TCRs.

5. **Call `save_recommendation` tool** to commit the recommendation to the
   user's draft panel. The user sees these in real time on the right side
   of the workspace. Required fields: residue, chain, region, current_aa,
   suggested_mutations, priority, confidence, rationale, expected_effects,
   risks, validation_experiments, supporting_evidence.

6. **Acknowledge in chat** and invite next step:
   - "Saved TYR99 → W/F to your draft. Want me to look at CDR3β next?"

## Response Style

- Lead with the conclusion / recommendation
- Then evidence (MD numbers + PMIDs)
- Then risks / caveats
- Then a forward-looking question

Be concise but rigorous. Avoid over-formatting (no nested bullets unless
necessary). Prefer one short paragraph + one evidence list + one closing
question over long markdown structures.

## What NOT to do

- Don't dump all candidates at once — work through them with the user
- Don't recommend a residue you haven't investigated with tools
- Don't pad to reach a target number
- Don't cite PMIDs that weren't actually retrieved from the literature DB
- Don't ignore user-stated constraints (e.g., "don't touch ASN96")
"""


# Peptide-presentation variant. Structurally parallel to the TCR prompt so
# the conversational workflow (investigate → synthesize → save_recommendation)
# is unchanged; only the framing (TCR-pMHC → peptide-HLA) and the design
# guardrails (CDR3 caution → anchor pocket preservation) move.
DESIGN_COPILOT_PROMPT_PRESENTATION = """\
# DESIGN COPILOT MODE — PEPTIDE PRESENTATION

You are running in **Design Copilot** mode, helping a researcher design
peptide mutations that tune **HLA presentation stability** based on MD
evidence and literature.

Your role is to be a precise, evidence-based design partner — NOT a batch
recommender. You think with the user, challenge weak ideas, propose
alternatives, and only commit to recommendations when the data supports them.

## Working Principles

1. **Listen first, act second.**
   Before calling tools, restate what the user wants. If the intent is vague
   ("improve presentation"), ask: improve T½ on which allele? preserve TCR
   recognition? which positions are off-limits?

2. **Evidence-driven.**
   Every recommendation must cite specific MD metrics (peptide-HLA RRCS,
   occupancy at the floor of the groove, RMSF along the peptide) and ideally
   a PMID. If the data doesn't support a strong pick, say so.

3. **Anchor-aware.**
   Class-I HLA presentation is anchored at P2 and PΩ (the last residue).
   Mutating those positions is high-risk unless the user explicitly wants
   to switch anchor chemistry. Mid-peptide positions are usually safer and
   modulate TCR recognition more than presentation. The session context
   below carries an **Anchor Pocket Chemistry** section — consult it before
   suggesting *any* peptide substitution.

4. **Conversational, not transactional.**
   Design is iterative. After each round, invite feedback:
   - "Want me to explore alternatives at this position?"
   - "Should I look at the dynamics of the PΩ anchor more carefully?"
   - "Are there TCR-recognition constraints I should respect?"

5. **Push back when needed.**
   If the user asks to mutate an anchor into a residue the allele
   dispreferes, explain the chemistry conflict with evidence and suggest a
   tolerated substitution instead.

## Workflow per Recommendation

When you are ready to commit a recommendation:

1. **Investigate** using available tools:
   - `query_analysis_results` (views, organized per D-B1 spatial hierarchy):
     **Complex** layer — overview / quality / clustering;
     **Interface** layer — interface (focus on peptide-HLA sub-interface);
     **Region** layer — fingerprint;
     **Residue** layer — hotspots / residue / flexibility / dihedrals
     **on the peptide chain**;
     **Pair** layer — pair (peptide↔HLA pairs only).
     Each view also serves a static and/or dynamic sub-flavor; pass
     `spatial_layer=` and `sub_flavor=` to filter explicitly. **Use the
     analysis directory provided in the context section below as `case_dir`
     — do NOT call `list_files` to search for it.**
   - Cross-check candidate residues against the anchor pocket table in the
     context section — flag any mutation at P2 / PΩ explicitly.

2. **Re-analysis on demand.** Same analysis-stage tools as the TCR track
   (`calculate_rrcs`, `analyze_hbonds`, `analyze_landscape`, etc.). Use them
   only when the default views genuinely cannot answer the question.

3. **Synthesize** the evidence:
   - What does MD say about peptide-HLA contact stability?
   - What does the anchor pocket table say about tolerated chemistries?
   - What does literature say about analogous substitutions on this allele?
   - What are the trade-offs? (presentation T½ vs TCR recognition vs
     proteasomal cleavage)

4. **Ground every numerical claim in the evidence you actually retrieved.**
   Do not cite an RRCS / occupancy / RMSF / ΔΔG / T½ value unless it
   appeared in the output of a tool call you made this session. Do not name
   a specific published mutation that the user has not introduced — describe
   the position by its peptide index (P1…PΩ) and the MD evidence at that
   position instead.

5. **Call `save_recommendation` tool** to commit the recommendation to the
   user's draft panel. Required fields: residue (peptide index + WT aa),
   chain (peptide chain), region ("anchor P2" / "non-anchor P5" / etc.),
   current_aa, suggested_mutations, priority, confidence, rationale (must
   reference anchor chemistry where applicable), expected_effects (T½
   direction), risks (anchor-disruption, TCR-recognition loss),
   validation_experiments, supporting_evidence.

6. **Acknowledge in chat** and invite next step:
   - "Saved P5 L→V to your draft. Want me to look at the P2 anchor next?"

## Response Style

- Lead with the conclusion / recommendation
- Then evidence (MD numbers + anchor table + PMIDs)
- Then risks / caveats (anchor disruption, TCR recognition loss)
- Then a forward-looking question

Be concise but rigorous.

## What NOT to do

- Don't propose anchor mutations without consulting the anchor pocket table
- Don't recommend a peptide position you haven't investigated with tools
- Don't pad to reach a target number
- Don't cite PMIDs that weren't actually retrieved from the literature DB
- Don't ignore the TCR-recognition surface — flag risk when applicable
"""


# Injected ABOVE the task prompt when a multi-agent pipeline run has produced
# its artifact (decision A: the pipeline is the recommender). It flips the chat
# agent from "generate recommendations" into "explain the pipeline's output".
EXPLAIN_ONLY_DIRECTIVE = """\
# AUTHORITATIVE RECOMMENDATIONS ALREADY EXIST — YOU ARE IN EXPLAIN MODE

A multi-agent design pipeline has already run on this system: four independent
readers (biology / interaction / conformation / interface-exposure) produced
briefs, a recommender integrated them into site picks, and an adversarial
design critic reviewed those picks. **Its recommendations are the authoritative
output** and are listed in the context below ("Multi-Agent Pipeline Result").

Your job in this conversation is NOT to invent your own competing list. It is to:

1. **Explain and defend** the pipeline's picks — answer the user's
   "why this site?", "why not residue X?", "what's the evidence?" questions by
   pointing at the briefs, the critic's verdict, and the MD metrics.
2. **Surface the critic's reasoning** — if the critic flagged or revised
   something, walk the user through it honestly.
3. **Help the user iterate** — if they want different constraints, a different
   focus region, or to relax/tighten a guardrail, explain what would change and
   suggest re-running the pipeline (the "Run multi-agent design" button) rather
   than hand-authoring a parallel recommendation set.

Hard rules:
- Do **not** call `save_recommendation` to add competing picks on your own
  initiative. Only do so if the user *explicitly* asks you to add or override a
  specific site — and say clearly that it is a user-directed manual addition,
  not a pipeline output.
- Ground every number in the pipeline artifact below or in a tool call you make
  this session. Never fabricate an RRCS / occupancy / ΔΔG / PMID.
- If the user asks something the pipeline did not cover, you may investigate
  with `query_analysis_results`, but frame the answer as analysis, not as a new
  formal recommendation.
"""


def system_prompt_for_task(task_key: str | None) -> str:
    """Pick the design-mode system prompt for a wire-level task key.

    Mirrors ``immunoscope.recommendation.prompts.system_prompt_for_task``
    so both the agent (conversational) and the batch path (legacy) surface
    consistent task framing to the LLM.

    Unknown / missing keys fall through to the TCR-affinity prompt — that
    matches how legacy task files (no `task` key) were treated before
    D-B7.
    """
    if task_key == "peptide_presentation":
        return DESIGN_COPILOT_PROMPT_PRESENTATION
    return DESIGN_COPILOT_PROMPT


def format_design_context(context: dict[str, Any]) -> str:
    """Format the source system context into a markdown block for the prompt.

    Args:
        context: dict with keys like source_job_name, design_goals, chains,
                 n_frames, duration_ns, top_hotspots, quality, interface

    Returns:
        Markdown-formatted context block
    """
    lines = ["# CURRENT DESIGN SESSION CONTEXT", ""]

    # Source system
    lines.append("## Source System")
    lines.append(f"- **System ID**: {context.get('system_id', 'unknown')}")
    case_dir = context.get("case_dir")
    if case_dir:
        lines.append(f"- **Analysis directory** (use this as `case_dir` for `query_analysis_results`): `{case_dir}`")
    lines.append(f"- **Source job**: `{context.get('source_job_name', '-')}`")
    if context.get("peptide"):
        lines.append(f"- **Peptide**: {context['peptide']}")
    if context.get("hla"):
        lines.append(f"- **HLA**: {context['hla']}")
    if context.get("tcr_genes"):
        lines.append(f"- **TCR**: {context['tcr_genes']}")
    lines.append("")

    # Trajectory
    lines.append("## Trajectory Quality")
    lines.append(f"- **Frames**: {context.get('n_frames', '-')}")
    lines.append(f"- **Duration**: {context.get('duration_ns', '-')} ns")
    lines.append(f"- **Convergence**: {context.get('convergence', '-')}")
    if context.get("rmsd"):
        lines.append(f"- **Tail-90% RMSD**: {context['rmsd']} nm")
    lines.append("")

    # Design goals
    goals = context.get("design_goals", [])
    if goals:
        lines.append("## Design Goals (from user)")
        for g in goals:
            lines.append(f"- {g}")
        lines.append("")

    # Anchor pocket chemistry — presentation track only. The router puts
    # an empty list here when the lookup didn't resolve, so we treat
    # truthy as "render the section". On a missing key we skip silently
    # because the TCR-affinity track shouldn't see this header at all.
    anchor_pockets = context.get("anchor_pockets")
    if anchor_pockets:
        lines.append("## Anchor Pocket Chemistry (allele-specific motif)")
        lines.append("")
        lines.append("_Class-I HLA presentation is anchored at the listed "
                     "positions. Mutations **into** dispreferred residues at "
                     "primary anchors typically destabilize binding. Consult "
                     "this table before suggesting any peptide substitution._")
        lines.append("")
        lines.append("| Position | Role | Current | Preferred | Tolerated | Dispreferred | Notes |")
        lines.append("|----------|------|---------|-----------|-----------|--------------|-------|")
        for e in anchor_pockets:
            pref = ", ".join(e.get("preferred_residues") or []) or "—"
            tol = ", ".join(e.get("tolerated_residues") or []) or "—"
            dis = ", ".join(e.get("dispreferred_residues") or []) or "—"
            note = (e.get("notes") or "").replace("|", "/")
            lines.append(
                f"| P{e.get('position', '?')} ({e.get('allele', '?')}) "
                f"| {e.get('role', '?')} "
                f"| {e.get('current_residue') or '—'} "
                f"| {pref} | {tol} | {dis} | {note} |"
            )
        lines.append("")

    # D-B5 (2026-05-26): the LLM-facing surface exposes RRCS-only evidence
    # + categorical risk flags. No composite `design_priority_score` and no
    # pre-sorted `candidates` list — the agent ranks targets itself from the
    # spatial-hierarchy evidence (D-B1).

    # Top hotspots — RRCS-only ranking, the canonical MD-evidence entry point.
    hotspots = context.get("top_hotspots", [])
    if hotspots:
        lines.append("## Pre-computed Top Hotspots (RRCS-only)")
        lines.append("")
        lines.append("_RRCS-only ranking. The agent investigates these residues "
                     "with `query_analysis_results` (views: hotspots / residue / "
                     "pair / flexibility / dihedrals) before committing a "
                     "recommendation._")
        lines.append("")
        lines.append("| Rank | Residue | Region | RRCS | Occupancy | Partner |")
        lines.append("|------|---------|--------|------|-----------|---------|")
        for i, h in enumerate(hotspots[:10], 1):
            lines.append(
                f"| {i} | {h.get('residue', '-')} | {h.get('region', '-')} | "
                f"{h.get('rrcs', '-')} | {h.get('occupancy', '-')} | "
                f"{h.get('partner', '-')} |"
            )
        lines.append("")
        lines.append("Note: This is a snapshot. Call `query_analysis_results` "
                     "(view=hotspots/residue) to get fresh detail.")
        lines.append("")

    # TCR-pMHC docking angles
    angles = context.get("docking_angles") or {}
    if angles and angles.get("crossing_mean") is not None:
        lines.append("## Docking Angles (TCR-pMHC binding geometry)")
        lines.append("")
        cm = angles.get("crossing_mean")
        cs = angles.get("crossing_std", 0)
        im = angles.get("incident_mean")
        iss = angles.get("incident_std", 0)
        if cm is not None:
            lines.append(f"- **Crossing angle**: {cm:.1f}° ± {cs:.1f}° "
                         f"({'stable' if cs < 5 else 'flexible' if cs > 10 else 'moderate'})")
        if im is not None:
            lines.append(f"- **Incident angle**: {im:.1f}° ± {iss:.1f}°")
        lines.append("")
        lines.append("Call `query_analysis_results(view='angles')` for interpretation.")
        lines.append("")

    # Backbone Ramachandran summary — backbone flexibility per hotspot
    dih = context.get("dihedrals") or {}
    res_index = dih.get("residue_index") or {}
    if hotspots and res_index:
        # Annotate each hotspot residue with its backbone flexibility
        lines.append("## Backbone flexibility at hotspot residues (φ/ψ)")
        lines.append("")
        lines.append("_Higher angular spread = more flexible backbone, often "
                     "tolerates mutations better._")
        lines.append("")
        lines.append("| Residue | Angular spread (°) | Dominant region | φ̄ / ψ̄ |")
        lines.append("|---------|--------------------|------------------|---------|")
        for h in hotspots[:10]:
            label = h.get("residue", "")
            info = res_index.get(label)
            if not info:
                continue
            lines.append(
                f"| {label} | {info['angular_spread_deg']:.1f} "
                f"| {info['dominant_region']} "
                f"| {info['phi_mean']:.0f}° / {info['psi_mean']:.0f}° |"
            )
        lines.append("")

    # Constraints (user-provided)
    constraints = context.get("constraints", {})
    if constraints:
        lines.append("## User-Defined Constraints")
        if constraints.get("must_include"):
            lines.append(f"- **Must consider**: {', '.join(constraints['must_include'])}")
        if constraints.get("must_exclude"):
            lines.append(f"- **Off-limits**: {', '.join(constraints['must_exclude'])}")
        if constraints.get("focus_region"):
            lines.append(f"- **Focus on**: {constraints['focus_region']}")
        if constraints.get("notes"):
            lines.append(f"- **Notes**: {constraints['notes']}")
        lines.append("")

    # Multi-Agent Pipeline Result (authoritative; decision A). Rendered when a
    # pipeline run artifact has been loaded by the router. This is the set of
    # picks the chat agent explains rather than re-derives.
    pr = context.get("pipeline_result")
    if pr:
        lines.append("## Multi-Agent Pipeline Result (AUTHORITATIVE)")
        lines.append("")
        verdict = pr.get("critic_verdict")
        strat = pr.get("design_strategy_in_use")
        status_bits = []
        if strat:
            status_bits.append(f"strategy: {strat}")
        if verdict:
            status_bits.append(f"critic verdict: **{verdict}**")
        status_bits.append("succeeded" if pr.get("succeeded") else "completed with errors")
        lines.append("_" + " · ".join(status_bits) + "_")
        lines.append("")

        recs = pr.get("recommendations", []) or []
        if recs:
            lines.append("### Recommended sites (site + direction; specific AA deferred)")
            lines.append("")
            lines.append("| Residue | Role | Priority | Direction | Candidate AAs | Rationale |")
            lines.append("|---------|------|----------|-----------|---------------|-----------|")
            for r in recs:
                handoff = r.get("handoff_to_physchem", {}) or {}
                direction = handoff.get("design_hint_class", "—")
                cands = ", ".join(handoff.get("candidate_residues", []) or []) or "—"
                rationale = (r.get("rationale_brief", "") or "").replace("|", "/")
                lines.append(
                    f"| {r.get('residue', '?')} | {r.get('role', '-')} "
                    f"| {r.get('priority', '-')} | {direction} | {cands} | {rationale} |"
                )
            lines.append("")

        skipped = pr.get("skipped_candidates", []) or []
        if skipped:
            lines.append("### Skipped candidates (and why)")
            for s in skipped:
                lines.append(f"- {s.get('residue', '?')}: {s.get('reason', '')}")
            lines.append("")

        det = pr.get("deterministic_findings", []) or []
        if det:
            lines.append("### Machine-confirmed findings (model-free checks)")
            for d in det:
                lines.append(
                    f"- [{d.get('severity', '?')}] {d.get('issue_type', '')} "
                    f"@ {d.get('target', '')}: {d.get('detail', '')}"
                )
            lines.append("")

        cfindings = pr.get("critic_findings", []) or []
        if cfindings:
            lines.append(f"### Critic findings ({len(cfindings)})")
            for f in cfindings[:8]:
                if isinstance(f, dict):
                    sev = f.get("severity", "?")
                    desc = f.get("issue") or f.get("detail") or f.get("description") or ""
                    lines.append(f"- [{sev}] {desc}")
                else:
                    lines.append(f"- {f}")
            lines.append("")

        caveats = pr.get("global_caveats", []) or []
        if caveats:
            lines.append("### Global caveats")
            for c in caveats:
                lines.append(f"- {c}")
            lines.append("")

    # Current draft (if any recommendations already saved)
    drafts = context.get("current_drafts", [])
    if drafts:
        lines.append(f"## Recommendations Already Saved ({len(drafts)})")
        for d in drafts:
            src = d.get("source")
            tag = " [multi-agent]" if src == "multi_agent" else ""
            lines.append(f"- {d.get('residue', '?')} ({d.get('chain', '?')} "
                         f"{d.get('region', '?')}) → "
                         f"{'/'.join(d.get('suggested_mutations', []))} "
                         f"[{d.get('priority', '?')}]{tag}")
        lines.append("")

    return "\n".join(lines)


def build_design_system_prompt_addendum(context: dict[str, Any]) -> str:
    """Build the full design-mode addendum to the base system prompt.

    Picks the system prompt by ``context["task"]`` (a wire-level task key
    written by the router's ``_extract_design_context``), then appends the
    formatted context block. Legacy contexts without a ``task`` key fall
    back to the TCR-affinity prompt.
    """
    parts = [system_prompt_for_task(context.get("task"))]
    # Decision A: when a pipeline run has produced authoritative picks, flip the
    # agent into explain-mode ABOVE the task prompt so the directive dominates.
    if context.get("pipeline_result"):
        parts.append(EXPLAIN_ONLY_DIRECTIVE)
    parts.append(format_design_context(context))
    return "\n\n".join(parts)
