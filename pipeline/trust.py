"""Composite per-trajectory trust tier (A / B / C).

The single headline "how much should I believe this trajectory's story" badge,
rolled up from the honesty signals already computed elsewhere:

  * dynamics-trust verdict   (essential_dynamics.dynamics_trust — equilibration,
                              PC1 cosine content, cross-replica RMSIP, CDR reliability)
  * BSA reconciliation       (struct_metrics.bsa_reconciliation — do the two
                              independent buried-surface estimates agree)
  * CDR-decomposition reliability
  * trajectory equilibration

Equilibration and CDR reliability contribute BOTH via the dynamics verdict and
again here — deliberate, so a severe QC problem pushes the composite to C rather
than being diluted inside a single sub-verdict. A missing input is never a
demerit (absence of evidence ≠ low trust).

NOTE on ``mode_soft``: recognition-mode ``soft`` axes live in
``meta.recognition_mode`` (written by recognition_mode.py, which runs AFTER
extract_analysis and READS analysis.json — a cycle). To keep this stage
dependency-free it defaults to None; callers in extract_analysis pass None so the
soft-axis term is inert on every build. (Wire it to the meta value only if you
accept a one-cycle staleness on a fresh build.)
"""
from __future__ import annotations


def trust_tier(dyn_verdict=None, bsa_recon=None, cdr_reliable=None,
               mode_soft=None) -> dict:
    """A/B/C composite trust tier + contributing reasons for one trajectory."""
    demerits: float = 0
    reasons: list[str] = []

    lvl = (dyn_verdict or {}).get("level")
    if lvl == "low":
        demerits += 2; reasons.append("low dynamics trust")
    elif lvl == "medium":
        demerits += 1; reasons.append("medium dynamics trust")

    if bsa_recon is not None and bsa_recon.get("consistent") is False:
        demerits += 1; reasons.append("BSA estimates disagree")

    if cdr_reliable is False:
        demerits += 1; reasons.append("CDR decomposition unreliable")

    n_soft = len(mode_soft or [])
    if n_soft:
        demerits += min(n_soft, 2) * 0.5
        reasons.append(f"{n_soft} recognition axis(es) near a tercile boundary")

    tier = "A" if demerits < 1 else "B" if demerits < 3 else "C"
    return {"tier": tier, "reasons": reasons}
