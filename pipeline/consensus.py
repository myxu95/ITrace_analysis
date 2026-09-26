"""Replica-consensus helpers — collapse a complex's three per-replica values into one
served value, instead of taking the primary replica and silently dropping
inter-replica disagreement.

This is the Tier-1 replica-consensus layer. It is imported by both collapse points —
aggregate.py (Explore distributions + reproducibility) and build_manifest.py
(Browse/Detail complex rows) — so the two surfaces report the same value, rather than
each picking the run3/primary replica independently.

**Quantitative descriptors use :func:`numeric_consensus`** (mean plus the measured
inter-replica spread). It replaced the boolean ``*_agree`` flags that the withdrawn
recognition-mode label axes used: a reported spread tells a user how far three
independent runs actually disagreed, where a flag only said whether they happened to
land in the same bin — so the flag read as "reproducible" for a complex whose replicas
differed by nearly the width of a bin, and as "not reproducible" for one that missed a
cut by a hair. See :mod:`pipeline.interface_descriptors`.

:func:`majority_label` remains for genuinely categorical, non-derived attributes. Do
not reintroduce a vote over binned continuous scalars — publish the scalars.

Nothing here re-runs analysis: it combines per-replica values already on disk."""
from __future__ import annotations

import math
from collections import Counter


def numeric_consensus(values):
    """Combine a complex's per-replica values for one quantitative descriptor.

    Returns ``{"mean", "sd", "min", "max", "range", "n"}`` over the non-null replica
    values, or ``None`` when none were observed. ``sd`` is the sample standard
    deviation (0.0 for a single replica); ``range`` is max - min, which is the more
    honest spread statistic at n=3 and is the one quoted in the descriptor docs.
    """
    vals = [float(v) for v in values if v is not None]
    if not vals:
        return None
    n = len(vals)
    mean = sum(vals) / n
    var = sum((v - mean) ** 2 for v in vals) / (n - 1) if n > 1 else 0.0
    return {
        "mean": round(mean, 4),
        "sd": round(math.sqrt(var), 4),
        "min": round(min(vals), 4),
        "max": round(max(vals), 4),
        "range": round(max(vals) - min(vals), 4),
        "n": n,
    }

# Equilibration severity order for a worst-case (most conservative) complex flag, so a
# complex is only reported 'ok' when *every* replica equilibrated.
EQ_RANK = {"ok": 0, "drifting": 1, "outlier": 2}


def majority_label(labels, primary=None, tie_break=None):
    """Majority vote over a complex's per-replica labels for one categorical axis.

    Returns ``(label, agree, n_observed)``:
      * ``label``      — the most common non-null label; on a tie prefer ``tie_break``,
        else ``primary`` if it is among the tied labels, else the lexicographically
        smallest tied label (deterministic, so rebuilds are stable).
      * ``agree``      — ``True`` when every observed replica shares one label,
        ``False`` when they disagree, ``None`` when fewer than two replicas produced a
        label (a single run can neither agree nor disagree).
      * ``n_observed`` — how many replicas produced a non-null label.
    """
    labs = [l for l in labels if l is not None]
    if not labs:
        return None, None, 0
    counts = Counter(labs)
    top = max(counts.values())
    tied = sorted(l for l, c in counts.items() if c == top)
    if len(tied) == 1:
        label = tied[0]                       # clear majority
    elif tie_break is not None:
        label = tie_break                     # tie -> hedged middle bin (straddle)
    elif primary in tied:
        label = primary                       # no hedge defined: keep the primary
    else:
        label = tied[0]                       # else deterministic (lexicographic)
    agree = (len(counts) == 1) if len(labs) >= 2 else None
    return label, agree, len(labs)


def majority_bool(values, primary=None, weights=None):
    """Majority vote over a complex's per-replica booleans (e.g. reversed_polarity).

    Returns ``(value, agree, n_observed)``. A strict majority wins; on an even split
    the decision falls back to the mean of ``weights`` (e.g. each replica's reversed
    fraction) > 0.5 when available, else the ``primary`` value. ``agree`` mirrors
    :func:`majority_label` (None for fewer than two observations)."""
    vals = [bool(v) for v in values if v is not None]
    if not vals:
        return None, None, 0
    n_true = sum(vals)
    if 2 * n_true != len(vals):
        value = n_true > len(vals) / 2
    elif weights:
        w = [x for x in weights if x is not None]
        value = (sum(w) / len(w) > 0.5) if w else bool(primary)
    else:
        value = bool(primary)
    agree = (len(set(vals)) == 1) if len(vals) >= 2 else None
    return value, agree, len(vals)


def qc_consensus(labels):
    """Worst-case equilibration consensus over a complex's replicas.

    Returns ``{"consensus", "reproducible", "counts", "n"}``. The consensus is the
    *worst* label across replicas (ok < drifting < outlier), so the served complex
    flag never hides a drifting/outlier run; ``reproducible`` is True only when every
    replica shares one label (None for fewer than two)."""
    labs = [l for l in labels if l]
    counts = {"ok": 0, "drifting": 0, "outlier": 0}
    for l in labs:
        if l in counts:
            counts[l] += 1
    if not labs:
        return {"consensus": None, "reproducible": None, "counts": counts, "n": 0}
    worst = max(labs, key=lambda l: EQ_RANK.get(l, 0))
    reproducible = (len(set(labs)) == 1) if len(labs) >= 2 else None
    return {"consensus": worst, "reproducible": reproducible, "counts": counts, "n": len(labs)}
