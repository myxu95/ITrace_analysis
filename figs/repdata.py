"""Shared loader for replicas.tsv, built by build_replicas.py.

Two readers, because the two figures that use this table count differently.

``library()`` is for Figure 3: the counting unit is the COMPLEX and the three
replicas are already pooled, so every row is one entry in the atlas and the
spread the figure draws is spread across complexes.

``replicas()`` is for Figure 4, where the replicas are the subject rather than
something to average away.  It returns the three per-replica values of one
quantity and drops the complexes whose replicas cannot answer the question:

  4prh, 5hhm and 7rk7 for the contact quantities only -- run3's chain mapping
  double-counts a chain, inflating its contact counts to 225-239 stable pairs
  where the rest of the library tops out at 73.  Their angle, Q, RMSF and RMSIP
  come from other code paths and are in line with their siblings, so those
  quantities keep all three complexes.  RMSF in particular is averaged over
  every residue rather than over a chain selection, so the doubled chain does
  not reach it: the three complexes' run-3 means sit within 0.2-0.6 A of their
  siblings while their stable-contact counts are inflated three- to fourfold.
"""
import csv

from figstyle import HERE

TSV = HERE / "replicas.tsv"

N_TOTAL = 245
QUANTITIES = {          # column stem -> (label, unit, decimals)
    "inc": ("Incident angle", "°", 1),
    "cross": ("Crossing angle", "°", 1),
    "q": ("Native-contact fraction", "", 3),
    "stable": ("Stable contact pairs", "", 0),
    "rmsip": ("Essential-subspace RMSIP", "", 3),
    "rmsf": ("Mean Cα RMSF", "Å", 2),
}


def _rows():
    return list(csv.DictReader(TSV.open(), delimiter="\t"))


def library():
    """-> [{pdb_id, inc_p5/p50/p95, cross_p5/p50/p95, flag}] , one per complex.

    Percentiles are over the pooled 3 x 1001 frames of the complex, so p5-p95
    is the range the angle normally occupies and p95-p5 is how far it travels.
    """
    out = []
    for r in _rows():
        d = {"pdb_id": r["pdb_id"], "flag": r["flag"]}
        for k in ("inc", "cross"):
            for p in (5, 50, 95):
                d[f"{k}_p{p}"] = float(r[f"{k}_p{p}"])
        out.append(d)
    assert len(out) == N_TOTAL, len(out)
    return out


# Which flags disqualify which quantity.  The exclusions differ by quantity
# because the defects damage different code paths: chain_selection corrupts only
# what is built on the interface chain selection.  dup_analysis is kept in the
# table because the builder still derives it, but no complex carries it since
# 6g9q's three analysis records were rewritten on 2026-09-14; if it ever fires
# again it costs a complex one independent replica and so ruins every agreement
# statistic, which is why it disqualifies every quantity.
# If a future figure draws buried surface area or any interaction count, those
# come off the same selection and take {"dup_analysis", "chain_selection"} too.
DROP = {"inc": {"dup_analysis"}, "cross": {"dup_analysis"},
        "q": {"dup_analysis"}, "rmsip": {"dup_analysis"},
        "rmsf": {"dup_analysis"},
        "stable": {"dup_analysis", "chain_selection"}}


def replicas(stem, drop_flagged=True):
    """-> [{pdb_id, values: [v1, v2, v3], mean, spread}] for one quantity.

    ``spread`` is max - min of the three replica values: the whole of the
    between-replica disagreement, not a summary of it.  ``drop_flagged``
    removes the complexes whose replicas cannot answer the question for THIS
    quantity -- see the module docstring for which, and why it differs by
    quantity rather than being one blanket exclusion list.
    """
    assert stem in QUANTITIES, stem
    out = []
    for r in _rows():
        if drop_flagged and r["flag"] in DROP[stem]:
            continue
        v = [float(r[f"{stem}_run{i}"]) for i in (1, 2, 3) if r.get(f"{stem}_run{i}")]
        if len(v) < 3:
            continue
        out.append({"pdb_id": r["pdb_id"], "values": v,
                    "mean": sum(v) / 3, "spread": max(v) - min(v)})
    return out


def angle_variation(drop_flagged=True):
    """-> [{pdb_id, band, spread, within_sd, ratio, flag}] per complex.

    Two measurements of how far the incident angle moves, both in degrees, so
    they can be put on one ruler.

      ``band``      p95 - p5 of the complex's pooled 3 x 1001 frames: the whole
                    angular interval its trajectories occupy.  This is the
                    quantity Figure 3d draws as a band, reused here as the
                    denominator rather than redrawn.
      ``spread``    max - min of the three replica MEANS: how far independent
                    runs disagree about where that interval sits.
      ``within_sd`` mean of the three replicas' own SD: the part of the motion
                    that happens inside one trajectory, with no between-run
                    contribution in it at all.

    ``ratio`` is spread / band.  The band is a total, so it contains the
    disagreement the ratio is measuring -- that is what makes it a decomposition
    and not an independent comparison, and it is why ``within_sd`` is returned
    alongside: the second yardstick has the between-run part removed, and the
    manuscript quotes both rather than the flattering one.
    """
    out = []
    for r in _rows():
        if drop_flagged and r["flag"] in DROP["inc"]:
            continue
        v = [float(r[f"inc_run{i}"]) for i in (1, 2, 3)]
        band = float(r["inc_p95"]) - float(r["inc_p5"])
        sd = sum(float(r[f"incsd_run{i}"]) for i in (1, 2, 3)) / 3
        out.append({"pdb_id": r["pdb_id"], "band": band,
                    "spread": max(v) - min(v), "within_sd": sd,
                    "ratio": (max(v) - min(v)) / band, "flag": r["flag"]})
    return out


def within_sd():
    """-> {pdb_id: mean of the three replicas' own incident-angle SD}.

    The yardstick Figure 4 measures the between-replica spread against: how far
    the angle moves inside a single 200 ns trajectory.
    """
    return {r["pdb_id"]: sum(float(r[f"incsd_run{i}"]) for i in (1, 2, 3)) / 3
            for r in _rows()}


def quant(values, p):
    """Linear-interpolated quantile; one definition across every figure."""
    v = sorted(values)
    k = (len(v) - 1) * p
    f = int(k)
    c = min(f + 1, len(v) - 1)
    return v[f] + (v[c] - v[f]) * (k - f)


def icc(groups):
    """One-way random-effects ICC(1,1), for [[v1, v2, v3], ...].

    Shrout and Fleiss (1979) model 1: each complex is a random subject and its
    replicas are interchangeable measurements of it, so the coefficient is
    (MSB - MSW) / (MSB + (k - 1) MSW).  One number for "does a single replica
    reproduce its own complex's value": 1 means the replicas of a complex are
    identical and all the variance in the library is between complexes; 0 means
    a replica says nothing about which complex it came from.

    This used to be between / (between + within) with an unscaled between term,
    which is a variance ratio and not an ICC -- it ran about 0.06 high for Q and
    0.005-0.01 high for the angles.  The name is a defined statistic, so the
    definition is the published one.
    """
    k = len(groups[0])
    n = len(groups)
    assert all(len(g) == k for g in groups), "ICC needs balanced groups"
    means = [sum(g) / k for g in groups]
    gm = sum(means) / n
    msb = k * sum((m - gm) ** 2 for m in means) / (n - 1)
    msw = sum(sum((v - m) ** 2 for v in g)
              for g, m in zip(groups, means)) / (n * (k - 1))
    return (msb - msw) / (msb + (k - 1) * msw)


def icc_k(groups):
    """ICC(1,k): the reliability of the MEAN of the k replicas.

    icc() asks how far one replica can be trusted on its own.  This asks how far
    the average of the three can, which is what someone who downloads a complete
    entry actually holds.  Spearman-Brown on icc(): k r / (1 + (k - 1) r).
    """
    k = len(groups[0])
    r = icc(groups)
    return k * r / (1 + (k - 1) * r)
