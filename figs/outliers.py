#!/usr/bin/env python
"""Which complexes stand out in the figures, on which axis, and whether it is real.

The eye picks detached points out of a ranked strip immediately.  This puts a
number on that impression, so the same complexes get named the same way in
every panel and in the text, and -- more importantly -- so that "looks odd" is
separated from "is wrong".  Nothing here is a published descriptor: these are
triage flags for deciding what to re-check and what to annotate.  Anything that
reaches print reaches it as the measured scalar, never as a rank of our own
library.

Criterion is the Tukey fence -- outside Q1 - k*IQR .. Q3 + k*IQR, k = 1.5 for
"out" and 3.0 for "far out".  It is computed on the library's own distribution,
which is what makes it agree with what a reader sees, and it is shape-free, so
the long right tail every spread quantity has does not manufacture flags the
way a mean +- n*SD rule would.

Two axes are deliberately scanned inside a stratum rather than pooled:

  bulge      binned by peptide length.  Pooled, the scan returns the five
             13-mers and seven 11-mers and calls them outliers, but bulge rises
             with length (Pearson r = +0.77 over 245 complexes) and the figure
             that draws bulge already bins by length, so pooling would flag the
             length distribution's tail instead of anything about the bulge.

  crossing   reverse-polarity complexes are held out.  A crossing angle past
             90 deg is a documented minority binding mode, not a defect, and
             leaving them in widens the fence enough to hide everything else.

What a flag MEANS is read off the COUPLING of the quantities, not off its
size.  For each complex the replica with the widest incident angle is tested
against its three partners: is it also the lowest-Fnat, lowest-BSA and
highest-RMSF replica?

  3/3   the angle opened because the interface let go -- a partial unbinding
        event, sampled and recorded, which is the point of running replicas.
  1-2/3 the angle moved without the interface degrading in step.  3utt run2 is
        the clean example: lowest Fnat of its three, but the HIGHEST buried
        surface area, i.e. the contacts rearranged rather than dissolved.
        1fo0 inverts it outright -- its widest-angle replica is its best-bonded
        one, so the angle is rocking on an intact footprint.

Neither reading is a defect.  A defect looks different and is tested
separately: ``build_replicas.defect_state`` re-derives the two known analysis
faults from the trajectories themselves (duplicated analysis blocks; a chain
mapping that drives contacts up while buried surface area ALSO rises, which is
physically impossible).  Re-run over all 245 complexes and 735 trajectories it
returns exactly the four complexes already diagnosed and nothing else, so no
complex on this page is flagged here because its data is broken.

    python outliers.py            report to stdout, outliers.tsv beside it
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import dyndata
import repdata
from repdata import quant

DATA = Path("/home/xmy/work/data/immunotrace/web_data_1000")
K_OUT, K_FAR = 1.5, 3.0
REVERSED = ("5sws", "5swz", "7jwi", "9gv7")

# already diagnosed, repair in flight -- listed so the tally is complete.
# 6g9q left this list on 2026-09-17: its three analysis records were rewritten
# on 2026-09-14 and now carry distinct Fnat, incident angle and RMSF.
KNOWN_DEFECT = {
    "4prh": "run3 chain mapping labels both TCR chains alpha",
    "5hhm": "run3 chain mapping labels both TCR chains alpha",
    "7rk7": "run3 chain mapping labels both TCR chains alpha",
}


def fences(vals, k):
    q1, q3 = quant(vals, 0.25), quant(vals, 0.75)
    return q1 - k * (q3 - q1), q3 + k * (q3 - q1)


def scan(name, series, side, fmt, strata=None):
    """series: {pdb_id: value}.  strata: {pdb_id: key} to fence within."""
    groups = defaultdict(dict)
    for pid, v in series.items():
        groups[strata[pid] if strata else None][pid] = v
    hits = []
    for key, sub in groups.items():
        if len(sub) < 8:                      # too few to fence honestly
            continue
        lo15, hi15 = fences(list(sub.values()), K_OUT)
        lo30, hi30 = fences(list(sub.values()), K_FAR)
        for pid, v in sub.items():
            deg = None
            if side in ("hi", "both"):
                deg = "far" if v > hi30 else "out" if v > hi15 else None
            if deg is None and side in ("lo", "both"):
                deg = "far" if v < lo30 else "out" if v < lo15 else None
            if deg:
                hits.append((pid, v, deg, key))
    hits.sort(key=lambda t: t[1] if side == "lo" else -t[1])
    vals = list(series.values())
    lo15, hi15 = fences(vals, K_OUT)
    return {"name": name, "n": len(vals), "side": side, "fmt": fmt,
            "strat": strata is not None, "med": quant(vals, 0.5),
            "lo15": lo15, "hi15": hi15, "hits": hits}


def traj_table():
    """-> {pdb_id: [{run, inc, fnat, fnat_min, bsa, rmsf, stable}, ...]}"""
    out = {}
    for d in sorted(DATA.glob("*_run*")):
        pid, run = d.name.rsplit("_run", 1)
        p = d / "analysis" / "analysis.json"
        if not p.exists():
            continue
        a = json.loads(p.read_text())
        out.setdefault(pid, []).append({
            "run": int(run),
            "inc": ((a.get("angle") or {}).get("incident_deg") or {}).get("mean"),
            "fnat": (a.get("fnat") or {}).get("fnat_mean"),
            "fnat_min": (a.get("fnat") or {}).get("fnat_min"),
            "bsa": ((a.get("bsa") or {}).get("buried_surface_area") or {}).get("mean"),
            "rmsf": (a.get("rmsf") or {}).get("mean_rmsf_angstrom"),
            "stable": (a.get("contact") or {}).get("n_stable_pairs"),
        })
    return out


def coupling(runs):
    """Does the loosest replica lose the most interface?  -> (n_agree, of 3)."""
    ok = [r for r in runs if None not in (r["inc"], r["fnat"], r["bsa"], r["rmsf"])]
    if len(ok) < 3:
        return None
    loose = max(ok, key=lambda r: r["inc"])
    return (sum([loose["fnat"] == min(r["fnat"] for r in ok),
                 loose["bsa"] == min(r["bsa"] for r in ok),
                 loose["rmsf"] == max(r["rmsf"] for r in ok)]), 3)


def main():
    lib = {r["pdb_id"]: r for r in repdata.library()}
    bulge_rows = dyndata.bulge()
    bulge = {r["pdb_id"]: r["bulge"] for r in bulge_rows}
    plen = {r["pdb_id"]: r["peptide_length"] for r in bulge_rows}
    wsd = repdata.within_sd()
    rep = {q: {r["pdb_id"]: r for r in repdata.replicas(q)}
           for q in ("inc", "q", "rmsip", "rmsf")}
    fwd = {p: r for p, r in lib.items() if p not in REVERSED}
    tt = traj_table()

    axes = [
        scan("Fig3 incident median", {p: r["inc_p50"] for p, r in lib.items()}, "both", "%.1f deg"),
        scan("Fig3 incident p5-p95 span", {p: r["inc_p95"] - r["inc_p5"] for p, r in lib.items()}, "hi", "%.1f deg"),
        scan("Fig3 crossing median", {p: r["cross_p50"] for p, r in fwd.items()}, "both", "%.1f deg"),
        scan("Fig3 crossing p5-p95 span", {p: r["cross_p95"] - r["cross_p5"] for p, r in fwd.items()}, "hi", "%.1f deg"),
        scan("Fig3 peptide bulge (per length)", bulge, "hi", "%.2f A", strata=plen),
        scan("Fig4 incident replica spread", {p: r["spread"] for p, r in rep["inc"].items()}, "hi", "%.1f deg"),
        scan("Fig4 incident spread / own SD", {p: r["spread"] / wsd[p] for p, r in rep["inc"].items() if wsd.get(p)}, "hi", "%.2f x"),
        scan("Fig4 Q replica spread", {p: r["spread"] for p, r in rep["q"].items()}, "hi", "%.3f"),
        scan("Fig4 RMSF replica spread", {p: r["spread"] for p, r in rep["rmsf"].items()}, "hi", "%.2f A"),
        scan("Fig4 RMSF mean", {p: r["mean"] for p, r in rep["rmsf"].items()}, "both", "%.2f A"),
        scan("Fig4 Q mean", {p: r["mean"] for p, r in rep["q"].items()}, "lo", "%.3f"),
        scan("Fig4 RMSIP mean", {p: r["mean"] for p, r in rep["rmsip"].items()}, "lo", "%.3f"),
    ]

    print("Outlier scan -- Tukey fences, k=1.5 (out), k=3.0 (far out)\n")
    tally = defaultdict(list)
    for a in axes:
        f = a["fmt"]
        bound = ("per stratum" if a["strat"] else
                 "< " + f % a["lo15"] if a["side"] == "lo" else
                 "> " + f % a["hi15"] if a["side"] == "hi" else
                 "outside " + f % a["lo15"] + " .. " + f % a["hi15"])
        print("%-34s n=%3d  median %-10s fence %-24s %2d flagged (%d far)"
              % (a["name"], a["n"], f % a["med"], bound, len(a["hits"]),
                 sum(1 for h in a["hits"] if h[2] == "far")))
        for pid, v, deg, key in a["hits"][:6]:
            print("        %-6s %-11s %-8s %s" % (pid, f % v,
                  "" if key is None else "(%s-mer)" % key,
                  "FAR OUT" if deg == "far" else ""))
        if len(a["hits"]) > 6:
            print("        ... %d more" % (len(a["hits"]) - 6))
        print()
        for pid, v, deg, key in a["hits"]:
            tally[pid].append((a["name"], v, deg, f))

    # ---- interface loss, from the trajectories rather than the summary table
    lost = sorted((p, r["run"], r["fnat"]) for p, rs in tt.items() for r in rs
                  if r["fnat_min"] is not None and r["fnat_min"] <= 1e-9)
    print("Interface fully broken at least once (fnat_min = 0):"
          "  %d trajectories in %d complexes of 735 / 245" % (len(lost), len({p for p, _, _ in lost})))
    for p, r, fa in lost:
        print("        %-6s run%d   fnat_mean %.3f" % (p, r, fa))
    print()

    print("=" * 78)
    print("Flagged on two or more axes -- the complexes a reader meets repeatedly.")
    print("`coupling` = of the widest-angle replica's three partners (lowest Fnat,")
    print("lowest BSA, highest RMSF), how many agree.  3/3 = the angle opened")
    print("because the interface let go; 0-2/3 = it moved on a footprint that")
    print("held.  Both are real motion -- neither is a defect.\n")
    multi = sorted(tally.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    for pid, hits in multi:
        if len(hits) < 2:
            continue
        c = coupling(tt.get(pid, []))
        far = sum(1 for h in hits if h[2] == "far")
        print("%-6s %d axes%s%s%s" % (
            pid, len(hits), "  %d far out" % far if far else "",
            "   coupling %d/%d" % c if c else "",
            "   *** %s" % KNOWN_DEFECT[pid] if pid in KNOWN_DEFECT else ""))
        for name, v, deg, f in hits:
            print("        %-34s %-11s %s" % (name, f % v, "FAR OUT" if deg == "far" else ""))
        print()

    n_multi = sum(1 for _, h in multi if len(h) > 1)
    print("%d complexes flagged on one axis, %d on two or more, %d never flagged."
          % (len(multi) - n_multi, n_multi, repdata.N_TOTAL - len(tally)))

    fam = defaultdict(list)
    for p in sorted(tally):
        fam[p[:3]].append(p)
    ser = {k: v for k, v in fam.items() if len(v) > 1}
    print("\n%d of the flagged complexes are not independent: they share a PDB" % sum(len(v) for v in ser.values()))
    print("deposition series, so one system can put several neighbouring points")
    print("into the same tail.")
    for k, v in sorted(ser.items()):
        print("        %s" % " ".join(v))

    print("\nReverse-polarity complexes, held out of the crossing scan because the")
    print("mode is real:")
    for pid in REVERSED:
        print("        %-6s crossing p50 %6.1f deg   incident p50 %5.1f deg"
              % (pid, lib[pid]["cross_p50"], lib[pid]["inc_p50"]))

    out = HERE / "outliers.tsv"
    with out.open("w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["pdb_id", "axis", "value", "degree", "n_axes", "coupling", "known_defect"])
        for pid, hits in multi:
            c = coupling(tt.get(pid, []))
            for name, v, deg, f in hits:
                wr.writerow([pid, name, "%.4f" % v, deg, len(hits),
                             "%d/3" % c[0] if c else "", KNOWN_DEFECT.get(pid, "")])
    print("\n%s  %d rows" % (out.name, sum(len(h) for _, h in multi)))


if __name__ == "__main__":
    main()
