#!/usr/bin/env python
"""Freeze the per-complex replica table that Figures 3 and 4 draw from.

Two different questions are answered from the same trajectories and they need
the numbers aggregated two different ways, so both live in one table:

  Figure 3 asks how the LIBRARY is spread -- what range of approach geometry
  the 245 complexes cover, and how far one complex's angle travels during a
  simulation.  For that the three replicas of a complex are pooled into a
  single 3 x 1001-frame sample and summarised by percentiles.  Those
  percentiles are read from angle_percentiles.tsv rather than computed here:
  analysis.json stores the angle as a 150-point plot downsample, so pooling it
  would describe 450 points as 3003.  Run build_angle_frames.py first, and
  again after any trajectory repair.

  Figure 4 asks whether the three replicas of one complex AGREE.  For that the
  replicas must stay apart, so each one's own mean is written separately and
  the figure reads their spread.

Percentiles, not min/max: over 3003 frames a single excursion sets a min or a
max, so the extremes describe one frame rather than the trajectory.  p5-p95 is
what "the angle normally sits here" means.

Two provenance flags travel with the rows so no figure can silently average a
known-bad number (see also the audit in the manuscript notes):

  dup_analysis   two replicas carrying bit-identical analysis blocks for every
                 trajectory-derived quantity while their trajectories differ, so
                 the complex has two independent analyses and not three and its
                 between-replica spread is deflated by construction.  6g9q was
                 the one case; its three analysis records were rewritten on
                 2026-09-14 and now differ, so no complex carries this flag.
                 defect_state() still tests for it on every build.
  chain_selection  4prh, 5hhm, 7rk7 -- run3's chain mapping labels both TCR
                 chains alpha, so the alpha-beta contacts inside the TCR are
                 counted as TCR-MHC interface contacts.  tcr_cdr.by_chain
                 records it plainly: beta n = 0, alpha_contribution = 1.0, and
                 the framework-to-MHC term carries the excess (n = 237 and 231
                 against 18 and 5 in their own siblings) while the CDR3 terms
                 stay ordinary.  Everything built on that selection inflates
                 together -- n_pairs 551-649 against a library maximum of 238,
                 n_stable_pairs 225-239 against 73, buried surface 5.2-5.7x,
                 hydrogen bonds 420-472 against 183, hydrophobic contacts
                 80-104 against 36.  Salt bridges are the exception and stay
                 normal, which is itself the tell: a protein has few salt
                 bridges inside itself but many hydrogen bonds.
                 The trajectories are sound.  Fnat (0.62-0.80), the approach
                 angles and the one-sided areas in bsa_decomposition all match
                 their siblings, so nothing physically dissociated; only the
                 chain-labelled analyses are wrong.  Those are bsa, contact and
                 interactions.  Angle, Fnat, RMSF and RMSIP come from other
                 code paths and stay.

Six further trajectories were checked and cleared.  Their buried surface is
1.4-1.7x what their decomposition implies, but their contact counts are
ordinary and they hold the five lowest Fnat values in the library: a loose
interface is where two solvent-accessibility engines disagree most, which is a
definitional difference and not the chain-mapping fault.  In the lowest Fnat
quintile 5 of 146 trajectories exceed 1.3x; in the highest, none of 147 do.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

DATA = Path("/home/xmy/work/data/immunotrace/web_data_1000")
HERE = Path(__file__).resolve().parent
COMP = HERE / "composition.tsv"
ANGLES = HERE / "angle_percentiles.tsv"
OUT = HERE / "replicas.tsv"

# 6g9q used to sit here: its run-1 and run-2 analysis records were byte-identical.
# The three records were rewritten on 2026-09-14 and now carry distinct Fnat,
# incident angle and RMSF, so the defect is gone and the complex is no longer
# excluded.  defect_state() re-derives this from the data on every build, so an
# entry left here after a repair would silently drop a good complex.
DUP_ANALYSIS: set[str] = set()
CHAIN_SELECTION = {"4prh", "5hhm", "7rk7"}


def defect_state(aj):
    """-> the defect names one complex's live analyses still show.

    The two exclusions below are written out by hand so a reader can see who is
    dropped and why, but a hand-written exclusion outlives the defect it names:
    once the pipeline is repaired, a stale entry here would silently drop a
    complex that is now perfectly good.  So the flags are re-derived from the
    data on every build and ``main`` complains if the two disagree in either
    direction.
    """
    out = set()
    # Duplicated analysis: two replicas carrying identical trajectory-derived
    # numbers.  Independent 200 ns runs never agree to the last decimal on all
    # three of these at once, so a collision means one block was copied.
    sig = [(a["fnat"]["fnat_mean"], a["angle"]["incident_deg"]["mean"],
            a["rmsf"]["mean_rmsf_angstrom"]) for a in aj]
    if len(sig) != len(set(sig)):
        out.add("dup_analysis")
    # Chain double-assignment: the pipeline's own reconciliation compares the
    # buried surface scalar against the decomposition and already catches it.
    for a in aj:
        rec = (a.get("bsa_decomposition") or {}).get("reconciliation") or {}
        if rec.get("consistent") is False:
            out.add("chain_selection")
    return out


def main():
    if not ANGLES.exists():
        raise SystemExit(
            f"{ANGLES.name} is missing -- run build_angle_frames.py first "
            f"(~20 min over 735 trajectories), and re-run it after any "
            f"trajectory repair lands in web_data_1000.")
    angles = {r["pdb_id"]: r for r in csv.DictReader(ANGLES.open(), delimiter="\t")}
    rows, missing, seen = [], [], {}
    for r in csv.DictReader(COMP.open(), delimiter="\t"):
        pid = r["pdb_id"]
        runs = sorted(DATA.glob(f"{pid}_run*/analysis/analysis.json"))
        if not runs:
            missing.append(pid)
            continue
        aj = [json.loads(p.read_text()) for p in runs]
        seen[pid] = defect_state(aj)
        if pid not in angles:
            raise RuntimeError(f"{pid}: not in {ANGLES.name}; regenerate it")

        # Pooled percentiles come from build_angle_frames.py, which recomputes the
        # angle at full frame resolution.  analysis.json only carries the 150-point
        # plot downsample, and the figure states a 3 x 1001-frame basis.
        ap = angles[pid]
        if int(ap["n_frames"]) != len(aj) * 1001:
            raise RuntimeError(
                f"{pid}: angle_percentiles.tsv has {ap['n_frames']} frames for "
                f"{len(aj)} replicas; expected {len(aj) * 1001}")
        row = {"pdb_id": pid, "n_replicas": len(aj)}
        for tag in ("inc", "cross"):
            for p in (5, 50, 95):
                row[f"{tag}_p{p}"] = ap[f"{tag}_p{p}"]
        # per-replica, kept apart: Figure 4 reads the spread of these
        for i, a in enumerate(aj, 1):
            row[f"inc_run{i}"] = f"{a['angle']['incident_deg']['mean']:.3f}"
            row[f"incsd_run{i}"] = f"{a['angle']['incident_deg']['std']:.3f}"
            row[f"cross_run{i}"] = f"{a['angle']['crossing_deg']['mean']:.3f}"
            row[f"q_run{i}"] = f"{a['fnat']['q_mean']:.4f}"
            row[f"stable_run{i}"] = str(a["contact"]["n_stable_pairs"])
            row[f"rmsip_run{i}"] = f"{a['essential_dynamics']['subspace_rmsip']:.4f}"
            row[f"rmsf_run{i}"] = f"{a['rmsf']['mean_rmsf_angstrom']:.4f}"
        row["flag"] = ("dup_analysis" if pid in DUP_ANALYSIS else
                       "chain_selection" if pid in CHAIN_SELECTION else "")
        rows.append(row)

    with OUT.open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)
    print(f"{OUT.name:16s} {len(rows)} complexes, "
          f"{sum(r['n_replicas'] for r in rows)} trajectories")
    short = [r["pdb_id"] for r in rows if r["n_replicas"] != 3]
    print(f"   flagged  dup_analysis {sorted(DUP_ANALYSIS)}   chain_selection {sorted(CHAIN_SELECTION)}")
    for name, declared in (("dup_analysis", DUP_ANALYSIS),
                           ("chain_selection", CHAIN_SELECTION)):
        found = {p for p, d in seen.items() if name in d}
        if found - declared:
            print(f"   !! {name} NOT DECLARED: {sorted(found - declared)}"
                  f"  -- new occurrence, add it before trusting this table")
        if declared - found:
            print(f"   !! {name} NO LONGER PRESENT: {sorted(declared - found)}"
                  f"  -- repaired upstream? drop it from the set and redraw")
    if short:
        print("   NOT THREE REPLICAS:", short)
    if missing:
        print("   NO TRAJECTORY:", missing)


if __name__ == "__main__":
    main()
