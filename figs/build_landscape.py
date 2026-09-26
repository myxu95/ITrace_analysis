#!/usr/bin/env python
"""Freeze the conformational-landscape tables Figure 4 draws from.

Every deposited trajectory carries three related analyses, all recomputed by the
pipeline from the stored frames and all served under ``analysis/analysis.json``:

  peptide_dpca    dihedral principal-component analysis of the peptide -- its
                  backbone phi/psi, the chi1/chi2 of every side chain that has
                  them, and the Calpha-Calpha distances, normalised so no block
                  dominates -- reduced to two components and rendered as an
                  80x80 free-energy surface, kcal/mol above the global minimum.

  tcr_cdr3_dpca   the same construction for the CDR3 alpha and beta loops
                  together, superposed on the Valpha/Vbeta framework so that
                  loop motion is measured relative to the domain carrying it
                  rather than to the box.

  coupled_states  the joint occupancy of peptide basin x CDR3 basin over the
                  1001 frames, with the normalised mutual information between
                  the two assignments and a verdict derived from it.

BASIN, NOT SUBSTATE.  Each dPCA block stores two different decompositions of the
same scores and they do not agree: ``basins`` are local minima of the free-energy
surface (kept if within 2.5 kcal/mol of the global minimum and holding at least
8% of the frames, at most four), while ``n_substates`` counts k-means clusters
of the raw two-dimensional scores.  For 1nam run1 that is 4 against 3.  The
coupling matrix is built by Voronoi assignment to the BASIN centres, so the
basin count is the one that indexes the matrix and the only one this table
carries; anything drawn or written from it must say basin.

The reliability flag is deliberately not carried either.  ``tcr_cdr3_dpca.
reliable`` gates the pipeline's per-CDR CONTACT decomposition, not the
conformational landscape -- the landscape is computed and stored for every
trajectory regardless -- so importing that flag here would silently rebrand 222
sound landscapes as unreliable ones.

Two outputs:

  landscape.tsv            one row per trajectory: the two basin counts, the two
                           PC1/PC2 variance fractions, the mutual information,
                           the verdict, and the entry's cross-replica agreement.
                           This is what the figure's legend reads its library
                           context off.

  landscape_exemplar.json  the plottable parts for the one drawn trajectory: the
                           two free-energy grids, their basin centres and
                           populations, and the joint-occupancy matrix.  Frozen
                           as JSON and not as a table because a grid is not a
                           table; the per-basin side-chain rotamer lists and the
                           raw feature matrices are dropped, since the figure
                           does not draw them and they are a hundredfold larger
                           than what it does.

A verdict of "n/a" means the matrix is degenerate: one side resolved fewer than
two basins, so there is nothing for the mutual information to be computed
between.  It is a property of that trajectory's landscape, not a failure, and
those rows stay in the table with their reason.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

DATA = Path("/home/xmy/work/data/immunotrace/web_data_1000")
HERE = Path(__file__).resolve().parent
OUT_TSV = HERE / "landscape.tsv"
OUT_EX = HERE / "landscape_exemplar.json"

EXEMPLAR = "3w0w_run2"
ZMAX = 6.0          # the pipeline clips the free-energy surface here
BASIN_ZMAX = 2.5    # peptide_dihedrals.py: a minimum this far above the global one
BASIN_MINPOP = 0.08


def keep_dpca(b):
    """The plottable part of a dPCA block: the grid, the basin centres, the axes."""
    return dict(
        pc1_var_frac=b["pc1_var_frac"],
        pc2_var_frac=b["pc2_var_frac"],
        fel={k: b["fel"][k] for k in ("x", "y", "z", "x_label", "y_label", "unit")},
        basins=[{k: v[k] for k in ("pc1", "pc2", "dG", "pop")} for v in b["basins"]],
    )


def main():
    rows, ex = [], None
    for p in sorted(DATA.glob("*_run*/analysis/analysis.json")):
        tid = p.parent.parent.name
        pid, run = tid.rsplit("_run", 1)
        d = json.loads(p.read_text())
        pep, cdr, cs = d["peptide_dpca"], d["tcr_cdr3_dpca"], d["coupled_states"]
        ra = cs.get("replica_agreement") or {}
        assert len(pep["basins"]) == cs["n_peptide"], tid
        assert len(cdr["basins"]) == cs["n_cdr3"], tid
        rows.append(dict(
            pdb_id=pid, run=int(run),
            pep_basins=len(pep["basins"]), cdr3_basins=len(cdr["basins"]),
            pep_pc1=f"{pep['pc1_var_frac']:.4f}", pep_pc2=f"{pep['pc2_var_frac']:.4f}",
            cdr3_pc1=f"{cdr['pc1_var_frac']:.4f}", cdr3_pc2=f"{cdr['pc2_var_frac']:.4f}",
            nmi="" if cs["nmi"] is None else f"{cs['nmi']:.4f}",
            coupling=cs["coupling"], coupling_reason=cs.get("reason") or "",
            replica_n=ra.get("n", ""), replica_n_coupled=ra.get("n_coupled", ""),
            replica_reproducible="yes" if ra.get("reproducible") else ""))
        if tid == EXEMPLAR:
            ex = dict(
                traj_id=tid,
                peptide=keep_dpca(pep), cdr3=keep_dpca(cdr),
                coupled={k: cs[k] for k in ("matrix", "peptide_pop", "cdr3_pop",
                                            "enrichment", "nmi", "coupling",
                                            "top_pair", "replica_agreement")},
                zmax=ZMAX, basin_zmax=BASIN_ZMAX, basin_min_pop=BASIN_MINPOP)

    rows.sort(key=lambda r: (r["pdb_id"], r["run"]))
    with OUT_TSV.open("w", newline="") as fh:
        wr = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)
    if ex is None:
        raise SystemExit(f"exemplar {EXEMPLAR} not found")
    OUT_EX.write_text(json.dumps(ex, separators=(",", ":")))

    n_pid = len({r["pdb_id"] for r in rows})
    print(f"{OUT_TSV.name}   {len(rows)} trajectories, {n_pid} complexes")
    print(f"{OUT_EX.name}   {EXEMPLAR}, "
          f"{ex['coupled']['n_peptide'] if 'n_peptide' in ex['coupled'] else len(ex['peptide']['basins'])}"
          f" x {len(ex['cdr3']['basins'])} basins, NMI {ex['coupled']['nmi']}, "
          f"{OUT_EX.stat().st_size / 1024:.0f} kB")
    from collections import Counter
    print("   coupling:", dict(Counter(r["coupling"] for r in rows)))
    print("   peptide basins:", dict(sorted(Counter(r["pep_basins"] for r in rows).items())))
    print("   CDR3 basins:   ", dict(sorted(Counter(r["cdr3_basins"] for r in rows).items())))
    miss = [r for r in rows if not r["replica_n"]]
    if miss:
        print(f"   !! {len(miss)} rows carry no replica_agreement block")


if __name__ == "__main__":
    main()
