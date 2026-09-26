#!/usr/bin/env python
"""Re-derive TCR V/J germline assignments for the 245 ITrace complexes.

Why this exists
---------------
ANARCI HMM-types the alpha chain of many classic alpha/beta TCRs as a *delta*
chain, because the TRAV and TRDV germline V-gene sets overlap (already noted in
pipeline/annotate_tcr.py).  What was NOT handled is the knock-on effect on the
J call: anarci.run_germline_assignment() searches
    all_germlines["J"][chain_type][species]
i.e. once the chain is typed "D" the J gene can ONLY be a TRDJ.  The result is
22/245 complexes carrying an alpha-chain J gene named TRDJ1..TRDJ4 with
germline identity 0.50-0.64, where the true TRAJ scores 0.93-1.00.

Here we renumber every TCR chain, then score the *same* IMGT state sequence
against BOTH the TRAJ and TRDJ germline sets (and TRBJ/TRGJ for the beta slot)
and keep the best hit over the union.  V genes are left as ANARCI called them
(TRAV/DV names are legitimate IMGT names).

Writes vj_genotypes.tsv (one row per complex).
"""
from __future__ import annotations
import glob, json, os, sys
from pathlib import Path

DATA = Path("/home/xmy/work/data/immunotrace/web_data_1000")
OUT = Path(__file__).resolve().parent / "vj_genotypes.tsv"

from anarci.anarci import anarci, all_germlines, get_identity

ALPHA_J = ("A", "D")   # TRAJ + TRDJ, alpha structural slot
BETA_J = ("B", "G")    # TRBJ + TRGJ, beta structural slot


def state_sequence(seq: str, numbering, start: int) -> str:
    """Rebuild ANARCI's 128-position match-state string from the numbering."""
    idx, pos_of = start, {}
    for (pos, ins), aa in numbering:
        if aa != "-":
            if ins == " ":
                pos_of[pos] = idx
            idx += 1
    return "".join(seq[pos_of[k]] if k in pos_of else "-" for k in range(1, 129))


def best_j(state_seq: str, species: str, chain_types) -> tuple[str, float]:
    hits = {}
    for ct in chain_types:
        for gene, germ in all_germlines["J"].get(ct, {}).get(species, {}).items():
            hits[gene] = get_identity(state_seq, germ)
    if not hits:
        return None, 0.0
    g = max(hits, key=hits.get)
    return g, hits[g]


def species_of(meta: dict) -> str:
    ents = (meta.get("rcsb") or {}).get("entities") or []
    orgs = {e.get("organism") for e in ents
            if e.get("role") in ("tcr_alpha", "tcr_beta") and e.get("organism")}
    if not orgs:
        tcr = {c for c, r in (meta.get("chain_roles") or {}).items()
               if r in ("TCRα", "TCRβ")}
        for e in ents:
            if set(e.get("chains") or []) & tcr and e.get("organism"):
                orgs.add(e["organism"])
    if not orgs:
        for e in ents:
            if "CELL RECEPTOR" in (e.get("description") or "").upper() and e.get("organism"):
                orgs.add(e["organism"])
    return "|".join(sorted(orgs)) or "?"


def main() -> None:
    metas = {}
    for f in sorted(DATA.glob("*/meta.json")):
        pdb = f.parent.name.split("_run")[0]
        metas.setdefault(pdb, json.loads(f.read_text()))

    cols = ["pdb", "species",
            "a_chain", "a_v", "a_v_id", "a_j_stored", "a_j_stored_id",
            "a_j_fixed", "a_j_fixed_id", "a_cdr3",
            "b_chain", "b_v", "b_v_id", "b_j_stored", "b_j_stored_id",
            "b_j_fixed", "b_j_fixed_id", "b_cdr3"]
    rows = []
    for n, (pdb, m) in enumerate(sorted(metas.items()), 1):
        roles = m.get("chain_roles") or {}
        seqs = {c["id"]: c["sequence"] for c in m.get("chains", [])}
        sp_txt = species_of(m)
        sp = "mouse" if sp_txt == "Mus musculus" else "human"
        rec = {"pdb": pdb, "species": sp_txt}
        ident = {}
        aj = DATA / f"{pdb}_run1" / "analysis" / "analysis.json"
        if aj.exists():
            ident = json.loads(aj.read_text())["identity"]["tcr_identity"]
        for slot, role, cts in (("a", "TCRα", ALPHA_J), ("b", "TCRβ", BETA_J)):
            ch = [c for c, r in roles.items() if r == role]
            meta_ch = (m.get("tcr") or {}).get("chains", {}).get(
                "alpha" if slot == "a" else "beta", {})
            rec[f"{slot}_chain"] = ch[0] if ch else ""
            rec[f"{slot}_v"] = meta_ch.get("v_gene") or ""
            rec[f"{slot}_v_id"] = ident.get(f"{'alpha' if slot=='a' else 'beta'}_v_identity") or ""
            rec[f"{slot}_j_stored"] = meta_ch.get("j_gene") or ""
            rec[f"{slot}_j_stored_id"] = ident.get(f"{'alpha' if slot=='a' else 'beta'}_j_identity") or ""
            rec[f"{slot}_cdr3"] = meta_ch.get("cdr3") or ""
            rec[f"{slot}_j_fixed"] = ""
            rec[f"{slot}_j_fixed_id"] = ""
            if not ch or ch[0] not in seqs:
                continue
            seq = seqs[ch[0]]
            num, det, _ = anarci([(pdb, seq)], scheme="imgt", database="ALL",
                                 allow={"A", "B", "G", "D"}, assign_germline=False,
                                 allowed_species=["human", "mouse"])
            if not num[0]:
                continue
            numbering, start, _end = num[0][0]
            ss = state_sequence(seq, numbering, start)
            g, i = best_j(ss, sp, cts)
            rec[f"{slot}_j_fixed"] = g or ""
            rec[f"{slot}_j_fixed_id"] = f"{i:.3f}"
        rows.append(rec)
        if n % 40 == 0:
            print(f"  ... {n}/{len(metas)}", file=sys.stderr)

    with OUT.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"wrote {OUT}  ({len(rows)} complexes)")


if __name__ == "__main__":
    main()
