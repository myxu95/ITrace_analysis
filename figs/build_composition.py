#!/usr/bin/env python
"""Build composition.tsv: one row per complex, for the composition figure.

Counting unit is the COMPLEX (245), not the trajectory (735).

Two fields need repair before they can be plotted:

  * MHC allele.  manifest.json carries a four-digit HLA allele for the 219
    human complexes but None for the 26 non-human ones, and per-complex
    rcsb.hla_allele is gene-level ("HLA-A") and unusable.  The mouse heavy
    chains are typed here by sequence identity against three references taken
    from complexes whose curated entity gene symbol is unambiguous, so the
    25 mouse complexes enter the figure as H-2Db / H-2Kb / H-2Dd instead of
    being silently dropped.  The single macaque complex (7byd) is not typed
    by sequence -- its entity gene is only "Mamu-B"; the allele is taken from
    the deposition title and the source paper (Int Immunol 2020;32:805), both
    of which name Mamu-B*05104.

  * Two allele labels name the wrong molecule.  The human labels come from a
    sequence match against the IPD-IMGT/HLA allele database, which names every
    variant ever seen in a person, so an engineered point mutant of a common
    allele can land on the name of a rare natural one.  2uwe is HLA-A*02:01
    T163A and reproduces HLA-A*02:354 exactly; 6q3s is HLA-A*02:01 Y84C/A139C
    and was labelled HLA-A*02:624, which is HLA-A*02:01 Y84C and nothing else,
    so that label is a nearest match that drops A139C.  Both depositors built
    on A*02:01 and say so.  AS_BUILT restores the allele the construct was
    built from; mhc_allele_matched keeps the label the match returned.

  * antigen.category mixes an organism term ("murine") into an origin axis.
    The value is kept verbatim; only the display name is disambiguated, in
    make_composition.py.
"""
from __future__ import annotations

import csv
import difflib
import glob
import json
from pathlib import Path

DATA = Path("/home/xmy/work/data/immunotrace/web_data_1000")
OUT = Path(__file__).resolve().parent / "composition.tsv"

# references whose curated rcsb entity gene symbol / description is unambiguous
MOUSE_REF = {"H-2Kb": "1nam", "H-2Db": "3pqy", "H-2Dd": "5ivx"}
MIN_ID = 0.95           # every mouse chain must clear this against its best ref

# {pdb: (allele the construct was built from, label the sequence match returned)}.
# Both sides were checked.  Against the unmodified A*02:01 chains of 1ao7 and
# 2bnq, 2uwe differs by T163A alone and 6q3s by the disulfide pair Y84C/A139C
# alone; the substitutions are in provenance/mhc_construct_audit.tsv.  Against
# IPD-IMGT/HLA 3.65, A*02:354 (HLA07678) is A*02:01 T163A, so 2uwe reproduces it
# exactly, and A*02:624 (HLA14898) is A*02:01 Y84C, so 6q3s -- which also carries
# A139C -- is not that allele at all.  Both database entries are typed over exons
# 2-3 only, mature residues 1-182, so 139 and 163 are inside the sequenced region
# and neither difference is an artefact of partial coverage; and both rest on one
# unconfirmed observation in one individual.
AS_BUILT = {
    "2uwe": ("HLA-A*02:01", "HLA-A*02:354"),
    "6q3s": ("HLA-A*02:01", "HLA-A*02:624"),
}


def complexes():
    rows = json.loads((DATA / "manifest.json").read_text())
    rows = rows["trajectories"] if isinstance(rows, dict) else rows
    out = {}
    for r in rows:
        out.setdefault(r.get("pdb_id") or r["traj_id"].split("_")[0], r)
    return out


def meta(pdb):
    return json.loads(Path(sorted(glob.glob(f"{DATA}/{pdb}_run*/meta.json"))[0]).read_text())


def mhc_chain_seq(pdb):
    m = meta(pdb)
    ids = [k for k, v in (m.get("chain_roles") or {}).items() if v == "MHC"]
    if not ids:
        return None
    for c in m.get("chains", []):
        if c["id"] == ids[0]:
            return c.get("sequence")
    return None


def type_mouse(pdb, refs):
    s = (mhc_chain_seq(pdb) or "")[:275]
    sc = {k: difflib.SequenceMatcher(None, s, v, autojunk=False).ratio()
          for k, v in refs.items()}
    best = max(sc, key=sc.get)
    return best, sc[best]


def main():
    comp = complexes()
    refs = {k: (mhc_chain_seq(p) or "")[:275] for k, p in MOUSE_REF.items()}

    rows, weak = [], []
    for pdb in sorted(comp):
        r = comp[pdb]
        host = r.get("host_species") or ""
        allele, src, ident = r.get("hla_allele"), "manifest", ""
        if allele:
            src = "manifest/" + (r.get("hla_source") or "?")
        elif host == "Mus musculus":
            allele, q = type_mouse(pdb, refs)
            src, ident = "seq_identity", f"{q:.3f}"
            if q < MIN_ID:
                weak.append((pdb, allele, q))
        elif host == "Macaca mulatta":
            allele, src = "Mamu-B*05104", "deposition_title"
        else:
            allele, src = "unassigned", "none"

        matched = ""
        if pdb in AS_BUILT:
            allele, matched = AS_BUILT[pdb]
            src = "as_built"

        locus = ("H-2" if allele.startswith("H-2")
                 else "Mamu" if allele.startswith("Mamu")
                 else "HLA-" + allele.split("*")[0].replace("HLA-", "") if allele.startswith("HLA-")
                 else "other")
        rows.append(dict(
            pdb_id=pdb,
            peptide_length=r.get("peptide_length") or "",
            peptide_seq=r.get("peptide_seq") or "",
            antigen_category=r.get("antigen_category") or "",
            antigen_organism=r.get("antigen_organism") or "",
            release_year=r.get("release_year") or "",
            host_species=host,
            mhc_allele=allele,
            mhc_locus=locus,
            mhc_source=src,
            mhc_identity=ident,
            mhc_allele_matched=matched,
        ))

    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    print(f"{OUT.name}: {len(rows)} complexes")
    for f in ("peptide_length", "antigen_category", "host_species", "mhc_locus", "mhc_source"):
        print(f"  {f:17s}", dict(sorted(Counter(r[f] for r in rows).items(), key=lambda t: -t[1])))
    miss = [f for f in ("peptide_length", "release_year", "antigen_category")
            if any(r[f] == "" for r in rows)]
    print("  missing fields:", miss or "none")
    print("  weak mouse calls:", weak or "none")


if __name__ == "__main__":
    main()
