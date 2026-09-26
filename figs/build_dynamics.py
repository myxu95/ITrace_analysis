#!/usr/bin/env python
"""Freeze the two per-complex dynamics tables Figure 3 draws from.

Counting unit is the COMPLEX.  A complex is simulated three times from the same
starting structure, so a replica is not an independent observation of the
library; every value written here is the mean over that complex's replicas, and
the spread that the figures show is spread ACROSS COMPLEXES.

Two outputs, because they have different shapes:

  dynamics_bulge.tsv   one row per complex -- peptide length and the
                       time-averaged peptide bulge height, analysis.json
                       geometry.bulge_height_angstrom.

  dynamics_pep9.tsv    one row per (complex, peptide position) for the 9-mers,
                       which are the only length numerous enough to read
                       position by position -- analysis.json
                       interface.peptide_table, fields tcr_contact (max contact
                       occupancy to any TCR residue), hla_contact (the same to
                       any MHC residue) and sasa_nm2 (time-averaged
                       Shrake-Rupley residue SASA).

Both are small text files, so the figures rebuild without the 23 GB release
mount and every number in Figure 3 stays auditable from the repository alone.
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

DATA = Path("/home/xmy/work/data/immunotrace/web_data_1000")
HERE = Path(__file__).resolve().parent
COMP = HERE / "composition.tsv"
OUT_BULGE = HERE / "dynamics_bulge.tsv"
OUT_PEP9 = HERE / "dynamics_pep9.tsv"

FOCUS_LENGTH = 9          # the position table is written for 9-mers only
FIELDS = ("tcr_contact", "hla_contact", "sasa_nm2")


def replicas(pdb_id):
    return sorted(DATA.glob(f"{pdb_id}_run*/analysis/analysis.json"))


def mean(v):
    return sum(v) / len(v)


def main():
    comp = list(csv.DictReader(COMP.open(), delimiter="\t"))
    bulge_rows, pep_rows, missing = [], [], []

    for r in comp:
        pid, L = r["pdb_id"], int(r["peptide_length"])
        runs = replicas(pid)
        if not runs:
            missing.append(pid)
            continue
        aj = [json.loads(p.read_text()) for p in runs]

        b = [a["geometry"]["bulge_height_angstrom"] for a in aj
             if (a.get("geometry") or {}).get("bulge_height_angstrom") is not None]
        bulge_rows.append({"pdb_id": pid, "peptide_length": L,
                           "n_replicas": len(runs), "n_bulge": len(b),
                           "bulge_angstrom": f"{mean(b):.3f}" if b else ""})

        if L != FOCUS_LENGTH:
            continue
        acc = defaultdict(lambda: defaultdict(list))
        for a in aj:
            for p in (a.get("interface") or {}).get("peptide_table") or []:
                for f in FIELDS:
                    if p.get(f) is not None:
                        acc[p["position"]][f].append(p[f])
        for pos in sorted(acc):
            row = {"pdb_id": pid, "position": pos,
                   "n_replicas": len(acc[pos][FIELDS[0]])}
            for f in FIELDS:
                v = acc[pos][f]
                row[f] = f"{mean(v):.4f}" if v else ""
            pep_rows.append(row)

    def write(path, rows):
        with path.open("w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            wr.writeheader()
            wr.writerows(rows)
        print(f"{path.name:22s} {len(rows)} rows")

    write(OUT_BULGE, bulge_rows)
    write(OUT_PEP9, pep_rows)

    n9 = sum(1 for r in bulge_rows if r["peptide_length"] == FOCUS_LENGTH)
    print(f"   complexes {len(bulge_rows)}   9-mers {n9}   "
          f"positions/9-mer {len(pep_rows) / n9:.2f}")
    if missing:
        print("   NO TRAJECTORY:", missing)
    short = [r["pdb_id"] for r in bulge_rows if r["n_bulge"] != r["n_replicas"]]
    if short:
        print("   BULGE MISSING IN SOME REPLICA:", short)


if __name__ == "__main__":
    main()
