"""Merge ANARCI TCR annotation (gene usage + CDR loops) into meta.json.

ANARCI is run separately on the server's `anarci` conda environment over all
TCR chain sequences (IMGT scheme, --restrict tr), producing one CSV per ANARCI
chain_type (A=TCRalpha, B=TCRbeta, G=TCRgamma, D=TCRdelta). Each row gives the
germline V/J genes and the IMGT-numbered residues; CDR1/2/3 are read from the
standard IMGT position ranges. This distils that into each trajectory's
meta.json under `tcr`.

Chain pairing caveat
--------------------
ANARCI frequently labels the alpha chain of a classic alpha/beta TCR as a delta
chain (chain_type "D"), because the TRAV and TRDV germline V-gene sets overlap.
We therefore assign chains by *pairing role*, not by the raw ANARCI letter:

    alpha-position  <- locus TRA or TRD   (delta == mis-called alpha, OR true delta)
    beta-position   <- locus TRB or TRG   (gamma pairs opposite delta)

A complex is only classified as a gamma/delta TCR when a genuine gamma (TRG)
chain is present. In this dataset every "delta" chain is paired with a beta
chain, i.e. they are all alpha/beta TCRs whose alpha chain ANARCI mis-typed as
delta -- so all are reported as `type: "ab"`, with the original locus preserved
under each chain for transparency.

Usage (after fetching tcr_anarci_*.csv from the server):
    python -m pipeline.annotate_tcr /tmp/tcr_anarci_A.csv /tmp/tcr_anarci_B.csv ...
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

from . import config

# IMGT CDR definitions (position ranges, inclusive).
CDR_RANGES = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr3": (105, 117)}

# ANARCI chain_type -> immunological locus.
LOCUS = {"A": "TRA", "B": "TRB", "G": "TRG", "D": "TRD"}
# Pairing position each locus occupies. delta (TRD) sits opposite beta/gamma in
# the same structural slot as alpha; gamma (TRG) sits opposite delta.
ALPHA_LOCI = {"TRA", "TRD"}
BETA_LOCI = {"TRB", "TRG"}

META_COLS = 13  # Id..j_identity; remaining columns are IMGT positions


def _pos_number(header: str) -> int:
    digits = "".join(c for c in header if c.isdigit())
    return int(digits) if digits else -1


def parse_anarci_csv(path: Path) -> list[dict]:
    rows = []
    with path.open(newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        pos_cols = header[META_COLS:]
        pos_num = [_pos_number(h) for h in pos_cols]
        for raw in reader:
            d = dict(zip(header, raw))
            rec = {
                "id": d.get("Id"),
                "chain_type": d.get("chain_type"),
                "v_gene": d.get("v_gene") or None,
                "j_gene": d.get("j_gene") or None,
            }
            pos_vals = raw[META_COLS:]
            for cdr, (lo, hi) in CDR_RANGES.items():
                seq = "".join(
                    pos_vals[i] for i in range(len(pos_vals))
                    if lo <= pos_num[i] <= hi and pos_vals[i] not in ("-", "", " ")
                )
                rec[cdr] = seq
            rows.append(rec)
    return rows


def _position(locus: str) -> str | None:
    if locus in ALPHA_LOCI:
        return "alpha"
    if locus in BETA_LOCI:
        return "beta"
    return None


def build_annotation(csv_paths: list[Path]) -> dict[str, dict]:
    """traj_id -> {alpha/beta: {structural_chain, locus, v_gene, j_gene, cdr1/2/3}}."""
    by_traj: dict[str, dict] = {}
    for p in csv_paths:
        for rec in parse_anarci_csv(p):
            ident = rec["id"] or ""
            if "__" not in ident:
                continue
            traj_id, chain = ident.rsplit("__", 1)
            locus = LOCUS.get((rec["chain_type"] or "").upper())
            position = _position(locus) if locus else None
            if not position:
                continue
            entry = {
                "structural_chain": chain,
                "locus": locus,
                "v_gene": rec["v_gene"],
                "j_gene": rec["j_gene"],
                "cdr1": rec["cdr1"],
                "cdr2": rec["cdr2"],
                "cdr3": rec["cdr3"],
            }
            slot = by_traj.setdefault(traj_id, {})
            if position in slot:
                # Prefer a canonical locus (TRA/TRB) over a borrowed one (TRD) if
                # two chains land on the same position; otherwise keep the first.
                prev = slot[position]
                prefer_new = prev["locus"] in {"TRD", "TRG"} and locus in {"TRA", "TRB"}
                if not prefer_new:
                    continue
            slot[position] = entry
    return by_traj


def classify_type(chains: dict) -> str:
    loci = {c["locus"] for c in chains.values()}
    return "gd" if (loci & {"TRG"}) else "ab"


def main(argv: list[str]) -> int:
    csv_paths = [Path(a) for a in argv if a.endswith(".csv")]
    if not csv_paths:
        print("usage: python -m pipeline.annotate_tcr tcr_anarci_A.csv ...", file=sys.stderr)
        return 2

    by_traj = build_annotation(csv_paths)

    written = 0
    for traj_id, chains in by_traj.items():
        meta_path = config.WEB_DATA / traj_id / config.OUT_META
        if not meta_path.exists():
            continue
        with meta_path.open() as fh:
            meta = json.load(fh)
        meta["tcr"] = {"type": classify_type(chains), "chains": chains}
        with meta_path.open("w") as fh:
            json.dump(meta, fh, indent=2)
        written += 1

    n_gd = sum(1 for c in by_traj.values() if classify_type(c) == "gd")
    print(f"Annotated {written} trajectories with TCR gene usage + CDR loops "
          f"({len(by_traj)} had ANARCI hits; {n_gd} genuine gamma/delta).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
