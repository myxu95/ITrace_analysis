"""Authoritative chain -> biological-role mapping — the single source of truth.

Every per-chain role decision on the site historically trusted whatever each
analysis stage guessed, and those guesses disagree: md_analysis's per-chain
`component` label is correct-but-non-canonical for ~63 complexes (the source PDB
simply doesn't number chains A=MHC…E=TCRβ) and genuinely WRONG for 6 run3
trajectories (e.g. 4prh_run3 calls the real MHC chain "TCRβ"), even though the
topology is byte-identical to run1/run2. `meta` however resolves the roles
correctly and consistently across replicas (peptide_chain + ANARCI TCR chains),
so it is the reliable source.

    role(meta) -> {chain_id: "peptide"|"β2m"|"MHC"|"TCRα"|"TCRβ"|"other"}

Run as a module to backfill `chain_roles` into every meta.json:

    IMMUNO_WEB_DATA=.../web_data python -m pipeline.chain_roles [--dry-run]
"""
from __future__ import annotations

import argparse
import json

from . import config


def role_map(meta: dict) -> dict:
    """{chain_id: role} for one trajectory, derived only from meta.

    Order of resolution: the explicitly-annotated chains first (TCR α/β from
    ANARCI, peptide from peptide_chain), then MHC = largest remaining chain and
    β2m = a ~99-residue remaining chain. Peptide falls back to the shortest
    unassigned chain when meta lacks peptide_chain (e.g. 7byd), since the antigen
    peptide is always the shortest chain.
    """
    nres = {c["id"]: c.get("n_residues", 0) for c in (meta.get("chains") or [])}
    roles: dict[str, str] = {}
    tc = (meta.get("tcr") or {}).get("chains") or {}
    ta = (tc.get("alpha") or {}).get("structural_chain")
    tb = (tc.get("beta") or {}).get("structural_chain")
    pep = meta.get("peptide_chain")
    if ta in nres: roles[ta] = "TCRα"
    if tb in nres: roles[tb] = "TCRβ"
    if pep in nres: roles[pep] = "peptide"
    # peptide fallback: the shortest still-unassigned chain
    if "peptide" not in roles.values():
        cand = sorted((c for c in nres if c not in roles), key=lambda c: nres[c])
        if cand:
            roles[cand[0]] = "peptide"
    # remaining chains: largest = MHC heavy chain; a β2m-sized one = β2m; rest = other
    rem = sorted((c for c in nres if c not in roles), key=lambda c: nres[c], reverse=True)
    for i, c in enumerate(rem):
        roles[c] = "MHC" if i == 0 else ("β2m" if 80 <= nres[c] <= 130 else "other")
    return roles


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    n = 0
    for meta_path in sorted(config.WEB_DATA.glob("*/" + config.OUT_META)):
        meta = json.loads(meta_path.read_text())
        meta["chain_roles"] = role_map(meta)
        if not args.dry_run:
            meta_path.write_text(json.dumps(meta, indent=2))
        n += 1
    print(f"chain_roles: {n} meta.json {'(dry-run)' if args.dry_run else 'written'}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
