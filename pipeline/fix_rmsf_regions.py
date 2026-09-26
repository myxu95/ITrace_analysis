"""Correct RMSF region labels on trajectories whose md_analysis per-chain role
resolution disagrees with the authoritative meta.chain_roles (the 6 run3:
1fo0/1nam/2ol3/4prh/5hhm/7rk7, where MHC↔TCRβ or TCRα↔TCRβ were swapped).

The topology is byte-identical across a complex's replicas, so the region of a
given (chain, resid) is a STRUCTURAL constant. We therefore repair a mislabelled
replica by borrowing each residue's region from a sibling replica whose
rmsf-derived roles DO match meta.chain_roles. Values (RMSF) are untouched; only
the region tag + per-region summaries are corrected.

    IMMUNO_WEB_DATA=.../web_data python -m pipeline.chain_roles      # first (meta)
    IMMUNO_WEB_DATA=.../web_data python -m pipeline.fix_rmsf_regions [--dry-run]
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict

from . import config
from .chain_roles import role_map

_MHC_REG = {"MHC floor", "α1 helix", "α2 helix"}


def _rmsf_role_of_chain(residues) -> dict:
    byc = defaultdict(Counter)
    for r in residues:
        byc[r["chain"]][r["region"]] += 1
    out = {}
    for c, cnt in byc.items():
        reg = cnt.most_common(1)[0][0]
        out[c] = ("peptide" if reg == "peptide" else "β2m" if reg == "β2m"
                  else "MHC" if reg in _MHC_REG
                  else "TCRα" if reg.endswith("α") else "TCRβ" if reg.endswith("β") else "other")
    return out


def _agrees(residues, roles) -> bool:
    rr = _rmsf_role_of_chain(residues)
    return all(roles.get(c) == r for c, r in rr.items())


def _rebuild_summary(residues):
    region_vals = defaultdict(list)
    for r in residues:
        region_vals[r["region"]].append(r["rmsf"])
    from .interface_metrics import _RMSF_REGION_ORDER
    regions = {reg: {"mean": round(sum(v) / len(v), 2), "max": round(max(v), 2), "n": len(v)}
               for reg, v in region_vals.items()}
    order = [reg for reg in _RMSF_REGION_ORDER if reg in regions]
    return regions, order


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    # group trajectories by complex
    by_pdb = defaultdict(list)
    for meta_path in sorted(config.WEB_DATA.glob("*/" + config.OUT_META)):
        tid = meta_path.parent.name
        by_pdb[tid.rsplit("_", 1)[0]].append(tid)

    fixed, skipped = [], []
    for pdb, tids in by_pdb.items():
        # load each replica's meta-roles + rmsf residues
        data = {}
        for tid in tids:
            an = config.WEB_DATA / tid / "analysis" / "analysis.json"
            mp = config.WEB_DATA / tid / config.OUT_META
            if not (an.exists() and mp.exists()):
                continue
            meta = json.loads(mp.read_text())
            a = json.loads(an.read_text())
            res = (a.get("rmsf_profile") or {}).get("residues")
            if res:
                data[tid] = {"an": an, "a": a, "roles": role_map(meta), "res": res}
        good = [t for t, d in data.items() if _agrees(d["res"], d["roles"])]
        bad = [t for t, d in data.items() if t not in good]
        if not bad:
            continue
        if not good:
            skipped += [(t, "no good sibling") for t in bad]
            continue
        ref = data[good[0]]
        ref_region = {(r["chain"], r["resid"]): r["region"] for r in ref["res"]}
        for t in bad:
            d = data[t]
            keys = [(r["chain"], r["resid"]) for r in d["res"]]
            if any(k not in ref_region for k in keys):
                skipped.append((t, "topology mismatch vs sibling"))
                continue
            for r in d["res"]:
                r["region"] = ref_region[(r["chain"], r["resid"])]
            regions, order = _rebuild_summary(d["res"])
            d["a"]["rmsf_profile"]["regions"] = regions
            d["a"]["rmsf_profile"]["order"] = order
            if not args.dry_run:
                d["an"].write_text(json.dumps(d["a"], indent=2), encoding="utf-8")
            fixed.append((t, f"borrowed from {good[0]}"))

    print(f"fix_rmsf_regions: {len(fixed)} fixed{' (dry-run)' if args.dry_run else ''}, {len(skipped)} skipped.")
    for t, why in fixed: print(f"  FIXED  {t}: {why}")
    for t, why in skipped: print(f"  SKIP   {t}: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
