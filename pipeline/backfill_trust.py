"""Backfill the §1 trust-gate fields into served analysis.json without a full re-curation.

``dynamics_trust``, ``bsa_decomposition.reconciliation`` and ``trust_tier`` are
pure derivations of data ALREADY in each served ``analysis.json`` (+ the
``meta.quality.equilibration`` flag). A full ``extract_analysis --source`` re-run
would rebuild the whole file from the md_analysis source CSVs (re-running
cdr_contacts / interface_metrics) — overkill, slow, and risky if the source root
is incomplete. This tool recomputes ONLY the three trust fields in place, reusing
the exact same functions extract_analysis calls (single source of truth), so
served data can never regress. Idempotent; re-run whenever a threshold changes.

    IMMUNO_WEB_DATA=.../web_data python -m pipeline.backfill_trust [--ids a,b] [--dry-run]
"""
from __future__ import annotations

import argparse
import json

from . import config
from . import essential_dynamics
from . import struct_metrics
from . import trust


def compute_trust_fields(analysis: dict, meta: dict) -> dict:
    """Return {dynamics_trust, trust_tier, bsa_reconciliation-or-None} for one trajectory."""
    ed = analysis.get("essential_dynamics") or {}
    cdr_reliable = (analysis.get("tcr_cdr") or {}).get("reliable")
    recon = struct_metrics.bsa_reconciliation(analysis.get("bsa"), analysis.get("bsa_decomposition"))
    dyn = essential_dynamics.dynamics_trust(
        pc1_cosine=ed.get("pc1_cosine_content"),
        rmsip=ed.get("subspace_rmsip"),
        cdr_reliable=cdr_reliable,
    )
    tier = trust.trust_tier(
        dyn_verdict=dyn,
        bsa_recon=recon,
        cdr_reliable=cdr_reliable,
        mode_soft=None,
    )
    return {"dynamics_trust": dyn, "trust_tier": tier, "reconciliation": recon}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", default=None, help="Comma-separated traj ids (default: all).")
    ap.add_argument("--dry-run", action="store_true", help="Compute + print, do not write.")
    args = ap.parse_args(argv)
    want = {i.strip() for i in args.ids.split(",")} if args.ids else None

    n_ok = n_skip = 0
    from collections import Counter
    tiers, levels = Counter(), Counter()
    for an_path in sorted(config.WEB_DATA.glob("*/analysis/analysis.json")):
        tid = an_path.parent.parent.name
        if want is not None and tid not in want:
            continue
        meta_path = config.WEB_DATA / tid / config.OUT_META
        if not meta_path.exists():
            n_skip += 1
            continue
        analysis = json.loads(an_path.read_text())
        meta = json.loads(meta_path.read_text())
        fields = compute_trust_fields(analysis, meta)
        analysis["dynamics_trust"] = fields["dynamics_trust"]
        analysis["trust_tier"] = fields["trust_tier"]
        if fields["reconciliation"] is not None and analysis.get("bsa_decomposition"):
            analysis["bsa_decomposition"]["reconciliation"] = fields["reconciliation"]
        tiers[fields["trust_tier"]["tier"]] += 1
        levels[fields["dynamics_trust"]["level"]] += 1
        if args.dry_run:
            print(f"{tid}: tier={fields['trust_tier']['tier']} "
                  f"dyn={fields['dynamics_trust']['level']} recon={fields['reconciliation']}")
        else:
            an_path.write_text(json.dumps(analysis, indent=2), encoding="utf-8")
        n_ok += 1
    print(f"\nBackfill trust: {n_ok} trajectories{' (dry-run)' if args.dry_run else ' written'}, {n_skip} skipped.")
    print(f"  trust_tier: {dict(tiers)}")
    print(f"  dynamics_trust: {dict(levels)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
