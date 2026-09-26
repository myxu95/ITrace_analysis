"""Re-merge peptide-HLA contact / SASA / anchor into served analysis.json from the sidecar.

The peptide_table's ``hla_contact`` / ``sasa_nm2`` / ``anchor`` columns come from
``struct_metrics`` (role/geometry-resolved, written to the ``struct_metrics.json``
sidecar). But later curation passes that rebuilt the ``interface`` block (via
``interface_metrics``) recreated the peptide_table WITHOUT re-applying that merge,
so ~507/735 served trajectories ended up with an all-null "HLA contact" column on
the Peptide recognition profile — the sidecar was correct, the merge was simply
lost. (See ``_merge_peptide_struct`` in extract_analysis: the full path merges, a
partial rebuild does not.)

This tool re-applies ONLY the peptide_table merge from each existing sidecar,
reusing ``struct_metrics.merge_peptide_table`` (single source of truth). It does
NOT touch ``geometry`` / ``bsa_decomposition`` — those survived the rebuild and
carry the trust ``reconciliation`` sub-key we must not clobber. Idempotent:
re-merging an already-correct table writes the identical values.

    IMMUNO_WEB_DATA=.../web_data python -m pipeline.backfill_peptide_hla [--ids a,b] [--dry-run]
"""
from __future__ import annotations

import argparse
import json

from . import config
from . import struct_metrics


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", default=None, help="Comma-separated traj ids (default: all).")
    ap.add_argument("--dry-run", action="store_true", help="Compute + report, do not write.")
    args = ap.parse_args(argv)
    want = {i.strip() for i in args.ids.split(",")} if args.ids else None

    n_written = n_skip = n_nosidecar = n_notable = 0
    fixed_null = []  # trajectories that had an all-null HLA column before the remerge
    for an_path in sorted(config.WEB_DATA.glob("*/analysis/analysis.json")):
        tid = an_path.parent.parent.name
        if want is not None and tid not in want:
            continue
        sidecar_path = an_path.parent / struct_metrics.SIDECAR
        if not sidecar_path.exists():
            n_nosidecar += 1
            continue
        analysis = json.loads(an_path.read_text())
        table = (analysis.get("interface") or {}).get("peptide_table")
        if not table:
            n_notable += 1
            continue
        was_null = all(r.get("hla_contact") is None for r in table)
        peptide_hla = (json.loads(sidecar_path.read_text()).get("peptide_hla")) or {}
        n_rows = struct_metrics.merge_peptide_table(table, peptide_hla)
        if n_rows == 0:
            n_skip += 1
            continue
        if was_null:
            fixed_null.append(tid)
        if not args.dry_run:
            an_path.write_text(json.dumps(analysis, indent=2), encoding="utf-8")
        n_written += 1

    tag = " (dry-run)" if args.dry_run else " written"
    print(f"\nBackfill peptide-HLA: {n_written} trajectories{tag}; "
          f"{len(fixed_null)} had an all-null HLA column and were repaired.")
    if n_nosidecar or n_notable or n_skip:
        print(f"  skipped: {n_nosidecar} no sidecar, {n_notable} no peptide_table, "
              f"{n_skip} no matching resids.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
