#!/usr/bin/env python3
"""Repair the run3 identity metadata in the served analysis.json files.

Two defects, both introduced by the `parallel3_preprocess` run3 analysis batch, in
which ANARCI was not available so md_analysis's IntelligentChainIdentifier fell back
to ordering the three non-peptide/non-b2m chains by length:

  A. chain roles: for the 6 complexes where the length order disagrees with ANARCI,
     `identity.complex_identity` names the wrong chains. 1fo0 / 1nam / 2ol3 have TCR
     alpha and beta transposed (both are ~112-118 aa single-domain constructs); the
     derived `tcr_identity.alpha_length` / `beta_length` are transposed with them.
     4prh / 5hhm / 7rk7 transposed the MHC heavy chain with TCR beta and are NOT
     touched here -- their contact stage consumed the bad mapping and they need the
     full recompute in ``fix_run3_chain_mapping.py``.

  B. germline identities: `tcr_identity.{alpha,beta}_{v,j}_identity` are all exactly
     0.0 in 236 run3 records. That is the empty-dict default of
     BiologicalIdentityAnnotator when the ANARCI germline call returns nothing -- a
     missing measurement written as a real zero. The TCR alpha/beta sequences are
     byte-identical across replicas, and run1/run2 agree on every one of the 236, so
     the measured values are taken from the sibling replica.

`tcr_identity.*_genotype_confidence` is deliberately NOT touched: it is hard-coded to
"high" for every record by ``extract_analysis.py::_sync_identity_tcr``, so an edit here
would be reverted by the next curate pass. If it should become honest, fix it there.

Edits are line-level so every other byte of the file (including float formatting) is
preserved. Nothing is written without --apply.

    python -m pipeline.fix_run3_identity_metadata            # dry run, prints manifest
    python -m pipeline.fix_run3_identity_metadata --apply    # back up, write, verify
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
from collections import OrderedDict
from pathlib import Path

DATA = Path("/home/xmy/work/data/immunotrace/web_data_1000")
BACKUP = Path("/data/work/Immuno-Dyn/traj_repair_stage/identity_backup_2026-09-15")
LOG = Path("/data/work/Immuno-Dyn/traj_repair_stage/identity_write_log_2026-09-15.tsv")

IDENT_KEYS = ("alpha_v_identity", "alpha_j_identity", "beta_v_identity", "beta_j_identity")
LEN_KEYS = ("alpha_length", "beta_length")
# complexes whose run3 contact stage consumed the bad mapping -- role repair excluded
CONTACT_BROKEN = {"4prh", "5hhm", "7rk7"}

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "MSE": "M",
}


class Abort(RuntimeError):
    """Any condition that would make a write a guess rather than a repair."""


def chain_sequences(pdb: Path) -> dict[str, str]:
    """{chain_id: one-letter sequence}; a residue is keyed by resSeq+iCode."""
    out: "OrderedDict[str, list]" = OrderedDict()
    last: dict[str, str] = {}
    with pdb.open() as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            cid, key, resn = line[21], line[22:27], line[17:20].strip()
            if last.get(cid) == key:
                continue
            last[cid] = key
            out.setdefault(cid, []).append(THREE2ONE.get(resn, "X"))
    return {k: "".join(v) for k, v in out.items()}


def tcr_chains(meta: dict) -> tuple[str, str]:
    chains = (meta.get("tcr") or {}).get("chains") or {}
    a = (chains.get("alpha") or {}).get("structural_chain")
    b = (chains.get("beta") or {}).get("structural_chain")
    if not a or not b:
        raise Abort("meta.tcr is missing a structural_chain")
    return a, b


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def set_scalar(text: str, key: str, value) -> str:
    """Replace the value of a unique ``"key": <scalar>`` line, keeping indentation."""
    pat = re.compile(r'^(\s*"%s"\s*:\s*)(.+?)(,?)$' % re.escape(key), re.M)
    hits = pat.findall(text)
    if len(hits) != 1:
        raise Abort(f"key {key!r} occurs {len(hits)} times, expected exactly 1")
    return pat.sub(lambda m: m.group(1) + json.dumps(value) + m.group(3), text, count=1)


def plan_one(cx: str) -> dict | None:
    """Return the edit plan for one complex's run3 record, or None if nothing to do."""
    rec = DATA / f"{cx}_run3"
    target = rec / "analysis" / "analysis.json"
    if not target.exists():
        return None
    doc = json.loads(target.read_text())
    ident = (doc.get("identity") or {}).get("tcr_identity") or {}
    if not all(ident.get(k) == 0.0 for k in IDENT_KEYS):
        return None

    meta3 = json.loads((rec / "meta.json").read_text())
    a3, b3 = tcr_chains(meta3)
    seq3 = chain_sequences(rec / "topology.pdb")

    # sibling whose TCR alpha+beta sequences are byte-identical to run3's
    sibling = None
    for r in ("run1", "run2"):
        sib = DATA / f"{cx}_{r}"
        sj = sib / "analysis" / "analysis.json"
        if not (sj.exists() and (sib / "topology.pdb").exists()):
            continue
        am, bm = tcr_chains(json.loads((sib / "meta.json").read_text()))
        seqs = chain_sequences(sib / "topology.pdb")
        if seqs.get(am) == seq3.get(a3) and seqs.get(bm) == seq3.get(b3):
            sid = (json.loads(sj.read_text()).get("identity") or {}).get("tcr_identity") or {}
            if all(sid.get(k) == 0.0 for k in IDENT_KEYS):
                continue                      # sibling has no measurement either
            sibling = (r, sid)
            break
    if sibling is None:
        raise Abort(f"{cx}: no sibling replica with matching TCR sequences and a measurement")
    sib_run, sid = sibling

    edits: dict[str, object] = {k: sid[k] for k in IDENT_KEYS}

    # --- role repair (defect A), excluding the contact-broken three ---
    ci = (doc.get("identity") or {}).get("complex_identity") or {}
    role_fixed = False
    if cx not in CONTACT_BROKEN and (ci.get("tcr_alpha_chain"), ci.get("tcr_beta_chain")) != (a3, b3):
        if {ci.get("tcr_alpha_chain"), ci.get("tcr_beta_chain")} != {a3, b3}:
            raise Abort(f"{cx}: complex_identity TCR chains {ci} are not a transposition of {a3}/{b3}")
        edits["tcr_alpha_chain"] = a3
        edits["tcr_beta_chain"] = b3
        edits["alpha_length"] = sid["alpha_length"]
        edits["beta_length"] = sid["beta_length"]
        role_fixed = True

    before = {k: (ident.get(k) if k in IDENT_KEYS + LEN_KEYS else ci.get(k)) for k in edits}
    return dict(cx=cx, path=target, sibling=sib_run, edits=edits, before=before,
                role_fixed=role_fixed, contact_broken=cx in CONTACT_BROKEN)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true", help="write (default: dry run)")
    args = ap.parse_args(argv)

    complexes = sorted({p.name[:-5] for p in DATA.glob("*_run1")})
    plans, problems = [], []
    for cx in complexes:
        try:
            p = plan_one(cx)
        except Abort as exc:
            problems.append(str(exc))
            continue
        if p:
            plans.append(p)

    if problems:
        print("ABORT -- unresolved records:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 2

    role = [p for p in plans if p["role_fixed"]]
    broken = [p for p in plans if p["contact_broken"]]
    print(f"records to edit: {len(plans)}")
    print(f"  germline identity backfill : {len(plans)}")
    print(f"  chain-role repair          : {len(role)}  ({', '.join(p['cx'] for p in role) or '-'})")
    print(f"  contact-broken, roles left alone: {len(broken)}  ({', '.join(p['cx'] for p in broken) or '-'})")
    print(f"  sibling used: run1 x{sum(p['sibling']=='run1' for p in plans)}, "
          f"run2 x{sum(p['sibling']=='run2' for p in plans)}")
    print()
    for p in role:
        print(f"  {p['cx']}_run3 role repair:")
        for k, v in p["edits"].items():
            if k in ("tcr_alpha_chain", "tcr_beta_chain") + LEN_KEYS:
                print(f"     {k}: {p['before'][k]!r} -> {v!r}")
    print()
    print("  germline sample (first 5):")
    for p in plans[:5]:
        print(f"     {p['cx']}_run3  " + "  ".join(f"{k}=0.0->{p['edits'][k]}" for k in IDENT_KEYS))

    if not args.apply:
        print("\ndry run -- nothing written")
        return 0

    if BACKUP.exists():
        print(f"ABORT -- backup dir already exists: {BACKUP}", file=sys.stderr)
        return 2
    BACKUP.mkdir(parents=True)
    rows = ["record\tfile\tmd5_before\tmd5_after\tfields"]
    for p in plans:
        src: Path = p["path"]
        dst = BACKUP / f"{p['cx']}_run3.analysis.json"
        shutil.copy2(src, dst)
        before_md5 = md5(src)
        text = src.read_text()
        for k, v in p["edits"].items():
            text = set_scalar(text, k, v)
        src.write_text(text)
        doc = json.loads(src.read_text())          # parses, or we restore
        ident = doc["identity"]["tcr_identity"]
        ci = doc["identity"]["complex_identity"]
        for k, v in p["edits"].items():
            got = ci.get(k) if k.startswith("tcr_") and k.endswith("_chain") else ident.get(k)
            if got != v:
                shutil.copy2(dst, src)
                print(f"ABORT -- {p['cx']}: {k} did not take ({got!r}); restored", file=sys.stderr)
                return 2
        rows.append(f"{p['cx']}_run3\t{src}\t{before_md5}\t{md5(src)}\t{','.join(p['edits'])}")
    LOG.write_text("\n".join(rows) + "\n")
    print(f"\nwrote {len(plans)} files")
    print(f"backup: {BACKUP}")
    print(f"log:    {LOG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
