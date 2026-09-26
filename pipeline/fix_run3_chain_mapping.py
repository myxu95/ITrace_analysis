"""One-off repair — force the canonical A–E chain partition for the three run3
trajectories whose ORIGINAL md_analysis analysis mislabeled both TCR chains as
alpha (beta n=0), inflating the ``bsa`` scalar / ``contact`` / ``interactions``
blocks and tripping the CDR ``single_chain_collapse`` reliability gate.

    Targets:  4prh_run3   5hhm_run3   7rk7_run3
    Defect:   see memory  live-data-defects-2026-08-30  (Defect 2)

Root cause (grounded, read-only): a chain-LABELING error, NOT a physical chain
departure. The trajectories are healthy (fnat / incident angle / bsa_decomposition
are all sibling-consistent). The three run3 SOURCE PDBs carry clean, distinct
chains A=MHC-alpha, B=b2m, C=peptide, D=TCR-alpha, E=TCR-beta (verified: chain
letters and residue counts match the served meta.json byte-for-byte). Forcing this
known partition into every stage reproduces the correct decomposition and bypasses
whatever the original ANARCI/heuristic pass mislabeled. The forced map is identical
to ``ChainIdentificationNode._fallback_chain_assignment()``.

    ============================================================================
    *** THIS SCRIPT NEVER TOUCHES THE LIVE web_data_1000 TREE. ***
    ============================================================================
    - ``recompute`` writes md_analysis stage output to a SEPARATE scratch root.
    - ``stage``     COPIES the three live analysis dirs into a scratch web_data
                    (read-only from live) so their preservation sidecars travel.
    - ``embed``     re-curates analysis.json INTO THE STAGING COPY ONLY, with a
                    hard guard that refuses to write to the live tree.
    - ``validate``  reads staging vs. live (both read-only) and checks the
                    acceptance / reverse-check gates.

    The live DB stays byte-for-byte frozen. The gated atomic swap + FL_P8 data
    rsync + dual-origin restart is a SEPARATE, later, human-authorized step and is
    deliberately NOT implemented here. Do not run any mode against the live path.

Usage (nothing runs unless a mode is given; none of these writes to live)::

    python -m pipeline.fix_run3_chain_mapping recompute            # -> scratch stage output
    python -m pipeline.fix_run3_chain_mapping stage                # copy 3 live dirs -> staging
    python -m pipeline.fix_run3_chain_mapping embed                # curate into staging copy
    python -m pipeline.fix_run3_chain_mapping validate             # gate check (staging vs live)
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from pipeline import config

# --------------------------------------------------------------------------- #
# Frozen scope + invariants
# --------------------------------------------------------------------------- #
TARGET_IDS = ["4prh_run3", "5hhm_run3", "7rk7_run3"]

# Canonical, verified partition for all three targets (== node fallback A-E).
CANONICAL_MAP = {
    "mhc_alpha": "A",
    "b2m": "B",
    "peptide": "C",
    "tcr_alpha": "D",
    "tcr_beta": "E",
}
EXPECTED_CHAINS = set(CANONICAL_MAP.values())  # {A,B,C,D,E}

# Stages curated by extract_analysis (inter_cluster is deliberately NOT ingested,
# see extract_analysis.build_analysis note; skipping it saves the heaviest stage).
DEFAULT_STAGES = "identity,contact,hbond,saltbridge,hydrophobic,pipi,cationpi,angle,bsa,rmsf"

# The tree we must never write to.  Identified by inode, not by path string:
# /home/xmy/work and /data/work are the SAME directory through a bind mount, and
# Path.resolve() does not follow bind mounts, so any guard that compares resolved
# path strings is bypassed simply by naming the other side of the mount.
LIVE_WEB_DATA = Path("/home/xmy/work/data/immunotrace/web_data_1000").resolve()

# Default scratch layout, all under a single self-describing parent, none served.
SCRATCH_ROOT = config.PROJECT_ROOT / "chainfix_run3"
DEFAULT_RECOMPUTE_OUT = SCRATCH_ROOT / "analysis_output"
DEFAULT_STAGING = SCRATCH_ROOT / "web_data_stage"


# --------------------------------------------------------------------------- #
# Guards
# --------------------------------------------------------------------------- #
def _ident(p: Path):
    """(st_dev, st_ino) for *p*, or None if it does not exist."""
    try:
        st = p.stat()
    except OSError:
        return None
    return (st.st_dev, st.st_ino)


def _assert_not_live(path: Path, what: str) -> Path:
    """Abort loudly if *path* is, is inside, or contains the live served tree.

    The comparison is on (st_dev, st_ino) rather than on the resolved path text.
    A string comparison is not sound here: the served tree is reachable both as
    /home/xmy/work/... and as /data/work/... through a bind mount, resolve()
    leaves both spellings intact, and so the string guard accepts the /data
    spelling of the very tree it exists to protect.
    """
    rp = path.resolve()
    live = _ident(LIVE_WEB_DATA)
    if live is None:
        sys.exit(
            f"REFUSED: cannot stat the live tree {LIVE_WEB_DATA}, so the guard "
            f"cannot prove {what} is outside it. Refusing to guess."
        )

    # rp is the live tree, or sits inside it.  Non-existent ancestors stat to
    # None and are skipped, so this works for a path about to be created.
    q = rp
    while True:
        if _ident(q) == live:
            sys.exit(
                f"REFUSED: {what} is (or is inside) the LIVE web_data_1000 tree.\n"
                f"    given: {path}\n    resolves to: {rp}\n"
                f"    matches: {LIVE_WEB_DATA} (dev,ino)={live}\n"
                f"This script must never write to the served DB."
            )
        if q == q.parent:
            break
        q = q.parent

    # rp contains the live tree (only meaningful if rp already exists).
    rid = _ident(rp)
    if rid is not None:
        a = LIVE_WEB_DATA
        while True:
            if _ident(a) == rid:
                sys.exit(
                    f"REFUSED: {what} CONTAINS the LIVE web_data_1000 tree.\n"
                    f"    given: {path}\n    resolves to: {rp}\n"
                    f"This script must never write to the served DB."
                )
            if a == a.parent:
                break
            a = a.parent
    return rp


def _assert_targets_only(ids: list[str]) -> None:
    stray = [i for i in ids if i not in TARGET_IDS]
    if stray:
        sys.exit(f"REFUSED: this repair only touches {TARGET_IDS}; got stray ids {stray}")


# --------------------------------------------------------------------------- #
# Forced chain identification (the A-forced core)
# --------------------------------------------------------------------------- #
def _install_forced_chain_identification() -> None:
    """Monkeypatch ChainIdentificationNode.execute so that, for our targets, it
    writes the verified canonical A-E partition instead of running ANARCI.

    Because run_analysis executes every stage as its own BatchExecutor pipeline
    (each re-instantiates ChainIdentificationNode and re-identifies chains),
    patching the CLASS method forces the mapping into EVERY stage uniformly --
    which is exactly what a per-stage, re-reading pipeline requires. Injecting
    into only the first stage would be useless.
    """
    from md_analysis.pipeline.nodes.topology.chain_identification_node import (
        ChainIdentificationNode,
    )

    orig_execute = ChainIdentificationNode.execute

    def forced_execute(self, context):
        pdb = str(getattr(context, "structure_pdb", "") or "")
        is_target = any(tid in pdb for tid in TARGET_IDS)
        if not is_target:
            # Scope leak guard: this driver only ever loads the 3 targets. If we
            # ever see anything else, fail rather than silently force a mapping.
            raise RuntimeError(
                f"forced_execute saw a non-target structure_pdb: {pdb!r}. "
                "This repair must only run over the 3 run3 targets."
            )

        # Self-validate the forced map against the ACTUAL source PDB at run time:
        # the canonical roles are only valid if the source really is clean A-E.
        chains = _pdb_chain_ids(Path(pdb))
        if chains != EXPECTED_CHAINS:
            raise RuntimeError(
                f"source PDB {pdb} has chains {sorted(chains)}, expected "
                f"{sorted(EXPECTED_CHAINS)}; refusing to force a canonical map "
                "onto a non-canonical layout."
            )

        context.metadata["chain_mapping"] = dict(CANONICAL_MAP)
        context.metadata["chain_identification_method"] = "forced_canonical_AE"
        return context

    ChainIdentificationNode.execute = forced_execute
    # keep a handle in case a caller wants to restore
    forced_execute._orig = orig_execute  # type: ignore[attr-defined]


def _pdb_chain_ids(pdb: Path) -> set[str]:
    ids: set[str] = set()
    with open(pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                ids.add(line[21])
    return ids


# --------------------------------------------------------------------------- #
# Mode: recompute
# --------------------------------------------------------------------------- #
def mode_recompute(out_root: Path, stages: list[str], ids: list[str]) -> int:
    """Run the forced md_analysis stages for the 3 targets into a scratch root.

    Runs SERIALLY, in THIS process (no ProcessPoolExecutor), so the class-level
    monkeypatch is guaranteed live for every stage. Only ~3 systems x ~10 stages.
    """
    _assert_targets_only(ids)
    out_root = _assert_not_live(out_root, "recompute --out")
    out_root.mkdir(parents=True, exist_ok=True)

    _install_forced_chain_identification()

    # Import AFTER the patch is installed. run_trajectory builds pipelines lazily
    # via STAGE_FACTORIES, so the patched class method is what every stage uses.
    from pipeline.run_analysis import build_task_list, run_trajectory

    tasks = build_task_list(ids, None)
    got = {t["task_id"] for t in tasks}
    missing = [i for i in ids if i not in got]
    if missing:
        sys.exit(f"source trajectories not found under {config.SOURCE_ROOT}: {missing}")

    print(f"[recompute] forced canonical A-E map = {CANONICAL_MAP}")
    print(f"[recompute] source root = {config.SOURCE_ROOT}")
    print(f"[recompute] output root = {out_root}")
    print(f"[recompute] stages      = {','.join(stages)}")
    summary = {}
    for task in tasks:
        tid = task["task_id"]
        t0 = time.time()
        print(f"\n[recompute] === {tid} ===")
        res = run_trajectory(task, stages, str(out_root))
        summary[tid] = res
        ok = sum(1 for s in res["stages"].values() if s.get("status") == "success")
        print(f"[recompute] {tid}: {ok}/{len(res['stages'])} stages ok "
              f"({time.time() - t0:.0f}s)")
        for st, info in res["stages"].items():
            if info.get("status") != "success":
                print(f"    ! {st}: {info.get('errors')}")

    (out_root / "recompute_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(f"\n[recompute] done -> {out_root} (live web_data_1000 untouched)")
    return 0


# --------------------------------------------------------------------------- #
# Mode: stage  (copy the 3 live dirs into a scratch web_data so sidecars travel)
# --------------------------------------------------------------------------- #
def mode_stage(staging: Path, ids: list[str]) -> int:
    _assert_targets_only(ids)
    staging = _assert_not_live(staging, "stage --staging")
    staging.mkdir(parents=True, exist_ok=True)

    live = LIVE_WEB_DATA
    for tid in ids:
        src = live / tid
        if not src.is_dir():
            sys.exit(f"live source dir missing: {src}")
        dst = staging / tid
        if dst.exists():
            shutil.rmtree(dst)
        # copy2 preserves the preservation sidecars (rmsd_regions/essential_dynamics/
        # struct_metrics/concerted_motion/dpca/coupled) that embed reads back.
        shutil.copytree(src, dst, copy_function=shutil.copy2)
        print(f"[stage] copied {src}  ->  {dst}")
    print(f"[stage] staging ready at {staging} (read-only pull from live)")
    return 0


# --------------------------------------------------------------------------- #
# Mode: embed  (re-curate analysis.json into the STAGING copy, never live)
# --------------------------------------------------------------------------- #
def mode_embed(recompute_out: Path, staging: Path, ids: list[str]) -> int:
    _assert_targets_only(ids)
    recompute_out = recompute_out.resolve()
    staging = _assert_not_live(staging, "embed --staging")

    if not recompute_out.is_dir():
        sys.exit(f"recompute output not found: {recompute_out} (run `recompute` first)")

    # Preflight: every target must be staged AND carry every preservation sidecar,
    # otherwise extract_analysis would DROP fnat/rmsip/ed/dpca/concerted from the
    # rebuilt analysis.json. This is the reverse-check safety net for embed.
    need_sidecars = [
        "rmsd_regions.json", "essential_dynamics.json", "struct_metrics.json",
        "concerted_motion.json", "peptide_dihedrals.json", "tcr_cdr3_dpca.json",
        "coupled_states.json",
    ]
    for tid in ids:
        adir = staging / tid / "analysis"
        if not adir.is_dir():
            sys.exit(f"not staged: {adir} (run `stage` first)")
        missing = [s for s in need_sidecars if not (adir / s).exists()]
        if missing:
            sys.exit(f"{tid}: staging is missing preservation sidecars {missing}; "
                     "embed would drop those blocks -- abort.")

    # extract_analysis writes to config.WEB_DATA (read at import from IMMUNO_WEB_DATA).
    # Run it as a subprocess with IMMUNO_WEB_DATA pinned to the STAGING copy so the
    # live tree is impossible to hit. Guard the env value too.
    _assert_not_live(staging, "embed IMMUNO_WEB_DATA")
    env = dict(os.environ)
    env["IMMUNO_WEB_DATA"] = str(staging)

    cmd = [
        sys.executable, "-m", "pipeline.extract_analysis",
        "--source", str(recompute_out),
        "--ids", ",".join(ids),
    ]
    print(f"[embed] IMMUNO_WEB_DATA = {staging}  (NOT live)")
    print(f"[embed] $ {' '.join(cmd)}")
    proc = subprocess.run(cmd, env=env, cwd=str(config.PROJECT_ROOT))
    if proc.returncode != 0:
        sys.exit(f"extract_analysis failed (rc={proc.returncode})")
    print(f"[embed] curated into staging {staging} (live web_data_1000 untouched)")
    print("[embed] next: `validate` before ANY talk of ingest.")
    return 0


# --------------------------------------------------------------------------- #
# Mode: validate  (acceptance gates + reverse-check; all reads read-only)
# --------------------------------------------------------------------------- #
def _load(path: Path):
    return json.loads(path.read_text()) if path.exists() else None


def _g(a, *ks):
    for k in ks:
        a = (a or {}).get(k) if isinstance(a, dict) else None
    return a


def _by_chain_beta_n(a):
    bc = _g(a, "tcr_cdr", "by_chain")
    if isinstance(bc, dict):
        b = bc.get("beta") or bc.get("B") or {}
        return b.get("n") if isinstance(b, dict) else None
    if isinstance(bc, list):
        for row in bc:
            if str(row.get("chain", "")).lower() in ("beta", "b", "tcr_beta"):
                return row.get("n")
    return None


def _fr_to_mhc_n(a):
    """Framework->MHC pair count. Real nesting: tcr_cdr.by_region.fr.mhc.n.

    This is the block that ballooned under the mislabel (4prh 18->237, 5hhm ->223,
    7rk7 ->231): TCR-internal alpha-beta contacts miscounted as framework-MHC once
    the beta chain was folded into the MHC-side selection.
    """
    return _g(a, "tcr_cdr", "by_region", "fr", "mhc", "n")


def _alpha_contrib(a):
    return _g(a, "tcr_cdr", "alpha_contribution")


def _bsa_scalar(a):
    return _g(a, "bsa", "buried_surface_area", "mean")


def _recon_consistent(a):
    """bsa_decomposition's own self-consistency verdict (one-sided BSA scalar vs
    the two-sided decomposition). True in 524/735; False in EXACTLY the 3 targets.
    """
    return _g(a, "bsa_decomposition", "reconciliation", "consistent")


def _vbeta(a):
    """Vbeta buried-surface contribution (health precondition; key is Greek 'Vβ').

    Normal/unchanged for the 3 targets (386.6/467.7/415.7 vs sibling ~424.6). It
    would collapse toward ~0 if the beta chain had truly DETACHED, so a normal Vbeta
    confirms the defect is a labeling artifact, not a physical departure.
    """
    br = _g(a, "bsa_decomposition", "by_region")
    if isinstance(br, dict):
        return br.get("Vβ", br.get("Vbeta"))
    return None


def _ann_non_cdr(a):
    return _g(a, "contact", "annotation", "n_non_cdr_contacts")


def _ann_cdr3(a):
    return _g(a, "contact", "annotation", "n_cdr3_contacts")


def _inter_n(a, fam):
    return _g(a, "interactions", fam, "n_total_pairs")


def _incident(a):
    return _g(a, "angle", "incident_deg", "mean")


def _rmsf_mean(a):
    return _g(a, "rmsf", "mean_rmsf_angstrom")


def _tcr_rmsf_mean(a):
    """TCR-side mean RMSF. The sharpest reverse-check for an injection error: if the
    forced chain_mapping landed on the wrong chains, the TCR-side fluctuation is what
    moves first. Served run3 values (1.852/2.633/2.453) are already correct (RMSF was
    never contaminated), so a correct recompute must reproduce them.
    """
    return _g(a, "rmsf", "tcr_mean_rmsf_angstrom")


def _fnat(a):
    return _g(a, "fnat", "fnat_mean")


def _rmsip(a):
    return _g(a, "essential_dynamics", "subspace_rmsip")


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


def mode_validate(staging: Path, ids: list[str]) -> int:
    _assert_targets_only(ids)
    staging = staging.resolve()
    live = LIVE_WEB_DATA
    all_pass = True

    for tid in ids:
        base = tid.rsplit("_run", 1)[0]
        cand = _load(staging / tid / "analysis" / "analysis.json")
        served = _load(live / tid / "analysis" / "analysis.json")   # read-only reference
        sibs = [
            _load(live / f"{base}_run{r}" / "analysis" / "analysis.json")
            for r in (1, 2)
        ]
        sibs = [s for s in sibs if s]

        print("\n" + "=" * 72)
        print(f"[validate] {tid}   (candidate=staging, references=live, read-only)")
        if cand is None:
            print("  ! candidate analysis.json missing -- run recompute+stage+embed")
            all_pass = False
            continue

        def gate(name, ok, detail):
            nonlocal all_pass
            all_pass = all_pass and ok
            print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")

        # ---- health precondition: the trajectory was never physically broken ----
        # Vbeta buried-surface stays normal (~380-470); it would collapse toward 0 if
        # the beta chain had truly detached. Confirms this is a labeling artifact.
        vb = _vbeta(cand)
        sib_vb = [x for x in (_vbeta(s) for s in sibs) if x is not None]
        gate("Vbeta buried-surface normal (>=100, not detached)",
             vb is not None and vb >= 100,
             f"cand Vβ={_fmt(vb)}  siblings={[_fmt(x) for x in sib_vb]}")

        # ---- acceptance gates (contaminated blocks must recover) ----
        beta_n = _by_chain_beta_n(cand)
        gate("beta chain recovers (n>0)", bool(beta_n and beta_n > 0),
             f"by_chain.beta.n = {_fmt(beta_n)}  (was 0 when contaminated)")

        ac = _alpha_contrib(cand)
        sib_ac = [x for x in (_alpha_contrib(s) for s in sibs) if x is not None]
        gate("alpha_contribution ~0.5", ac is not None and 0.30 <= ac <= 0.70,
             f"cand={_fmt(ac)}  siblings={[_fmt(x) for x in sib_ac]}")

        fr = _fr_to_mhc_n(cand)
        gate("by_region.fr.mhc.n back to two-digit (<100)",
             fr is not None and fr < 100,
             f"cand fr.mhc.n = {_fmt(fr)}  (was 237/223/231 when contaminated)")

        # Sharpest annotation-level gates: the mislabel dumped non-CDR framework
        # contacts (545/543/632, DB-normal <50) and starved the real CDR3 count
        # (17/8/17, siblings ~40-70).
        n_noncdr = _ann_non_cdr(cand)
        gate("contact.annotation.n_non_cdr_contacts < 50",
             n_noncdr is not None and n_noncdr < 50,
             f"cand={_fmt(n_noncdr)}  (was 545/543/632 when contaminated)")

        n_cdr3 = _ann_cdr3(cand)
        gate("contact.annotation.n_cdr3_contacts in ~[40,70]",
             n_cdr3 is not None and 40 <= n_cdr3 <= 70,
             f"cand={_fmt(n_cdr3)}  (was 17/8/17 when contaminated)")

        npair = _g(cand, "contact", "n_pairs")
        gate("contact.n_pairs in ~[90,160]", npair is not None and 90 <= npair <= 160,
             f"cand={_fmt(npair)}  (siblings 98-134)")

        nst = _g(cand, "contact", "n_stable_pairs")
        gate("contact.n_stable_pairs in ~[40,60]", nst is not None and 40 <= nst <= 60,
             f"cand={_fmt(nst)}  (siblings 47-54)")

        hb = _inter_n(cand, "hbond")
        gate("hbond n_total in ~[60,110]", hb is not None and 60 <= hb <= 110,
             f"cand={_fmt(hb)}  (target 70-100; was 420-472)")

        hp = _inter_n(cand, "hydrophobic")
        gate("hydrophobic n_total in ~[4,20]", hp is not None and 4 <= hp <= 20,
             f"cand={_fmt(hp)}  (target 6-14; was 80-104)")

        # Primary BSA gate: the pipeline's own self-consistency verdict must flip
        # back to True (it is False in exactly the 3 targets). Stronger than a
        # hand-rolled ratio band and asserted `is True` (block is absent/None in
        # ~211 siblings incl. 4prh_run2 due to version drift).
        gate("bsa_decomposition.reconciliation.consistent is True",
             _recon_consistent(cand) is True,
             f"cand consistent = {_fmt(_recon_consistent(cand))}  (was False)")

        # ---- reverse-check: these must NOT move vs the currently-served run3 ----
        def unchanged(name, fn, tol):
            c, s = fn(cand), fn(served)
            if c is None or s is None:
                gate(f"[reverse] {name} present in both", False,
                     f"cand={_fmt(c)} served={_fmt(s)}")
                return
            d = abs(c - s)
            gate(f"[reverse] {name} unchanged (|d|<={tol})", d <= tol,
                 f"cand={_fmt(c)} served={_fmt(s)} |d|={_fmt(d)}")

        unchanged("fnat", _fnat, 1e-6)             # from preserved sidecar -> exact
        unchanged("rmsip", _rmsip, 1e-6)           # from preserved sidecar -> exact
        unchanged("incident angle", _incident, 1.0)   # recomputed -> allow 1 deg
        unchanged("rmsf mean", _rmsf_mean, 0.02)      # recomputed -> allow 0.02 A
        unchanged("tcr rmsf mean", _tcr_rmsf_mean, 0.05)  # canary: moves first if mapping mis-injected

    print("\n" + "=" * 72)
    print(f"[validate] OVERALL: {'ALL GATES PASS' if all_pass else 'GATES FAILED'}")
    if all_pass:
        print("[validate] Candidate is sibling-consistent AND reverse-check clean.")
        print("[validate] Ingest into live web_data_1000 remains BLOCKED pending")
        print("[validate] explicit sign-off from the manuscript authors AND the")
        print("[validate] repo owner. This script does not perform ingest.")
    return 0 if all_pass else 2


# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Force canonical A-E chain mapping for 3 run3 targets "
                    "(NEVER writes live web_data_1000).")
    sub = p.add_subparsers(dest="mode", required=True)

    pr = sub.add_parser("recompute", help="forced md_analysis stages -> scratch root")
    pr.add_argument("--out", type=Path, default=DEFAULT_RECOMPUTE_OUT)
    pr.add_argument("--stages", default=DEFAULT_STAGES)
    pr.add_argument("--ids", default=",".join(TARGET_IDS))

    ps = sub.add_parser("stage", help="copy 3 live analysis dirs -> staging web_data")
    ps.add_argument("--staging", type=Path, default=DEFAULT_STAGING)
    ps.add_argument("--ids", default=",".join(TARGET_IDS))

    pe = sub.add_parser("embed", help="curate analysis.json INTO staging (not live)")
    pe.add_argument("--out", type=Path, default=DEFAULT_RECOMPUTE_OUT)
    pe.add_argument("--staging", type=Path, default=DEFAULT_STAGING)
    pe.add_argument("--ids", default=",".join(TARGET_IDS))

    pv = sub.add_parser("validate", help="acceptance + reverse-check gates")
    pv.add_argument("--staging", type=Path, default=DEFAULT_STAGING)
    pv.add_argument("--ids", default=",".join(TARGET_IDS))

    args = p.parse_args(argv)
    ids = [s.strip() for s in args.ids.split(",") if s.strip()]

    if args.mode == "recompute":
        stages = [s.strip() for s in args.stages.split(",") if s.strip()]
        return mode_recompute(args.out, stages, ids)
    if args.mode == "stage":
        return mode_stage(args.staging, ids)
    if args.mode == "embed":
        return mode_embed(args.out, args.staging, ids)
    if args.mode == "validate":
        return mode_validate(args.staging, ids)
    p.error(f"unknown mode {args.mode}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
