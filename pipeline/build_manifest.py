"""Aggregate every per-trajectory meta.json into a single manifest.json.

The manifest is the compact array the list page loads in one request: one row
per trajectory with just the fields needed for the searchable/sortable table.

Usage:
    python -m pipeline.build_manifest
"""
from __future__ import annotations

import json
import re
import sys

from . import config
from . import antigen
from . import consensus

_PRECISE_HLA = re.compile(r"^HLA-[A-Z]+\*\d+:\d+")

# Structures where RCSB is only locus-coarse AND ImmunoScope's typer mis-scored the
# MHC chain (it matched short truncated `N` null alleles at ~55-67% identity). The
# correct 2-field allele was resolved by directly aligning the MHC heavy chain to the
# IMGT/HLA protein database (identity in comment), and cross-checks against the bound
# peptide. Highest priority so these known-bad cases are always corrected.
_MANUAL_HLA = {
    "4prh": "HLA-B*35:08",   # 100% identity to B*35:08
    "5hhm": "HLA-A*02:01",   # 100%; peptide GILGLVFTL (flu M1) is the canonical A*02:01 epitope
    "7rk7": "HLA-A*02:01",   # 97%; A*02:01-restricted peptide YMDGTMSQV
}

# Mouse (H-2Db) structures that use human beta-2-microglobulin, so the b2m-based
# host heuristic mislabels them as human. The MHC heavy chain is the definitive
# species marker (RCSB titles say "H2-Db"; the chain aligns ~75% to human HLA,
# i.e. homologous mouse MHC, not an HLA allele).
_MANUAL_HOST = {pdb: "Mus musculus" for pdb in ("7jwi", "7n4k", "7n5c", "7n5p", "7na5")}


def _round(x, n=3):
    return round(x, n) if isinstance(x, (int, float)) else None


def _host_species(meta: dict):
    """Host species of the MHC, used to split human from non-human complexes.
    Beta-2 microglobulin is the most reliable host marker (always the host's),
    then the MHC heavy chain, then the TCR."""
    rcsb = meta.get("rcsb") or {}
    by_role = {}
    for e in rcsb.get("entities") or []:
        role = e.get("role")
        org = e.get("organism") or (e.get("organisms") or [None])[0]
        if role and org and role not in by_role:
            by_role[role] = org
    for role in ("b2m", "mhc", "tcr_alpha", "tcr_beta", "tcr"):
        if by_role.get(role):
            return by_role[role]
    return None


def row_from_meta(meta: dict) -> dict:
    q = meta.get("quality") or {}
    rcsb = meta.get("rcsb") or {}
    tcr = meta.get("tcr") or {}
    tcr_chains = tcr.get("chains") or {}
    desc = meta.get("interface_descriptors") or {}
    ag = meta.get("antigen") or antigen.extract(meta)
    return {
        "traj_id": meta["traj_id"],
        "pdb_id": meta["pdb_id"],
        "title": rcsb.get("title"),
        "peptide_seq": meta.get("peptide_seq"),
        "peptide_length": meta.get("peptide_length"),
        "antigen_name": ag.get("name"),
        "antigen_organism": ag.get("organism"),
        "antigen_category": ag.get("category"),
        "hla_allele": rcsb.get("hla_allele"),
        "organisms": rcsb.get("organisms") or [],
        "tcr_genes": rcsb.get("tcr_genes") or [],
        "tcr_type": tcr.get("type"),
        "trav": (tcr_chains.get("alpha") or {}).get("v_gene"),
        "trbv": (tcr_chains.get("beta") or {}).get("v_gene"),
        "peptide_recognition_ratio": desc.get("peptide_recognition_ratio"),
        "hla_restriction_ratio": desc.get("hla_restriction_ratio"),
        "alpha_contribution": desc.get("alpha_contribution"),
        "cdr_decomposition_reliable": desc.get("cdr_decomposition_reliable"),
        "resolution": rcsb.get("resolution"),
        "release_year": rcsb.get("release_year"),
        "citation": rcsb.get("citation"),
        "n_frames_full": meta.get("n_frames_full"),
        "n_frames_web": meta.get("n_frames_web"),
        "n_frames_view": meta.get("n_frames_view"),
        "duration_ns": meta.get("duration_ns"),
        "n_atoms": meta.get("n_atoms"),
        "n_chains": meta.get("n_chains"),
        "has_rmsd_png": meta.get("has_rmsd_png", False),
    }


def _two_field(allele: str):
    """'HLA-A*02:01:03' -> ('HLA-A', 2, 1); None if not an allele string."""
    m = re.match(r"^(HLA-[A-Z]+)\*(\d+):(\d+)", allele or "")
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None


def _struct_hla(traj_id: str):
    """Sequence-typed HLA allele from the MHC heavy chain (ImmunoScope's IMGT/HLA
    match), resolved to the 2-field (protein) level — the deeper fields are
    synonymous/intronic and not determinable from the structure.

    Only used when the match is high-confidence (so non-human MHC, e.g. mouse H-2,
    which only match HLA weakly, are not mislabelled). Among equally-good tied
    candidates the canonical (lowest-numbered) 2-field allele is chosen, because
    the structure cannot distinguish synonymous sub-alleles and the common
    reference allele is the low-numbered one (HLA-A*02:01, not HLA-A*02:665)."""
    path = config.WEB_DATA / traj_id / "analysis" / "analysis.json"
    if not path.exists():
        return None
    try:
        with path.open() as fh:
            a = json.load(fh)
    except (ValueError, OSError):
        return None
    h = (a.get("identity") or {}).get("hla_identity") or {}
    if (h.get("confidence") or "").lower() != "high":
        return None
    if (h.get("identity") or 0) < 0.98 or (h.get("coverage") or 0) < 0.8:
        return None
    cands = h.get("top_candidates") or [{"allele": h.get("best_candidate_allele"),
                                         "identity": h.get("identity") or 1.0}]
    top_id = max((c.get("identity", 0) for c in cands), default=0)
    twos = [t for c in cands if c.get("identity", 0) >= top_id - 1e-9
            for t in [_two_field(c.get("allele"))] if t]
    if not twos:
        return None
    locus, grp, prot = min(twos)
    return f"{locus}*{grp:02d}:{prot:02d}"


# Fields that describe the complex (shared across replicas of the same PDB) vs
# fields that vary per replica (per-run trajectory).
_SHARED_FIELDS = ("title", "peptide_seq", "peptide_length", "hla_allele", "hla_predicted",
                  "hla_source", "organisms", "tcr_genes", "tcr_type", "trav", "trbv",
                  "resolution", "release_year", "citation", "host_species", "is_human",
                  "antigen_name", "antigen_organism", "antigen_category")
# NOTE: rmsd_mean/std/max/variation are deliberately absent -- RMSD-drift
# descriptors are computed internally but never served (do-not-surface-QC rule).
_REPLICA_FIELDS = ("traj_id", "duration_ns",
                   "n_frames_web", "n_frames_view", "n_frames_full",
                   "n_atoms", "n_chains", "has_rmsd_png",
                   "peptide_recognition_ratio", "hla_restriction_ratio",
                   "alpha_contribution", "cdr_decomposition_reliable")


def _run_label(traj_id: str, pdb_id: str) -> str:
    """`1g6r_run3` (pdb `1g6r`) -> `run3`."""
    return traj_id[len(pdb_id):].lstrip("_") or "run"


def _primary_replica(reps: list[dict]) -> dict:
    """Representative replica for the complex row: prefer run3, else the longest run."""
    run3 = [r for r in reps if r["traj_id"].endswith("_run3")]
    if run3:
        return run3[0]
    return max(reps, key=lambda r: (r.get("duration_ns") or 0))


def build_complexes(rows: list[dict]) -> list[dict]:
    """Group trajectory rows by pdb_id into one entry per complex with a replica list.

    Replica identity stays encoded in each traj_id (no persisted run field); the
    complex row copies shared metadata from a primary replica (run3 preferred) and
    lists every replica with its per-run varying fields.
    """
    by_pdb: dict[str, list[dict]] = {}
    for r in rows:
        by_pdb.setdefault(r["pdb_id"], []).append(r)
    complexes = []
    for pdb_id, reps in sorted(by_pdb.items()):
        reps = sorted(reps, key=lambda r: r["traj_id"])
        prim = _primary_replica(reps)
        entry = {"pdb_id": pdb_id, "primary_traj_id": prim["traj_id"], "n_replicas": len(reps)}
        for k in _SHARED_FIELDS:
            entry[k] = prim.get(k)
        # Replica-invariant representatives surfaced for the browse table / sorting.
        for k in ("duration_ns",):
            entry[k] = prim.get(k)
        # Quantitative interface descriptors: the complex-level value is the mean over
        # replicas, published together with the MEASURED inter-replica spread. This
        # replaced three tercile-binned labels and their boolean *_agree flags — a bin
        # cut out of this dataset's own distribution described a complex's rank inside
        # our collection, not the molecule (see pipeline.interface_descriptors).
        for key in ("peptide_recognition_ratio", "hla_restriction_ratio", "alpha_contribution"):
            stat = consensus.numeric_consensus([r.get(key) for r in reps])
            entry[key] = stat["mean"] if stat else None
            entry[f"{key}_spread"] = stat
        entry["cdr_decomposition_reliable"] = all(
            r.get("cdr_decomposition_reliable") is not False for r in reps)
        entry["replicas"] = [
            {"run_label": _run_label(r["traj_id"], pdb_id),
             **{k: r.get(k) for k in _REPLICA_FIELDS}}
            for r in reps
        ]
        complexes.append(entry)
    return complexes


def main(argv: list[str]) -> int:
    metas = sorted(config.WEB_DATA.glob("*/" + config.OUT_META))
    rows = []
    for mp in metas:
        with mp.open() as fh:
            meta = json.load(fh)
        row = row_from_meta(meta)
        host = _MANUAL_HOST.get(row["pdb_id"]) or _host_species(meta)
        row["host_species"] = host
        row["is_human"] = (host == "Homo sapiens")
        # Resolve the HLA allele to the most precise reliable value:
        #   precise crystallographic (HLA-A*02:01...)  -> keep (source rcsb)
        #   coarse/missing crystallographic + high-conf sequence typing -> upgrade
        #   coarse crystallographic, no confident typing -> keep coarse
        rcsb_hla = row["hla_allele"]
        struct_hla = _struct_hla(row["traj_id"])
        manual_hla = _MANUAL_HLA.get(row["pdb_id"])
        if manual_hla:
            row["hla_allele"] = manual_hla
            row["hla_predicted"] = True
            row["hla_source"] = "structure"
        elif rcsb_hla and _PRECISE_HLA.match(rcsb_hla):
            row["hla_source"] = "rcsb"
        elif struct_hla:
            row["hla_allele"] = struct_hla
            row["hla_predicted"] = True
            row["hla_source"] = "structure"
        elif rcsb_hla:
            row["hla_source"] = "rcsb-coarse"
        else:
            row["hla_source"] = None
        # persist the resolved allele into meta.json so the detail page (which reads
        # meta, not the manifest) shows the precise genotype too (idempotent)
        if meta.get("hla_resolved") != row["hla_allele"] or meta.get("hla_source") != row["hla_source"]:
            meta["hla_resolved"] = row["hla_allele"]
            meta["hla_source"] = row["hla_source"]
            mp.write_text(json.dumps(meta, indent=2))
        rows.append(row)
    rows.sort(key=lambda r: r["traj_id"])
    complexes = build_complexes(rows)

    # Report the stride the DELIVERED files actually carry, not config.STRIDE.
    # config.STRIDE only governs how extract.py reads the source; the published
    # traj.xtc has been full-resolution (1001 frames, 1:1 with the source) since
    # the 2026-08-25 dedup repair, so the old constant would have written a
    # stale "5" into a manifest reviewers download. Null if the library is mixed.
    obs = {r.get("n_frames_full") // r["n_frames_web"]
           for r in rows if r.get("n_frames_full") and r.get("n_frames_web")}
    manifest = {
        "count": len(rows),
        "n_complexes": len(complexes),
        "stride": obs.pop() if len(obs) == 1 else None,
        "trajectories": rows,
        "complexes": complexes,
    }
    with config.MANIFEST.open("w") as fh:
        json.dump(manifest, fh, indent=1)
    print(f"Wrote {config.MANIFEST} with {len(rows)} trajectories / {len(complexes)} complexes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
