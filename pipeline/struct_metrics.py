"""Trajectory-based structural metrics (peptide-HLA contact, SASA, bulge, Rg).

The contact-CSV layer only sees pHLA-vs-TCR pairs, so peptide-in-groove (anchor)
engagement, peptide solvent exposure, bulge, and whole-complex Rg / SASA were all
missing. Without them the peptide table reads as if low-TCR-contact anchor
positions were "unimportant", when in fact they are buried HLA anchors.

This module is **purely trajectory-derived** (``topology.pdb`` + ``traj.xtc``),
like ``concerted_motion``: it globs ``WEB_DATA`` and takes no ``--source``. It
deliberately does NOT rely on the md_analysis rmsf CSV for chain / groove
identification — those ``component`` / ``mhc_subregion`` labels key on chain ids
that are not consistent across the PDB and mis-resolved the MHC heavy chain for
several complexes (e.g. 7rk7 / 5hhm / 4prh got an impossible 18-27 A bulge and
zero N-terminal peptide-HLA contact). Instead the MHC heavy chain and groove
platform are resolved by **role / geometry** (the same robust approach as
``docking_angle``): MHC = largest chain that is not peptide / TCRalpha / TCRbeta,
and the groove plane is the SVD of its alpha1/alpha2 platform Calpha.

Per peptide residue: max HLA-contact occupancy + per-residue SASA (low SASA +
low TCR contact = a buried anchor). Per system: peptide bulge height (max peptide
Calpha above the MHC groove plane), radius of gyration, and total SASA.

Writes a sidecar ``analysis/struct_metrics.json`` AND patches ``analysis.json``
directly (geometry block + peptide_table hla_contact / sasa_nm2 / anchor), so it
does not depend on a later ``extract_analysis`` pass (which needs the now-absent
md_analysis source root).

    IMMUNO_WEB_DATA=.../immuno-dyn python -m pipeline.struct_metrics [--ids a,b]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import mdtraj as md

from . import config
from .docking_angle import _chain_ca, _mhc_chain, _MANUAL_PEPTIDE_CHAIN, MHC_PLATFORM

CONTACT_CUTOFF_NM = 0.45     # 4.5 A closest-heavy peptide-HLA contact
SASA_STRIDE = 10             # frames; SASA is the expensive term
BULGE_STRIDE = 5             # frames; bulge plane fit is cheaper, sample every 5th
ANCHOR_TCR_MAX = 0.20        # an anchor is ignored by the TCR (< this occupancy)
ANCHOR_SASA_MAX_NM2 = 0.20   # an anchor is buried (< this mean per-residue SASA)
SIDECAR = "struct_metrics.json"


def _chain_residues(top, letter):
    return [r.index for ch in top.chains if ch.chain_id == letter for r in ch.residues]


# Canonical class-I alpha1 / alpha2 groove-helix spans (heavy-chain resSeq). The
# bulge is the peptide rise above the GROOVE-HELIX RIM, not the whole platform
# mid-plane (which sits ~10 A below the rim and inflates every bulge), so the
# plane is fit on these two helices.
HELIX1 = (50, 85)
HELIX2 = (138, 178)
HELIX1_WIDE = (40, 92)
HELIX2_WIDE = (128, 185)


def _ca_in_spans(top, chain_id, spans):
    out = []
    for ch in top.chains:
        if ch.chain_id != chain_id:
            continue
        for r in ch.residues:
            if any(lo <= r.resSeq <= hi for lo, hi in spans):
                ca = next((a.index for a in r.atoms if a.name == "CA"), None)
                if ca is not None:
                    out.append(ca)
    return out


def _groove_helix_ca(top, mhc_chain, platform_ca):
    """Calpha of the alpha1/alpha2 groove helices (the rim flanking the peptide).

    Canonical class-I resSeq spans first; widen, then fall back to the platform
    only if numbering is non-canonical (never needed for the standard class-I
    set, where this reliably finds ~77 helix Calpha)."""
    helix = _ca_in_spans(top, mhc_chain, (HELIX1, HELIX2))
    if len(helix) < 20:
        helix = _ca_in_spans(top, mhc_chain, (HELIX1_WIDE, HELIX2_WIDE))
    return helix if len(helix) >= 4 else platform_ca


def bsa_decomposition(sub, sasa_complex, pep_chain, mhc_chain, a_ch, b_ch) -> dict | None:
    """Interface BSA footprint split by partner / TCR chain.

    Per residue ΔSASA = SASA(in isolated partner) − SASA(in complex); summed by
    region gives each component's contribution to the buried interface. pMHC side
    (peptide / MHC / β2m) and TCR side (Vα / Vβ). nm² → Å². Reuses the strided
    complex SASA already computed; two extra isolated-partner SASA passes.
    """
    top = sub.topology
    tcr_chains = {c for c in (a_ch, b_ch) if c}
    if not tcr_chains:
        return None
    pmhc_res, tcr_res = [], []
    for r in top.residues:
        (tcr_res if r.chain.chain_id in tcr_chains else pmhc_res).append(r.index)
    if not pmhc_res or not tcr_res:
        return None
    pmhc_res.sort()
    tcr_res.sort()
    pmhc_atoms = np.array([a.index for ri in pmhc_res for a in top.residue(ri).atoms], dtype=np.int32)
    tcr_atoms = np.array([a.index for ri in tcr_res for a in top.residue(ri).atoms], dtype=np.int32)
    iso_pmhc = md.shrake_rupley(sub.atom_slice(pmhc_atoms), mode="residue").mean(axis=0)
    iso_tcr = md.shrake_rupley(sub.atom_slice(tcr_atoms), mode="residue").mean(axis=0)
    comp = sasa_complex.mean(axis=0)
    by: dict[str, float] = {}
    for i, ri in enumerate(pmhc_res):
        d = max(0.0, float(iso_pmhc[i] - comp[ri])) * 100.0
        ch = top.residue(ri).chain.chain_id
        key = "peptide" if ch == pep_chain else ("MHC" if ch == mhc_chain else "β2m")
        by[key] = by.get(key, 0.0) + d
    for i, ri in enumerate(tcr_res):
        d = max(0.0, float(iso_tcr[i] - comp[ri])) * 100.0
        ch = top.residue(ri).chain.chain_id
        key = "Vα" if ch == a_ch else ("Vβ" if ch == b_ch else "TCR-other")
        by[key] = by.get(key, 0.0) + d
    by = {k: round(v, 1) for k, v in by.items() if v > 0.5}
    if not by:
        return None
    pmhc_total = round(sum(v for k, v in by.items() if k in ("peptide", "MHC", "β2m")), 1)
    tcr_total = round(sum(v for k, v in by.items() if k in ("Vα", "Vβ", "TCR-other")), 1)
    return {"by_region": by, "pmhc_total": pmhc_total, "tcr_total": tcr_total,
            "total": round(pmhc_total + tcr_total, 1)}


def bsa_reconciliation(bsa_immuno, bsa_decomp) -> dict | None:
    """Cross-check the two independent buried-surface-area estimates.

    md_analysis's ``bsa.buried_surface_area.mean`` is a one-sided ΔSASA; the
    struct decomposition ``total`` sums BOTH partners of the same interface, so
    ``total/2`` is the comparable one-sided figure. A large gap flags an artifact
    in one of the two pipelines. Separately, ``pmhc_total ≈ tcr_total`` by
    construction (same interface counted from each side), so a big asymmetry is a
    second sanity flag. Thresholds PINNED to the 735-traj distribution (2026-07-01):
    half-vs-immuno median 8.7%, only 3/735 exceed 20%; pmhc/tcr asymmetry median
    3.2%, max 14.2% (→ 12% flags the extreme tail). Returns None if either missing.
    """
    im = ((bsa_immuno or {}).get("buried_surface_area") or {}).get("mean")
    total = (bsa_decomp or {}).get("total")
    pmhc = (bsa_decomp or {}).get("pmhc_total")
    tcr = (bsa_decomp or {}).get("tcr_total")
    out: dict = {}
    if im is not None and total is not None:
        one_sided = total / 2.0
        denom = max(one_sided, im, 1e-6)
        out["decomp_one_sided"] = round(one_sided, 1)
        out["immuno_mean"] = round(im, 1)
        out["decomp_half_vs_immuno_pct"] = round(100.0 * abs(one_sided - im) / denom, 1)
        out["consistent"] = out["decomp_half_vs_immuno_pct"] <= 20.0
    if pmhc is not None and tcr is not None:
        m = max(pmhc, tcr, 1e-6)
        out["pmhc_tcr_asymmetry_pct"] = round(100.0 * abs(pmhc - tcr) / m, 1)
        out["pmhc_tcr_balanced"] = out["pmhc_tcr_asymmetry_pct"] <= 12.0
    return out or None


def compute(traj_id: str, meta: dict) -> dict | None:
    """All trajectory-derived structural metrics for one complex (WEB_DATA only)."""
    wd = config.WEB_DATA / traj_id
    top_path, traj_path = wd / config.OUT_TOPOLOGY, wd / config.OUT_TRAJ
    if not top_path.exists() or not traj_path.exists():
        return None

    pep_chain = meta.get("peptide_chain") or _MANUAL_PEPTIDE_CHAIN.get(meta.get("pdb_id"))
    tc = (meta.get("tcr") or {}).get("chains") or {}
    a_ch = (tc.get("alpha") or {}).get("structural_chain")
    b_ch = (tc.get("beta") or {}).get("structural_chain")
    if not pep_chain:
        return None
    # MHC heavy chain by role/geometry (largest non-peptide / non-TCR chain).
    mhc_chain = _mhc_chain(meta, pep_chain, a_ch, b_ch)
    if not mhc_chain:
        return None

    t = md.load(str(traj_path), top=str(top_path))
    if t.n_frames < 2:
        return None
    top = t.topology
    pep = _chain_residues(top, pep_chain)
    mhc = _chain_residues(top, mhc_chain)
    if not pep or not mhc:
        return None
    pep_resid = {ri: top.residue(ri).resSeq for ri in pep}

    # peptide-HLA contact occupancy (closest-heavy), max per peptide residue
    pairs = [(pr, mr) for pr in pep for mr in mhc]
    d, _ = md.compute_contacts(t, pairs, scheme="closest-heavy")  # (frames, npairs) nm
    occ = (d < CONTACT_CUTOFF_NM).mean(axis=0)
    hla_occ: dict[int, float] = {}
    for (pr, _mr), o in zip(pairs, occ):
        hla_occ[pr] = max(hla_occ.get(pr, 0.0), float(o))

    # per-residue SASA (nm^2), strided; total SASA is the per-frame residue sum
    sub = t[::SASA_STRIDE]
    sasa = md.shrake_rupley(sub, mode="residue")  # (nframes, nres)
    pep_sasa = {ri: round(float(sasa[:, ri].mean()), 3) for ri in pep}
    total_sasa = sasa.sum(axis=1)
    bsa = bsa_decomposition(sub, sasa, pep_chain, mhc_chain, a_ch, b_ch)

    # peptide bulge: max peptide Calpha height above the alpha1/alpha2 GROOVE-HELIX
    # rim plane (oriented toward the peptide). MHC chain is now role/geometry
    # resolved, so the helices come from the correct chain.
    platform_ca = [i for i in _chain_ca(top, mhc_chain) if i is not None][:MHC_PLATFORM]
    groove_ca = _groove_helix_ca(top, mhc_chain, platform_ca)
    pep_ca = [i for i in _chain_ca(top, pep_chain) if i is not None]
    bulge = None
    if len(groove_ca) >= 4 and pep_ca:
        gca = np.array(groove_ca)
        pca = np.array(pep_ca)
        heights = []
        for f in range(0, t.n_frames, BULGE_STRIDE):
            G = t.xyz[f, gca]
            c = G.mean(0)
            _, _, vt = np.linalg.svd(G - c)
            n = vt[2]
            pc = t.xyz[f, pca] - c
            if (pc.mean(0) @ n) < 0:
                n = -n
            heights.append(float((pc @ n).max()))
        bulge = round(float(np.mean(heights)) * 10.0, 2)  # nm -> A

    rg = md.compute_rg(t)

    return {
        "peptide_hla": {
            str(pep_resid[ri]): {"hla_contact": round(hla_occ.get(ri, 0.0), 3),
                                 "sasa_nm2": pep_sasa[ri]}
            for ri in pep
        },
        "geometry": {
            "bulge_height_angstrom": bulge,
            "rg_mean_nm": round(float(rg.mean()), 3),
            "rg_std_nm": round(float(rg.std()), 4),
            "total_sasa_mean_nm2": round(float(total_sasa.mean()), 1),
        },
        "bsa_decomposition": bsa,
    }


def merge_peptide_table(table, peptide_hla) -> int:
    """Merge per-residue HLA-contact / SASA / anchor into a peptide_table in place.

    An **anchor** engages the HLA but not the TCR — buried in an HLA pocket (the
    class-I P2 and C-terminal anchors). "Buried" is a low absolute per-residue
    SASA (< ANCHOR_SASA_MAX_NM2; the saturating max-to-any HLA-contact occupancy
    cannot resolve burial); "ignored by the TCR" is tcr_contact < ANCHOR_TCR_MAX.
    An absolute (not below-median) cutoff is used so a peptide with no truly
    buried residue gets no spurious anchor, and so both the B-pocket (P2) and the
    F-pocket (PΩ) anchors are caught even when one is far more buried than the other.

    The peptide_hla dict is keyed by ``str(resSeq)`` (matches peptide_table.resid).
    Returns the number of rows updated. This is the single merge rule shared by the
    full curation path (``_patch_analysis``) and the ``backfill_peptide_hla`` remerge.
    """
    if not table or not peptide_hla:
        return 0
    n = 0
    for r in table:
        ph = peptide_hla.get(str(r.get("resid")))
        if not ph:
            continue
        r["hla_contact"] = ph.get("hla_contact")
        r["sasa_nm2"] = ph.get("sasa_nm2")
        buried = ph.get("sasa_nm2") is not None and ph["sasa_nm2"] < ANCHOR_SASA_MAX_NM2
        r["anchor"] = bool(buried and (r.get("tcr_contact") or 0) < ANCHOR_TCR_MAX)
        n += 1
    return n


def _patch_analysis(analysis: dict, result: dict) -> None:
    """Merge geometry + BSA + per-residue HLA-contact / SASA / anchor into analysis.json."""
    analysis["geometry"] = result["geometry"]
    if result.get("bsa_decomposition"):
        analysis["bsa_decomposition"] = result["bsa_decomposition"]
    table = (analysis.get("interface") or {}).get("peptide_table")
    merge_peptide_table(table, result["peptide_hla"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", default=None, help="Comma-separated traj ids (default: all in WEB_DATA).")
    args = ap.parse_args(argv)
    only = set(args.ids.split(",")) if args.ids else None

    web = config.WEB_DATA
    n_ok = n_skip = 0
    for meta_path in sorted(web.glob("*/" + config.OUT_META)):
        traj_id = meta_path.parent.name
        if only and traj_id not in only:
            continue
        an_path = meta_path.parent / "analysis" / "analysis.json"
        meta = json.loads(meta_path.read_text())
        try:
            result = compute(traj_id, meta)
        except Exception as exc:  # never let one bad trajectory abort the batch
            print(f"[warn] {traj_id}: {exc}")
            result = None
        if result is None:
            n_skip += 1
            continue
        an_dir = meta_path.parent / "analysis"
        an_dir.mkdir(parents=True, exist_ok=True)
        (an_dir / SIDECAR).write_text(json.dumps(result))
        if an_path.exists():
            analysis = json.loads(an_path.read_text())
            _patch_analysis(analysis, result)
            an_path.write_text(json.dumps(analysis))
        n_ok += 1
    print(f"Struct metrics: {n_ok} written, {n_skip} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
