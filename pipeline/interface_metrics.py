"""Interface recognition metrics derived from the per-residue analysis output.

Implements several layers of the bravo analysis plan that share the same raw
inputs (the contact / RMSF / interaction stage CSVs), so they are computed in
one pass per complex:

  * contact entropy        -- how concentrated the TCR-pHLA interface is (3.1)
  * peptide position table -- per-peptide-residue recognition profile  (2.2)
  * interface hotspots     -- multi-evidence hotspot scoring           (6.2)

Each is independent and returns ``None`` when its inputs are missing, so the
module degrades gracefully on incomplete systems. The result is embedded into
``analysis.json`` under ``interface`` by ``extract_analysis``.

Run standalone to (re)compute + patch analysis.json for all systems:
    IMMUNO_WEB_DATA=.../immuno-dyn \
        python -m pipeline.interface_metrics --source /home/xmy/immuno_analysis_run3
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from . import config
from . import cdr_contacts  # reuse topology parsing + CDR residue mapping

CONTACT_FREQ_REL = "analysis/contacts/residue_contact_frequencies.csv"

# Per-residue-pair interaction files for interaction-type diversity (hotspots).
INTERACTION_FILES = {
    "hbond": "hbond/{id}/analysis/interactions/hydrogen_bonds/residue_pair_hbonds.csv",
    "saltbridge": "saltbridge/{id}/analysis/interactions/salt_bridges/residue_pair_salt_bridges.csv",
    "hydrophobic": "hydrophobic/{id}/analysis/interactions/hydrophobic_contacts/residue_pair_hydrophobic_contacts.csv",
    "pipi": "pipi/{id}/analysis/interactions/pi_interactions/residue_pair_pi_pi.csv",
    "cationpi": "cationpi/{id}/analysis/interactions/cation_pi_interactions/residue_pair_cation_pi.csv",
}
RMSF_REL = "rmsf/{id}/analysis/rmsf/residue_rmsf.csv"
RMSF_STABLE_CUTOFF = 3.0  # Angstrom; RMSF at/above this contributes zero stability
GROOVE_SUBREGIONS = {"alpha1_helix", "alpha2_helix"}


def _read_contact_pairs(source_root: Path, traj_id: str) -> list[dict] | None:
    p = source_root / "contact" / traj_id / CONTACT_FREQ_REL
    if not p.exists():
        return None
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh))


def contact_entropy(rows: list[dict]) -> dict | None:
    """Shannon entropy of the TCR-pHLA contact occupancy distribution (bravo 3.1).

    Low entropy  -> recognition concentrated on a few residue pairs (hotspots).
    High entropy -> contacts spread diffusely over the interface.
    Reported both in nats and normalised to [0,1] (divided by ln N) so complexes
    with different interface sizes are comparable.
    """
    occ = [float(r["contact_frequency"]) for r in rows
           if r.get("contact_frequency") not in (None, "") and float(r["contact_frequency"]) > 0]
    total = sum(occ)
    if total <= 0 or len(occ) < 2:
        return None
    ps = [o / total for o in occ]
    h = -sum(p * math.log(p) for p in ps)
    h_max = math.log(len(occ))
    ps_sorted = sorted(ps, reverse=True)
    return {
        "n_pairs": len(occ),
        "entropy_nats": round(h, 3),
        "entropy_normalized": round(h / h_max, 3) if h_max > 0 else None,
        "top1_share": round(ps_sorted[0], 3),
        "top3_share": round(sum(ps_sorted[:3]), 3),
    }


def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _tcr_cdr3_resids(traj_id: str, meta: dict, source_root: Path) -> dict[str, set]:
    """{structural_chain_letter: set(cdr3 resids)} for the alpha/beta TCR chains."""
    tcr = meta.get("tcr")
    if not tcr or not tcr.get("chains"):
        return {}
    summ = source_root / "contact" / traj_id / cdr_contacts.SUMMARY_REL
    if not summ.exists():
        return {}
    topo = Path(json.loads(summ.read_text()).get("topology", ""))
    if not topo.exists():
        return {}
    chains = cdr_contacts.parse_topology_chains(topo)
    out: dict[str, set] = {}
    for _pos, c in tcr["chains"].items():
        letter = c.get("structural_chain")
        if letter in chains:
            cmap = cdr_contacts.cdr_resids_for_chain(chains[letter], c)
            out[letter] = cmap.get("cdr3", set())
    return out


def peptide_position_table(traj_id: str, meta: dict, source_root: Path,
                           contact_rows: list[dict]) -> list[dict] | None:
    """Per-peptide-residue recognition + dynamics profile (bravo 2.2).

    For each peptide position: max TCR contact occupancy, max occupancy to a
    CDR3 loop residue, number of distinct TCR partner residues, backbone RMSF,
    and max H-bond occupancy to the TCR. Anchor (peptide-HLA) contacts are not
    in the pHLA-vs-TCR contact output, so they are left for the trajectory pass.
    """
    pep = meta.get("peptide_chain")
    if not pep:
        return None
    cdr3 = _tcr_cdr3_resids(traj_id, meta, source_root)

    per: dict[int, dict] = {}
    for r in contact_rows:
        if r.get("chain_id_1") != pep:
            continue
        try:
            resid = int(r["resid_1"])
            occ = float(r["contact_frequency"])
        except (ValueError, KeyError):
            continue
        d = per.setdefault(resid, {"resname": r.get("resname_1"), "tcr_max": 0.0,
                                   "partners": set(), "cdr3_max": 0.0})
        d["tcr_max"] = max(d["tcr_max"], occ)
        d["partners"].add((r.get("chain_id_2"), r.get("resid_2")))
        tchain = r.get("chain_id_2")
        try:
            tresid = int(r["resid_2"])
        except (ValueError, KeyError):
            tresid = None
        if tresid is not None and tresid in cdr3.get(tchain, ()):
            d["cdr3_max"] = max(d["cdr3_max"], occ)

    # per-residue backbone RMSF (canonical peptide residue set)
    rmsf: dict[int, float] = {}
    resname: dict[int, str] = {}
    for r in _read_csv(source_root / "rmsf" / traj_id / "analysis/rmsf/residue_rmsf.csv"):
        if r.get("chain_id") == pep:
            try:
                rid = int(r["resid"])
                rmsf[rid] = float(r["rmsf_angstrom"])
                resname[rid] = r.get("resname") or resname.get(rid, "")
            except (ValueError, KeyError):
                pass

    # max H-bond occupancy to TCR per peptide residue (peptide is on the pHLA side)
    hb: dict[int, float] = {}
    for r in _read_csv(source_root / "hbond" / traj_id
                       / "analysis/interactions/hydrogen_bonds/residue_pair_hbonds.csv"):
        for a in ("1", "2"):
            if r.get(f"chain_id_{a}") == pep:
                try:
                    rid = int(r[f"resid_{a}"])
                    hb[rid] = max(hb.get(rid, 0.0), float(r["contact_frequency"]))
                except (ValueError, KeyError):
                    pass

    all_resids = sorted(set(per) | set(rmsf))
    if not all_resids:
        return None
    table = []
    for i, rid in enumerate(all_resids):
        d = per.get(rid, {})
        table.append({
            "position": i + 1,
            "resid": rid,
            "resname": resname.get(rid) or d.get("resname") or "",
            "tcr_contact": round(d.get("tcr_max", 0.0), 3),
            "cdr3_contact": round(d.get("cdr3_max", 0.0), 3),
            "n_tcr_partners": len(d.get("partners", ())),
            "rmsf": round(rmsf[rid], 2) if rid in rmsf else None,
            "hbond": round(hb.get(rid, 0.0), 3),
        })
    return table


def _cdr_by_chain(traj_id: str, meta: dict, source_root: Path) -> dict[str, dict]:
    """{structural_chain_letter: {cdr1/2/3: set(resid)}}."""
    out: dict[str, dict] = {}
    tcr = meta.get("tcr") or {}
    summ = source_root / "contact" / traj_id / cdr_contacts.SUMMARY_REL
    if not summ.exists():
        return out
    topo = Path(json.loads(summ.read_text()).get("topology", ""))
    if not topo.exists():
        return out
    chains = cdr_contacts.parse_topology_chains(topo)
    for c in (tcr.get("chains") or {}).values():
        letter = c.get("structural_chain")
        if letter in chains:
            cmap = cdr_contacts.cdr_resids_for_chain(chains[letter], c)
            out[letter] = {k: cmap.get(k, set()) for k in ("cdr1", "cdr2", "cdr3")}
    return out


def hotspots(traj_id: str, meta: dict, source_root: Path, contact_rows: list[dict]) -> list[dict] | None:
    """Multi-evidence interface hotspot scoring (bravo 6.2).

    For every interface residue, combine: contact occupancy (engagement), low
    RMSF (conformational stability), and number of distinct interaction types
    (contact / H-bond / salt bridge / hydrophobic / pi-pi / cation-pi). Residues
    that are persistently engaged, rigid and multi-modal score highest. Each is
    categorised (peptide / HLA-groove / CDR3 / CDR1-2 / TCR-FR) so the list maps
    onto recognition roles and feeds FEP-candidate selection.
    """
    # per-residue identity + RMSF + structural role
    info: dict = {}
    for r in _read_csv(source_root / RMSF_REL.format(id=traj_id)):
        try:
            key = (r["chain_id"], int(r["resid"]))
            info[key] = {
                "resname": r.get("resname"),
                "rmsf": float(r["rmsf_angstrom"]) if r.get("rmsf_angstrom") else None,
                "component": r.get("component"),
                "mhc_subregion": r.get("mhc_subregion"),
            }
        except (ValueError, KeyError):
            pass

    # contact engagement (both sides of every pair)
    eng: dict = {}
    for r in contact_rows:
        try:
            occ = float(r["contact_frequency"])
        except (ValueError, KeyError):
            continue
        for s in ("1", "2"):
            try:
                key = (r[f"chain_id_{s}"], int(r[f"resid_{s}"]))
            except (ValueError, KeyError):
                continue
            d = eng.setdefault(key, {"occ": 0.0, "partners": set()})
            d["occ"] = max(d["occ"], occ)
            other = "2" if s == "1" else "1"
            d["partners"].add((r.get(f"chain_id_{other}"), r.get(f"resid_{other}")))

    # interaction-type membership per residue
    types: dict = {k: {"contact"} for k in eng}
    for tname, rel in INTERACTION_FILES.items():
        for r in _read_csv(source_root / rel.format(id=traj_id)):
            for s in ("1", "2"):
                try:
                    key = (r[f"chain_id_{s}"], int(r[f"resid_{s}"]))
                except (ValueError, KeyError):
                    continue
                types.setdefault(key, set()).add(tname)

    cdr_by_chain = _cdr_by_chain(traj_id, meta, source_root)

    rows = []
    for key, e in eng.items():
        if e["occ"] <= 0:
            continue
        ch, rid = key
        m = info.get(key, {})
        comp = m.get("component")
        rmsf = m.get("rmsf")
        n_types = len(types.get(key, {"contact"}))

        if comp == "peptide":
            side, region, category = "peptide", "peptide", "peptide"
        elif comp == "HLA_alpha":
            groove = m.get("mhc_subregion") in GROOVE_SUBREGIONS
            side, region, category = "hla", ("groove" if groove else "non_groove"), "HLA-groove" if groove else "HLA"
        elif comp in ("TCR_alpha", "TCR_beta"):
            side = "tcr"
            cdrs = cdr_by_chain.get(ch, {})
            region = "fr"
            for ck in ("cdr3", "cdr2", "cdr1"):
                if rid in cdrs.get(ck, ()):
                    region = ck
                    break
            greek = "α" if comp == "TCR_alpha" else "β"
            category = {"cdr3": f"CDR3{greek}", "cdr2": f"CDR1/2{greek}",
                        "cdr1": f"CDR1/2{greek}", "fr": f"TCR-FR{greek}"}[region]
        elif comp == "beta2m":
            side, region, category = "b2m", "-", "β2m"
        else:
            side, region, category = "other", "-", "other"

        occ_score = e["occ"]
        stab_score = max(0.0, 1 - rmsf / RMSF_STABLE_CUTOFF) if rmsf is not None else 0.5
        div_score = min(n_types, 5) / 5.0
        score = 0.4 * occ_score + 0.3 * stab_score + 0.3 * div_score
        rows.append({
            "chain": ch, "resid": rid, "resname": m.get("resname"),
            "side": side, "region": region, "category": category,
            "occupancy": round(e["occ"], 3),
            "rmsf": round(rmsf, 2) if rmsf is not None else None,
            "n_interaction_types": n_types, "n_partners": len(e["partners"]),
            "score": round(score, 3),
        })

    hs = [h for h in rows if h["score"] >= 0.5 and h["occupancy"] >= 0.5]
    hs.sort(key=lambda h: h["score"], reverse=True)
    return hs[:20] or None


# Residue-specific primary substitution for FEP-candidate generation. Each leads
# with a chemically informative point mutation (remove a specific H-bond, vary
# size, reverse charge), so candidates are not all collapsed onto alanine; the
# generic alanine scan is appended as a secondary, knockout-style control.
_PRIMARY_MUT = {
    "TYR": "→Phe (remove hydroxyl H-bond)",
    "PHE": "→Tyr/Trp (aromatic ring variation)",
    "TRP": "→Phe (reduce indole bulk)",
    "HIS": "→Asn/Gln (remove charge & aromaticity)",
    "ARG": "→Lys (shorten) or →Glu (charge reversal)",
    "LYS": "→Arg (extend) or →Glu (charge reversal)",
    "ASP": "→Asn (neutralize) or →Lys (charge reversal)",
    "GLU": "→Gln (neutralize) or →Arg (charge reversal)",
    "LEU": "→Ile/Val (size variation)",
    "ILE": "→Leu/Val (size variation)",
    "VAL": "→Ile (enlarge) or →Ala (truncate)",
    "MET": "→Leu/Ile (remove sulfur)",
    "ALA": "→Ser (add polarity) or →Val (add bulk)",
    "SER": "→Thr (β-branch) or →Ala (remove –OH)",
    "THR": "→Ser (remove methyl) or →Val (remove –OH)",
    "ASN": "→Asp (H-bond tuning) or →Gln (extend)",
    "GLN": "→Glu (H-bond tuning) or →Asn (shorten)",
    "CYS": "→Ser (remove thiol)",
    "GLY": "→Ala (add methyl; backbone-sensitive)",
    "PRO": "→Ala (restore backbone H-bond; interpret cautiously)",
}


def _mutation_suggestion(resname: str) -> str:
    r = (resname or "").upper()
    primary = _PRIMARY_MUT.get(r)
    if not primary:
        return "→Ala (alanine scan)"
    # When the residue-specific suggestion already centres on Ala (Gly/Pro), don't
    # tack a redundant alanine scan on the end.
    if "→Ala" in primary:
        return primary
    return f"{primary}; →Ala (knockout scan)"


def fep_candidates(hs: list[dict], meta: dict, pep_pos: dict | None = None) -> list[dict] | None:
    """FEP-ready mutation candidates derived from the top interface hotspots (7.2).

    pep_pos maps a peptide residue id to its sequence position (from the peptide
    table); when absent, the peptide resid is shown as the position fallback.
    """
    if not hs:
        return None
    pep_pos = pep_pos or {}

    out = []
    for h in hs:
        side, region = h["side"], h["region"]
        if side == "peptide":
            pos = pep_pos.get(h["resid"], h["resid"])
            label = f"Peptide P{pos} {h['resname']}"
            ctype = "peptide TCR-facing (escape candidate)"
            risk = "may abrogate TCR recognition; check HLA-binding if near an anchor"
            validation = "FEP + SPR; T-cell activation assay"
        elif side == "tcr" and region == "cdr3":
            label = f"{h['category']} {h['resname']}{h['resid']}"
            ctype = "CDR3 affinity / specificity"
            risk = "alters TCR affinity; may shift specificity"
            validation = "FEP + SPR; alanine scan"
        elif side == "hla":
            label = f"HLA {h['resname']}{h['resid']}"
            ctype = "HLA contact / restriction"
            risk = "may affect peptide presentation or allele restriction"
            validation = "FEP; alanine scan"
        else:
            continue
        out.append({
            "residue": label,
            "candidate_type": ctype,
            "evidence": f"occupancy {h['occupancy']:.2f}, {h['n_interaction_types']} interaction types,"
                        f" RMSF {h['rmsf']:.2f} Å" if h["rmsf"] is not None
                        else f"occupancy {h['occupancy']:.2f}, {h['n_interaction_types']} interaction types",
            "suggested_mutation": _mutation_suggestion(h["resname"]),
            "risk": risk,
            "recommended_validation": validation,
            "score": h["score"],
        })
        if len(out) >= 8:
            break
    return out or None


# ---- per-pair H-bond / salt-bridge persistence tables (bravo §4.x) -------------
# Surfaces the specific ionic clips / anchor H-bonds that are stable, split by
# immunological role: anchor (peptide-MHC), specificity (peptide-TCR), restriction
# (MHC-TCR). The per-pair CSVs already carry occupancy + segment statistics.
_ROLE_OF_COMP = {"peptide": "peptide", "HLA_alpha": "mhc", "beta2m": "b2m",
                 "TCR_alpha": "tcr", "TCR_beta": "tcr"}
_PAIR_FILES = {
    "hbonds": "hbond/{id}/analysis/interactions/hydrogen_bonds/residue_pair_hbonds.csv",
    "saltbridges": "saltbridge/{id}/analysis/interactions/salt_bridges/residue_pair_salt_bridges.csv",
}


def _chain_components(traj_id: str, source_root: Path) -> dict[str, str]:
    """{chain_id: component} (HLA_alpha / beta2m / peptide / TCR_alpha / TCR_beta)
    from the per-residue RMSF CSV — the robust, role-labelled chain map."""
    comp: dict[str, str] = {}
    for r in _read_csv(source_root / RMSF_REL.format(id=traj_id)):
        ch, c = r.get("chain_id"), r.get("component")
        if ch and c and ch not in comp:
            comp[ch] = c
    return comp


def _pair_group(role_a: str, role_b: str) -> str | None:
    s = {role_a, role_b}
    if s == {"peptide", "mhc"}:
        return "anchor"          # peptide ↔ MHC groove
    if s == {"peptide", "tcr"}:
        return "specificity"     # peptide ↔ TCR (antigen readout)
    if s == {"mhc", "tcr"}:
        return "restriction"     # MHC ↔ TCR (germline restriction)
    return None                  # intra-component / β2m → not an interface clip


def pairwise_interactions(traj_id: str, source_root: Path, top: int = 12) -> dict | None:
    """Top inter-component H-bond / salt-bridge pairs with occupancy + persistence.

    longest_frac = max_consecutive_frames / total_frames (longest unbroken stretch);
    n_segments = how fragmented the contact is (1 = one stable clip, high = flickering).
    """
    chain_comp = _chain_components(traj_id, source_root)
    if not chain_comp:
        return None
    role = lambda ch: _ROLE_OF_COMP.get(chain_comp.get(ch), "other")  # noqa: E731
    out: dict = {}
    for kind, rel in _PAIR_FILES.items():
        items = []
        for r in _read_csv(source_root / rel.format(id=traj_id)):
            try:
                occ = float(r["contact_frequency"])
            except (ValueError, KeyError):
                continue
            if occ <= 0:
                continue
            c1, c2 = r.get("chain_id_1"), r.get("chain_id_2")
            grp = _pair_group(role(c1), role(c2))
            if grp is None:
                continue
            try:
                tot = float(r.get("total_frames") or 0)
                maxc = float(r.get("max_consecutive_frames") or 0)
            except ValueError:
                tot = maxc = 0.0
            lab1 = r.get("residue_label_1") or f"{r.get('resname_1', '')}{r.get('resid_1', '')}"
            lab2 = r.get("residue_label_2") or f"{r.get('resname_2', '')}{r.get('resid_2', '')}"
            items.append({
                "group": grp,
                "a": lab1, "a_role": role(c1), "b": lab2, "b_role": role(c2),
                "occupancy": round(occ, 3),
                "longest_frac": round(maxc / tot, 3) if tot else None,
                "n_segments": int(r.get("n_segments") or 0),
            })
        items.sort(key=lambda x: x["occupancy"], reverse=True)
        if items:
            out[kind] = items[:top]
    return out or None


# ---- region-keyed RMSF flexibility profile (bravo §5.x) ------------------------
# Promotes the per-residue RMSF (already region-labelled in residue_rmsf.csv) into a
# queryable flexibility fingerprint: which loops are rigid (specific) vs plastic.
_RMSF_REGION_ORDER = ["peptide", "CDR1α", "CDR2α", "CDR3α", "CDR1β", "CDR2β", "CDR3β",
                      "α1 helix", "α2 helix", "MHC floor", "β2m", "TCR-FRα", "TCR-FRβ"]


def _rmsf_region(r: dict) -> str:
    comp = r.get("component")
    if comp == "peptide":
        return "peptide"
    if comp == "beta2m":
        return "β2m"
    if comp == "HLA_alpha":
        return {"alpha1_helix": "α1 helix", "alpha2_helix": "α2 helix"}.get(
            r.get("mhc_subregion"), "MHC floor")
    if comp in ("TCR_alpha", "TCR_beta"):
        g = "α" if comp == "TCR_alpha" else "β"
        tr = r.get("tcr_region")
        return f"{tr}{g}" if tr in ("CDR1", "CDR2", "CDR3") else f"TCR-FR{g}"
    return "other"


def rmsf_profile(traj_id: str, source_root: Path) -> dict | None:
    """Per-residue Cα RMSF labelled by structural region + per-region summaries."""
    rows = _read_csv(source_root / RMSF_REL.format(id=traj_id))
    if not rows:
        return None
    residues = []
    region_vals: dict[str, list] = {}
    for r in rows:
        try:
            rmsf = float(r["rmsf_angstrom"])
        except (ValueError, KeyError):
            continue
        reg = _rmsf_region(r)
        residues.append({"chain": r.get("chain_id"), "resid": int(r["resid"]),
                         "resname": r.get("resname"), "region": reg, "rmsf": round(rmsf, 2)})
        region_vals.setdefault(reg, []).append(rmsf)
    if not residues:
        return None
    regions = {reg: {"mean": round(sum(v) / len(v), 2), "max": round(max(v), 2), "n": len(v)}
               for reg, v in region_vals.items()}
    order = [reg for reg in _RMSF_REGION_ORDER if reg in regions]
    return {"residues": residues, "regions": regions, "order": order}


def compute(traj_id: str, meta: dict, source_root: Path) -> dict | None:
    """All interface metrics for one complex. Returns None if no contact data."""
    rows = _read_contact_pairs(source_root, traj_id)
    if not rows:
        return None
    out: dict = {}
    ent = contact_entropy(rows)
    if ent:
        out["contact_entropy"] = ent
    pep_table = peptide_position_table(traj_id, meta, source_root, rows)
    if pep_table:
        out["peptide_table"] = pep_table
    hs = hotspots(traj_id, meta, source_root, rows)
    if hs:
        out["hotspots"] = hs
        pep_pos = {r["resid"]: r["position"] for r in (pep_table or [])}
        fep = fep_candidates(hs, meta, pep_pos)
        if fep:
            out["fep_candidates"] = fep
    pw = pairwise_interactions(traj_id, source_root)
    if pw:
        out["pairwise"] = pw
    return out or None


def _iter_traj_dirs(source_root: Path):
    contact_root = source_root / "contact"
    if not contact_root.is_dir():
        return
    for d in sorted(contact_root.iterdir()):
        if d.is_dir():
            yield d.name


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", required=True, type=Path)
    args = ap.parse_args(argv)

    n_ok = n_skip = 0
    for traj_id in _iter_traj_dirs(args.source):
        meta_path = config.WEB_DATA / traj_id / config.OUT_META
        if not meta_path.exists():
            n_skip += 1
            continue
        meta = json.loads(meta_path.read_text())
        result = compute(traj_id, meta, args.source)
        rprof = rmsf_profile(traj_id, args.source)
        if result is None and rprof is None:
            n_skip += 1
            continue
        an_json = config.WEB_DATA / traj_id / "analysis" / "analysis.json"
        if an_json.exists():
            data = json.loads(an_json.read_text())
            if result is not None:
                data["interface"] = result
            if rprof is not None:
                data["rmsf_profile"] = rprof
            an_json.write_text(json.dumps(data))
        n_ok += 1

    print(f"Interface metrics: {n_ok} written, {n_skip} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
