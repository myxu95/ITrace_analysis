"""Re-decompose TCR-pHLA contacts by *correct* ANARCI CDR boundaries.

Why this exists
---------------
The upstream immunoscope contact stage annotates each TCR contact residue as
CDR1/2/3 or framework using its own regex-based CDR detection. That regex keys
on the first Cys...aromatic motif and mislocates the loops -- e.g. for 1ao7 it
labels CDR1 residues as "CDR3" and dumps the genuine CDR3 (and CDR2) contacts
into "non_cdr". It also never finds CDR1/CDR2 at all (counts are always 0). Any
CDR-contribution metric built on those labels is therefore wrong.

This module re-derives the loop boundaries from the ANARCI/IMGT annotation that
`annotate_tcr` wrote into `meta.tcr` (run that first), maps each CDR's amino-acid
sequence onto residue numbers in the trajectory topology, and re-classifies the
per-residue contact frequencies the contact stage already computed. From that it
reports, per complex:

  * peptide-recognition ratio  = CDR3-peptide contacts / all TCR-pHLA contacts
  * HLA-restriction ratio      = (CDR1+CDR2)-HLA contacts / all TCR-pHLA contacts
  * a region x partner breakdown (counts and occupancy-weighted)
  * alpha vs beta chain contribution to the interface

Both an occupancy-weighted form (primary; a persistent contact counts more than
a transient one) and a raw residue-pair-count form are produced.

Usage (after annotate_tcr, reads the contact stage output):
    IMMUNO_WEB_DATA=.../immuno-dyn-run3 \
        python -m pipeline.cdr_contacts --source /home/xmy/immuno_analysis_run3
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from . import config

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
    # common protonation / modified variants seen in CHARMM-processed PDBs
    "HSD": "H", "HSE": "H", "HSP": "H", "HID": "H", "HIE": "H", "HIP": "H",
    "CYX": "C", "CYM": "C", "MSE": "M", "SEC": "U", "PYL": "O",
}
CDR_KEYS = ("cdr1", "cdr2", "cdr3")
CONTACT_REL = "analysis/contacts/residue_contact_frequencies.csv"
SUMMARY_REL = "analysis/contacts/residue_contact_summary.json"

# A real TCR-pHLA interface has at most a few hundred contacting residue PAIRS
# (clean dataset max ~300). Far above that, the upstream contact stage glitched
# and the occupancy_total denominator is bogus, diluting every ratio (e.g.
# 6g9q_run2 = 5105 pairs -> pep_recog 4.7% vs the clean run3 replica's 26.5%).
MAX_TCR_PHLA_PAIRS = 500


def is_unreliable(tcr: dict | None) -> str | None:
    """Reason the CDR decomposition should NOT be trusted, else None.

    - "degenerate": occupancy exists but no CDR loop was matched (a parse/numbering
      failure), so peptide/HLA ratios are 0 by artifact, not biology (e.g. 6vma).
    - "pair_count_outlier": implausibly many TCR-pHLA residue pairs -> the upstream
      contact stage over-counted and the ratio denominator is meaningless.
    - "single_chain_collapse": one TCR chain contributes ~zero occupancy while the
      other takes everything (e.g. the beta chain missing entirely from the contact
      output). A docked TCR always engages with both chains, so alpha=1.0/beta=0.0 is
      a chain-labeling / incomplete-contact artefact, not biology (e.g. 7rk7/5hhm/4prh
      run3, whose contact output lost the beta chain).
    """
    if not tcr:
        return None
    occ_total = tcr.get("occupancy_total") or 0
    by_region = tcr.get("by_region") or {}
    if occ_total > 0 and by_region:
        cdr_occ = sum((by_region.get(k) or {}).get(p, {}).get("occ", 0.0)
                      for k in CDR_KEYS for p in ("peptide", "mhc"))
        if cdr_occ == 0.0:
            return "degenerate"
    if (tcr.get("n_tcr_phla_pairs") or 0) > MAX_TCR_PHLA_PAIRS:
        return "pair_count_outlier"
    a_c, b_c = tcr.get("alpha_contribution"), tcr.get("beta_contribution")
    if occ_total > 0 and a_c is not None and b_c is not None and min(a_c, b_c) < 0.02:
        return "single_chain_collapse"
    return None


def parse_topology_chains(pdb_path: Path) -> dict[str, list[tuple[int, str]]]:
    """chain_id -> ordered list of (resid, one_letter), one entry per residue."""
    chains: dict[str, list[tuple[int, str]]] = {}
    seen: dict[str, set[int]] = {}
    with pdb_path.open() as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21]
            try:
                resid = int(line[22:26])
            except ValueError:
                continue
            if resid in seen.setdefault(chain, set()):
                continue
            seen[chain].add(resid)
            resname = line[17:20].strip()
            chains.setdefault(chain, []).append((resid, THREE_TO_ONE.get(resname, "X")))
    return chains


def cdr_resids_for_chain(residues: list[tuple[int, str]], cdr_seqs: dict) -> dict:
    """Map each CDR amino-acid sequence onto residue numbers via substring match.

    Returns {cdr_key: set(resid)} plus "_unmatched": [cdr_keys not located].
    """
    seq = "".join(aa for _, aa in residues)
    resids = [r for r, _ in residues]
    out: dict[str, set[int]] = {k: set() for k in CDR_KEYS}
    unmatched = []
    for key in CDR_KEYS:
        sub = (cdr_seqs.get(key) or "").strip()
        if not sub:
            continue
        idx = seq.find(sub)
        if idx < 0:
            unmatched.append(key)
            continue
        out[key] = set(resids[idx:idx + len(sub)])
    out["_unmatched"] = unmatched
    return out


def _region_of(resid: int, cdr_map: dict) -> str:
    for key in CDR_KEYS:
        if resid in cdr_map.get(key, ()):
            return key
    return "fr"


def _blank_bucket() -> dict:
    return {"peptide": {"n": 0, "occ": 0.0}, "mhc": {"n": 0, "occ": 0.0}}


def decompose(traj_id: str, meta: dict, source_root: Path) -> dict | None:
    """Compute the CDR contact decomposition for one trajectory.

    Returns None when the prerequisites are missing (no TCR annotation, no
    contact output, or no topology).
    """
    tcr = meta.get("tcr")
    if not tcr or not tcr.get("chains"):
        return None
    contact_csv = source_root / "contact" / traj_id / CONTACT_REL
    if not contact_csv.exists():
        return None

    # Topology path is recorded by the contact stage; fall back to source layout.
    topo = None
    summ = source_root / "contact" / traj_id / SUMMARY_REL
    if summ.exists():
        topo = Path(json.loads(summ.read_text()).get("topology", ""))
    if not topo or not topo.exists():
        return None
    chains = parse_topology_chains(topo)

    # structural chain letter -> ("alpha"/"beta", cdr_resid_map)
    chain_roles: dict[str, tuple[str, dict]] = {}
    cdr_resids_dump: dict[str, dict] = {}
    unmatched_any = []
    for position, c in tcr["chains"].items():
        letter = c.get("structural_chain")
        if not letter or letter not in chains:
            continue
        cmap = cdr_resids_for_chain(chains[letter], c)
        unmatched_any += [f"{position}:{k}" for k in cmap.pop("_unmatched", [])]
        chain_roles[letter] = (position, cmap)
        cdr_resids_dump[position] = {k: sorted(v) for k, v in cmap.items()}

    if not chain_roles:
        return None

    peptide_chain = meta.get("peptide_chain")

    by_region = {k: _blank_bucket() for k in (*CDR_KEYS, "fr")}
    by_chain = {"alpha": {"n": 0, "occ": 0.0}, "beta": {"n": 0, "occ": 0.0}}
    n_total = 0
    occ_total = 0.0

    with contact_csv.open(newline="") as fh:
        for row in csv.DictReader(fh):
            # side 1 = pHLA, side 2 = TCR (selection_1 = pHLA, selection_2 = TCR)
            tcr_chain = row["chain_id_2"]
            role = chain_roles.get(tcr_chain)
            if role is None:
                continue
            position, cmap = role
            try:
                tcr_resid = int(row["resid_2"])
                occ = float(row["contact_frequency"])
            except (ValueError, KeyError):
                continue
            partner = "peptide" if row["chain_id_1"] == peptide_chain else "mhc"
            region = _region_of(tcr_resid, cmap)

            by_region[region][partner]["n"] += 1
            by_region[region][partner]["occ"] += occ
            by_chain[position]["n"] += 1
            by_chain[position]["occ"] += occ
            n_total += 1
            occ_total += occ

    if n_total == 0:
        return None

    def ratio(num_occ: float) -> float:
        return round(num_occ / occ_total, 4) if occ_total else 0.0

    def ratio_n(num_n: int) -> float:
        return round(num_n / n_total, 4) if n_total else 0.0

    cdr3_pep_occ = by_region["cdr3"]["peptide"]["occ"]
    cdr3_pep_n = by_region["cdr3"]["peptide"]["n"]
    cdr12_hla_occ = by_region["cdr1"]["mhc"]["occ"] + by_region["cdr2"]["mhc"]["occ"]
    cdr12_hla_n = by_region["cdr1"]["mhc"]["n"] + by_region["cdr2"]["mhc"]["n"]

    result = {
        "method": "anarci_imgt",
        "n_tcr_phla_pairs": n_total,
        "occupancy_total": round(occ_total, 3),
        "by_region": {
            k: {p: {"n": v[p]["n"], "occ": round(v[p]["occ"], 3)} for p in ("peptide", "mhc")}
            for k, v in by_region.items()
        },
        "by_chain": {k: {"n": v["n"], "occ": round(v["occ"], 3)} for k, v in by_chain.items()},
        "peptide_recognition_ratio": ratio(cdr3_pep_occ),
        "hla_restriction_ratio": ratio(cdr12_hla_occ),
        "peptide_recognition_ratio_paircount": ratio_n(cdr3_pep_n),
        "hla_restriction_ratio_paircount": ratio_n(cdr12_hla_n),
        "alpha_contribution": ratio(by_chain["alpha"]["occ"]),
        "beta_contribution": ratio(by_chain["beta"]["occ"]),
        "cdr_resids": cdr_resids_dump,
        "unmatched_cdrs": unmatched_any,
    }
    reason = is_unreliable(result)
    result["reliable"] = reason is None
    if reason:
        result["unreliable_reason"] = reason
    return result


def _iter_traj_dirs(source_root: Path):
    contact_root = source_root / "contact"
    if not contact_root.is_dir():
        return
    for d in sorted(contact_root.iterdir()):
        if d.is_dir():
            yield d.name


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", required=True, type=Path,
                    help="analysis output root (the run_analysis --output dir)")
    ap.add_argument("--patch-analysis", action="store_true", default=True,
                    help="write result into web_data/<id>/analysis/analysis.json")
    args = ap.parse_args(argv)

    n_ok = n_skip = n_unmatched = 0
    for traj_id in _iter_traj_dirs(args.source):
        meta_path = config.WEB_DATA / traj_id / config.OUT_META
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        result = decompose(traj_id, meta, args.source)
        if result is None:
            n_skip += 1
            continue
        if result["unmatched_cdrs"]:
            n_unmatched += 1
        # sidecar (always) + patch analysis.json (if present)
        an_dir = config.WEB_DATA / traj_id / "analysis"
        an_dir.mkdir(parents=True, exist_ok=True)
        (an_dir / "tcr_cdr.json").write_text(json.dumps(result, indent=2))
        an_json = an_dir / "analysis.json"
        if an_json.exists():
            data = json.loads(an_json.read_text())
            data["tcr_cdr"] = result
            an_json.write_text(json.dumps(data))
        n_ok += 1

    print(f"CDR contact decomposition: {n_ok} written, {n_skip} skipped, "
          f"{n_unmatched} had >=1 unmatched CDR.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
