"""Enrich trajectories with annotations fetched from the RCSB PDB REST API.

For every unique 4-char PDB id we fetch the entry + its polymer entities, then
derive: title, resolution, method, release year, primary citation, and a
per-entity table (description, organism, genes, chains). From the entity
descriptions we heuristically assign chain roles (MHC heavy / b2m / peptide /
TCR alpha / TCR beta) and extract a best-effort HLA allele and antigen name.

Results are cached in web_data/rcsb_cache.json (keyed by pdb_id) so re-runs are
cheap and resumable, then merged into each trajectory's meta.json.

Usage:
    python -m pipeline.enrich_rcsb            # fetch (cached) + merge into meta
    python -m pipeline.enrich_rcsb --refresh  # ignore cache, re-fetch all
"""
from __future__ import annotations

import json
import re
import sys
import time

import requests
from tqdm import tqdm

from . import config


def _get(url: str, retries: int = 3) -> dict | None:
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=config.RCSB_TIMEOUT)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.json()
        except requests.RequestException:
            if attempt == retries - 1:
                return None
            time.sleep(1.5 * (attempt + 1))
    return None


def normalize_allele(text: str) -> str | None:
    """`HLA-A 0201` / `HLA-A*0201` -> `HLA-A*02:01` (best effort)."""
    if not text:
        return None
    m = re.search(r"HLA[- ]?([A-Z]+)\d?\*?\s*(\d{2})(\d{2})?", text, re.I)
    if not m:
        # Fallback: only accept a real gene-letter allele like "HLA-A" / "HLA-B7",
        # never generic phrasing such as "HLA CLASS I".
        m2 = re.search(r"\bHLA[- ]([A-EG])\b\*?\s*(\d{2})?(\d{2})?", text, re.I)
        if not m2:
            return None
        gene, g1, g2 = m2.group(1), m2.group(2), m2.group(3)
        allele = f"HLA-{gene.upper()}"
        if g1:
            allele += f"*{g1}" + (f":{g2}" if g2 else "")
        return allele
    gene, g1, g2 = m.group(1), m.group(2), m.group(3)
    allele = f"HLA-{gene.upper()}*{g1}"
    if g2:
        allele += f":{g2}"
    return allele


def classify_entity(desc: str, n_res: int | None) -> str:
    """Map a polymer-entity description to a coarse role."""
    d = (desc or "").lower()
    if "beta-2" in d or "b2m" in d or "microglobulin" in d:
        return "b2m"
    if any(k in d for k in ("histocompatibility", "hla", "mhc", "h-2", "class i")):
        return "mhc"
    if "alpha" in d and ("t-cell" in d or "t cell" in d or "tcr" in d or "receptor" in d):
        return "tcr_alpha"
    if "beta" in d and ("t-cell" in d or "t cell" in d or "tcr" in d or "receptor" in d):
        return "tcr_beta"
    if any(k in d for k in ("t-cell receptor", "t cell receptor", "tcr")):
        return "tcr"
    if "peptide" in d or (n_res is not None and n_res <= 30):
        return "peptide"
    return "other"


def fetch_pdb(pdb_id: str) -> dict:
    """Fetch and distill RCSB annotation for one PDB id."""
    entry = _get(config.RCSB_ENTRY_URL.format(pdb_id=pdb_id))
    if entry is None:
        return {"pdb_id": pdb_id, "found": False}

    cit = entry.get("rcsb_primary_citation", {}) or {}
    rel = entry.get("rcsb_accession_info", {}).get("initial_release_date", "") or ""
    res = entry.get("rcsb_entry_info", {}).get("resolution_combined") or []
    info = {
        "pdb_id": pdb_id,
        "found": True,
        "title": (entry.get("struct", {}) or {}).get("title"),
        "method": (entry.get("exptl", [{}]) or [{}])[0].get("method"),
        "resolution": round(res[0], 2) if res else None,  # 2 dp (Å); RCSB sends float artifacts
        "release_year": rel[:4] if rel else None,
        "citation": {
            "journal": cit.get("journal_abbrev"),
            "year": cit.get("year"),
            "doi": cit.get("pdbx_database_id_DOI"),
            "title": cit.get("title"),
            "pubmed": cit.get("pdbx_database_id_PubMed"),
        },
        "entities": [],
    }

    entity_ids = (
        entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids")
        or []
    )
    hla_allele = None
    organisms: set[str] = set()
    tcr_genes: list[str] = []

    for eid in entity_ids:
        pe = _get(config.RCSB_POLYMER_URL.format(pdb_id=pdb_id, entity_id=eid))
        if pe is None:
            continue
        desc = (pe.get("rcsb_polymer_entity", {}) or {}).get("pdbx_description")
        src = pe.get("rcsb_entity_source_organism", []) or [{}]
        organism = src[0].get("ncbi_scientific_name")
        genes = [g.get("value") for g in (src[0].get("rcsb_gene_name") or [])]
        chains = (
            pe.get("rcsb_polymer_entity_container_identifiers", {}).get(
                "auth_asym_ids"
            )
            or []
        )
        n_res = (pe.get("entity_poly", {}) or {}).get(
            "rcsb_sample_sequence_length"
        )
        role = classify_entity(desc, n_res)

        if organism:
            organisms.add(organism)
        if role == "mhc":
            hla_allele = hla_allele or normalize_allele(desc) or normalize_allele(
                " ".join(genes)
            )
        if role in ("tcr_alpha", "tcr_beta", "tcr"):
            tcr_genes.extend([g for g in genes if g])

        info["entities"].append({
            "entity_id": eid,
            "description": desc,
            "organism": organism,
            "genes": genes,
            "chains": chains,
            "role": role,
            "length": n_res,
        })

    info["hla_allele"] = hla_allele
    info["organisms"] = sorted(organisms)
    info["tcr_genes"] = sorted(set(tcr_genes))
    # Antigen / what the structure is about: use the title as the human label.
    info["antigen"] = info["title"]
    return info


def load_cache() -> dict:
    if config.RCSB_CACHE.exists():
        with config.RCSB_CACHE.open() as fh:
            return json.load(fh)
    return {}


def save_cache(cache: dict) -> None:
    config.WEB_DATA.mkdir(parents=True, exist_ok=True)
    tmp = config.RCSB_CACHE.with_suffix(".json.tmp")
    with tmp.open("w") as fh:
        json.dump(cache, fh, indent=2)
    tmp.replace(config.RCSB_CACHE)


def main(argv: list[str]) -> int:
    refresh = "--refresh" in argv
    cache = {} if refresh else load_cache()

    # Collect unique pdb ids from existing meta.json files (extract.py output).
    metas = sorted(config.WEB_DATA.glob("*/" + config.OUT_META))
    pdb_to_trajs: dict[str, list] = {}
    for mp in metas:
        with mp.open() as fh:
            meta = json.load(fh)
        pdb_to_trajs.setdefault(meta["pdb_id"], []).append(mp)

    print(f"{len(pdb_to_trajs)} unique PDB ids across {len(metas)} trajectories")

    # Fetch (with caching).
    for pdb_id in tqdm(sorted(pdb_to_trajs), unit="pdb"):
        if pdb_id not in cache:
            cache[pdb_id] = fetch_pdb(pdb_id)
            save_cache(cache)  # incremental, resumable

    # Merge annotation back into each meta.json.
    for pdb_id, meta_paths in pdb_to_trajs.items():
        rcsb = cache.get(pdb_id)
        for mp in meta_paths:
            with mp.open() as fh:
                meta = json.load(fh)
            meta["rcsb"] = rcsb
            with mp.open("w") as fh:
                json.dump(meta, fh, indent=2)

    n_found = sum(1 for v in cache.values() if v.get("found"))
    print(f"Done. {n_found}/{len(cache)} PDB entries found and merged.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
