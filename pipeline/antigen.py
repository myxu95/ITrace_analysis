"""Antigen-peptide annotation: name, source organism, and a coarse category.

The peptide presented in each pHLA-TCR complex is the antigen. RCSB already gives
us, per polymer entity, a free-text ``description`` and source ``organism``
(fetched by ``enrich_rcsb`` into ``meta.rcsb.entities``); the peptide entity
carries the antigen's name (e.g. "TAX PEPTIDE", "Epstein Barr Virus peptide") and
the organism it derives from (e.g. "Influenza A virus", "Homo sapiens"). This
module surfaces that as a compact ``meta.antigen`` block and derives a coarse
**category** (viral / bacterial / tumor-self / murine / synthetic / other) for
browsing and filtering.

It is metadata-only (reads ``meta.rcsb`` already on disk) — no network. Run after
``enrich_rcsb`` and before ``build_manifest``:

    IMMUNO_WEB_DATA=.../immuno-dyn python -m pipeline.antigen
"""
from __future__ import annotations

import json

from . import config

# Keyword -> category. Order matters: viral/bacterial are checked before the
# human/mouse hosts, because many viral organism names contain "human"
# ("Human immunodeficiency virus 1", "human gammaherpesvirus 4", ...).
# These run against BOTH the source organism AND the RCSB description name, because
# RCSB leaves the source organism empty for almost every short (synthetic/processed)
# peptide — the antigen identity then lives only in the description (e.g. "Epstein
# Barr Virus peptide", "TAX PEPTIDE", "HCMV pp65 fragment").
_VIRAL = ("virus", "viral", "hiv", "influenza", "coronavirus", "sars-cov", "sars ",
          "herpesvirus", "epstein", "barr", "hepatitis", "cytomegalovirus", "retrovirus",
          "lentivirus", "arenavirus", "papillomavirus", "orthohepacivirus",
          "choriomeningitis", "lcmv", "htlv", "tax peptide", "tax(", "tax ", "hcv ", "ebv ",
          "ebv,", "ebv-", "ebna", "cmv ", "hcmv", "pp65", "gpul40", "ul40", "bzlf1",
          "nuclear antigen", "nucleocapsid", "trans-activator")
_BACTERIAL = ("mycobacterium", "bacteri", "tuberculosis", "salmonella", "listeria",
              "escherichia", "streptococc", "staphylococc", "borrelia", "yersinia",
              "klebsiella", "pseudomonas", "clostrid")
# Tumour / self antigens named in the description (no source organism on RCSB).
_TUMOR = ("ny-eso", "cancer/testis", "cancer testis", "mart-1", "mart1", "melan",
          "melanoma", "gp100", "mage", "wt1", "telomerase", "hud", "tumor", "tumour")
# Designed / mimic peptides named in the description -> synthetic.
_SYNTH = ("mimotope", "synthetic", "construct")

# Display order / known labels (also the browse filter's option order).
CATEGORIES = ("viral", "bacterial", "tumor/self", "murine", "synthetic", "other")


def _host_category(host: str | None) -> str | None:
    """Last-resort bucket from the MHC host organism, for peptides RCSB never
    annotated (no source organism, no informative name). A peptide presented on a
    mouse MHC is most likely a murine self/allo peptide; on a human MHC, self/tumour."""
    h = (host or "").strip().lower()
    if "homo sapiens" in h:
        return "tumor/self"
    if "mus musculus" in h or "mouse" in h:
        return "murine"
    return None


def category(organism: str | None, name: str | None = None,
             host: str | None = None) -> str | None:
    """Coarse antigen-source category.

    Precedence: the description ``name`` wins first (its keyword is the strongest
    signal for these short peptides, and must beat a misleading peptide organism
    like "Homo sapiens" on an "EBV peptide"); then the peptide source organism;
    then a fallback to the MHC ``host`` organism so unannotated peptides still bucket."""
    n = (name or "").strip().lower()
    if n:
        if any(k in n for k in _VIRAL):
            return "viral"
        if any(k in n for k in _BACTERIAL):
            return "bacterial"
        if any(k in n for k in _TUMOR):
            return "tumor/self"
        if any(k in n for k in _SYNTH):
            return "synthetic"
        if "self" in n:
            return "tumor/self"
    o = (organism or "").strip().lower()
    if o:
        if "synthetic" in o or "construct" in o:
            return "synthetic"
        if any(k in o for k in _VIRAL):
            return "viral"
        if any(k in o for k in _BACTERIAL):
            return "bacterial"
        if "homo sapiens" in o:
            return "tumor/self"
        if "mus musculus" in o or "mouse" in o:
            return "murine"
        return "other"
    # No organism and no informative name: fall back to the MHC host organism.
    return _host_category(host)


def _peptide_entity(meta: dict) -> dict | None:
    for e in ((meta.get("rcsb") or {}).get("entities") or []):
        if e.get("role") == "peptide":
            return e
    return None


def _mhc_host(meta: dict) -> str | None:
    """Source organism of the MHC entity — the host whose cells present the peptide.
    Used as the last-resort category fallback for peptides RCSB never annotated."""
    for e in ((meta.get("rcsb") or {}).get("entities") or []):
        if e.get("role") == "mhc" and e.get("organism"):
            return e.get("organism")
    return None


def extract(meta: dict) -> dict:
    """``{name, organism, category}`` for one complex's antigen peptide.

    ``name`` is dropped when it is just the peptide sequence again (RCSB sometimes
    sets the description to the sequence itself), since the sequence is already shown.
    """
    pe = _peptide_entity(meta) or {}
    name = (pe.get("description") or "").strip() or None
    organism = (pe.get("organism") or "").strip() or None
    seq = (meta.get("peptide_seq") or "").strip()
    if name and seq and name.replace(" ", "").upper() == seq.upper():
        name = None
    return {"name": name, "organism": organism,
            "category": category(organism, name, _mhc_host(meta))}


def main(argv=None) -> int:
    web = config.WEB_DATA
    n = n_cat = 0
    for meta_path in sorted(web.glob("*/" + config.OUT_META)):
        meta = json.loads(meta_path.read_text())
        ag = extract(meta)
        meta["antigen"] = ag
        meta_path.write_text(json.dumps(meta, indent=2))
        n += 1
        if ag["category"]:
            n_cat += 1
    print(f"Antigen annotation: {n} complexes, {n_cat} with a source category.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
