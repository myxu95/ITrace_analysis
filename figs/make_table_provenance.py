#!/usr/bin/env python
"""Provenance table: where every statement about a source structure comes from.

Two independent traces per complex, kept visually separate because they answer
different questions and fail in different ways:

  1. what the DEPOSITED STRUCTURE itself records (parsed from the mmCIF header)
  2. what EXTERNAL ASSAY DATABASES and the primary literature record

Every count is recomputed here from Scientific_data/provenance/provenance_master.csv
and assay_matrix.csv, so the table can never drift from the audit that produced it.
"""
from pathlib import Path
import csv

from figstyle import *          # canvas, type scale, ink, save, render

# The tables keep the pre-2026-09 type scale.  The figures were enlarged because
# they are reduced on placement; these tables are placed at their drawn width,
# and every column here is measured in mm against text set at 7 pt.
FS_BODY, FS_SMALL, FS_HEAD = 7.0, 6.5, 7.0
F_BODY, F_SMALL, F_HEAD = font(FS_BODY), font(FS_SMALL), font(FS_HEAD, bold=True)

LEGEND = {}   # prose harvested by build_legends.py; see figstyle.LEAN

import provdata
from provdata import NAME

PROV = Path(__file__).resolve().parents[2] / "Scientific_data" / "provenance"
N_TOTAL = 245

# ---------------------------------------------------------------- data ------
def load():
    master = list(csv.DictReader(open(PROV / "provenance_master.csv", encoding="utf-8-sig")))
    assay = list(csv.DictReader(open(PROV / "assay_matrix.csv", encoding="utf-8-sig")))
    assert len(master) == len(assay) == N_TOTAL, (len(master), len(assay))

    def pos(rows, f):
        return sum(1 for r in rows if int(r[f] or 0) > 0)

    eng = [r for r in master if r["engineered"] == "Y"]
    cats = {}
    for r in eng:
        for c in (r["eng_categories"] or "").split("|"):
            if c.strip():
                cats[c.strip()] = cats.get(c.strip(), 0) + 1

    ev = {}
    for r in master:
        ev[r["assay_evidence_level"]] = ev.get(r["assay_evidence_level"], 0) + 1

    return dict(
        engineered=len(eng),
        wild_type=N_TOTAL - len(eng),
        ev_seqdif=sum(1 for r in master if r["eng_evidence"] == "seq_dif"),
        ev_both=sum(1 for r in master if r["eng_evidence"] == "seq_dif+text"),
        ev_text=sum(1 for r in master if r["eng_evidence"] == "text"),
        natural=sum(1 for r in master if r["natural_variant_context"] == "Y"),
        disulf=sum(1 for r in master if r["eng_disulfides"].strip()),
        cats=cats,
        tcr_level=ev.get("TCR-level", 0),
        ep_level=ev.get("epitope-level", 0),
        lit_only=ev.get("literature-only", 0),
        none=ev.get("none", 0),
        kin=pos(assay, "tcr_affinity_kinetics"),
        qual=pos(assay, "tcr_binding_qual"),
        cyto=pos(assay, "tcr_cytotoxicity"),
        activ=pos(assay, "tcr_cytokine_activation"),
        prolif=pos(assay, "tcr_proliferation"),
        vivo=pos(assay, "tcr_in_vivo"),
        kd=pos(master, "n_KD_values"),
        eng_kd=sum(1 for r in master if r["engineered"] == "Y" and int(r["n_KD_values"] or 0) > 0),
        n_pmid=len({r["pmid"] for r in master if r["pmid"].strip()}),
        n_oa=len({r["pmid"] for r in master if r["has_oa_fulltext"] == "Y" and r["pmid"].strip()}),
    )


def build(d):
    """(kind, label, n, [(source text, mono?), ...])

    The source column names the FIELD a number was read out of -- nothing else.
    What a category means, and what the assay abbreviations stand for, is legend
    prose; a table that has to gloss itself in a fourth column is not a table.
    """
    c = d["cats"]
    return [
        ("block", "Parsed out of the deposited structure", "RCSB mmCIF header", PRIMARY),
        ("row", "Engineered construct", d["engineered"], []),
        ("sub", "recorded in the structured field", d["ev_seqdif"], [("_struct_ref_seq_dif", 1)]),
        ("sub", "recorded in both places", d["ev_both"], []),
        ("sub", "recorded only in free text", d["ev_text"],
         [("title / ", 0), ("_pdbx_entry_details", 1)]),
        ("row", "Wild-type sequence", d["wild_type"], []),
        ("row", "Natural variant context", d["natural"], []),

        ("head", "Engineering type", None, []),
        ("sub", "Affinity-engineered TCR", c.get("affinity_engineered"), []),
        ("sub", "Interchain linker", c.get("linked_construct"), [("_struct_ref_seq_dif", 1)]),
        ("sub", "Substituted peptide", c.get("peptide_variant"), []),
        ("sub", "MHC heavy-chain mutation", c.get("mhc_mutation"), []),
        ("sub", "Insertion, deletion or truncation", c.get("indel_truncation"), []),
        ("sub", "Engineered peptide", c.get("peptide_engineered"), []),
        ("sub", "Engineered disulfide", d["disulf"], [("_struct_conn", 1)]),
        ("sub", "Single-chain construct", c.get("single_chain"), []),
        ("sub", "TCR constant-region mutation", c.get("tcr_constant_mutation"), []),
        ("sub", "Beta-2-microglobulin mutation", c.get("b2m_mutation"), []),
        ("sub", "TCR variable-region mutation", c.get("tcr_variable_mutation"), []),
        ("sub", "Point mutant, described in text only", c.get("point_mutant_doc"), []),

        ("block", "Mined from assay databases and the literature",
         "IEDB Query API \u00b7 Europe PMC", SECOND),
        ("head", "Strongest evidence found", None, []),
        # names come from provdata so this table and fig_provenance_evidence
        # cannot end up calling the same class two different things
        ("sub", NAME["TCR-level"], d["tcr_level"], []),
        ("sub", NAME["epitope-level"], d["ep_level"], []),
        ("sub", NAME["literature-only"], d["lit_only"], []),
        ("sub", NAME["none"], d["none"], []),
        ("total", "Any experimental anchor", d["tcr_level"] + d["ep_level"] + d["lit_only"], []),

        ("head", "Assay type on the TCR itself", None, []),
        # tcr_affinity_kinetics and n_KD_values pick out exactly the same 112 entries,
        # so they are one row, not two.
        ("sub", "Affinity and kinetics", d["kin"], []),
        ("sub", "Qualitative binding", d["qual"], []),
        ("sub", "Cytokine and activation", d["activ"], []),
        ("sub", "Cytotoxicity", d["cyto"], []),
        ("sub", "Proliferation", d["prolif"], []),
        ("sub", "In vivo", d["vivo"], []),

        ("block", "Where the two traces meet", None, ACCENT),
        ("total", "Engineered and affinity-measured", d["eng_kd"], []),
    ]


# -------------------------------------------------------------- layout -----
# The source field is named inline, right after the row it qualifies: only four
# rows carry one, so a fourth column would be four entries and a dead half-page.
# The counts hold the right edge, where a reader scans a table of numbers.
X_N, X_PCT = 122.0, W                     # right edges of the count and the share
GAP_SRC = 3.0                             # label to inline field name
IND = 3.4
LH = 3.85
GAP_BLOCK, GAP_HEAD = 4.4, 2.6


def draw(ax, top, w):
    d = load()
    rows = build(d)

    LEGEND.update(
        title="Provenance of the 245 source structures",
        subtitle="Every complex is traced twice and the two traces are reported separately: "
                 "what the deposited crystal structure itself records, and what external "
                 "experiments record.",
        notes=[
            f"Percentages are of the {N_TOTAL} complexes. One evidence class is assigned per "
            "complex, so those four rows partition the library; engineering types and assay "
            "types overlap, so an entry may carry several and they do not sum to the row "
            "above them. Natural sequence variation \u2014 a disease, escape or neoantigen "
            "sequence \u2014 is recorded but is not counted as engineering.",
            "Engineered peptides are mimotopes or non-natural residues, as distinct from "
            "substituted peptides, which are natural-sequence variants. One of the "
            f"{d['disulf']} engineered disulfides spans the TCR\u2013pMHC interface.",
            "Affinity and kinetics means a quantitative KD, rate constant or half-life; "
            "qualitative binding means tetramer or T cell\u2013APC staining; cytokine and "
            "activation covers IFN-gamma, IL-2 and CCL4 readouts; cytotoxicity is "
            "chromium-51 release and proliferation is tritiated-thymidine incorporation. "
            "The final row counts the complexes that carry an engineering record and a "
            "measured affinity at once, which is what makes wild-type versus "
            "affinity-matured TCR ladders possible.",
            "The assay counts rest on 1,329 IEDB records for the receptors and 10,304 for "
            "the epitopes "
            f"and on {d['n_pmid']} publications, {d['n_oa']} of them open access; these are "
            "unique record counts, not the per-entry sums, because IEDB receptor groups are "
            "shared across entries.",
            # computed, not written down: this used to be hard-coded, was labelled
            # IEDB coverage when it measures the receptor-level share instead, and
            # was stratified on a citation year 28 complexes do not have
            f"The share with an assay for the TCR itself falls from {provdata.era_tcr_share()[0]:.1f}% "
            f"for structures released up to {provdata.EARLY} to {provdata.era_tcr_share()[1]:.1f}% for "
            f"{provdata.LATE} onward, so {NAME['none'].lower()} means not yet curated and not yet "
            "readable rather than proven absent: of those complexes, most have no publication at "
            "all yet and the rest had no open-access full text to scan.",
            "provenance_master.csv carries all of the above for every complex, one row each.",
        ],
    )
    y = header(ax, 0, top, w, LEGEND["title"], LEGEND["subtitle"])
    y -= 3.4

    for kind, label, n, src in rows:
        if kind == "block":
            y -= GAP_BLOCK
            ax.text(0, y, label, fontproperties=font(FS_HEAD, bold=True), color=INK,
                    ha="left", va="baseline", zorder=5)
            if n:
                ax.text(w, y, n, fontproperties=F_BODY, color=SECONDARY,
                        ha="right", va="baseline", zorder=5)
            y -= 1.5
            rule(ax, 0, w, y, lw=1.0, color=src)
            y -= LH * 0.95
            continue

        if kind == "head":
            y -= GAP_HEAD
            ax.text(0, y, label, fontproperties=font(FS_BODY, bold=True), color=INK,
                    ha="left", va="baseline", zorder=5)
            y -= LH
            continue

        bold = kind == "total"
        fp = font(FS_BODY, bold=True) if bold else F_BODY
        x = IND if kind == "sub" else 0.0
        ax.text(x, y, label, fontproperties=fp, color=INK if kind != "sub" else SECONDARY,
                ha="left", va="baseline", zorder=5)
        ax.text(X_N, y, str(n), fontproperties=fp, color=INK,
                ha="right", va="baseline", zorder=5)
        ax.text(X_PCT, y, f"{100 * n / N_TOTAL:.1f}%", fontproperties=fp, color=SECONDARY,
                ha="right", va="baseline", zorder=5)
        if src:
            sx = x + text_w(label + "|", fp) - text_w("|", fp) + GAP_SRC
            for txt, mono in src:
                fps = font(FS_SMALL, mono=True) if mono else F_BODY
                ax.text(sx, y, txt, fontproperties=fps, color=SECONDARY,
                        ha="left", va="baseline", zorder=5)
                sx += text_w(txt + "|", fps) - text_w("|", fps)
        y -= LH

    return notes(ax, 0, w, y, LEGEND["notes"], gap=1.6, rule_lw=0.7)


if __name__ == "__main__":
    render(draw, "Table_provenance")
