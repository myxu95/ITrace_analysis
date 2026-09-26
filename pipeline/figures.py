"""NAR publication main figures (bravo §10), generated from the aggregated data.

Six cross-system figures rendered with matplotlib (needs the `imscope` env, which
has matplotlib). Each function reads aggregate.json / manifest.json / per-system
analysis.json under the default web-data path and saves a PNG into <web_data>/figures/.

    /home/xmy/miniforge3/envs/imscope/bin/python -m pipeline.figures
"""
from __future__ import annotations

import os

from . import config


def fig1_overview(out_dir):
    import json
    import os
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from collections import Counter

    # ---- palette ----
    PEP = "#B5838D"
    TCR = "#6D8A96"
    HLA = "#A08E7B"
    STABLE = "#84A59D"
    FLEX = "#CB997E"
    NEUTRAL = "#A3A7AB"

    data_dir = "/home/xmy/work/data/immunotrace/web_data_1000"
    agg = json.load(open(os.path.join(data_dir, "aggregate.json")))
    man = json.load(open(os.path.join(data_dir, "manifest.json")))
    trajs = man["trajectories"]

    n_complexes = agg["n"]
    # Panels a & b describe the HLA dataset, so they are computed over HUMAN
    # complexes only (HLA is human MHC) and PER COMPLEX (not per trajectory), so
    # both totals equal the number of human complexes.
    human = [c for c in man["complexes"] if c.get("is_human")]
    n_human = len(human)

    # ---- panel a: peptide-length distribution (human complexes) ----
    plen = Counter(c["peptide_length"] for c in human if c.get("peptide_length"))
    lengths = sorted(plen.keys(), key=lambda x: int(x))
    lcounts = [plen[k] for k in lengths]
    n_lengths = len(lengths)

    # ---- panel b: HLA alleles collapsed to gene (locus) level (human complexes) ----
    def hla_gene(a):
        if not a:
            return None
        return a.split("*")[0]

    genes = Counter()
    for c in human:
        g = hla_gene(c.get("hla_allele"))
        if g:
            genes[g] += 1
    distinct_alleles = len({c.get("hla_allele") for c in human if c.get("hla_allele")})
    top_genes = genes.most_common(10)
    glabels = [g for g, _ in top_genes][::-1]
    gvals = [c for _, c in top_genes][::-1]

    # ---- panel c: host species ----
    host = agg["global"]["host"]
    host_items = sorted(host.items(), key=lambda kv: -kv[1])
    hlabels = [k for k, _ in host_items]
    hvals = [v for _, v in host_items]

    # ---- figure ----
    fig = plt.figure(figsize=(9.5, 7.0))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.32,
                          left=0.085, right=0.965, top=0.93, bottom=0.09)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, 0])
    axd = fig.add_subplot(gs[1, 1])

    def clean(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # --- (a) peptide length ---
    xa = np.arange(len(lengths))
    barsa = axa.bar(xa, lcounts, color=PEP, edgecolor="white", linewidth=0.7, width=0.74)
    axa.set_xticks(xa)
    axa.set_xticklabels(lengths)
    axa.set_xlabel("Peptide length (residues)", fontsize=10)
    axa.set_ylabel("Complexes", fontsize=10)
    axa.set_title(f"Peptide-length distribution  ·  {sum(lcounts)} human complexes", fontsize=9.5)
    axa.set_ylim(0, max(lcounts) * 1.16)
    axa.yaxis.grid(True, alpha=0.25)
    axa.set_axisbelow(True)
    for rect, v in zip(barsa, lcounts):
        axa.text(rect.get_x() + rect.get_width() / 2, v + max(lcounts) * 0.02,
                 str(v), ha="center", va="bottom", fontsize=8.5, color="#3a3a3a")
    clean(axa)

    # --- (b) HLA genes ---
    yb = np.arange(len(glabels))
    barsb = axb.barh(yb, gvals, color=HLA, edgecolor="white", linewidth=0.7, height=0.72)
    axb.set_yticks(yb)
    axb.set_yticklabels(glabels)
    axb.set_xlabel("Complexes", fontsize=10)
    axb.set_title(f"HLA gene usage  ·  {sum(gvals)} human complexes", fontsize=9.5)
    axb.set_xlim(0, max(gvals) * 1.16)
    axb.xaxis.grid(True, alpha=0.25)
    axb.set_axisbelow(True)
    for rect, v in zip(barsb, gvals):
        axb.text(v + max(gvals) * 0.015, rect.get_y() + rect.get_height() / 2,
                 str(v), ha="left", va="center", fontsize=8.5, color="#3a3a3a")
    clean(axb)

    # --- (c) host species (donut) ---
    palette_c = [TCR, FLEX, NEUTRAL, STABLE]
    cols = [palette_c[i % len(palette_c)] for i in range(len(hlabels))]
    wedges, _ = axc.pie(hvals, colors=cols, startangle=90,
                        counterclock=False,
                        wedgeprops=dict(width=0.42, edgecolor="white", linewidth=1.4))
    axc.set_aspect("equal")
    axc.set_title("Host species", fontsize=10.5)
    total_host = sum(hvals)
    short = {"Homo sapiens": "Homo sapiens", "Mus musculus": "Mus musculus",
             "Macaca mulatta": "Macaca mulatta"}
    leg_labels = ["{}  ({}, {:.0f}%)".format(short.get(l, l), v, 100.0 * v / total_host)
                  for l, v in zip(hlabels, hvals)]
    axc.legend(wedges, leg_labels, loc="center", frameon=False, fontsize=8.0,
               bbox_to_anchor=(0.5, -0.18), handlelength=1.0, handletextpad=0.5)
    axc.text(0, 0, "{}\nsystems".format(total_host), ha="center", va="center",
             fontsize=11, fontweight="bold", color="#3a3a3a")

    # --- (d) summary big-number panel ---
    axd.axis("off")
    n_human = host.get("Homo sapiens", 0)
    pct_human = 100.0 * n_human / n_complexes

    metrics = [
        (str(n_complexes), "pHLA-TCR complexes", PEP),
        (str(len(genes)) + " / " + str(distinct_alleles), "HLA genes / distinct alleles", HLA),
        (str(n_lengths), "peptide lengths (8-13 aa)", TCR),
        ("200 ns", "per-system MD trajectory", STABLE),
    ]
    positions = [(0.02, 0.74), (0.52, 0.74), (0.02, 0.30), (0.52, 0.30)]
    for (big, lab, col), (x, y) in zip(metrics, positions):
        axd.text(x, y, big, transform=axd.transAxes, fontsize=21,
                 fontweight="bold", color=col, ha="left", va="center")
        axd.text(x, y - 0.135, lab, transform=axd.transAxes, fontsize=8.6,
                 color="#4a4a4a", ha="left", va="center")
    _scope = ("Force field CHARMM36m  |  explicit-solvent MD  |  "
              "{:.0f}% human".format(pct_human))
    axd.text(0.02, 0.04, _scope,
             transform=axd.transAxes, fontsize=7.6, color="#6a6a6a",
             ha="left", va="center", style="italic")
    axd.set_title("Resource scope", fontsize=10.5, loc="left")
    axd.add_patch(plt.Rectangle((0.0, 0.0), 1.0, 0.92, transform=axd.transAxes,
                                fill=False, edgecolor=NEUTRAL, linewidth=1.0, alpha=0.55))

    # --- panel letters ---
    for ax, lab in [(axa, "a"), (axb, "b"), (axc, "c"), (axd, "d")]:
        ax.annotate(lab, xy=(0.0, 1.0), xycoords="axes fraction",
                    xytext=(-38, 16), textcoords="offset points",
                    fontsize=13, fontweight="bold", va="top", ha="left")

    fig.suptitle("Figure 1 — ImmunoTrace database overview",
                 fontsize=12.5, fontweight="bold", y=0.985)

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "fig1_overview.png")
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path

def fig2_qc(out_dir):
    import json, glob, os
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ---- Morandi palette ----
    C_TCR = "#6D8A96"
    C_STABLE = "#84A59D"   # ok / stable
    C_FLEX = "#CB997E"     # drifting / dynamic
    C_NEUTRAL = "#A3A7AB"
    C_OUTLIER = "#9B5C5C"  # darker red-mauve for outliers

    DATA = "/home/xmy/work/data/immunotrace/web_data_1000"

    # ---- (a) equilibration breakdown from aggregate ----
    eq_order = ["ok", "drifting", "outlier"]
    eq_labels = ["OK", "Drifting", "Outlier"]
    eq_colors = [C_STABLE, C_FLEX, C_OUTLIER]

    # ---- (b) tail-90% mean RMSD distribution (nm -> Angstrom) ----
    rmsd_vals, rmsd_eq = [], []
    for f in sorted(glob.glob(os.path.join(DATA, "*/meta.json"))):
        try:
            q = json.load(open(f)).get("quality", {})
            v = q.get("tail90_mean_rmsd_nm")
            if v is not None:
                rmsd_vals.append(v * 10.0)
                rmsd_eq.append(q.get("equilibration", "ok"))
        except Exception:
            pass
    rmsd_vals = np.array(rmsd_vals)
    rmsd_eq = np.array(rmsd_eq)

    # Panel (a) counts equilibration PER TRAJECTORY (like panels b/c), so all three
    # QC panels share n = number of trajectories rather than mixing a per-complex bar.
    eq_counts = [int((rmsd_eq == k).sum()) for k in eq_order]
    eq_total = len(rmsd_eq)

    # ---- (c) mean RMSF distribution (Angstrom) ----
    rmsf_vals = []
    for f in sorted(glob.glob(os.path.join(DATA, "*/analysis/analysis.json"))):
        try:
            v = json.load(open(f)).get("rmsf", {}).get("mean_rmsf_angstrom")
            if v is not None:
                rmsf_vals.append(v)
        except Exception:
            pass
    rmsf_vals = np.array(rmsf_vals)

    # =====================================================================
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11,
                         "axes.labelsize": 10, "font.family": "sans-serif"})
    fig = plt.figure(figsize=(9.5, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0],
                          hspace=0.42, wspace=0.30)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    def clean(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    def panel_letter(ax, letter):
        ax.text(-0.16, 1.06, letter, transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="top", ha="left")

    # ---------- (a) equilibration ----------
    xpos = np.arange(len(eq_order))
    ax_a.bar(xpos, eq_counts, color=eq_colors, width=0.62,
             edgecolor="white", linewidth=0.8)
    for x, c in zip(xpos, eq_counts):
        pct = 100.0 * c / eq_total if eq_total else 0
        ax_a.text(x, c + max(eq_counts) * 0.02, f"{c}\n({pct:.0f}%)",
                  ha="center", va="bottom", fontsize=8.5, linespacing=1.0)
    ax_a.set_xticks(xpos)
    ax_a.set_xticklabels(eq_labels)
    ax_a.set_ylabel("Trajectories")
    ax_a.set_title("Equilibration status", fontsize=10.5)
    ax_a.set_ylim(0, max(eq_counts) * 1.22)
    ax_a.grid(axis="y", alpha=0.25)
    ax_a.set_axisbelow(True)
    clean(ax_a)
    panel_letter(ax_a, "a")
    ax_a.text(0.97, 0.92, f"n = {eq_total}", transform=ax_a.transAxes,
              ha="right", va="top", fontsize=8.5, color=C_NEUTRAL)

    # ---------- (b) tail-90% mean RMSD ----------
    bins = np.linspace(rmsd_vals.min(), rmsd_vals.max(), 28)
    ok_mask = rmsd_eq == "ok"
    drift_mask = rmsd_eq == "drifting"
    out_mask = rmsd_eq == "outlier"
    ax_b.hist([rmsd_vals[ok_mask], rmsd_vals[drift_mask], rmsd_vals[out_mask]],
              bins=bins, stacked=True,
              color=[C_STABLE, C_FLEX, C_OUTLIER],
              edgecolor="white", linewidth=0.3,
              label=["OK", "Drifting", "Outlier"])
    med = np.median(rmsd_vals)
    ax_b.axvline(med, color=C_NEUTRAL, ls="--", lw=1.2)
    ax_b.text(med, ax_b.get_ylim()[1] * 0.70, f" median\n {med:.1f} Å",
              ha="left", va="top", fontsize=8.2, color="#555555",
              linespacing=1.05)
    n_out = int(out_mask.sum())
    ax_b.annotate(f"{n_out} outliers\n> 6 Å",
                  xy=(rmsd_vals[out_mask].mean(), 1.5),
                  xytext=(rmsd_vals.max() * 0.62, max(eq_counts) * 0.10 + 4),
                  fontsize=8.2, color=C_OUTLIER, ha="center",
                  arrowprops=dict(arrowstyle="->", color=C_OUTLIER, lw=1.0))
    ax_b.set_xlabel("Tail-90% mean backbone RMSD (Å)")
    ax_b.set_ylabel("Trajectories")
    ax_b.set_title("Backbone stability", fontsize=10.5)
    ax_b.legend(frameon=False, fontsize=8, loc="upper right",
                handlelength=1.1, borderpad=0.2)
    ax_b.grid(axis="y", alpha=0.25)
    ax_b.set_axisbelow(True)
    clean(ax_b)
    panel_letter(ax_b, "b")

    # ---------- (c) mean RMSF ----------
    bins_c = np.linspace(rmsf_vals.min(), rmsf_vals.max(), 36)
    ax_c.hist(rmsf_vals, bins=bins_c, color=C_TCR,
              edgecolor="white", linewidth=0.3, alpha=0.92)
    med_f = np.median(rmsf_vals)
    mean_f = np.mean(rmsf_vals)
    ax_c.axvline(med_f, color=C_NEUTRAL, ls="--", lw=1.2)
    ax_c.text(med_f, ax_c.get_ylim()[1] * 0.95,
              f" median {med_f:.2f} Å   (mean {mean_f:.2f} Å)",
              ha="left", va="top", fontsize=8.6, color="#555555")
    n_flex = int((rmsf_vals > 4.0).sum())
    ax_c.axvspan(4.0, rmsf_vals.max(), color=C_FLEX, alpha=0.12)
    ax_c.text(4.1, ax_c.get_ylim()[1] * 0.55,
              f"flexible tail\n{n_flex} systems > 4 Å",
              fontsize=8.2, color=C_FLEX, va="top")
    ax_c.set_xlabel("Per-system mean Cα RMSF (Å)")
    ax_c.set_ylabel("Trajectories")
    ax_c.set_title("Residue-level flexibility", fontsize=10.5)
    ax_c.grid(axis="y", alpha=0.25)
    ax_c.set_axisbelow(True)
    clean(ax_c)
    panel_letter(ax_c, "c")
    ax_c.text(0.99, 0.92, f"n = {len(rmsf_vals)} systems",
              transform=ax_c.transAxes, ha="right", va="top",
              fontsize=8.5, color=C_NEUTRAL)

    fig.suptitle("Figure 2 — Trajectory quality control",
                 fontsize=12.5, fontweight="bold", y=1.005)

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "fig2_qc.png")
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path

def fig3_docking(out_dir):
    import json, glob, os
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    # ---- Morandi palette ----
    PEP = "#B5838D"; TCR = "#6D8A96"; HLA = "#A08E7B"
    STABLE = "#84A59D"; FLEX = "#CB997E"; NEUTRAL = "#A3A7AB"

    DATA = "/home/xmy/work/data/immunotrace/web_data_1000"
    agg = json.load(open(os.path.join(DATA, "aggregate.json")))
    g = agg["global"]

    # ---- per-COMPLEX crossing vs incident (panel b); replicas averaged so the
    # scatter counts each complex once, consistent with the per-complex panel a ----
    man = json.load(open(os.path.join(DATA, "manifest.json")))
    cr, inc, rv = [], [], []
    for cx in man.get("complexes", []):
        cs, ins, revs = [], [], []
        for r in cx.get("replicas", []):
            ap = os.path.join(DATA, r["traj_id"], "analysis", "analysis.json")
            if not os.path.exists(ap):
                continue
            a = json.load(open(ap)).get("angle", {})
            c = a.get("crossing_deg", {}).get("mean")
            i = a.get("incident_deg", {}).get("mean")
            if c is not None:
                cs.append(c)
            if i is not None:
                ins.append(i)
            revs.append(bool(a.get("reversed_polarity")))
        if cs and ins:
            cm = float(np.mean(cs))
            reversed_complex = sum(revs) > len(revs) / 2 if revs else False
            cr.append(cm); inc.append(float(np.mean(ins)))
            rv.append(reversed_complex)
    cr = np.array(cr); inc = np.array(inc)
    rv = np.array(rv, dtype=bool)

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11,
                         "axes.labelsize": 10, "font.family": "sans-serif"})
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.0))
    (axA, axB), (axC, axD) = axes

    def despine(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    def panel_label(ax, letter):
        ax.text(-0.13, 1.06, letter, transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="top", ha="left")

    # =================== (a) crossing-angle histogram ===================
    vals = np.array(g["crossing_values"], dtype=float)
    rev_vals = np.array(g.get("crossing_values_reversed", []), dtype=float)
    bins = np.arange(0, 185, 5)   # directed crossing spans 0-180 (>90 = reverse polarity)
    axA.hist(vals, bins=bins, color=TCR, edgecolor="white", linewidth=0.6,
             label="forward polarity")
    axA.set_xlim(0, 180)
    if rev_vals.size:   # overlay the (few) reverse-polarity systems distinctly
        axA.hist(rev_vals, bins=bins, color=FLEX, edgecolor="white", linewidth=0.6,
                 label=f"reverse-polarity (n={rev_vals.size})", zorder=4)
    med = g["crossing"]["median"]
    axA.axvline(med, color=PEP, lw=1.8, ls="--", label=f"median {med:.1f}°")
    axA.set_xlabel("TCR–pMHC crossing angle (°, directed 0–180; >90° = reverse polarity)")
    axA.set_ylabel("number of systems")
    axA.set_title(f"Crossing-angle distribution (n = {len(vals)})")
    axA.legend(frameon=False, fontsize=8.5, loc="upper right")
    axA.grid(axis="y", alpha=0.25)
    despine(axA); panel_label(axA, "a")

    # =================== (b) crossing vs incident scatter ===================
    fwd = ~rv
    axB.scatter(cr[fwd], inc[fwd], s=20, c=TCR, alpha=0.8,
                edgecolors="white", linewidths=0.3,
                label=f"forward polarity (n={int(fwd.sum())})")
    axB.scatter(cr[rv], inc[rv], s=30, c=FLEX, alpha=0.9, marker="D",
                edgecolors="white", linewidths=0.4,
                label=f"reverse-polarity (n={int(rv.sum())})")
    axB.set_xlim(0, 180)
    axB.set_xlabel("crossing angle (°, directed 0–180)")
    axB.set_ylabel("incident angle (°)")
    axB.set_title("Docking geometry: crossing vs incident")
    axB.legend(frameon=False, fontsize=8, loc="upper right")
    axB.grid(alpha=0.2)
    despine(axB); panel_label(axB, "b")

    # =================== shared bar helper (c, d) ===================
    def grouped_bar(ax, dmap, order, color, xlabel):
        order = [k for k in order if k in dmap]
        meds = [dmap[k]["median"] for k in order]
        p25 = [dmap[k]["p25"] for k in order]
        p75 = [dmap[k]["p75"] for k in order]
        ns = [dmap[k]["n"] for k in order]
        lo = [max(0.0, m - a) for m, a in zip(meds, p25)]
        hi = [b - m for m, b in zip(meds, p75)]
        x = np.arange(len(order))
        ax.bar(x, meds, color=color, edgecolor="white", width=0.66, zorder=2)
        ax.errorbar(x, meds, yerr=[lo, hi], fmt="none", ecolor="#444444",
                    elinewidth=1.1, capsize=4, capthick=1.1, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels(order, rotation=0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("crossing angle (°)")
        ymax = max(p75) * 1.16
        ax.set_ylim(0, ymax)
        for xi, m, nn in zip(x, meds, ns):
            ax.text(xi, ymax * 0.035, f"n={nn}", ha="center", va="bottom",
                    fontsize=7.5, color="white" if m > ymax * 0.12 else "#333333")
        ax.grid(axis="y", alpha=0.25)
        despine(ax)

    # =================== (c) crossing by HLA group ===================
    hla_map = agg["crossing_by_hla"]
    hla_order = sorted(hla_map.keys(), key=lambda k: -hla_map[k]["n"])
    grouped_bar(axC, hla_map, hla_order, HLA, "HLA locus")
    axC.set_title("Crossing angle by HLA restriction")
    axC.tick_params(axis="x", labelsize=8.5)
    panel_label(axC, "c")

    # =================== (d) crossing by peptide length ===================
    len_map = agg["crossing_by_length"]
    len_order = sorted(len_map.keys(), key=lambda k: int(k))
    grouped_bar(axD, len_map, len_order, PEP, "peptide length (residues)")
    axD.set_title("Crossing angle by peptide length")
    panel_label(axD, "d")

    fig.suptitle("Figure 3 — TCR docking-geometry landscape",
                 fontsize=12.5, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])

    out = os.path.join(out_dir, "fig3_docking.png")
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out

def fig4_recognition(out_dir):
    import json, os
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ---- Morandi palette ----
    C_PEP = "#B5838D"   # peptide
    C_TCR = "#6D8A96"   # tcr / cdr3
    C_HLA = "#A08E7B"   # hla
    C_STABLE = "#84A59D"
    C_FLEX = "#CB997E"
    C_NEU = "#A3A7AB"

    agg = "/home/xmy/work/data/immunotrace/web_data_1000/aggregate.json"
    with open(agg) as fh:
        d = json.load(fh)

    # ---------- figure scaffold ----------
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11,
                         "axes.labelsize": 10, "figure.facecolor": "white",
                         "axes.facecolor": "white"})
    fig = plt.figure(figsize=(9.5, 6.6), dpi=200)
    # wspace is wide because panel (b) carries a colorbar on its right edge and
    # panel (c) carries long y tick labels on its left edge -- they collide below ~0.5.
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0],
                          hspace=0.42, wspace=0.52)
    axa = fig.add_subplot(gs[0, :])   # panel a spans top row
    axb = fig.add_subplot(gs[1, 0])
    axc = fig.add_subplot(gs[1, 1])

    def strip(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    def panel_label(ax, txt, dx=-0.06, dy=1.04):
        ax.text(dx, dy, txt, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="top", ha="right")

    # ================= panel (a) =================
    cp = d["conserved_peptide"]
    sm = cp.get("sasa_mean", {})
    pos = [p for p in cp["positions"] if str(p) in cp["tcr_mean"] and str(p) in sm]
    tcr = np.array([cp["tcr_mean"][str(p)] for p in pos])
    sasa = np.array([sm[str(p)] * 100.0 for p in pos])  # nm^2 -> A^2
    x = np.arange(len(pos))
    axa.bar(x, tcr, 0.6, color=C_TCR, label="TCR contact",
            edgecolor="white", linewidth=0.6)
    axa.set_xticks(x)
    axa.set_xticklabels(["P%d" % p for p in pos])
    axa.set_xlabel("Peptide position (9-mer)")
    axa.set_ylabel("Mean TCR contact occupancy")
    axa.set_ylim(0, 1.30)
    axa.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axa.set_title("Conserved peptide-position recognition map", pad=6)
    axa.grid(axis="y", alpha=0.25)
    # Buried SASA on a twin axis (low = buried in an HLA pocket = anchor),
    # inverted so buried reads upward. Replaces the saturated max-to-any
    # HLA-contact axis, which was ~0.9 at every position and hid the anchors.
    axa2 = axa.twinx()
    axa2.plot(x, sasa, "-o", color=C_HLA, lw=1.8, ms=5, label="buried SASA")
    axa2.set_ylabel("← buried   SASA (Å²)   exposed →", fontsize=8)
    axa2.invert_yaxis()
    axa2.spines["top"].set_visible(False)
    # mark canonical anchor positions (P2 + P-Omega) at their SASA minima
    for xi, p, s in zip(x, pos, sasa):
        if p in (2, len(pos)):
            axa2.annotate("anchor", xy=(xi, s), xytext=(0, -14),
                          textcoords="offset points", ha="center", va="top",
                          fontsize=7.5, color=C_HLA, annotation_clip=False,
                          arrowprops=dict(arrowstyle="-|>", color=C_HLA, lw=1.0))
    h1, l1 = axa.get_legend_handles_labels()
    h2, l2 = axa2.get_legend_handles_labels()
    # legend centred at the top: the P2 / PΩ "anchor" callouts sit at the left and
    # right edges (buried-SASA minima), so a centred legend clears both. The raised
    # ylim (1.30) opens a clear band above the bars for it.
    axa.legend(h1 + h2, l1 + l2, frameon=False, loc="upper center", ncol=2,
               bbox_to_anchor=(0.5, 1.0), fontsize=8.5)
    # annotate central recognition bulge region
    axa.axvspan(2.5, 7.5, color=C_PEP, alpha=0.07, zorder=0)
    axa.text(5.0, 0.04, "TCR-read central bulge", ha="center",
             va="bottom", fontsize=8, color=C_PEP, style="italic")
    strip(axa)
    panel_label(axa, "a", dx=-0.035)

    # ================= panel (b) =================
    g = d["global"]
    pr = np.array(g["pep_recog_values"], dtype=float)   # CDR3 -> peptide
    hr = np.array(g["hla_restr_values"], dtype=float)    # CDR1/2 -> HLA
    n = len(pr)
    hb = axb.hexbin(pr, hr, gridsize=20, cmap="BuPu", mincnt=1,
                    linewidths=0.2, edgecolors="white")
    cb = fig.colorbar(hb, ax=axb, fraction=0.038, pad=0.02)
    cb.set_label("systems", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    # reference diagonal (equal engagement)
    lim = max(pr.max(), hr.max()) * 1.05
    axb.plot([0, lim], [0, lim], ls="--", lw=1.0, color=C_NEU, zorder=3)
    # median crosshair
    axb.axvline(np.median(pr), color=C_TCR, lw=1.0, ls=":", alpha=0.8)
    axb.axhline(np.median(hr), color=C_HLA, lw=1.0, ls=":", alpha=0.8)
    axb.set_xlim(0, lim)
    axb.set_ylim(0, lim)
    axb.set_xlabel("CDR3 → peptide recognition")
    axb.set_ylabel("CDR1/2 → HLA restriction")
    axb.set_title("Peptide vs HLA engagement per system", pad=6)
    axb.text(0.97, 0.06, "n = %d" % n, transform=axb.transAxes,
             ha="right", va="bottom", fontsize=8, color=C_NEU)
    strip(axb)
    panel_label(axb, "b")

    # ================= panel (c) =================
    # This panel used to bin systems by "recognition mode" and plot the peptide-
    # recognition ratio per bin — but the bins were terciles of that very ratio, so it
    # plotted a variable against cuts of itself. What replaces it is the quantity a
    # reader of a data descriptor actually needs: how far the three independent
    # replicas of one complex disagreed on each descriptor.
    g4 = d["global"]
    rows = [("CDR3 → peptide", "pep_recog", C_PEP),
            ("CDR1/2 → HLA", "hla_restr", C_HLA),
            ("α-chain share", "alpha_contrib", C_TCR)]
    rows = [(lab, np.asarray(g4.get(f"{k}_replica_range_values") or [], dtype=float), col)
            for lab, k, col in rows]
    rows = [r for r in rows if r[1].size]
    ypos = np.arange(len(rows))[::-1]
    rng_j = np.random.default_rng(0)
    for y, (lab, vals, colr) in zip(ypos, rows):
        bp = axc.boxplot(vals, positions=[y], vert=False, widths=0.46, whis=(5, 95),
                         showfliers=False, patch_artist=True, zorder=2)
        bp["boxes"][0].set(facecolor=colr, alpha=0.35, edgecolor=colr, linewidth=1.1)
        for part in ("whiskers", "caps"):
            for a in bp[part]:
                a.set(color=colr, linewidth=1.0)
        bp["medians"][0].set(color="#333333", linewidth=1.6)
        axc.scatter(vals, y + rng_j.normal(0, 0.075, vals.size), s=7, color=colr,
                    alpha=0.45, linewidths=0, zorder=3)
        axc.text(1.005, y, "med %.3f" % float(np.median(vals)), transform=
                 axc.get_yaxis_transform(), va="center", fontsize=7.5, color="#4d4d4d")
    # the tolerance the reproducibility flags are cut at (aggregate.REPLICA_TOL)
    tol = ((g4.get("reproducibility") or {}).get("pep_recog") or {}).get("tolerance")
    if tol:
        axc.axvline(tol, ls="--", lw=1.0, color=C_NEU, zorder=1)
        axc.text(tol, len(rows) - 0.42, " tolerance %.2f" % tol, fontsize=7.5,
                 color=C_NEU, va="top", ha="left")
    axc.set_yticks(ypos)
    axc.set_yticklabels([r[0] for r in rows], fontsize=9)
    axc.set_xlabel("Inter-replica range (3 independent runs)")
    axc.set_xlim(left=0)
    axc.set_ylim(-0.6, len(rows) - 0.4)
    axc.set_title("Replica reproducibility of the descriptors", pad=6)
    axc.grid(axis="x", alpha=0.25)
    axc.text(0.985, 0.04, "n = %d complexes" % len(rows[0][1]), transform=axc.transAxes,
             ha="right", va="bottom", fontsize=8, color=C_NEU)
    strip(axc)
    panel_label(axc, "c")

    fig.suptitle("Figure 4 — Interface recognition map", fontsize=12.5, fontweight="bold", y=0.985)
    fig.subplots_adjust(top=0.91)

    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "fig4_recognition.png")
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out

# fig5_states (interface-state counts / dominant-state occupancy / "rigid-lock <->
# dynamic continuum") was DELETED, not merely unregistered. Every panel was built on
# interface_clustering, whose cluster count comes from cutting an average-linkage
# dendrogram at a fixed absolute t=0.35 on a distance matrix that has been normalised
# by its own maximum element — so one outlier frame pair in the 201-frame subsample
# moves the answer. Measured on the 90 complexes whose three replicas were untouched
# by the 2026-08-25 frame repair, the dominant-state population had a median
# inter-replica range of 29.6 percentage points (34.6% of the library range), with
# 22/90 exceeding half the library range. The descriptor is not reproducible enough to
# publish, so the figure has no source data. Do not restore it from the stale PNG:
# delete web_data_1000/figures/fig5_states.png on BOTH origins (see the fig2_qc
# precedent -- /data is a blanket static mount, so an unregistered figure stays
# publicly fetchable until the file itself is removed).


def fig6_design(out_dir):
    import json, glob, os
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    # ---- Morandi palette ----
    C_PEP = "#B5838D"   # peptide
    C_TCR = "#6D8A96"   # tcr / cdr3
    C_HLA = "#A08E7B"   # hla
    C_STABLE = "#84A59D"
    C_FLEX = "#CB997E"
    C_NEU = "#A3A7AB"

    DATA = "/home/xmy/work/data/immunotrace/web_data_1000"

    # ---- gather hotspots across the library ----
    side_scores = {"peptide": [], "tcr": [], "hla": []}
    cat_counts = {}
    files = sorted(glob.glob(os.path.join(DATA, "*", "analysis", "analysis.json")))
    n_sys = 0
    for f in files:
        try:
            d = json.load(open(f))
        except Exception:
            continue
        hot = d.get("interface", {}).get("hotspots", [])
        if not hot:
            continue
        n_sys += 1
        for h in hot:
            s = h.get("side")
            sc = h.get("score")
            if s in side_scores and sc is not None:
                side_scores[s].append(float(sc))
            c = h.get("category")
            if c is not None:
                cat_counts[c] = cat_counts.get(c, 0) + 1

    # ---- 1ao7_run3 FEP candidates for the worked example ----
    ex_id = "1ao7_run3"
    ex = json.load(open(os.path.join(DATA, ex_id, "analysis", "analysis.json")))
    feps = ex.get("interface", {}).get("fep_candidates", [])[:6]

    # ================= FIGURE =================
    fig = plt.figure(figsize=(10.2, 7.3), dpi=200)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05],
                          width_ratios=[1.15, 1.0],
                          hspace=0.42, wspace=0.34)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    side_color = {"peptide": C_PEP, "tcr": C_TCR, "hla": C_HLA}
    side_label = {"peptide": "Peptide", "tcr": "TCR (CDR)", "hla": "HLA"}

    # ---------- Panel a: hotspot score distribution split by side (violin) ----------
    order = ["peptide", "tcr", "hla"]
    data_a = [np.array(side_scores[s]) for s in order]
    positions = np.arange(1, len(order) + 1)
    parts = ax_a.violinplot(data_a, positions=positions, showextrema=False,
                            widths=0.78)
    for pc, s in zip(parts["bodies"], order):
        pc.set_facecolor(side_color[s])
        pc.set_edgecolor("white")
        pc.set_alpha(0.85)
        pc.set_linewidth(0.8)
    # overlay median + IQR
    for i, s in enumerate(order):
        v = data_a[i]
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        x = positions[i]
        ax_a.vlines(x, q1, q3, color="#3a3a3a", lw=4, zorder=3)
        ax_a.scatter([x], [med], color="white", edgecolor="#3a3a3a",
                     s=22, zorder=4, linewidth=0.9)
    ax_a.set_xticks(positions)
    ax_a.set_xticklabels([f"{side_label[s]}\n(n={len(side_scores[s])})" for s in order],
                         fontsize=8.5)
    ax_a.set_ylabel("Interface hotspot score", fontsize=10)
    ax_a.set_title("Hotspot strength by interface side", fontsize=10, pad=6)
    ax_a.set_ylim(0.48, 0.95)
    ax_a.grid(axis="y", alpha=0.25)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)

    # ---------- Panel b: hotspot category counts across library (bar) ----------
    cat_side = {
        "peptide": C_PEP,
        "CDR3α": C_TCR, "CDR3β": C_TCR,
        "CDR1/2α": C_TCR, "CDR1/2β": C_TCR,
        "TCR-FRα": C_TCR, "TCR-FRβ": C_TCR,
        "HLA-groove": C_HLA, "HLA": C_HLA,
        "β2m": C_NEU,
    }
    items = sorted(cat_counts.items(), key=lambda kv: kv[1], reverse=True)
    labels_b = [k for k, _ in items]
    vals_b = [v for _, v in items]
    colors_b = [cat_side.get(k, C_NEU) for k in labels_b]
    ypos = np.arange(len(labels_b))[::-1]
    ax_b.barh(ypos, vals_b, color=colors_b, edgecolor="white", linewidth=0.6, height=0.74)
    ax_b.set_yticks(ypos)
    ax_b.set_yticklabels(labels_b, fontsize=8.2)
    for y, v in zip(ypos, vals_b):
        ax_b.text(v + max(vals_b) * 0.015, y, str(v), va="center",
                  fontsize=7.6, color="#3a3a3a")
    ax_b.set_xlabel("Hotspot residues across library", fontsize=9.5)
    ax_b.set_title("Hotspot category counts", fontsize=10, pad=6)
    ax_b.set_xlim(0, max(vals_b) * 1.14)
    ax_b.grid(axis="x", alpha=0.25)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    legend_b = [
        Patch(facecolor=C_PEP, label="Peptide"),
        Patch(facecolor=C_TCR, label="TCR"),
        Patch(facecolor=C_HLA, label="HLA"),
        Patch(facecolor=C_NEU, label="β2m"),
    ]
    ax_b.legend(handles=legend_b, fontsize=7.6, frameon=False,
                loc="lower right", handlelength=1.0, labelspacing=0.3)

    # ---------- Panel c: worked example FEP candidate table (1ao7_run3) ----------
    ax_c.axis("off")
    col_titles = ["Rank", "Residue", "Candidate type", "Score",
                  "Suggested mutation", "Validation"]
    col_x = [0.02, 0.06, 0.245, 0.45, 0.525, 0.86]
    col_w_align = ["center", "left", "left", "center", "left", "left"]

    n_rows = len(feps)
    y_top = 0.86
    y_bot = 0.06
    row_h = (y_top - y_bot) / (n_rows + 1)

    ax_c.add_patch(plt.Rectangle((0.0, y_top - row_h * 0.5), 1.0, row_h,
                                 transform=ax_c.transAxes, facecolor="#4a5560",
                                 edgecolor="none", zorder=1))
    for x, t, al in zip(col_x, col_titles, col_w_align):
        ax_c.text(x, y_top, t, transform=ax_c.transAxes, ha=al, va="center",
                  fontsize=8.6, fontweight="bold", color="white", zorder=2)

    type_short = {
        "peptide TCR-facing (escape candidate)": "peptide escape",
        "CDR3 affinity / specificity": "CDR3 affinity",
        "HLA contact / restriction": "HLA restriction",
    }
    side_of_type = {
        "peptide escape": C_PEP,
        "CDR3 affinity": C_TCR,
        "HLA restriction": C_HLA,
    }

    for i, fp in enumerate(feps):
        y = y_top - row_h * (i + 1)
        if i % 2 == 0:
            ax_c.add_patch(plt.Rectangle((0.0, y - row_h * 0.5), 1.0, row_h,
                                         transform=ax_c.transAxes,
                                         facecolor="#f2efec", edgecolor="none", zorder=0))
        ct = type_short.get(fp.get("candidate_type", ""), fp.get("candidate_type", ""))
        dot = side_of_type.get(ct, C_NEU)
        res = fp.get("residue", "")
        sc = fp.get("score", "")
        mut = fp.get("suggested_mutation", "")
        mut_short = mut.split(";")[0].strip() if mut else ""
        if len(mut_short) > 48:                 # keep within the column, never overrun Validation
            mut_short = mut_short[:47].rstrip() + "…"
        val = fp.get("recommended_validation", "")
        val_short = val.split(";")[0].strip() if val else ""

        ax_c.text(col_x[0], y, str(i + 1), transform=ax_c.transAxes, ha="center",
                  va="center", fontsize=8.0, color="#3a3a3a")
        ax_c.text(col_x[1], y, res, transform=ax_c.transAxes, ha="left",
                  va="center", fontsize=8.0, color="#222", fontweight="bold")
        ax_c.scatter([col_x[2] - 0.018], [y], transform=ax_c.transAxes,
                     s=42, color=dot, edgecolor="white", linewidth=0.6, zorder=3,
                     clip_on=False)
        ax_c.text(col_x[2], y, ct, transform=ax_c.transAxes, ha="left",
                  va="center", fontsize=7.8, color="#3a3a3a")
        ax_c.text(col_x[3], y, f"{sc:.3f}", transform=ax_c.transAxes, ha="center",
                  va="center", fontsize=8.0, color="#3a3a3a")
        ax_c.text(col_x[4], y, mut_short, transform=ax_c.transAxes, ha="left",
                  va="center", fontsize=7.5, color="#3a3a3a")
        ax_c.text(col_x[5], y, val_short, transform=ax_c.transAxes, ha="left",
                  va="center", fontsize=7.6, color="#3a3a3a")

    ax_c.set_title(
        f"Worked example: top FEP / mutation candidates for {ex_id}  "
        f"(HLA-A*02:01 / NLVPMVATV-like, well-characterised complex)",
        fontsize=9.6, pad=8, loc="left")
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)

    # ---- panel letters ----
    for ax, lab in [(ax_a, "a"), (ax_b, "b"), (ax_c, "c")]:
        ax.text(-0.02 if ax is not ax_c else -0.01, 1.06, lab,
                transform=ax.transAxes, fontsize=14, fontweight="bold",
                va="bottom", ha="right")

    fig.suptitle("Figure 6 — Interface hotspots & mutation-design utility",
                 fontsize=12.5, y=0.995, fontweight="bold")

    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "fig6_design.png")
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


# fig2_qc (equilibration/RMSD-QC panel) is intentionally NOT generated or served —
# equilibration diagnostics are kept internal-only.
ALL_FIGURES = [fig1_overview, fig3_docking, fig4_recognition, fig6_design]


def main(argv=None) -> int:
    out_dir = str(config.WEB_DATA / "figures")
    os.makedirs(out_dir, exist_ok=True)
    for fn in ALL_FIGURES:
        path = fn(out_dir)
        print(f"  wrote {path}")
    print(f"Figures: {len(ALL_FIGURES)} written to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
