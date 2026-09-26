"""Concerted-motion / allosteric-coupling analysis (per-complex, charlie section 1).

Whole-complex Cα dynamic cross-correlation (DCCM), a persistent-contact network
weighted by |corr|, Louvain communities + betweenness hubs, collective-mode
descriptors, and a peptide<->CDR3 cross-interface coupling scalar. Everything is
computed on the 200-frame web trajectory; the module is purely trajectory-derived
(topology.pdb + traj.xtc) so it globs WEB_DATA and takes no --source.

It writes a sidecar analysis/concerted_motion.json (re-embedded by
extract_analysis) AND patches analysis.json["concerted_motion"], and emits two
PNGs into web_data/<id>/analysis/: dccm_heatmap.png, concerted_network.png.

The JSON also carries ``communities``: per-module residue blocks (chain + resSeq
ranges, size-ranked) so the detail page can color the live 3D structure by dynamic
module (MolViewSpec) and show where each concerted module sits — see
frontend/js/detail.js (colorByModules / buildModuleMVS).

Run (immuno-web has mdtraj + networkx 3.6.1 + matplotlib):

    IMMUNO_WEB_DATA=/home/xmy/work/data/immunotrace/web_data_1000 \\
        python -m pipeline.concerted_motion [--ids a,b]
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import config

SIDECAR = "concerted_motion.json"
ANALYSIS_KEY = "concerted_motion"

CONTACT_CUTOFF_NM = 1.0      # 10 Angstrom on Cα (coarse contact proxy)
CONTACT_OCCUPANCY = 0.75     # persistent if the pair is within cutoff in >=75% of frames
MIN_CORR = 1e-6              # drop near-zero-correlation edges before -log
LOUVAIN_SEED = 42
LAYOUT_SEED = 7
HEATMAP_PNG = "dccm_heatmap.png"
NETWORK_PNG = "concerted_network.png"

ROLE_ORDER = ("MHC", "PEP", "TCRa", "TCRb")   # node-block order in the DCCM

# Morandi house palette (figures.py / detail.js C object). TCRβ gets a lighter
# teal than TCRα so the two TCR variable domains read distinctly in the network.
COL = {"PEP": "#B5838D", "TCRa": "#6D8A96", "TCRb": "#9DB4BD", "MHC": "#A08E7B"}

_AA3TO1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
           "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
           "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
           "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def _chain_ca(top, chid):
    """Ordered (ca_atom_index, resSeq, resName) for a chain, by chain id."""
    out = []
    for ch in top.chains:
        if ch.chain_id != chid:
            continue
        for r in ch.residues:
            ca = next((a.index for a in r.atoms if a.name == "CA"), None)
            if ca is not None:
                out.append((ca, int(r.resSeq), r.name))
    return out


def _mhc_chain(meta, pep_ch, a_ch, b_ch):
    """MHC heavy chain = the largest chain that is not peptide / TCRα / TCRβ."""
    cand = [c for c in (meta.get("chains") or []) if c["id"] not in (pep_ch, a_ch, b_ch)]
    return max(cand, key=lambda c: c["n_residues"])["id"] if cand else None


def _resolve_cdr3_nodes(node_chain, node_resname, a_ch, b_ch, cdr3a, cdr3b):
    """Node indices whose 1-letter residues form the CDR3 subsequence in each TCR chain.

    Exact subsequence match of the meta CDR3 string within the chain's ordered
    residue list. Returns the union over the α and β CDR3 loops (may be empty)."""
    out = []
    for chid, cdr3 in ((a_ch, cdr3a), (b_ch, cdr3b)):
        if not cdr3:
            continue
        idx = [i for i, c in enumerate(node_chain) if c == chid]
        seq = "".join(_AA3TO1.get(node_resname[i], "X") for i in idx)
        pos = seq.find(cdr3)
        if pos >= 0:
            out.extend(idx[pos:pos + len(cdr3)])
    return out


def _ranges(seqs):
    """Compress a residue-number list into contiguous [beg, end] runs."""
    out = []
    for s in sorted(set(int(x) for x in seqs)):
        if out and s == out[-1][1] + 1:
            out[-1][1] = s
        else:
            out.append([s, s])
    return out


def _community_blocks(comm_label, node_chain, node_resseq, node_role):
    """Per-Louvain-community residue blocks for 3D module coloring.

    Returns a list (sorted by community id, which is already size-rank: 0 = largest)
    of ``{id, n, roles: {role: count}, ranges: {chain: [[beg_resSeq, end_resSeq], …]}}``.
    The per-chain resSeq ranges map directly onto MolViewSpec component selectors.
    """
    from collections import defaultdict
    by_comm = defaultdict(lambda: {"residues": defaultdict(list),
                                   "roles": defaultdict(int), "n": 0})
    for i, k in enumerate(comm_label):
        k = int(k)
        if k < 0:
            continue
        c = by_comm[k]
        c["residues"][node_chain[i]].append(int(node_resseq[i]))
        c["roles"][node_role[i]] += 1
        c["n"] += 1
    out = []
    for k in sorted(by_comm):
        c = by_comm[k]
        out.append({
            "id": k,
            "n": c["n"],
            "roles": dict(c["roles"]),
            "ranges": {ch: _ranges(seqs) for ch, seqs in sorted(c["residues"].items())},
        })
    return out


def compute(traj_id: str, meta: dict) -> dict | None:
    """Concerted-motion summary for one trajectory, or None on missing inputs.

    Pure except for writing the two PNGs into the trajectory's analysis dir (the
    figure data is in hand here, so we render rather than recompute downstream)."""
    import mdtraj as md

    pep_ch = meta.get("peptide_chain")
    tc = (meta.get("tcr") or {}).get("chains") or {}
    a_ch = (tc.get("alpha") or {}).get("structural_chain")
    b_ch = (tc.get("beta") or {}).get("structural_chain")
    if not (pep_ch and a_ch and b_ch):
        return None
    mhc_ch = _mhc_chain(meta, pep_ch, a_ch, b_ch)
    if not mhc_ch:
        return None

    wd = config.WEB_DATA / traj_id
    top_path, traj_path = wd / config.OUT_TOPOLOGY, wd / config.OUT_TRAJ
    if not top_path.exists() or not traj_path.exists():
        return None
    t = md.load(str(traj_path), top=str(top_path))
    if t.n_frames < 2:
        return None
    top = t.topology

    # --- ordered whole-complex Cα node set: MHC, PEP, TCRa, TCRb (b2m excluded) ---
    role_chain = {"MHC": mhc_ch, "PEP": pep_ch, "TCRa": a_ch, "TCRb": b_ch}
    ca_idx, node_role, node_chain, node_resseq, node_resname = [], [], [], [], []
    block_sizes = {}
    for role in ROLE_ORDER:
        n = 0
        for (ca, resseq, resname) in _chain_ca(top, role_chain[role]):
            ca_idx.append(ca); node_role.append(role); node_chain.append(role_chain[role])
            node_resseq.append(resseq); node_resname.append(resname); n += 1
        block_sizes[role] = n
    N = len(ca_idx)
    if N < 10 or block_sizes["PEP"] < 2:
        return None

    # --- DCCM after Cα superposition (exactly symmetric, in [-1, 1]) ---
    sub = t.atom_slice(ca_idx)
    sub.superpose(sub, frame=0)
    xyz = sub.xyz.astype(np.float64)                  # (F, N, 3) nm
    disp = xyz - xyz.mean(axis=0, keepdims=True)
    F = disp.shape[0]
    cov = np.einsum("fix,fjx->ij", disp, disp) / F
    msf = np.diag(cov).copy()
    denom = np.sqrt(np.outer(msf, msf)); denom[denom == 0] = 1.0
    dccm = np.clip(cov / denom, -1.0, 1.0)

    # --- memory-bounded persistent-contact occupancy (NxN accumulator) ---
    occ = np.zeros((N, N), dtype=np.float32)
    cut2 = CONTACT_CUTOFF_NM ** 2
    for f in range(F):
        d = xyz[f][:, None, :] - xyz[f][None, :, :]
        occ += (np.einsum("ijx,ijx->ij", d, d) <= cut2)
    occ /= F
    contact = occ >= CONTACT_OCCUPANCY
    np.fill_diagonal(contact, False)

    # --- network: persistent contacts weighted by |corr|, length = -log|corr| ---
    G = nx.Graph()
    G.add_nodes_from(range(N))
    iu, ju = np.where(np.triu(contact, k=1))
    for i, j in zip(iu.tolist(), ju.tolist()):
        c = abs(float(dccm[i, j]))
        if c < MIN_CORR:
            continue
        G.add_edge(i, j, weight=c, length=-float(np.log(c)))
    n_edges = G.number_of_edges()

    comms = nx.community.louvain_communities(G, weight="weight", seed=LOUVAIN_SEED)
    modularity = float(nx.community.modularity(G, comms, weight="weight"))
    btw = nx.betweenness_centrality(G, weight="length", normalized=True)
    comm_label = np.full(N, -1, dtype=int)
    for k, cset in enumerate(sorted(comms, key=len, reverse=True)):
        for i in cset:
            comm_label[i] = k

    # --- collective-mode descriptors from the Cα covariance spectrum ---
    evals = np.linalg.eigvalsh(cov)
    evals = np.clip(evals[::-1], 0, None)
    tot = float(evals.sum()) or 1.0
    frac = evals / tot
    n_modes_90 = int(np.searchsorted(np.cumsum(frac), 0.90) + 1)
    p = frac[frac > 0]
    collectivity = float(np.exp(-np.sum(p * np.log(p))) / N)   # normalized participation ratio

    # --- cross-interface coupling: peptide <-> CDR3 mean |corr| ---
    pep_nodes = [i for i in range(N) if node_role[i] == "PEP"]
    cdr3_nodes = _resolve_cdr3_nodes(
        node_chain, node_resname, a_ch, b_ch,
        (tc.get("alpha") or {}).get("cdr3"), (tc.get("beta") or {}).get("cdr3"))
    pep_tcr_coupling = None
    if pep_nodes and cdr3_nodes:
        pep_tcr_coupling = float(np.abs(dccm[np.ix_(pep_nodes, cdr3_nodes)]).mean())

    # whole-complex off-diagonal mean |corr|
    iu2, ju2 = np.triu_indices(N, k=1)
    mean_correlation = float(np.abs(dccm)[iu2, ju2].mean())

    # --- top betweenness hubs (label = ROLE:RESNAME+RESSEQ) ---
    def label(i):
        return f"{node_role[i]}:{node_resname[i]}{node_resseq[i]}"
    hubs = sorted(range(N), key=lambda i: btw[i], reverse=True)[:10]
    hub_list = [{"residue": label(i), "betweenness": round(float(btw[i]), 4)} for i in hubs]

    # --- render the two PNGs into the analysis dir ---
    an_dir = wd / "analysis"
    an_dir.mkdir(parents=True, exist_ok=True)
    _plot_dccm(dccm, block_sizes, an_dir / HEATMAP_PNG)
    _plot_network(G, node_role, comm_label, set(pep_nodes), an_dir / NETWORK_PNG)

    return {
        "n_residues": int(N),
        "block_sizes": block_sizes,
        "mean_correlation": round(mean_correlation, 4),
        "peptide_tcr_coupling": (round(pep_tcr_coupling, 4) if pep_tcr_coupling is not None else None),
        "collectivity": round(collectivity, 4),
        "n_modes_90pct": n_modes_90,
        "n_communities": len(comms),
        "modularity": round(modularity, 4),
        "n_contact_edges": int(n_edges),
        "top_hubs": hub_list,
        # per-residue community membership (chain + resSeq ranges) for 3D module
        # coloring on the detail page; size-ranked (module 0 = largest).
        "communities": _community_blocks(comm_label, node_chain, node_resseq, node_role),
        "figures": [HEATMAP_PNG, NETWORK_PNG],
        "method": "ca_dccm_persistent_contact_louvain",
    }


def _style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _plot_dccm(dccm, block_sizes, out_path):
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(dccm, cmap="RdBu_r", vmin=-1.0, vmax=1.0,
                   interpolation="nearest", origin="upper")
    bounds, acc = [], 0
    for role in ROLE_ORDER:
        n = block_sizes.get(role, 0)
        if n == 0:
            continue
        bounds.append((role, acc, acc + n)); acc += n
    for _, s, _e in bounds[1:]:
        ax.axhline(s - 0.5, color="white", lw=0.8)
        ax.axvline(s - 0.5, color="white", lw=0.8)
    ticks = [(s + e - 1) / 2 for _, s, e in bounds]
    labels = [r for r, _, _ in bounds]
    ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=9)
    ax.set_yticks(ticks); ax.set_yticklabels(labels, fontsize=9)
    ax.set_title("Cα dynamic cross-correlation (DCCM)", fontsize=10.5, pad=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("correlation", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _plot_network(G, node_role, comm_label, pep_set, out_path):
    from matplotlib.lines import Line2D
    pos = nx.spring_layout(G, weight="weight", seed=LAYOUT_SEED)
    fig, ax = plt.subplots(figsize=(6.0, 5.4))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#D9D6D1", width=0.4, alpha=0.6)
    colors = [COL[node_role[i]] for i in G.nodes()]
    sizes = [40 if i in pep_set else 14 for i in G.nodes()]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=colors, node_size=sizes,
                           linewidths=0.4, edgecolors="white")
    ax.set_axis_off()
    n_comm = len(set(int(c) for c in comm_label if c >= 0))
    ax.set_title(f"Persistent-contact network · {n_comm} dynamic communities",
                 fontsize=10.5, pad=8)
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=7,
                      markerfacecolor=COL[r], markeredgecolor="white", label=r)
               for r in ("MHC", "PEP", "TCRa", "TCRb")]
    ax.legend(handles=handles, loc="lower right", fontsize=8, frameon=False)
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", default=None, help="Comma-separated traj ids (default: all).")
    args = ap.parse_args(argv)
    only = set(args.ids.split(",")) if args.ids else None

    web = config.WEB_DATA
    n_ok = n_skip = 0
    for meta_path in sorted(web.glob("*/" + config.OUT_META)):
        tid = meta_path.parent.name
        if only and tid not in only:
            continue
        an_dir = meta_path.parent / "analysis"
        an_path = an_dir / "analysis.json"
        if not an_path.exists():
            n_skip += 1
            continue
        meta = json.loads(meta_path.read_text())
        try:
            result = compute(tid, meta)
        except Exception as exc:           # never let one bad trajectory abort the batch
            print(f"[warn] {tid}: {exc}")
            result = None
        if result is None:
            n_skip += 1
            continue
        (an_dir / SIDECAR).write_text(json.dumps(result))         # sidecar (re-embedded later)
        analysis = json.loads(an_path.read_text())                # idempotent patch
        analysis[ANALYSIS_KEY] = result
        for png in (HEATMAP_PNG, NETWORK_PNG):
            if (an_dir / png).exists() and png not in analysis.get("figures", []):
                analysis.setdefault("figures", []).append(png)
        an_path.write_text(json.dumps(analysis, indent=2))
        n_ok += 1

    print(f"Concerted motion computed for {n_ok} trajectories, {n_skip} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
