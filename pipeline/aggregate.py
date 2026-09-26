"""Cross-system aggregation — the database's comparative layer (bravo §5, figs 1 & 3).

Per-complex analysis answers "what does this complex do"; this answers "how do
complexes compare". It joins every analysis.json with the manifest and writes a
single ``aggregate.json`` (served at /api/aggregate, consumed by the Explore
page and the publication-figure module):

  * global distributions   — crossing/incident angle, recognition ratios, bulge,
    contact entropy, peptide length, recognition modes
  * grouped comparisons    — crossing angle by HLA group and by peptide length,
    recognition ratios by recognition mode and by TRBV family
  * conserved peptide map   — mean TCR / HLA contact by peptide position

All crossing angles are kept (no filtering). The crossing angle is the DIRECTED angle
in 0-180 deg (``docking_angle``): forward dockings keep their acute 0-90 value,
reverse-polarity dockings are reflected past 90 deg (180 - acute) so they land in the
upper half toward 180. Docking POLARITY is reported explicitly (``reversed_polarity`` +
split ``crossing_values_forward``/``_reversed`` + ``reversed_polarity_n``).

    IMMUNO_WEB_DATA=.../immuno-dyn python -m pipeline.aggregate
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from . import config
from . import cdr_contacts
from . import consensus

MIN_GROUP = 5          # keep groups with at least this many systems
TOP_GROUPS = 8         # cap displayed HLA / TRBV groups
DOMINANT_LENGTH = 9    # peptide length for the conserved-position map

# When a complex has several replicas, numeric metrics are averaged across them and
# replica-invariant categoricals are taken from the primary replica. Docking polarity
# VARIES across replicas, so it is NOT taken from the primary — it is voted across
# them (see the collapse loop). The three quantitative interface descriptors keep
# their measured per-replica spread alongside the mean (`replica_stats`), which is
# what the withdrawn recognition-mode label axes used to hide behind a boolean flag.
_NUMERIC = ("crossing", "incident", "pep_recog", "hla_restr", "alpha_contrib", "bulge", "entropy",
            "fnat", "cosine", "rmsip", "bsa_total")
_CATEG = ("pdb_id", "hla_group", "peptide_length", "trbv", "host")

# Extra raw per-complex features the layered "map of the database" needs beyond
# _NUMERIC. Collapsed across replicas by mean, exactly like _NUMERIC. All ENCODING
# (sin/cos of angle, logit of ratios, per-residue counts, RMSF ratios) happens later
# in _database_map so replica collapse stays a simple mean of raw quantities.
_MAP_EXTRA = ("rg_mean_nm", "total_sasa_mean_nm2",
              "cm_mean_corr", "cm_pep_tcr_coupling", "cm_modularity", "cm_collectivity",
              "cm_n_communities", "cm_n_contact_edges", "cm_n_residues",
              "rmsf_iface", "rmsf_cdr3", "rmsf_fr", "rmsf_peptide", "rmsf_groove",
              "ed_pc1_var_frac")

# Interface-side RMSF regions (contact-bearing); framework / groove used for ratios.
_RMSF_IFACE = ("peptide", "α1 helix", "α2 helix",
               "CDR1α", "CDR2α", "CDR3α", "CDR1β", "CDR2β", "CDR3β")


def _rmsf_region_mean(regions: dict, names) -> float | None:
    """Mean of region-mean Cα RMSF over `names` (None if none are present)."""
    vals = [(regions.get(n) or {}).get("mean") for n in names]
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else None

# Inter-replica reproducibility tolerances: a metric is "reproducible" for a complex
# when the spread (max-min) across its replicas is <= the tolerance. Defaults are ~p90
# of the observed inter-replica |Δ| over the two-replica complexes (see
# docs/research/2026-06-16_replica_aggregation.md). The dominant-state population is
# handled separately (tercile agreement), as it is the least reproducible metric.
REPLICA_TOL = {"crossing": 10.0, "incident": 8.0,        # degrees
               "pep_recog": 0.10, "hla_restr": 0.10, "alpha_contrib": 0.10,  # ratios
               "bulge": 1.5,                              # Angstrom
               "fnat": 0.15}                              # fraction of native contacts
_REPLICA_METRICS = tuple(REPLICA_TOL)


def _hla_group(allele: str | None) -> str | None:
    """Coarsen to the HLA locus ('HLA-A*02:01' -> 'HLA-A').

    Most manifest alleles are only annotated to the locus, so grouping at the
    allele-group level ('HLA-A*02') would split the same locus into a large
    bare 'HLA-A' bucket plus a few finely-annotated ones — inconsistent and
    confusing. Locus-level keeps the groups consistent and comparable.
    """
    if not allele:
        return None
    return allele.split("*")[0].strip() or None


def _stat(values: list[float]) -> dict:
    """Summary stats with proper linearly-interpolated quantiles (numpy), not the
    biased nearest-rank s[int(p*n)] this used to use (which mis-placed medians and
    the IQR error bars on the grouped charts, worst for small groups)."""
    s = [v for v in values if v is not None]
    if not s:
        return {"n": 0}
    a = np.asarray(s, dtype=float)
    p25, med, p75 = (float(q) for q in np.percentile(a, [25, 50, 75], method="linear"))
    return {"n": len(s), "mean": round(float(a.mean()), 3), "median": round(med, 3),
            "p25": round(p25, 3), "p75": round(p75, 3),
            "min": round(float(a.min()), 3), "max": round(float(a.max()), 3)}


def _replica_stat(values: list, tol: float | None) -> dict | None:
    """n-agnostic per-metric replica summary (works for 1, 2, 3+ replicas).

    spread = half-range (max-min)/2 over the replicas; reproducible = whether the
    full range is within `tol` (None for a single replica — reproducibility is
    undefined without a second run). `parallel1`/`run1` will simply add a value."""
    vals = [float(v) for v in values if v is not None]
    n = len(vals)
    if n == 0:
        return None
    rng = (max(vals) - min(vals)) if n >= 2 else 0.0
    reproducible = None if (n < 2 or tol is None) else bool(rng <= tol)
    return {"mean": round(sum(vals) / n, 4), "values": [round(v, 4) for v in vals],
            "n": n, "spread": round(rng / 2.0, 4), "reproducible": reproducible}


# ---------------------------------------------------------------------------
# Layered "map of the database" — 2-block Multiple Factor Analysis.
#
# The 11-static-feature z-scored PCA it replaces had several rigour problems:
# directed angles z-scored as if linear, extensive size features leaking onto PC1,
# three collinear interface-engagement columns triple-weighting one concept, and
# NaN -> column-mean imputation. It also captured only static geometry.
#
# The rebuild (see docs/research/2026-07-06_layered_embedding.py for the empirical
# study) organises features into TWO blocks and normalises each so neither can
# dominate a served axis:
#   Block G — geometry / quasi-static interface (L1): docking angles (sin/cos of the
#     directed crossing + incident tilt), bulge, Rg, buried-surface FRACTION, contact
#     entropy, and the three recognition ratios (logit-encoded).
#   Block D — dynamics: L2 contact-network (DCCM mean corr, peptide-TCR coupling,
#     modularity, communities/edges per residue, collectivity) + L3 flexibility
#     (interface RMSF, CDR3-vs-framework & peptide-vs-groove mobility ratios, fnat,
#     PC1 variance fraction).
# Empirically the two blocks are more separable than the raw three layers (block
# canonical corr ~0.77 vs ~0.9). Per-block classic-MFA normalisation (divide by the
# block's first singular value) then a joint PCA gives PC1 balanced ~55/45 and keeps
# no near-degenerate block from taking an axis. Angles are encoded on the circle,
# every extensive count is intensified (per-residue / fraction), and a complex whose
# whole contact-network block is absent (7byd) is FLAGGED (network_missing), never
# imputed from geometric neighbours.
# ---------------------------------------------------------------------------
_MAP_SEED = 42


def _map_feature_matrix(recs: list[dict]):
    """Build the encoded 2-block feature matrix. Returns (names, blocks, layers, M)."""
    def num(r, k):
        v = r.get(k)
        return float(v) if v is not None else np.nan

    def logit(v, lo=0.02, hi=0.98):
        if v is None:
            return np.nan
        v = min(max(float(v), lo), hi)
        return float(np.log(v / (1.0 - v)))

    def logratio(a, b):
        if a is None or b is None or not b:
            return np.nan
        try:
            return float(np.log(float(a) / float(b)))
        except (ValueError, ZeroDivisionError):
            return np.nan

    def per_res(a, n):
        return float(a) / float(n) if (a is not None and n) else np.nan

    def circ(r, fn):
        c = r.get("crossing")
        return float(fn(np.radians(c))) if c is not None else np.nan

    def bsa_frac(r):
        b, s = r.get("bsa_total"), r.get("total_sasa_mean_nm2")
        return float(b) / float(s) if (b and s) else np.nan

    # (name, block, layer, encoder)
    spec = [
        ("crossing_sin", "G", "L1", lambda r: circ(r, np.sin)),
        ("crossing_cos", "G", "L1", lambda r: circ(r, np.cos)),
        ("incident", "G", "L1", lambda r: num(r, "incident")),
        ("bulge", "G", "L1", lambda r: num(r, "bulge")),
        ("rg", "G", "L1", lambda r: num(r, "rg_mean_nm")),
        ("bsa_frac", "G", "L1", bsa_frac),
        ("entropy", "G", "L1", lambda r: num(r, "entropy")),
        ("pep_recog", "G", "L1", lambda r: logit(r.get("pep_recog"))),
        ("hla_restr", "G", "L1", lambda r: logit(r.get("hla_restr"))),
        ("alpha_contrib", "G", "L1", lambda r: logit(r.get("alpha_contrib"))),
        ("dccm_corr", "D", "L2", lambda r: num(r, "cm_mean_corr")),
        ("pep_tcr_coupling", "D", "L2", lambda r: num(r, "cm_pep_tcr_coupling")),
        ("modularity", "D", "L2", lambda r: num(r, "cm_modularity")),
        ("communities_per_res", "D", "L2", lambda r: per_res(r.get("cm_n_communities"), r.get("cm_n_residues"))),
        ("edges_per_res", "D", "L2", lambda r: per_res(r.get("cm_n_contact_edges"), r.get("cm_n_residues"))),
        ("collectivity", "D", "L2", lambda r: num(r, "cm_collectivity")),
        ("rmsf_iface", "D", "L3", lambda r: num(r, "rmsf_iface")),
        ("cdr3_vs_fr", "D", "L3", lambda r: logratio(r.get("rmsf_cdr3"), r.get("rmsf_fr"))),
        ("pep_vs_groove", "D", "L3", lambda r: logratio(r.get("rmsf_peptide"), r.get("rmsf_groove"))),
        ("fnat", "D", "L3", lambda r: num(r, "fnat")),
        ("pc1_var_frac", "D", "L3", lambda r: num(r, "ed_pc1_var_frac")),
    ]
    names = [s[0] for s in spec]
    blocks = [s[1] for s in spec]
    layers = [s[2] for s in spec]
    M = np.array([[s[3](r) for s in spec] for r in recs], dtype=float)
    return names, blocks, layers, M


def _robust_scale(X: np.ndarray) -> np.ndarray:
    """Median / IQR scaling, then median-impute residual NaN (post-scale median ~ 0)."""
    med = np.nanmedian(X, axis=0)
    iqr = np.nanpercentile(X, 75, axis=0) - np.nanpercentile(X, 25, axis=0)
    iqr[iqr < 1e-9] = 1.0
    Z = (X - med) / iqr
    Z[np.isnan(Z)] = 0.0
    return Z


def _horn_k(Z: np.ndarray, n_iter: int = 200, pct: int = 95) -> int:
    """Horn parallel analysis: # real eigenvalues exceeding the column-permuted null
    (deterministic — fixed seed). Data-driven per-block component count."""
    n, p = Z.shape
    Zc = Z - Z.mean(0)
    real = np.linalg.svd(Zc, compute_uv=False) ** 2
    rng = np.random.default_rng(_MAP_SEED)
    perm = np.zeros((n_iter, min(n, p)))
    for t in range(n_iter):
        Zs = np.column_stack([rng.permutation(Zc[:, j]) for j in range(p)])
        perm[t] = (np.linalg.svd(Zs, compute_uv=False) ** 2)[:min(n, p)]
    thr = np.percentile(perm, pct, axis=0)
    return max(1, min(int(np.sum(real[:len(thr)] > thr)), p))


def _block_pca(X: np.ndarray):
    """RobustScale -> PCA. Returns (scores n×p, loadings p×p over scaled feats, eig, k)."""
    Z = _robust_scale(X)
    Zc = Z - Z.mean(0)
    U, S, Vt = np.linalg.svd(Zc, full_matrices=False)
    return U * S, Vt, S ** 2, _horn_k(Z)


def _database_map(recs: list[dict]) -> dict | None:
    """2D layered 2-block MFA embedding — the navigable "map of the database".
    Block G (geometry / L1) and Block D (dynamics / L2+L3) are each PCA-reduced to a
    Horn-selected rank, classic-MFA normalised (÷ first singular value), concatenated
    and jointly PCA'd. Returns per-complex points, raw-feature PC loadings (layer-
    tagged), and the block-contribution table."""
    if len(recs) < 4:
        return None
    names, blocks, layers, M = _map_feature_matrix(recs)
    gcol = [i for i, b in enumerate(blocks) if b == "G"]
    dcol = [i for i, b in enumerate(blocks) if b == "D"]
    net_cols = [i for i, ly in enumerate(layers) if ly == "L2"]
    # complexes whose whole contact-network block is absent (e.g. 7byd) — flag, never
    # fabricate a dynamics position from geometric neighbours.
    net_missing = {i for i in range(len(recs)) if all(np.isnan(M[i, j]) for j in net_cols)}

    sG, vtG, eigG, kG = _block_pca(M[:, gcol])
    sD, vtD, eigD, kD = _block_pca(M[:, dcol])
    Gn = sG[:, :kG] / np.sqrt(eigG[0])          # classic MFA: block top-axis inertia = 1
    Dn = sD[:, :kD] / np.sqrt(eigD[0])
    block_of_col = (["G"] * kG) + (["D"] * kD)
    Jc = np.hstack([Gn, Dn])
    Jc = Jc - Jc.mean(0)
    Uj, Sj, Vtj = np.linalg.svd(Jc, full_matrices=False)
    scores = Uj[:, :2] * Sj[:2]
    vf = Sj ** 2 / np.sum(Sj ** 2)

    # Interpretable loadings: correlation of each (robust-scaled) raw feature with the
    # joint PC score — a standard biplot loading in [-1, 1], readable in terms of the
    # real features rather than opaque block components.
    Zall = _robust_scale(M)                     # per-column scaled, order == names

    def corr_loadings(pc_scores):
        s = pc_scores - pc_scores.mean()
        out = np.zeros(len(names))
        for j in range(len(names)):
            f = Zall[:, j] - Zall[:, j].mean()
            denom = f.std() * s.std()
            out[j] = float(np.dot(f, s) / (len(s) * denom)) if denom > 1e-12 else 0.0
        return out

    rl = [corr_loadings(scores[:, pc]) for pc in range(2)]
    # deterministic axis orientation: pin the largest-|loading| feature positive.
    for pc in range(2):
        if rl[pc][int(np.argmax(np.abs(rl[pc])))] < 0:
            scores[:, pc] *= -1
            rl[pc] = -rl[pc]
    contrib = {}
    for pc in range(2):
        w = Vtj[pc] ** 2
        contrib[f"pc{pc+1}"] = {b: round(float(np.sum([w[i] for i, bb in enumerate(block_of_col)
                                                       if bb == b])), 3) for b in ("G", "D")}

    # t-SNE on the joint block-PC space — the "spread" view (axes non-quantitative).
    tsne = None
    perp = max(5, min(30, (len(recs) - 1) // 3))
    try:
        from sklearn.manifold import TSNE
        tsne = TSNE(n_components=2, perplexity=perp, init="pca",
                    learning_rate="auto", random_state=_MAP_SEED).fit_transform(Jc)
    except Exception:
        tsne = None

    pts = []
    for i, r in enumerate(recs):
        pt = {
            "id": r.get("id"), "pdb_id": r.get("pdb_id"),
            "x": round(float(scores[i, 0]), 3), "y": round(float(scores[i, 1]), 3),
            "hla_group": r.get("hla_group"), "peptide_length": r.get("peptide_length"),
            "pep_recog": r.get("pep_recog"), "hla_restr": r.get("hla_restr"),
            "alpha_contrib": r.get("alpha_contrib"),
        }
        if i in net_missing:
            pt["network_missing"] = True
        if tsne is not None:
            pt["tx"] = round(float(tsne[i, 0]), 2)
            pt["ty"] = round(float(tsne[i, 1]), 2)
        pts.append(pt)

    def loadings(pc):
        return [[f"{layers[j]}·{names[j]}", round(float(rl[pc][j]), 2)]
                for j in np.argsort(-np.abs(rl[pc]))[:5]]

    return {"points": pts, "n_features": len(names), "features": names,
            "method": "layered-2block-mfa", "block_kG": int(kG), "block_kD": int(kD),
            "block_contribution": contrib, "n_network_missing": len(net_missing),
            "pc1_var": round(float(vf[0]), 3), "pc2_var": round(float(vf[1]), 3),
            "pc1_loadings": loadings(0), "pc2_loadings": loadings(1),
            "has_tsne": tsne is not None, "tsne_perplexity": perp if tsne is not None else None}


def build() -> dict:
    web = config.WEB_DATA
    manifest = json.loads((web / "manifest.json").read_text())
    rows = manifest if isinstance(manifest, list) else manifest.get("trajectories", [])
    by_id = {r["traj_id"]: r for r in rows}

    # per-trajectory scalar record (one per analysis.json)
    traj_recs = []
    for an_path in sorted(web.glob("*/analysis/analysis.json")):
        tid = an_path.parent.parent.name
        a = json.loads(an_path.read_text())
        row = by_id.get(tid, {})
        ang = a.get("angle") or {}
        tcr = a.get("tcr_cdr") or {}
        geo = a.get("geometry") or {}
        ent = (a.get("interface") or {}).get("contact_entropy") or {}
        cm = a.get("concerted_motion") or {}          # L2 contact network (absent for 7byd)
        ed = a.get("essential_dynamics") or {}
        regs = ((a.get("rmsf_profile") or {}).get("regions")) or {}
        # Drop ratios from an unreliable CDR decomposition (degenerate, e.g. 6vma;
        # or a contact-stage pair-count artifact, e.g. 6g9q_run2 whose inflated
        # denominator would otherwise dilute the per-complex replica average).
        degen = cdr_contacts.is_unreliable(tcr) is not None
        traj_recs.append({
            "id": tid,
            "pdb_id": tid.split("_", 1)[0],
            "hla_group": _hla_group(row.get("hla_allele")),
            "peptide_length": row.get("peptide_length"),
            "trbv": (row.get("trbv") or "").split("*")[0] or None,
            "host": row.get("host_species"),
            "duration_ns": row.get("duration_ns"),
            "crossing": (ang.get("crossing_deg") or {}).get("mean"),
            "incident": (ang.get("incident_deg") or {}).get("mean"),
            "reversed_polarity": ang.get("reversed_polarity") if ang.get("crossing_deg") else None,
            "reversed_polarity_fraction": ang.get("reversed_polarity_fraction") if ang.get("crossing_deg") else None,
            "pep_recog": None if degen else tcr.get("peptide_recognition_ratio"),
            "hla_restr": None if degen else tcr.get("hla_restriction_ratio"),
            "alpha_contrib": None if degen else tcr.get("alpha_contribution"),
            "bulge": geo.get("bulge_height_angstrom"),
            "entropy": ent.get("entropy_normalized"),
            "fnat": (a.get("fnat") or {}).get("fnat_mean"),
            "cosine": (a.get("essential_dynamics") or {}).get("pc1_cosine_content"),
            "rmsip": (a.get("essential_dynamics") or {}).get("subspace_rmsip"),
            "bsa_total": (a.get("bsa_decomposition") or {}).get("total"),
            # --- layered-map extras (raw; encoded later in _database_map) ---
            "rg_mean_nm": geo.get("rg_mean_nm"),
            "total_sasa_mean_nm2": geo.get("total_sasa_mean_nm2"),
            "cm_mean_corr": cm.get("mean_correlation"),
            "cm_pep_tcr_coupling": cm.get("peptide_tcr_coupling"),
            "cm_modularity": cm.get("modularity"),
            "cm_collectivity": cm.get("collectivity"),
            "cm_n_communities": cm.get("n_communities"),
            "cm_n_contact_edges": cm.get("n_contact_edges"),
            "cm_n_residues": cm.get("n_residues"),
            "rmsf_iface": _rmsf_region_mean(regs, _RMSF_IFACE),
            "rmsf_cdr3": _rmsf_region_mean(regs, ("CDR3α", "CDR3β")),
            "rmsf_fr": _rmsf_region_mean(regs, ("TCR-FRα", "TCR-FRβ")),
            "rmsf_peptide": (regs.get("peptide") or {}).get("mean"),
            "rmsf_groove": _rmsf_region_mean(regs, ("α1 helix", "α2 helix")),
            "ed_pc1_var_frac": ed.get("pc1_variance_fraction"),
            "pep_table": (a.get("interface") or {}).get("peptide_table") or [],
        })

    # Collapse replicas: one record per complex (pdb_id). Average numeric metrics
    # over the replicas that have them; take categoricals from the primary (run3).
    by_pdb = defaultdict(list)
    for r in traj_recs:
        by_pdb[r["pdb_id"]].append(r)

    def _primary(group):
        # Same rule as build_manifest._primary_replica: prefer run3, else the
        # longest run — so the catalogue and Explore pick the same representative.
        return next((r for r in group if r["id"].endswith("_run3")),
                    max(group, key=lambda r: (r.get("duration_ns") or 0)))

    recs = []
    pos_tcr = defaultdict(list)   # position -> [tcr_contact] for DOMINANT_LENGTH peptides
    pos_hla = defaultdict(list)
    pos_sasa = defaultdict(list)  # position -> [sasa_nm2] (buried = anchor; the HLA-burial axis)
    for group in by_pdb.values():
        prim = _primary(group)
        rec = {k: prim.get(k) for k in _CATEG}
        rec["id"] = prim["id"]
        for k in _NUMERIC + _MAP_EXTRA:
            vals = [r[k] for r in group if r.get(k) is not None]
            rec[k] = round(sum(vals) / len(vals), 4) if vals else None
        # Docking polarity: consensus over the replicas' booleans (tie -> mean reversed
        # fraction), not the primary alone.
        pol, pol_agree, _ = consensus.majority_bool(
            [r.get("reversed_polarity") for r in group], primary=prim.get("reversed_polarity"),
            weights=[r.get("reversed_polarity_fraction") for r in group])
        rec["reversed_polarity"] = pol
        rec["reversed_polarity_reproducible"] = pol_agree
        # Inter-replica reproducibility: keep the mean (above) but also expose the
        # per-metric spread + a reproducibility flag instead of silently discarding
        # the agreement signal (docs/research/2026-06-16_replica_aggregation.md).
        rec["n_replicas"] = len(group)
        rstats = {}
        for k in _REPLICA_METRICS:
            st = _replica_stat([r.get(k) for r in group], REPLICA_TOL[k])
            if st is not None:
                rstats[k] = st
        rec["replica_stats"] = rstats
        recs.append(rec)
        # conserved-position map uses the primary replica only (no replica double-weighting)
        if prim.get("peptide_length") == DOMINANT_LENGTH:
            for p in prim["pep_table"]:
                if p.get("tcr_contact") is not None:
                    pos_tcr[p["position"]].append(p["tcr_contact"])
                if p.get("hla_contact") is not None:
                    pos_hla[p["position"]].append(p["hla_contact"])
                if p.get("sasa_nm2") is not None:
                    pos_sasa[p["position"]].append(p["sasa_nm2"])

    def col(key):
        return [r[key] for r in recs if r.get(key) is not None]

    def grouped(group_key, value_key, top=None, min_n=MIN_GROUP):
        g = defaultdict(list)
        for r in recs:
            if r.get(group_key) is not None and r.get(value_key) is not None:
                g[r[group_key]].append(r[value_key])
        items = [(k, _stat(v)) for k, v in g.items() if len(v) >= min_n]
        items.sort(key=lambda kv: kv[1]["n"], reverse=True)
        if top:
            items = items[:top]
        return {k: st for k, st in items}

    def counts(key):
        return dict(Counter(r[key] for r in recs if r.get(key) is not None).most_common())

    def repro_summary():
        """Dataset-level reproducibility for the Explore disclosure: replication
        coverage + per-metric fraction of replicated complexes whose replicas agree."""
        out = {"n_complexes": len(recs),
               "n_replicated": sum(1 for r in recs if (r.get("n_replicas") or 1) >= 2)}
        for k in _REPLICA_METRICS:
            flags = [r["replica_stats"][k]["reproducible"] for r in recs
                     if (r.get("replica_stats", {}).get(k) or {}).get("reproducible") is not None]
            if flags:
                out[k] = {"n": len(flags), "frac_reproducible": round(sum(flags) / len(flags), 3),
                          "tolerance": REPLICA_TOL[k]}
        # Docking polarity is a genuine binary geometry (forward vs reverse), not a bin
        # cut out of a continuous scalar, so replica agreement is meaningful for it.
        flags = [r["reversed_polarity_reproducible"] for r in recs
                 if r.get("reversed_polarity_reproducible") is not None]
        if flags:
            out["reversed_polarity"] = {"n": len(flags),
                                        "frac_reproducible": round(sum(flags) / len(flags), 3)}
        return out

    return {
        "n": len(recs),
        "n_trajectories": len(rows),
        "global": {
            "crossing": _stat(col("crossing")), "crossing_values": col("crossing"),
            # directed crossing angles split by docking polarity (disjoint): the
            # forward set sits in 0-90, the (few) reverse-polarity dockings in 90-180,
            # so the histogram shows them distinctly in the upper half.
            "crossing_values_forward": [r["crossing"] for r in recs
                                        if r.get("crossing") is not None and not r.get("reversed_polarity")],
            "crossing_values_reversed": [r["crossing"] for r in recs
                                         if r.get("crossing") is not None and r.get("reversed_polarity")],
            "reversed_polarity_n": sum(1 for r in recs if r.get("reversed_polarity") is True),
            "incident": _stat(col("incident")),
            # Per-complex inter-replica RANGE for each descriptor (2 x half-range).
            # Published alongside the values themselves so a reader can see how far the
            # three independent runs actually disagreed, rather than being handed a
            # boolean that only says whether they landed in the same bin.
            **{f"{k}_replica_range_values":
               [round(2.0 * r["replica_stats"][k]["spread"], 4) for r in recs
                if (r.get("replica_stats", {}).get(k) or {}).get("n", 0) >= 2]
               for k in _REPLICA_METRICS},
            "pep_recog": _stat(col("pep_recog")), "pep_recog_values": col("pep_recog"),
            "hla_restr": _stat(col("hla_restr")), "hla_restr_values": col("hla_restr"),
            "bulge": _stat(col("bulge")), "bulge_values": col("bulge"),
            "entropy": _stat(col("entropy")),
            "fnat": _stat(col("fnat")), "fnat_values": col("fnat"),
            "cosine": _stat(col("cosine")), "cosine_values": col("cosine"),
            "rmsip": _stat(col("rmsip")), "rmsip_values": col("rmsip"),
            "bsa_total": _stat(col("bsa_total")), "bsa_total_values": col("bsa_total"),
            "peptide_length": counts("peptide_length"),
            "alpha_contrib": _stat(col("alpha_contrib")), "alpha_contrib_values": col("alpha_contrib"),
            "host": counts("host"),
            "reproducibility": repro_summary(),
        },
        "map": _database_map(recs),
        "crossing_by_hla": grouped("hla_group", "crossing", top=TOP_GROUPS),
        "crossing_by_length": grouped("peptide_length", "crossing"),
        "bulge_by_length": grouped("peptide_length", "bulge"),
        "pep_recog_by_trbv": grouped("trbv", "pep_recog", top=TOP_GROUPS),
        "conserved_peptide": {
            "length": DOMINANT_LENGTH,
            "n": max((len(v) for v in pos_tcr.values()), default=0),
            "positions": sorted(pos_tcr),
            "tcr_mean": {p: round(sum(v) / len(v), 3) for p, v in sorted(pos_tcr.items())},
            "hla_mean": {p: round(sum(v) / len(v), 3) for p, v in sorted(pos_hla.items())},
            # mean buried SASA (nm^2) per position — the HLA-burial axis on the
            # conserved chart (low SASA = buried in the groove = anchor). Replaces
            # the saturated max-to-any hla_contact, recovering the P2/POmega signal.
            "sasa_mean": {p: round(sum(v) / len(v), 4) for p, v in sorted(pos_sasa.items())},
        },
    }


def main(argv=None) -> int:
    agg = build()
    (config.WEB_DATA / "aggregate.json").write_text(json.dumps(agg))
    g = agg["global"]
    print(f"Aggregate: {agg['n']} systems. crossing n={g['crossing']['n']} "
          f"median={g['crossing'].get('median')}; HLA groups={len(agg['crossing_by_hla'])}; "
          f"lengths={list(agg['crossing_by_length'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
