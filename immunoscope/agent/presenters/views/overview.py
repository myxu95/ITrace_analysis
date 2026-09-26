"""Overview view - system TL;DR."""

import csv
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from ..base import Presenter, RenderContext, PresenterError
from ..loaders import read_json
from .. import register_view


_PAIR_LABEL_RE = re.compile(r"([A-Za-z]):(\d+)")


@register_view("overview")
class OverviewPresenter(Presenter):
    """
    System overview with TL;DR.

    Provides high-level summary of system metadata, completed modules,
    and key findings. Intended as the first call when Agent lands on
    a case directory.
    """

    view_name = "overview"
    # D-B1: surfaces Complex-level RMSD / COM-distance scalars (static)
    # alongside TICA / FEL / cluster summaries (dynamic) for the system TL;DR.
    spatial_layer = "complex"
    sub_flavors_served = ("static", "dynamic")
    max_length = 2000

    def render(self, context: RenderContext) -> str:
        """Render overview."""
        locator = context.locator
        case_id = locator.get_case_id()

        lines = [f"# Overview: {case_id}\n"]

        # System metadata
        lines.append(self._render_system_metadata(locator))

        # Module status
        lines.append(self._render_modules(locator))

        # TL;DR
        lines.append(self._render_tldr(locator))

        # Cross-modal highlights (residues surfacing across multiple analyses)
        cross_modal = self._render_cross_modal(locator)
        if cross_modal:
            lines.append(cross_modal)

        # Next steps
        lines.append(self._render_next_steps())

        output = "\n".join(lines)
        return self._cap_length(output, "Use specific views for details")

    def _render_system_metadata(self, locator) -> str:
        """Render system metadata section."""
        # Try to load biological identity
        identity_path = locator.find_file(
            "identity",
            "analysis/identity/biological_identity.json"
        )
        if not identity_path:
            identity_path = locator.find_file(
                "identity",
                "biological_identity.json"
            )

        identity = read_json(identity_path)

        if not identity:
            return "## System\n\n_Metadata not available_\n"

        lines = ["## System\n"]

        # Peptide (field: peptide_identity)
        peptide = identity.get("peptide_identity", identity.get("peptide", {}))
        if isinstance(peptide, dict):
            seq = peptide.get("sequence", "")
            length = peptide.get("length", len(seq) if seq else 0)
            if seq:
                lines.append(f"- **Peptide**: {seq} ({length} aa)")
            else:
                lines.append(f"- **Peptide**: N/A")
        else:
            lines.append(f"- **Peptide**: {peptide}")

        # HLA (field: hla_identity)
        hla = identity.get("hla_identity", identity.get("hla", {}))
        if isinstance(hla, dict):
            allele = hla.get("best_candidate_allele", hla.get("allele", "N/A"))
            locus = hla.get("best_locus", hla.get("locus", ""))
            lines.append(f"- **HLA**: {allele} (locus {locus})")
        else:
            lines.append(f"- **HLA**: {hla}")

        # TCR (field: tcr_identity)
        tcr = identity.get("tcr_identity", identity.get("tcr", {}))
        if isinstance(tcr, dict):
            alpha_v = tcr.get("alpha_v_gene", "N/A")
            beta_v = tcr.get("beta_v_gene", "N/A")
            lines.append(f"- **TCR**: {alpha_v} / {beta_v}")

            cdr3_alpha = tcr.get("cdr3_alpha_sequence", tcr.get("cdr3_alpha", "N/A"))
            cdr3_beta = tcr.get("cdr3_beta_sequence", tcr.get("cdr3_beta", "N/A"))
            lines.append(f"- **CDR3α**: {cdr3_alpha}")
            lines.append(f"- **CDR3β**: {cdr3_beta}")

        return "\n".join(lines) + "\n"

    def _render_modules(self, locator) -> str:
        """Render module status section."""
        # Load run_summary.json
        summary_path = locator.case_dir / "run_summary.json"
        summary = read_json(summary_path)

        if not summary:
            return "## Modules\n\n_Run summary not available_\n"

        lines = ["## Modules\n"]

        module_results = summary.get("module_results", {})
        for module, result in module_results.items():
            status = result.get("status", "unknown")
            icon = "✅" if status == "completed" else "❌"

            # Get one-line headline for each module
            headline = self._get_module_headline(locator, module, result)
            lines.append(f"{icon} **{module}**: {status}{headline}")

        return "\n".join(lines) + "\n"

    def _get_module_headline(
        self,
        locator,
        module: str,
        result: Dict[str, Any]
    ) -> str:
        """Get one-line headline for a module."""
        # Try to extract key metric from module's summary
        root = result.get("root")
        if not root:
            return ""

        root_path = Path(root)

        if module == "rrcs":
            summary_path = root_path / "analysis/interactions/rrcs/rrcs_summary.json"
            summary = read_json(summary_path)
            if summary:
                n_nonzero = summary.get("n_nonzero_pairs", 0)
                top_pairs = summary.get("top_pairs", [])
                if top_pairs:
                    top_rrcs = top_pairs[0].get("mean_rrcs", 0)
                    return f" — {n_nonzero} nonzero pairs, top RRCS {top_rrcs:.2f}"

        elif module == "bsa":
            summary_path = root_path / "analysis/interface/interface_summary.json"
            summary = read_json(summary_path)
            if summary:
                # BSA is nested: buried_surface_area.mean
                bsa_data = summary.get("buried_surface_area", {})
                if isinstance(bsa_data, dict):
                    mean_bsa = bsa_data.get("mean", 0)
                else:
                    mean_bsa = summary.get("mean_bsa", 0)
                return f" — mean BSA {mean_bsa:.1f} Å²"

        elif module == "rmsf":
            summary_path = root_path / "analysis/rmsf/rmsf_summary.json"
            summary = read_json(summary_path)
            if summary:
                # Try both field naming conventions
                tcr_mean = summary.get("tcr_mean_rmsf_angstrom",
                                       summary.get("tcr_mean_rmsf", 0))
                return f" — TCR mean {tcr_mean:.2f} Å"

        elif module == "inter_cluster":
            summary_path = root_path / "analysis/conformation/interface_clustering/interface_clustering_summary.json"
            summary = read_json(summary_path)
            if summary:
                n_clusters = summary.get("n_clusters", 0)
                dominant_frac = summary.get("dominant_cluster_fraction", 0)
                return f" — {n_clusters} clusters, dominant {dominant_frac*100:.0f}%"

        return ""

    def _render_tldr(self, locator) -> str:
        """Render TL;DR section with key findings."""
        lines = ["## TL;DR\n"]

        # Try to extract top hotspot
        rrcs_root = locator.get_module_root("rrcs")
        if rrcs_root:
            summary_path = rrcs_root / "analysis/interactions/rrcs/rrcs_summary.json"
            summary = read_json(summary_path)
            if summary:
                top_pairs = summary.get("top_pairs", [])
                if top_pairs:
                    top = top_pairs[0]
                    residues = top.get("residues", "")
                    mean_rrcs = top.get("mean_rrcs", 0)
                    occupancy = top.get("nonzero_fraction", 0)
                    lines.append(
                        f"- **Strongest hotspot**: {residues}, "
                        f"mean RRCS {mean_rrcs:.2f}, occupancy {occupancy*100:.0f}%"
                    )

        # Try to extract dominant cluster
        cluster_root = locator.get_module_root("inter_cluster")
        if cluster_root:
            summary_path = cluster_root / "analysis/conformation/interface_clustering/interface_clustering_summary.json"
            summary = read_json(summary_path)
            if summary:
                n_clusters = summary.get("n_clusters", 0)
                dominant_frac = summary.get("dominant_cluster_fraction", 0)
                lines.append(
                    f"- **Conformational state**: {n_clusters} clusters, "
                    f"dominant basin {dominant_frac*100:.0f}% of frames"
                )

        if len(lines) == 1:
            lines.append("_Key findings not available_")

        return "\n".join(lines) + "\n"

    def _render_cross_modal(self, locator) -> str:
        """Surface residues appearing across multiple analyses (≥2 of RRCS / RMSF / dominant cluster)."""
        rrcs_top, rrcs_meta = self._top_rrcs_residues(locator, top_n=20)
        rmsf_top, rmsf_meta = self._top_rmsf_residues(locator, top_n=20)
        cluster_top = self._dominant_cluster_residues(locator)

        # Need at least two signals to do cross-modal at all
        signals_present = sum(1 for s in (rrcs_top, rmsf_top, cluster_top) if s)
        if signals_present < 2:
            return ""

        # Score each residue by how many analyses surface it
        scored: List[Tuple[str, int, Dict[str, Any]]] = []
        candidates: Set[str] = rrcs_top | rmsf_top | cluster_top
        for residue in candidates:
            evidence: List[str] = []
            count = 0
            if residue in rrcs_top:
                count += 1
                meta = rrcs_meta.get(residue, {})
                evidence.append(
                    f"RRCS rank {meta.get('rank', '?')} "
                    f"(mean {meta.get('mean_rrcs', 0):.2f})"
                )
            if residue in rmsf_top:
                count += 1
                meta = rmsf_meta.get(residue, {})
                evidence.append(
                    f"RMSF {meta.get('rmsf', 0):.2f} Å "
                    f"(rank {meta.get('rank', '?')})"
                )
            if residue in cluster_top:
                count += 1
                evidence.append("in dominant cluster signature")
            if count >= 2:
                scored.append((residue, count, {"evidence": evidence}))

        if not scored:
            return ""

        # Sort by count desc, then prefer ones that hit RRCS+RMSF (action-relevant)
        scored.sort(key=lambda t: (-t[1], t[0]))
        scored = scored[:5]

        lines = ["## Cross-modal highlights\n"]
        lines.append(
            "_Residues surfacing across multiple analyses "
            "(≥ 2 of: top RRCS / top RMSF / dominant cluster signature). "
            "These are the strongest cross-validated targets._\n"
        )
        for residue, _count, payload in scored:
            evidence_str = "; ".join(payload["evidence"])
            lines.append(f"- **{residue}**: {evidence_str}")
        return "\n".join(lines) + "\n"

    @staticmethod
    def _top_rrcs_residues(locator, top_n: int = 20) -> Tuple[Set[str], Dict[str, Dict[str, Any]]]:
        """Return set of TCR residue labels from the top-N RRCS pairs, plus per-residue metadata."""
        rrcs_root = locator.get_module_root("rrcs")
        if not rrcs_root:
            return set(), {}
        path = rrcs_root / "analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv"
        if not path.exists():
            return set(), {}

        result: Set[str] = set()
        meta: Dict[str, Dict[str, Any]] = {}
        try:
            with path.open(encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
        except OSError:
            return set(), {}

        def _rrcs(row: Dict[str, str]) -> float:
            try:
                return float(row.get("mean_rrcs", 0) or 0)
            except (TypeError, ValueError):
                return 0.0

        rows.sort(key=_rrcs, reverse=True)
        rank = 0
        for row in rows:
            label = (row.get("tcr_residue_label") or "").strip()
            if not label or label in result:
                continue
            rank += 1
            result.add(label)
            meta[label] = {"rank": rank, "mean_rrcs": _rrcs(row)}
            if rank >= top_n:
                break
        return result, meta

    @staticmethod
    def _top_rmsf_residues(locator, top_n: int = 20) -> Tuple[Set[str], Dict[str, Dict[str, Any]]]:
        """Return set of TCR residue labels from the top-N RMSF rows, plus metadata."""
        rmsf_root = locator.get_module_root("rmsf")
        if not rmsf_root:
            return set(), {}
        path = rmsf_root / "analysis/rmsf/residue_rmsf.csv"
        if not path.exists():
            return set(), {}

        try:
            with path.open(encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
        except OSError:
            return set(), {}

        def _rmsf(row: Dict[str, str]) -> float:
            try:
                return float(row.get("rmsf_angstrom", 0) or 0)
            except (TypeError, ValueError):
                return 0.0

        # Limit to TCR side (have tcr_chain or tcr_region) for design-relevance
        tcr_rows = [r for r in rows if (r.get("tcr_chain") or r.get("tcr_region"))]
        tcr_rows.sort(key=_rmsf, reverse=True)

        result: Set[str] = set()
        meta: Dict[str, Dict[str, Any]] = {}
        rank = 0
        for row in tcr_rows:
            resname = (row.get("resname") or "").strip().upper()
            resid = (row.get("resid") or "").strip()
            if not resname or not resid:
                continue
            label = f"{resname}{resid}"
            if label in result:
                continue
            rank += 1
            result.add(label)
            meta[label] = {"rank": rank, "rmsf": _rmsf(row)}
            if rank >= top_n:
                break
        return result, meta

    @staticmethod
    def _dominant_cluster_residues(locator) -> Set[str]:
        """Return residue labels from the dominant cluster's signature pairs.

        We pick the cluster with the highest population_percent from summary_table.csv,
        then read its signature pair_labels from cluster_specific_signatures.csv.
        Pair labels look like 'E:98__A:155[region]' — we parse the TCR-side
        chain:resid and convert to a label using residue_rmsf.csv for the resname.
        Degrades silently if any input is missing.
        """
        cluster_root = locator.get_module_root("inter_cluster")
        if not cluster_root:
            return set()
        base = cluster_root / "analysis/conformation/interface_clustering"
        summary_path = base / "summary_table.csv"
        signatures_path = base / "cluster_specific_signatures.csv"
        if not summary_path.exists() or not signatures_path.exists():
            return set()

        try:
            with summary_path.open(encoding="utf-8") as f:
                summary_rows = list(csv.DictReader(f))
        except OSError:
            return set()

        if not summary_rows:
            return set()

        def _pop(row: Dict[str, str]) -> float:
            try:
                return float(row.get("population_percent", 0) or 0)
            except (TypeError, ValueError):
                return 0.0

        summary_rows.sort(key=_pop, reverse=True)
        dominant_id = (summary_rows[0].get("cluster_id") or "").strip()
        if not dominant_id:
            return set()

        # Map chain:resid → resname via residue_rmsf.csv (best available source)
        chain_resid_to_name: Dict[Tuple[str, str], str] = {}
        rmsf_root = locator.get_module_root("rmsf")
        if rmsf_root:
            rmsf_csv = rmsf_root / "analysis/rmsf/residue_rmsf.csv"
            if rmsf_csv.exists():
                try:
                    with rmsf_csv.open(encoding="utf-8") as f:
                        for row in csv.DictReader(f):
                            chain = (row.get("chain_id") or "").strip()
                            resid = (row.get("resid") or "").strip()
                            resname = (row.get("resname") or "").strip().upper()
                            if chain and resid and resname:
                                chain_resid_to_name.setdefault((chain, resid), resname)
                except OSError:
                    pass

        result: Set[str] = set()
        try:
            with signatures_path.open(encoding="utf-8") as f:
                for row in csv.DictReader(f):
                    if (row.get("cluster_id") or "").strip() != dominant_id:
                        continue
                    pair_label = row.get("pair_label", "") or ""
                    for chain, resid in _PAIR_LABEL_RE.findall(pair_label):
                        resname = chain_resid_to_name.get((chain, resid))
                        if resname:
                            result.add(f"{resname}{resid}")
                        else:
                            result.add(f"{chain}{resid}")
        except OSError:
            return set()
        return result

    def _render_next_steps(self) -> str:
        """Render next steps section."""
        lines = [
            "## Next views\n",
            "- `hotspots` — ranked RRCS contact pairs",
            "- `interface` — BSA statistics and composition",
            "- `flexibility` — RMSF by region",
            "- `clustering` — conformational states",
            "- `quality` — trajectory quality metrics",
            "- `residue` — single-residue drill-down (requires filter: residue)",
            "- `pair` — single-pair dynamics (requires filters: residue1, residue2)",
        ]
        return "\n".join(lines)
