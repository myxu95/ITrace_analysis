"""
分析数据格式化器

从ImmunoScope分析目录提取关键信息，格式化为LLM可理解的统一JSON格式
"""
import json
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional
import logging

from immunoscope.analysis.features import CaseLocator, FeatureSet, compute_feature

from .anchor_pocket import lookup_anchor_pockets
from .schemas import AnalysisDataSummary
from .task_spec import TaskSpec

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class AnalysisDataFormatter:
    """分析数据格式化器"""

    def __init__(self, analysis_dir: str):
        """
        初始化格式化器

        Args:
            analysis_dir: 分析结果目录（如 output/5c0a_run2_full_analysis）
        """
        self.analysis_dir = Path(analysis_dir)
        if not self.analysis_dir.exists():
            raise FileNotFoundError(f"分析目录不存在: {analysis_dir}")

        # Lazy-built {residue_label: rmsf_angstrom} cache, populated on first
        # _get_residue_rmsf call. None until first access.
        self._rmsf_lookup: Optional[Dict[str, float]] = None

        logger.info(f"初始化分析数据格式化器: {analysis_dir}")

    def format(
        self,
        max_candidates: int = 20,
        task: Optional[TaskSpec] = None,
    ) -> AnalysisDataSummary:
        """
        格式化分析数据

        Args:
            max_candidates: 最大候选残基数量
            task: which (subject, objective, target_subinterface) the
                recommendation is for. The candidate filter routes on this:
                TCR-affinity keeps the legacy CDR/non_cdr/framework filter;
                peptide-presentation keeps peptide-side residues with HLA
                or beta2m partners. Defaults to TCR-affinity.

        Returns:
            AnalysisDataSummary对象
        """
        if task is None:
            task = TaskSpec.tcr_affinity()

        logger.info("开始格式化分析数据...")
        logger.info(f"任务: {task.slug}")

        # 1. 提取系统元数据
        metadata = self._extract_metadata()

        # 1a. 任务为呈递时，附加 anchor pocket 元数据
        # The LLM needs to know which peptide positions are anchors *before*
        # it picks a candidate to mutate; without this, P2/PΩ mutations get
        # proposed as freely as P5 mutations.
        if task.subject == "peptide" and task.objective == "presentation":
            metadata["anchor_pockets"] = self._anchor_pocket_summary(metadata)

        # 2. 提取分析摘要
        analysis_summary = {
            "quality": self._extract_quality(),
            "interface": self._extract_interface(),
            "hotspots": self._extract_hotspots(),
            "flexibility": self._extract_flexibility(),
            "conformation": self._extract_conformation()
        }

        # 3. 提取候选残基
        candidate_residues = self._extract_candidate_residues(max_candidates, task=task)

        # 4. 构建摘要对象
        summary = AnalysisDataSummary(
            system_id=self.analysis_dir.name,
            metadata=metadata,
            analysis_summary=analysis_summary,
            candidate_residues=candidate_residues
        )

        # 5. 检查大小
        tokens = summary.estimate_tokens()
        logger.info(f"格式化完成，预估token数: {tokens}")

        if tokens > 12000:  # 留一些余量
            logger.warning(f"数据较大({tokens} tokens)，可能需要精简")

        return summary

    @staticmethod
    def _anchor_pocket_summary(metadata: Dict) -> List[Dict]:
        """Build a JSON-serializable anchor pocket summary from metadata.

        Reads the HLA allele and peptide sequence from the metadata dict
        and looks them up against the curated anchor table. Returns an
        empty list when neither side identifies an entry — the prompt
        builder then surfaces "no anchor info available" rather than
        fabricating chemistry.

        We use the peptide *sequence length* to drive PΩ position. If the
        peptide sequence is missing, we fall through to ``None`` and let
        the lookup return [] (rather than guessing 9-mer).
        """
        hla_allele = metadata.get("hla_allele") or metadata.get("hla")
        peptide_seq = (metadata.get("peptide") or "").strip()
        peptide_length = len(peptide_seq) if peptide_seq else None

        entries = lookup_anchor_pockets(hla_allele, peptide_length)
        if not entries:
            return []

        # Serialize — dataclass → plain dict, preserving the residue chemistry
        # tuples as lists so the JSON encoder can handle them.
        return [
            {
                "position": e.position,
                "allele": e.allele,
                "peptide_length": e.peptide_length,
                "role": e.role,
                "preferred_residues": list(e.preferred_residues),
                "tolerated_residues": list(e.tolerated_residues),
                "dispreferred_residues": list(e.dispreferred_residues),
                "notes": e.notes,
                "source": e.source,
                # When the peptide sequence is known, also include the current
                # residue at this anchor position so the LLM doesn't have to
                # re-derive it. peptide_seq is 1-indexed in the literature;
                # convert to Python's 0-indexed slice.
                "current_residue": (
                    peptide_seq[e.position - 1]
                    if peptide_seq and 1 <= e.position <= len(peptide_seq)
                    else None
                ),
            }
            for e in entries
        ]

    def _extract_metadata(self) -> Dict:
        """提取系统元数据"""
        identity_file = self.analysis_dir / "analysis" / "identity" / "analysis" / "identity" / "biological_identity.json"

        if not identity_file.exists():
            logger.warning(f"身份文件不存在: {identity_file}")
            return {}

        with open(identity_file, 'r') as f:
            identity = json.load(f)

        # 提取关键信息
        peptide_info = identity.get("peptide_identity", {})
        tcr_info = identity.get("tcr_identity", {})
        hla_info = identity.get("hla_identity", {})

        return {
            "peptide": peptide_info.get("sequence", ""),
            "hla": hla_info.get("best_locus", ""),
            "hla_allele": hla_info.get("best_candidate_allele", ""),
            "tcr_alpha_v": tcr_info.get("alpha_v_gene", ""),
            "tcr_beta_v": tcr_info.get("beta_v_gene", ""),
            "cdr3_alpha": tcr_info.get("cdr3_alpha_sequence", ""),
            "cdr3_beta": tcr_info.get("cdr3_beta_sequence", "")
        }

    def _extract_quality(self) -> Dict:
        """提取质量指标"""
        quality_file = self.analysis_dir / "analysis" / "quality" / "quality_report.json"

        if not quality_file.exists():
            logger.warning(f"质量文件不存在: {quality_file}")
            return {}

        with open(quality_file, 'r') as f:
            quality = json.load(f)

        rmsd_metrics = quality.get("rmsd_metrics", {})

        return {
            "is_converged": rmsd_metrics.get("is_converged", False),
            "rmsd_mean": round(rmsd_metrics.get("mean_rmsd", 0), 3),
            "rmsd_std": round(rmsd_metrics.get("std_rmsd", 0), 3),
            "overall_grade": quality.get("overall_grade", "N/A"),
            "is_qualified": quality.get("is_qualified", False)
        }

    def _extract_interface(self) -> Dict:
        """提取界面统计"""
        bsa_file = self.analysis_dir / "analysis" / "bsa" / "analysis" / "interface" / "interface_summary.json"

        if not bsa_file.exists():
            logger.warning(f"BSA文件不存在: {bsa_file}")
            return {}

        with open(bsa_file, 'r') as f:
            bsa = json.load(f)

        return {
            "bsa_mean": round(bsa.get("mean_bsa", 0), 1),
            "bsa_std": round(bsa.get("std_bsa", 0), 1),
            "n_contact_pairs": bsa.get("n_contact_pairs", 0)
        }

    def _extract_hotspots(self) -> Dict:
        """提取热点残基信息"""
        # 使用annotated CSV而不是summary JSON，因为CSV有完整的字段
        rrcs_csv = self.analysis_dir / "analysis" / "rrcs" / "analysis" / "interactions" / "rrcs" / "annotated_rrcs_pair_summary.csv"

        if not rrcs_csv.exists():
            logger.warning(f"RRCS CSV文件不存在: {rrcs_csv}")
            return {"top_pairs": [], "top_regions": []}

        try:
            df = pd.read_csv(rrcs_csv)
        except Exception:
            # 广义捕获：CSV 损坏不应阻塞整个报告生成，但 traceback 要保留以便排查
            logger.exception("读取RRCS CSV失败，返回空 hotspots")
            return {"top_pairs": [], "top_regions": []}

        # 提取top pairs（按mean_rrcs排序）
        top_pairs = []
        if not df.empty and "mean_rrcs" in df.columns:
            df_sorted = df.sort_values("mean_rrcs", ascending=False).head(10)

            for _, row in df_sorted.iterrows():
                top_pairs.append({
                    "tcr_residue": row.get("tcr_residue_label", ""),
                    "partner_residue": row.get("partner_residue_label", ""),
                    "partner_component": row.get("partner_component", ""),
                    "mean_rrcs": round(row.get("mean_rrcs", 0), 2),
                    "nonzero_fraction": round(row.get("rrcs_nonzero_fraction", 0), 2),
                    "tcr_region": row.get("tcr_region_detailed", row.get("tcr_region", ""))
                })

        # 提取top regions（从region summary读取）
        region_summary_file = self.analysis_dir / "analysis" / "rrcs" / "analysis" / "interactions" / "rrcs" / "region_interaction_summary.json"
        top_regions = []

        if region_summary_file.exists():
            try:
                with open(region_summary_file, 'r') as f:
                    region_data = json.load(f)

                for region in region_data.get("top_regions", [])[:5]:
                    top_regions.append({
                        "tcr_region": region.get("tcr_region", ""),
                        "partner_component": region.get("partner_component", ""),
                        "mean_rrcs_sum": round(region.get("mean_rrcs_sum", 0), 2),
                        "n_pairs": region.get("n_pairs", 0)
                    })
            except Exception:
                logger.exception("读取region summary失败（hotspots top_regions 留空）")

        return {
            "top_pairs": top_pairs,
            "top_regions": top_regions
        }

    def _extract_flexibility(self) -> Dict:
        """提取柔性统计"""
        rmsf_file = self.analysis_dir / "analysis" / "rmsf" / "analysis" / "rmsf" / "rmsf_summary.json"

        if not rmsf_file.exists():
            logger.warning(f"RMSF文件不存在: {rmsf_file}")
            return {}

        with open(rmsf_file, 'r') as f:
            rmsf = json.load(f)

        # 提取高柔性残基（从residue_rmsf.csv读取，因为summary JSON没有这个字段）
        high_flex = []
        residue_rmsf_csv = self.analysis_dir / "analysis" / "rmsf" / "analysis" / "rmsf" / "residue_rmsf.csv"

        if residue_rmsf_csv.exists():
            try:
                df = pd.read_csv(residue_rmsf_csv)
                if not df.empty and "rmsf_angstrom" in df.columns:
                    # 按RMSF降序排序，取top 5
                    df_sorted = df.sort_values("rmsf_angstrom", ascending=False).head(5)

                    for _, row in df_sorted.iterrows():
                        # 构建残基标签
                        resname = row.get("resname", "")
                        resid = row.get("resid", "")
                        residue_label = f"{resname}{resid}" if resname and resid else ""

                        # 获取region（优先使用tcr_region_detailed）
                        region = row.get("tcr_region_detailed") or row.get("tcr_region") or row.get("region_group", "")

                        high_flex.append({
                            "residue": residue_label,
                            "rmsf": round(row.get("rmsf_angstrom", 0), 2),
                            "region": region
                        })
            except Exception:
                logger.exception("读取residue RMSF CSV失败（high_flexibility_residues 留空）")

        return {
            "tcr_rmsf_mean": round(rmsf.get("tcr_mean_rmsf", 0), 2),
            "high_flexibility_residues": high_flex
        }

    def _extract_conformation(self) -> Dict:
        """提取构象信息"""
        cluster_file = self.analysis_dir / "analysis" / "inter_cluster" / "analysis" / "conformation" / "interface_clustering" / "interface_clustering_summary.json"

        if not cluster_file.exists():
            logger.warning(f"聚类文件不存在: {cluster_file}")
            return {}

        with open(cluster_file, 'r') as f:
            cluster = json.load(f)

        return {
            "n_clusters": cluster.get("n_clusters", 0),
            "dominant_population": round(cluster.get("dominant_cluster_population", 0), 3)
        }

    def _extract_candidate_residues(
        self,
        max_candidates: int,
        task: Optional[TaskSpec] = None,
    ) -> List[Dict]:
        """提取候选残基列表。

        优先使用 feature layer（带 chemistry/risk 注释，按 rrcs_contribution 排序）；
        失败或无候选时回退到 RRCS pair summary CSV。区分两种"没有结果"以便排错：
        - feature_candidates is None: feature layer 真正失败（已 logger.exception 记录 traceback）
        - feature_candidates == []   : feature layer 工作正常但所有候选都被任务过滤剔除了

        `task` controls which residues survive the filter. The legacy
        TCR-affinity path keeps only CDR/framework residues on the TCR side;
        the peptide-presentation path (L3, 2026-05-27) keeps peptide-side
        residues whose partner is HLA or beta2m.
        """
        if task is None:
            task = TaskSpec.tcr_affinity()

        feature_candidates = self._extract_feature_candidate_residues(max_candidates, task=task)
        if feature_candidates:
            logger.info(f"基于feature layer提取了 {len(feature_candidates)} 个候选残基")
            return feature_candidates

        if feature_candidates is None:
            logger.warning("feature layer 提取失败（见上方 traceback），回退到 RRCS pair summary CSV")
        else:
            logger.info(
                "feature layer 未产出符合任务 %s 的候选（全被过滤），回退到 RRCS pair summary CSV",
                task.slug,
            )

        # 从RRCS详细数据中提取
        rrcs_csv = self.analysis_dir / "analysis" / "rrcs" / "analysis" / "interactions" / "rrcs" / "annotated_rrcs_pair_summary.csv"

        if not rrcs_csv.exists():
            logger.warning(f"RRCS详细文件不存在: {rrcs_csv}")
            return []

        df = pd.read_csv(rrcs_csv)

        # 按mean_rrcs排序
        df = df.sort_values('mean_rrcs', ascending=False)

        # 提取候选残基。Subject 取决于任务：
        # - TCR-affinity: 用 tcr_residue_label/tcr_chain/tcr_region（旧路径）
        # - peptide-presentation: 用 peptide-side 残基（L3 才接入；暂时退化为不返回结果，
        #   让 caller 看到 fallback 也没有候选，而不是误用 TCR 字段）
        if task.subject == "tcr":
            subject_label_col = "tcr_residue_label"
            subject_chain_col = "tcr_chain"
            subject_region_col = "tcr_region"
        else:
            logger.warning(
                "RRCS CSV fallback 未实现 subject=%s 的列提取（peptide-side L3+ TODO），返回空列表",
                task.subject,
            )
            return []

        candidates = []
        seen_residues = set()

        for _, row in df.iterrows():
            subject_residue = row.get(subject_label_col, '')

            # 去重
            if subject_residue in seen_residues:
                continue

            seen_residues.add(subject_residue)

            candidates.append({
                "residue": subject_residue,
                "chain": row.get(subject_chain_col, ''),
                "region": row.get(subject_region_col, ''),
                # `rrcs_mean`: average RRCS over frames for this residue's
                # strongest contact pair. `occupancy`: fraction of frames
                # in which that pair had RRCS > 0. Both come from the
                # interaction-pair CSV.
                "rrcs_mean": round(row.get('mean_rrcs', 0), 2),
                "occupancy": round(row.get('rrcs_nonzero_fraction', 0), 2),
                "partner": f"{row.get('partner_residue_label', '')}@{row.get('partner_component', '')}",
                "rmsf": self._get_residue_rmsf(subject_residue)
            })

            if len(candidates) >= max_candidates:
                break

        logger.info(f"提取了 {len(candidates)} 个候选残基")
        return candidates

    def _extract_feature_candidate_residues(
        self,
        max_candidates: int,
        task: Optional[TaskSpec] = None,
    ) -> Optional[List[Dict]]:
        """Extract mutation candidates from the design feature layer.

        Per D-B5 (2026-05-26) the composite ``design_priority_score`` is no
        longer surfaced to the LLM; the LLM is responsible for contextual
        reasoning over the spatial hierarchy. Candidates are returned sorted
        by ``rrcs_contribution`` (single deterministic signal, matching the
        Pure-RRCS baseline), with risk and chemistry annotations preserved
        as categorical evidence rather than score components.

        Return contract: ``None`` if the feature layer raised (real failure —
        caller should fall back); an empty list if the feature layer
        succeeded but produced no candidates matching the task filter (also
        triggers fallback, but distinguishes the two cases for debugging).

        ``task`` selects the candidate filter (TCR-side vs peptide-side).
        """
        if task is None:
            task = TaskSpec.tcr_affinity()

        try:
            locator = CaseLocator(self.analysis_dir)
            features = FeatureSet(case_id=locator.get_case_id())
            for name in (
                "contact_count",
                "rrcs_contribution",
                "chemistry_tags",
                "risk_flag",
                "persistence_profile",
            ):
                compute_feature(name, locator, features)
        except Exception:
            # Broad catch is intentional: we want to fall back rather than
            # crash the whole report when the feature layer has a bug or
            # missing data. Use logger.exception so the traceback is kept —
            # the previous version dropped it via f-string formatting and
            # made silent failures look like "no candidates".
            logger.exception("feature layer候选残基提取失败，将回退到RRCS排序")
            return None

        candidates: List[Dict] = []
        for item in features.top_by("rrcs_contribution", n=max_candidates * 3):
            if not self._is_design_candidate_for_task(item, task):
                continue
            candidates.append({
                "residue": item.residue.label(),
                "chain": item.residue.chain,
                "position": item.residue.resid,
                "current_aa": item.residue.resname,
                "region": item.region or "",
                # `rrcs_mean`: per-residue aggregated RRCS contribution.
                # `occupancy`: average contact persistence over frames
                # (per ResidueFeatures.persistence_mean = "avg occupancy
                # over frames"). Same semantic as the CSV fallback's
                # `occupancy` so downstream consumers can treat them
                # uniformly.
                "rrcs_mean": round(item.rrcs_contribution or 0.0, 2),
                "occupancy": round(item.persistence_mean or 0.0, 2),
                "persistence_profile": item.persistence_profile or "",
                "contact_count": item.contact_count or 0,
                "partner": item.partner_chemistry or "",
                "chemistry_tags": list(item.chemistry_tags[:8]),
                "risk_flags": list(item.risk_flags[:8]),
                "rmsf": self._get_residue_rmsf(item.residue.label()),
            })
            if len(candidates) >= max_candidates:
                break

        return candidates

    @staticmethod
    def _is_design_candidate_for_task(item, task: TaskSpec) -> bool:
        """Route the candidate filter on the TaskSpec subject."""
        if task.subject == "tcr":
            return AnalysisDataFormatter._is_tcr_design_candidate(item)
        if task.subject == "peptide":
            return AnalysisDataFormatter._is_peptide_design_candidate(item)
        # `hla` subject is reserved for a future HLA-engineering track;
        # _WIRED_COMBINATIONS would already have rejected unwired subjects
        # at TaskSpec construction, so this branch should be unreachable
        # in practice — leaving the explicit guard so a future subject
        # added to TaskSpec without wiring a filter trips an obvious error.
        raise NotImplementedError(
            f"No candidate filter for task.subject={task.subject!r}; "
            "wire one in AnalysisDataFormatter before enabling this task."
        )

    @staticmethod
    def _is_tcr_design_candidate(item) -> bool:
        """Return True when a feature row represents a TCR-side design site.

        Tightened on 2026-05-27 with the L2 TaskSpec unlock: the original
        filter accepted ANY item whose partner was HLA/peptide-not-TCR via
        the second branch, which silently swept peptide-side residues into
        the TCR-affinity candidate pool now that the upstream feature set
        contains both sides. We now require the region label to be either
        TCR-y or non-{peptide,hla,beta2m}, so peptide-side rows are
        excluded explicitly rather than relying on the FeatureSet being
        TCR-only.
        """
        region = (item.region or "").lower()
        partner = (item.partner_chemistry or "").lower()
        if "cdr" in region or "non_cdr" in region or "framework" in region:
            return True
        # Reject anything whose region clearly marks a non-TCR component.
        if any(tag in region for tag in ("peptide", "hla", "mhc", "beta2m", "b2m")):
            return False
        if ("hla" in partner or "peptide" in partner) and "tcr" not in partner:
            return True
        return False

    @staticmethod
    def _is_peptide_design_candidate(item) -> bool:
        """Return True when a feature row represents a peptide-side design site.

        For the peptide-presentation track the candidates are peptide residues
        whose principal contact partner is on the HLA / beta2m groove side.
        We match on the same fields used by the TCR filter — `region` (which
        carries the ``peptide`` component label per ComplexResidueSemanticAnnotator)
        and `partner_chemistry` (which encodes the partner's component label).

        Excludes anything with a TCR partner — that is a peptide-TCR contact,
        relevant to TCR-affinity not to presentation.
        """
        region = (item.region or "").lower()
        partner = (item.partner_chemistry or "").lower()
        is_peptide_side = "peptide" in region
        if not is_peptide_side:
            return False
        if "tcr" in partner:
            return False
        return ("hla" in partner) or ("mhc" in partner) or ("b2m" in partner) or ("beta2m" in partner)

    def _get_residue_rmsf(self, residue: str) -> float:
        """Return per-residue RMSF (Å) by residue label, or 0.0 if unknown.

        Loads ``residue_rmsf.csv`` lazily and caches a {label: rmsf} dict.
        Residue label format follows the rest of this module: ``RESNAME +
        RESID`` (e.g. ``ASP92``), matching how candidate residues are
        labeled both in the CSV fallback branch and in the feature-layer
        branch (``ResidueKey.label()``).
        """
        if self._rmsf_lookup is None:
            self._rmsf_lookup = self._build_rmsf_lookup()
        if not residue:
            return 0.0
        return self._rmsf_lookup.get(residue, 0.0)

    def _build_rmsf_lookup(self) -> Dict[str, float]:
        """Read residue_rmsf.csv into a {residue_label: rmsf_angstrom} dict.

        Indexes each residue under both label conventions used in this
        module so the lookup works for both candidate code paths:
        - ``"{RESNAME}{RESID}"``        (CSV-fallback branch, e.g. ASP92)
        - ``"{RESNAME}-{CHAIN}{RESID}"`` (feature-layer branch via
          ``ResidueKey.label()``, e.g. ASP-D92)

        Returns an empty dict (not a sentinel) when the file is missing or
        malformed; callers fall back to 0.0 via _get_residue_rmsf.
        """
        residue_rmsf_csv = (
            self.analysis_dir
            / "analysis" / "rmsf" / "analysis" / "rmsf" / "residue_rmsf.csv"
        )
        if not residue_rmsf_csv.exists():
            logger.warning(f"residue_rmsf.csv 不存在，RMSF 查询返回 0.0: {residue_rmsf_csv}")
            return {}
        try:
            df = pd.read_csv(residue_rmsf_csv)
        except Exception:
            logger.exception("读取 residue_rmsf.csv 失败，RMSF 查询返回 0.0")
            return {}
        if df.empty or "rmsf_angstrom" not in df.columns:
            return {}
        lookup: Dict[str, float] = {}
        for _, row in df.iterrows():
            resname = str(row.get("resname", "") or "")
            resid = row.get("resid", "")
            chain = str(row.get("chain_id", "") or "")
            try:
                resid_str = str(int(resid))
            except (TypeError, ValueError):
                resid_str = str(resid)
            if not resname or not resid_str:
                continue
            try:
                rmsf_value = float(row["rmsf_angstrom"])
            except (TypeError, ValueError):
                continue
            # Per-atom rows in the CSV can collide on the chain-less key
            # (same RESNAME+RESID in different chains). Take the first
            # write — RMSF differs slightly per atom anyway; CA backbone
            # is conventionally first in MDAnalysis dumps.
            lookup.setdefault(f"{resname}{resid_str}", rmsf_value)
            if chain:
                lookup.setdefault(f"{resname}-{chain}{resid_str}", rmsf_value)
        return lookup


# 使用示例
if __name__ == "__main__":
    formatter = AnalysisDataFormatter("output/5c0a_run2_full_analysis")
    summary = formatter.format()

    print(f"系统ID: {summary.system_id}")
    print(f"候选残基数: {len(summary.candidate_residues)}")
    print(f"预估tokens: {summary.estimate_tokens()}")

    # 保存为JSON
    with open("analysis_summary.json", "w") as f:
        f.write(summary.to_json_str())

    print("已保存到 analysis_summary.json")
