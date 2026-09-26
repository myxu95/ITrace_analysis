"""
推荐系统数据结构定义
"""
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from datetime import datetime


# ============================================================================
# 可解释性相关数据结构
# ============================================================================

class MetricScore(BaseModel):
    """指标评分"""
    value: float = Field(..., description="数值")
    percentile: float = Field(..., description="百分位（0-100）")
    interpretation: str = Field(..., description="解释（high/medium/low）")


class SelectionRationale(BaseModel):
    """残基选择理由"""
    rank: int = Field(..., description="在候选中的排名")
    total_candidates: int = Field(..., description="总候选数")

    key_metrics: Dict[str, MetricScore] = Field(
        ...,
        description="关键指标（rrcs, occupancy）—— D-B5 后不再包含复合 mutability 分数"
    )

    comparison_to_next: str = Field(
        ...,
        description="与下一名的对比（如：RRCS 高 2.2 倍）"
    )

    why_this_residue: str = Field(
        ...,
        description="选择这个残基的核心原因（1-2 句话）"
    )


class AminoAcidProperties(BaseModel):
    """氨基酸性质"""
    code: str = Field(..., description="单字母代码")
    name: str = Field(..., description="全名")
    type: str = Field(..., description="类型（aromatic/charged/polar/hydrophobic）")
    polarity: str = Field(..., description="极性（polar/nonpolar）")
    charge: str = Field(..., description="电荷（positive/negative/neutral）")
    size: str = Field(..., description="大小（small/medium/large）")


class MutationOption(BaseModel):
    """突变选项"""
    aa_code: str = Field(..., description="氨基酸代码")
    properties: AminoAcidProperties = Field(..., description="氨基酸性质")
    advantage: str = Field(..., description="优势（如：π-π stacking + H-bond）")
    expected_effect: str = Field(..., description="预期效果")
    confidence: str = Field(..., description="置信度（high/medium/low）")


class MutationRationale(BaseModel):
    """突变选择理由"""
    current_aa_properties: AminoAcidProperties = Field(
        ...,
        description="当前氨基酸性质"
    )

    suggested_mutations: Dict[str, MutationOption] = Field(
        ...,
        description="建议的突变及其优势"
    )

    rejected_options: Dict[str, str] = Field(
        default_factory=dict,
        description="被拒绝的突变选项及原因"
    )

    design_strategy: str = Field(
        ...,
        description="设计策略（如：增强芳香相互作用）"
    )


class InteractionPartner(BaseModel):
    """相互作用伙伴"""
    partner_residue: str = Field(..., description="伙伴残基（如：HLA-A*02:01 TYR99）")
    partner_chain: str = Field(..., description="链（HLA_alpha/HLA_beta/peptide）")
    interaction_type: str = Field(..., description="相互作用类型")

    distance: float = Field(..., description="平均距离（Å）")
    distance_std: float = Field(..., description="距离标准差（Å）")

    frequency: float = Field(..., description="接触频率（0-1）")
    persistence: float = Field(..., description="持续性（0-1）")

    strength: str = Field(..., description="强度（strong/moderate/weak）")


class DynamicsEvidence(BaseModel):
    """动力学证据"""
    rmsf: float = Field(..., description="RMSF 值（Å）")
    rmsf_interpretation: str = Field(..., description="柔性解释")

    contact_persistence: float = Field(..., description="接触持续性（0-1）")
    persistence_interpretation: str = Field(..., description="持续性解释")

    conformational_states: int = Field(..., description="构象状态数")
    dominant_state_population: float = Field(..., description="主导构象占比")


class CandidateComparison(BaseModel):
    """候选残基对比"""
    candidate_residue: str = Field(..., description="候选残基")
    rank: int = Field(..., description="排名")
    rrcs: float = Field(..., description="RRCS 值")
    occupancy: float = Field(..., description="接触持续度（0-1）")
    why_not_selected: str = Field(..., description="为什么没被选择")


class MutationComparison(BaseModel):
    """突变选项对比"""
    mutation: str = Field(..., description="突变（如：E→K）")
    score: float = Field(..., description="评分（0-1）")
    pros: List[str] = Field(..., description="优点")
    cons: List[str] = Field(..., description="缺点")
    selected: bool = Field(..., description="是否被选择")


class ComparisonAnalysis(BaseModel):
    """对比分析"""
    vs_other_candidates: List[CandidateComparison] = Field(
        default_factory=list,
        description="与其他候选的对比"
    )

    vs_other_mutations: List[MutationComparison] = Field(
        default_factory=list,
        description="与其他突变选项的对比"
    )

    summary: str = Field(..., description="对比总结")


# ============================================================================
# 原有数据结构（扩展）
# ============================================================================


class MutationRecommendation(BaseModel):
    """单个突变推荐"""

    residue: str = Field(..., description="残基标签（如ASP92）")
    chain: str = Field(..., description="链标识（TCR_alpha/TCR_beta）")
    position: int = Field(..., description="残基位置")
    region: str = Field(..., description="CDR区域（CDR1/CDR2/CDR3/Framework）")
    current_aa: str = Field(..., description="当前氨基酸（单字母）")
    suggested_mutations: List[str] = Field(..., description="建议的突变氨基酸列表")
    priority: str = Field(..., description="优先级（high/medium/low）")

    rationale: str = Field(..., description="设计理由（基于MD数据和文献）")

    expected_effects: Dict[str, str] = Field(
        default_factory=dict,
        description="预期效果（affinity/specificity/stability/expression）"
    )

    risks: List[str] = Field(
        default_factory=list,
        description="潜在风险列表"
    )

    validation_experiments: List[str] = Field(
        default_factory=list,
        description="建议的验证实验"
    )

    supporting_evidence: List[Dict] = Field(
        default_factory=list,
        description="支持证据（MD数据 + 文献引用）"
    )

    confidence: str = Field(..., description="置信度（high/medium/low）")

    # ========================================================================
    # 新增：可解释性字段
    # ========================================================================

    selection_rationale: Optional[SelectionRationale] = Field(
        None,
        description="选择理由：为什么选择这个残基"
    )

    mutation_rationale: Optional[MutationRationale] = Field(
        None,
        description="突变理由：为什么选择这些突变"
    )

    interaction_partners: List[InteractionPartner] = Field(
        default_factory=list,
        description="相互作用伙伴详情"
    )

    dynamics_evidence: Optional[DynamicsEvidence] = Field(
        None,
        description="动力学证据"
    )

    comparison_analysis: Optional[ComparisonAnalysis] = Field(
        None,
        description="对比分析"
    )


class RecommendationReport(BaseModel):
    """完整的推荐报告"""

    system_id: str = Field(..., description="系统标识")

    design_goals: List[str] = Field(..., description="设计目标列表")

    recommendations: List[MutationRecommendation] = Field(
        ...,
        description="突变推荐列表（8-12个）"
    )

    design_strategy: str = Field(..., description="整体设计策略说明")

    alternative_approaches: List[str] = Field(
        default_factory=list,
        description="替代方案"
    )

    experimental_roadmap: List[str] = Field(
        default_factory=list,
        description="实验路线图"
    )

    confidence: str = Field(..., description="整体置信度")

    metadata: Dict = Field(
        default_factory=dict,
        description="元数据（生成时间、模型版本等）"
    )

    class Config:
        json_schema_extra = {
            "example": {
                "system_id": "5c0a_run2",
                "design_goals": ["affinity_enhancement"],
                "recommendations": [
                    {
                        "residue": "ASP92",
                        "chain": "TCR_beta",
                        "position": 92,
                        "region": "CDR3",
                        "current_aa": "D",
                        "suggested_mutations": ["K", "R"],
                        "priority": "high",
                        "rationale": "高RRCS值(4.39)且与HLA K66形成强相互作用",
                        "expected_effects": {
                            "affinity": "增强2-3倍",
                            "specificity": "保持",
                            "stability": "轻微提升"
                        },
                        "risks": ["可能增加免疫原性"],
                        "validation_experiments": ["SPR测定", "细胞表达检测"],
                        "supporting_evidence": [
                            {
                                "type": "md_data",
                                "data": "RRCS=4.39, occupancy=0.90"
                            },
                            {
                                "type": "literature",
                                "pmid": "12345678",
                                "title": "CDR3 mutation enhances affinity"
                            }
                        ],
                        "confidence": "high"
                    }
                ],
                "design_strategy": "靶向CDR3β高RRCS残基，引入带电相互作用",
                "confidence": "high",
                "metadata": {
                    "generated_at": "2026-05-08T23:45:00",
                    "model": "deepseek-v4-pro",
                    "version": "1.0"
                }
            }
        }


class AnalysisDataSummary(BaseModel):
    """分析数据摘要（传递给LLM的格式）"""

    system_id: str

    metadata: Dict = Field(
        default_factory=dict,
        description="系统元数据（肽段、HLA、TCR基因型等）"
    )

    analysis_summary: Dict = Field(
        default_factory=dict,
        description="分析结果摘要（quality/interface/hotspots/flexibility/conformation）"
    )

    candidate_residues: List[Dict] = Field(
        default_factory=list,
        description="候选残基列表（按RRCS排序，top 10-20）"
    )

    def to_json_str(self, indent: int = 2) -> str:
        """转换为JSON字符串"""
        return self.model_dump_json(indent=indent)

    def estimate_tokens(self) -> int:
        """估算token数量（粗略）"""
        json_str = self.to_json_str()
        return len(json_str) // 4  # 粗略估算：4字符≈1token
