"""
可解释性报告生成器

生成人类可读的 Markdown 格式可解释性报告
"""
from typing import List, Dict
from pathlib import Path
from .schemas import (
    MutationRecommendation,
    RecommendationReport,
    SelectionRationale,
    MutationRationale,
    InteractionPartner,
    DynamicsEvidence,
    ComparisonAnalysis
)


class ExplainabilityReportGenerator:
    """生成人类可读的可解释性报告"""

    def generate_markdown_report(
        self,
        report: RecommendationReport,
        output_file: str
    ) -> None:
        """
        生成 Markdown 格式的可解释性报告

        Args:
            report: 推荐报告
            output_file: 输出文件路径
        """
        md_content = self._generate_header(report)

        for i, rec in enumerate(report.recommendations, 1):
            md_content += self._generate_recommendation_section(rec, i)

        md_content += self._generate_summary(report)

        # 保存文件
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(md_content)

    def _generate_header(self, report: RecommendationReport) -> str:
        """生成报告头部"""
        md = f"# 突变推荐可解释性报告\n\n"
        md += f"**系统**: {report.system_id}\n\n"
        md += f"**设计目标**: {', '.join(report.design_goals)}\n\n"
        md += f"**推荐数量**: {len(report.recommendations)}\n\n"
        md += f"**整体置信度**: {report.confidence}\n\n"
        md += f"---\n\n"

        return md

    def _generate_recommendation_section(
        self,
        rec: MutationRecommendation,
        index: int
    ) -> str:
        """生成单个推荐的详细解释"""

        md = f"\n## 推荐 {index}: {rec.residue} ({rec.region})\n\n"

        # 基本信息
        md += f"**位置**: {rec.chain} 链第 {rec.position} 位\n\n"
        md += f"**当前氨基酸**: {rec.current_aa}\n\n"
        md += f"**建议突变**: {' / '.join(rec.suggested_mutations)}\n\n"
        md += f"**优先级**: {rec.priority} | **置信度**: {rec.confidence}\n\n"
        md += f"---\n\n"

        # 1. 选择理由
        if rec.selection_rationale:
            md += "### 🎯 为什么选择这个残基？\n\n"
            md += self._format_selection_rationale(rec.selection_rationale)
            md += "\n"

        # 2. 突变理由
        if rec.mutation_rationale:
            md += "### 🧬 为什么选择这些突变？\n\n"
            md += self._format_mutation_rationale(rec.mutation_rationale)
            md += "\n"

        # 3. 相互作用伙伴
        if rec.interaction_partners:
            md += "### 🤝 相互作用伙伴\n\n"
            md += self._format_interaction_partners(rec.interaction_partners)
            md += "\n"

        # 4. 动力学证据
        if rec.dynamics_evidence:
            md += "### 📊 动力学证据\n\n"
            md += self._format_dynamics_evidence(rec.dynamics_evidence)
            md += "\n"

        # 5. 对比分析
        if rec.comparison_analysis:
            md += "### ⚖️ 对比分析\n\n"
            md += self._format_comparison_analysis(rec.comparison_analysis)
            md += "\n"

        # 6. 预期效果和风险
        md += "### 💡 预期效果\n\n"
        md += self._format_expected_effects(rec.expected_effects)
        md += "\n"

        md += "### ⚠️ 风险评估\n\n"
        md += self._format_risks(rec.risks)
        md += "\n"

        # 7. 验证实验
        if rec.validation_experiments:
            md += "### 🔬 建议的验证实验\n\n"
            for exp in rec.validation_experiments:
                md += f"- {exp}\n"
            md += "\n"

        md += "---\n\n"

        return md

    def _format_selection_rationale(self, rationale: SelectionRationale) -> str:
        """格式化选择理由"""
        md = f"**排名**: {rationale.rank} / {rationale.total_candidates}\n\n"

        md += "**关键指标**:\n\n"
        md += "| 指标 | 数值 | 百分位 | 解释 |\n"
        md += "|------|------|--------|------|\n"

        for metric, score in rationale.key_metrics.items():
            md += f"| {metric.upper()} | {score.value:.2f} | "
            md += f"{score.percentile:.0f}% | {score.interpretation} |\n"

        md += f"\n**对比**: {rationale.comparison_to_next}\n\n"
        md += f"**核心原因**: {rationale.why_this_residue}\n"

        return md

    def _format_mutation_rationale(self, rationale: MutationRationale) -> str:
        """格式化突变理由"""
        md = f"**当前氨基酸**: {rationale.current_aa_properties.name} "
        md += f"({rationale.current_aa_properties.code})\n\n"

        md += "| 性质 | 值 |\n"
        md += "|------|----|\n"
        md += f"| 类型 | {rationale.current_aa_properties.type} |\n"
        md += f"| 极性 | {rationale.current_aa_properties.polarity} |\n"
        md += f"| 电荷 | {rationale.current_aa_properties.charge} |\n"
        md += f"| 大小 | {rationale.current_aa_properties.size} |\n\n"

        md += "**建议突变**:\n\n"
        for aa, option in rationale.suggested_mutations.items():
            md += f"**{aa} ({option.properties.name})**:\n"
            md += f"- **优势**: {option.advantage}\n"
            md += f"- **预期效果**: {option.expected_effect}\n"
            md += f"- **置信度**: {option.confidence}\n\n"

        if rationale.rejected_options:
            md += "**被拒绝的选项**:\n\n"
            for aa, reason in rationale.rejected_options.items():
                md += f"- **{aa}**: {reason}\n"
            md += "\n"

        md += f"**设计策略**: {rationale.design_strategy}\n"

        return md

    def _format_interaction_partners(self, partners: List[InteractionPartner]) -> str:
        """格式化相互作用伙伴"""
        if not partners:
            return "无详细数据\n"

        md = "| 伙伴残基 | 链 | 类型 | 距离 (Å) | 频率 | 持续性 | 强度 |\n"
        md += "|---------|-----|------|---------|------|--------|------|\n"

        for p in partners:
            md += f"| {p.partner_residue} | {p.partner_chain} | {p.interaction_type} | "
            md += f"{p.distance:.1f} ± {p.distance_std:.1f} | "
            md += f"{p.frequency:.0%} | {p.persistence:.2f} | {p.strength} |\n"

        return md

    def _format_dynamics_evidence(self, evidence: DynamicsEvidence) -> str:
        """格式化动力学证据"""
        md = f"- **RMSF**: {evidence.rmsf:.2f} Å ({evidence.rmsf_interpretation})\n"
        md += f"- **接触持续性**: {evidence.contact_persistence:.2f} "
        md += f"({evidence.persistence_interpretation})\n"
        md += f"- **构象状态**: {evidence.conformational_states} 个\n"
        md += f"- **主导构象占比**: {evidence.dominant_state_population:.0%}\n"

        return md

    def _format_comparison_analysis(self, analysis: ComparisonAnalysis) -> str:
        """格式化对比分析"""
        md = "**与其他候选残基对比**:\n\n"

        if analysis.vs_other_candidates:
            md += "| 残基 | 排名 | RRCS | 接触持续度 | 为什么没被选择 |\n"
            md += "|------|------|------|-----------|---------------|\n"

            for comp in analysis.vs_other_candidates[:3]:
                md += f"| {comp.candidate_residue} | {comp.rank} | "
                md += f"{comp.rrcs:.2f} | {comp.occupancy:.2f} | "
                md += f"{comp.why_not_selected} |\n"
        else:
            md += "无对比数据\n"

        md += "\n**与其他突变选项对比**:\n\n"

        if analysis.vs_other_mutations:
            for comp in analysis.vs_other_mutations:
                status = "✅ 已选择" if comp.selected else "❌ 未选择"
                md += f"**{comp.mutation}** ({status}, 评分: {comp.score:.2f})\n"
                md += f"  - 优点: {', '.join(comp.pros)}\n"
                md += f"  - 缺点: {', '.join(comp.cons)}\n\n"
        else:
            md += "无对比数据\n"

        md += f"\n**总结**: {analysis.summary}\n"

        return md

    def _format_expected_effects(self, effects: Dict[str, str]) -> str:
        """格式化预期效果"""
        if not effects:
            return "无详细数据\n"

        md = ""
        for key, value in effects.items():
            md += f"- **{key.capitalize()}**: {value}\n"

        return md

    def _format_risks(self, risks: List[str]) -> str:
        """格式化风险"""
        if not risks:
            return "无明显风险\n"

        md = ""
        for risk in risks:
            md += f"- {risk}\n"

        return md

    def _generate_summary(self, report: RecommendationReport) -> str:
        """生成报告总结"""
        md = "\n## 📋 整体设计策略\n\n"
        md += f"{report.design_strategy}\n\n"

        if report.alternative_approaches:
            md += "### 替代方案\n\n"
            for i, approach in enumerate(report.alternative_approaches, 1):
                md += f"{i}. {approach}\n"
            md += "\n"

        if report.experimental_roadmap:
            md += "### 实验路线图\n\n"
            for i, step in enumerate(report.experimental_roadmap, 1):
                md += f"{i}. {step}\n"
            md += "\n"

        md += "---\n\n"
        md += f"*报告生成时间: {report.metadata.get('generated_at', 'N/A')}*\n"
        model_label = report.metadata.get('model', 'N/A')
        provider_label = report.metadata.get('provider')
        if provider_label and model_label != 'N/A':
            md += f"*LLM 模型: {provider_label}/{model_label}*\n"
        else:
            md += f"*LLM 模型: {model_label}*\n"

        return md
