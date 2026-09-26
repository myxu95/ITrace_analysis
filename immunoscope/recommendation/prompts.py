"""
Prompt模板设计

为LLM提供System Prompt和User Prompt模板

The legacy `SYSTEM_PROMPT` is kept as the TCR-affinity wording (back-compat
for any external caller importing it directly). New code should call
``system_prompt_for_task(task)`` so the prompt matches the active TaskSpec.
"""
from typing import Dict, List, Optional
import json

from .task_spec import TaskSpec


# System Prompt - 定义AI角色、原则和输出格式（TCR-affinity 路径）
SYSTEM_PROMPT = """You are ImmunoScope Mutation Advisor, an expert system for TCR-pMHC mutation design.

Your role:
1. Analyze MD simulation data to identify key interaction sites
2. Combine literature knowledge to recommend rational mutation strategies
3. Assess potential risks and benefits of mutations
4. Provide clear design rationale and experimental validation suggestions

Design principles:
- Prioritize residues with high RRCS values and high observation frequency
- Consider impact on stability, expression, and immunogenicity
- Avoid disrupting TCR backbone structure
- Prefer conservative mutations (single point or few mutations)
- Base recommendations on both MD data AND literature evidence

Output format:
- Return ONLY valid JSON matching the RecommendationReport schema
- Number of recommendations should be driven by evidence quality, not a fixed target
  (typically 3-10; prefer fewer high-confidence picks over many low-confidence ones)
- Each recommendation must cite specific MD metrics and literature PMIDs
- Provide confidence levels (high/medium/low) based on evidence strength

Critical requirements:
- Recommended residues MUST be from the candidate list provided
- Suggested amino acids MUST be standard 20 amino acids (single letter code)
- All PMIDs must be from the provided literature list
- Explain WHY each mutation is recommended (mechanism-based rationale)

EXPLAINABILITY REQUIREMENTS (CRITICAL):
For each recommendation, you MUST provide detailed explainability fields:

1. selection_rationale: WHY this residue was selected over others
   - Include rank, percentile scores for key metrics
   - Quantitative comparison to next candidate (e.g., "RRCS 2.2x higher")
   - Core reason in 1-2 sentences

2. mutation_rationale: WHY these specific mutations were chosen
   - Current amino acid properties (type, polarity, charge, size)
   - Each suggested mutation's advantage and expected effect
   - At least 2 rejected options with reasons (e.g., "K: positive charge may repel HLA")
   - Overall design strategy

3. interaction_partners: Specific interaction details
   - List 3-5 key partners with distance, frequency, persistence
   - Interaction types (H-bond, hydrophobic, aromatic, salt bridge)

4. dynamics_evidence: MD dynamics support
   - RMSF value and interpretation
   - Contact persistence and interpretation
   - Conformational states

5. comparison_analysis: Comparative justification
   - Why other top candidates were NOT selected
   - Why other mutation options were NOT chosen
   - Summary of the comparison

These explainability fields are MANDATORY. Without them, the recommendation is incomplete.
"""


# System Prompt - 呈递任务路径（peptide × presentation × peptide_hla）
# Replaces the TCR-affinity framing in three places:
#   1. role description (TCR-pMHC → peptide-HLA)
#   2. design principles (avoid TCR backbone → preserve anchor pockets)
#   3. evidence vocabulary (CDR3 → anchor P2/PΩ)
# Kept structurally parallel to SYSTEM_PROMPT so the prompt parser and the
# downstream schema enforcer keep working unchanged.
SYSTEM_PROMPT_PRESENTATION = """You are ImmunoScope Presentation Advisor, an expert system for peptide-HLA presentation design.

Your role:
1. Analyze MD simulation data to identify peptide residues that drive HLA-groove anchoring
2. Combine literature knowledge to recommend peptide mutations that improve presentation stability
3. Assess potential risks (TCR-recognition loss, immunogenicity shift, anchor disruption)
4. Provide clear design rationale and experimental validation suggestions

Design principles:
- Prioritize peptide residues with high RRCS contribution to the HLA groove and persistent contacts
- Preserve anchor pocket chemistry (commonly P2 and PΩ for class-I MHC) unless explicitly relaxing it
- Avoid mutations that disrupt the TCR-facing surface unless the goal is to also retune TCR recognition
- Prefer conservative substitutions at anchor positions (hydrophobic→hydrophobic, etc.)
- Base recommendations on both MD data AND literature evidence

Output format:
- Return ONLY valid JSON matching the RecommendationReport schema
- Number of recommendations should be driven by evidence quality, not a fixed target
  (typically 3-10; prefer fewer high-confidence picks over many low-confidence ones)
- Each recommendation must cite specific MD metrics and literature PMIDs
- Provide confidence levels (high/medium/low) based on evidence strength

Critical requirements:
- Recommended residues MUST be from the candidate list provided (peptide-side)
- Suggested amino acids MUST be standard 20 amino acids (single letter code)
- All PMIDs must be from the provided literature list
- Explain WHY each mutation is recommended (mechanism-based rationale, anchor-aware)

EXPLAINABILITY REQUIREMENTS (CRITICAL):
For each recommendation, you MUST provide detailed explainability fields:

1. selection_rationale: WHY this peptide residue was selected over others
   - Include rank, percentile scores for key metrics
   - Quantitative comparison to next candidate (e.g., "RRCS 2.2x higher")
   - Whether the residue is at an anchor pocket (P2 / PΩ) and the implication
   - Core reason in 1-2 sentences

2. mutation_rationale: WHY these specific mutations were chosen
   - Current amino acid properties (type, polarity, charge, size)
   - Each suggested mutation's advantage and expected effect on groove fit
   - At least 2 rejected options with reasons (e.g., "G: too small for anchor pocket")
   - Overall design strategy

3. interaction_partners: Specific HLA-groove contacts
   - List 3-5 key partners with distance, frequency, persistence
   - Interaction types (H-bond, hydrophobic, aromatic, salt bridge)

4. dynamics_evidence: MD dynamics support
   - RMSF value and interpretation
   - Contact persistence and interpretation
   - Conformational states

5. comparison_analysis: Comparative justification
   - Why other top candidates were NOT selected
   - Why other mutation options were NOT chosen
   - Summary of the comparison

These explainability fields are MANDATORY. Without them, the recommendation is incomplete.
"""


def system_prompt_for_task(task: Optional[TaskSpec] = None) -> str:
    """Return the system prompt matching the active TaskSpec.

    Routing the system prompt on the task lets the LLM see consistent
    framing for the candidate filter and the wording. Falls back to the
    TCR-affinity prompt for the default / legacy path so existing call
    sites keep working without changes.
    """
    if task is None:
        task = TaskSpec.tcr_affinity()
    if task.subject == "peptide" and task.objective == "presentation":
        return SYSTEM_PROMPT_PRESENTATION
    return SYSTEM_PROMPT


def build_user_prompt(
    analysis_data: Dict,
    literature: List[Dict],
    design_goals: List[str],
    n_recommendations: int = 12,
    task: Optional[TaskSpec] = None,
) -> str:
    """
    构建User Prompt

    Args:
        analysis_data: 分析数据摘要（来自AnalysisDataFormatter）
        literature: 检索到的文献列表
        design_goals: 设计目标列表
        n_recommendations: 推荐数量
        task: TaskSpec controlling the task-framing line in the prompt
            header (TCR-pMHC affinity vs. peptide-HLA presentation). Defaults
            to the legacy TCR-affinity wording.

    Returns:
        格式化的User Prompt
    """
    if task is None:
        task = TaskSpec.tcr_affinity()

    # 提取关键信息
    metadata = analysis_data.get("metadata", {})
    analysis_summary = analysis_data.get("analysis_summary", {})
    candidate_residues = analysis_data.get("candidate_residues", [])

    # 格式化设计目标
    goals_text = ", ".join(design_goals)

    # 格式化系统信息
    system_info = f"""
Peptide: {metadata.get('peptide', 'N/A')}
HLA: {metadata.get('hla', 'N/A')} ({metadata.get('hla_allele', 'N/A')})
TCR Alpha V: {metadata.get('tcr_alpha_v', 'N/A')}
TCR Beta V: {metadata.get('tcr_beta_v', 'N/A')}
CDR3α: {metadata.get('cdr3_alpha', 'N/A')}
CDR3β: {metadata.get('cdr3_beta', 'N/A')}
"""

    # 格式化分析摘要
    quality = analysis_summary.get("quality", {})
    interface = analysis_summary.get("interface", {})
    conformation = analysis_summary.get("conformation", {})

    analysis_text = f"""
Quality:
  - Converged: {quality.get('is_converged', False)}
  - RMSD: {quality.get('rmsd_mean', 0):.3f} ± {quality.get('rmsd_std', 0):.3f} nm
  - Grade: {quality.get('overall_grade', 'N/A')}

Interface:
  - BSA: {interface.get('bsa_mean', 0):.1f} ± {interface.get('bsa_std', 0):.1f} Å²
  - Contact pairs: {interface.get('n_contact_pairs', 0)}

Conformation:
  - Clusters: {conformation.get('n_clusters', 0)}
  - Dominant population: {conformation.get('dominant_population', 0):.1%}
"""

    # 格式化候选残基（只显示top 10）
    # D-B5 (2026-05-26): composite priority score removed from LLM-facing
    # surface. Chemistry tags and risk flags are surfaced as categorical
    # evidence; ranking is by RRCS contribution alone.
    candidate_lines = []
    for i, res in enumerate(candidate_residues[:10]):
        design_bits = []
        if res.get("chemistry_tags"):
            design_bits.append("chemistry=" + ",".join(res["chemistry_tags"][:4]))
        if res.get("risk_flags"):
            design_bits.append("risks=" + ",".join(res["risk_flags"][:4]))
        design_text = "; " + "; ".join(design_bits) if design_bits else ""
        candidate_lines.append(
            f"{i+1}. {res['residue']:8s} ({res.get('region', ''):10s}, {res.get('chain', ''):5s}) - "
            f"RRCS={res.get('rrcs_mean', 0):5.2f}, occupancy={res.get('occupancy', 0):.2f}, "
            f"partner={res.get('partner', '')}{design_text}"
        )
    candidates_text = "\n".join(candidate_lines)

    # 格式化文献知识
    literature_text = "\n\n".join([
        f"[{i+1}] PMID: {paper['pmid']} ({paper['year']})\n"
        f"Title: {paper['title']}\n"
        f"Abstract: {paper['abstract'][:400]}...\n"
        f"Relevance: {paper['relevance_score']}"
        for i, paper in enumerate(literature[:5])  # 只显示top 5
    ])

    # 构建完整prompt — 第一行的任务描述按 TaskSpec 路由
    if task.subject == "peptide" and task.objective == "presentation":
        task_header_line = (
            f"Design peptide mutation recommendations for this peptide-HLA "
            f"system to achieve: {goals_text}"
        )
    else:
        task_header_line = (
            f"Design mutation recommendations for this TCR-pMHC system to "
            f"achieve: {goals_text}"
        )

    # 呈递任务专属：把 anchor pocket chemistry 渲染到 prompt 顶部，
    # 让 LLM 在挑候选之前先看到哪些位置是 anchor。
    # 这块对 TCR-affinity 任务不渲染，所以原有 prompt 完全向后兼容。
    anchor_pocket_section = ""
    if task.subject == "peptide" and task.objective == "presentation":
        anchor_entries = metadata.get("anchor_pockets") or []
        if anchor_entries:
            anchor_pocket_section = "\n\n# Anchor Pocket Chemistry (allele-specific motif knowledge)\n"
            for entry in anchor_entries:
                pref = ",".join(entry.get("preferred_residues") or []) or "—"
                tol = ",".join(entry.get("tolerated_residues") or []) or "—"
                disp = ",".join(entry.get("dispreferred_residues") or []) or "—"
                current = entry.get("current_residue") or "?"
                anchor_pocket_section += (
                    f"- P{entry['position']} ({entry['role']}): "
                    f"current={current}, "
                    f"preferred={{{pref}}}, tolerated={{{tol}}}, "
                    f"dispreferred={{{disp}}}\n"
                    f"  notes: {entry['notes']}\n"
                    f"  source: {entry['source']}\n"
                )
        else:
            # Tell the LLM explicitly that anchor info is unavailable —
            # silence here would be indistinguishable from "no anchors exist"
            anchor_pocket_section = (
                "\n\n# Anchor Pocket Chemistry\n"
                "- No allele-specific anchor data available for this system. "
                "Treat P2 and the C-terminal residue as candidate anchors "
                "regardless; mutations there carry presentation-loss risk.\n"
            )

    prompt = f"""# Task
{task_header_line}{anchor_pocket_section}

# System Information
{system_info.strip()}

# MD Analysis Summary
{analysis_text.strip()}

# Candidate Residues (ranked by RRCS contribution; chemistry + risk flags are categorical annotations, not weights)
{candidates_text}

# Relevant Literature Knowledge
{literature_text}

# Requirements
1. Recommend mutation sites with strong supporting evidence (up to {n_recommendations} sites maximum)
   - Quality over quantity: only include sites where MD data AND literature provide meaningful evidence
   - If evidence is limited or weak, return fewer recommendations (3-4 strong picks > 8 weak ones)
   - Do NOT pad with low-confidence recommendations to reach a target number
2. For each site, suggest 5-6 candidate amino acids consistent with the design direction
3. Explain design rationale based on MD data AND literature
4. Assess potential risks (stability, immunogenicity, expression)
5. Suggest validation experiments (SPR, cell assays, etc.)
6. Provide confidence level (high/medium/low) based on evidence strength

# Output Format
Return ONLY a valid JSON object matching this schema:

{{
  "system_id": "string",
  "design_goals": ["string"],
  "recommendations": [
    {{
      "residue": "string (e.g., ASP92)",
      "chain": "string (alpha/beta)",
      "position": integer,
      "region": "string (CDR1/CDR2/CDR3/Framework)",
      "current_aa": "string (single letter)",
      "suggested_mutations": ["string (single letters)"],
      "priority": "string (high/medium/low)",
      "rationale": "string (explain WHY based on MD + literature)",
      "expected_effects": {{
        "affinity": "string",
        "specificity": "string",
        "stability": "string"
      }},
      "risks": ["string"],
      "validation_experiments": ["string"],
      "supporting_evidence": [
        {{
          "type": "md_data",
          "data": "string (cite specific RRCS, occupancy, etc.)"
        }},
        {{
          "type": "literature",
          "pmid": "string",
          "title": "string",
          "relevance": "string"
        }}
      ],
      "confidence": "string (high/medium/low)",

      "selection_rationale": {{
        "rank": integer,
        "total_candidates": integer,
        "key_metrics": {{
          "rrcs": {{"value": float, "percentile": float, "interpretation": "string"}},
          "occupancy": {{"value": float, "percentile": float, "interpretation": "string"}}
        }},
        "comparison_to_next": "string (e.g., RRCS 2.2x higher than rank 2)",
        "why_this_residue": "string (1-2 sentences)"
      }},

      "mutation_rationale": {{
        "current_aa_properties": {{
          "code": "string",
          "name": "string",
          "type": "string (aromatic/charged/polar/hydrophobic)",
          "polarity": "string (polar/nonpolar)",
          "charge": "string (positive/negative/neutral)",
          "size": "string (small/medium/large)"
        }},
        "suggested_mutations": {{
          "Y": {{
            "aa_code": "Y",
            "properties": {{"code": "Y", "name": "Tyrosine", "type": "aromatic", "polarity": "polar", "charge": "neutral", "size": "large"}},
            "advantage": "string (e.g., pi-pi stacking + H-bond)",
            "expected_effect": "string",
            "confidence": "string"
          }}
        }},
        "rejected_options": {{
          "K": "string (reason why rejected)",
          "A": "string (reason why rejected)"
        }},
        "design_strategy": "string"
      }},

      "interaction_partners": [
        {{
          "partner_residue": "string (e.g., HLA-A*02:01 TYR99)",
          "partner_chain": "string (HLA_alpha/HLA_beta/peptide)",
          "interaction_type": "string (H-bond/hydrophobic/aromatic/salt-bridge)",
          "distance": float,
          "distance_std": float,
          "frequency": float,
          "persistence": float,
          "strength": "string (strong/moderate/weak)"
        }}
      ],

      "dynamics_evidence": {{
        "rmsf": float,
        "rmsf_interpretation": "string (high/moderate/low flexibility)",
        "contact_persistence": float,
        "persistence_interpretation": "string",
        "conformational_states": integer,
        "dominant_state_population": float
      }},

      "comparison_analysis": {{
        "vs_other_candidates": [
          {{
            "candidate_residue": "string",
            "rank": integer,
            "rrcs": float,
            "occupancy": float,
            "why_not_selected": "string"
          }}
        ],
        "vs_other_mutations": [
          {{
            "mutation": "string (e.g., E→K)",
            "score": float,
            "pros": ["string"],
            "cons": ["string"],
            "selected": boolean
          }}
        ],
        "summary": "string"
      }}
    }}
  ],
  "design_strategy": "string (overall strategy explanation)",
  "alternative_approaches": ["string"],
  "experimental_roadmap": ["string"],
  "confidence": "string (high/medium/low)",
  "metadata": {{
    "model": "string",
    "version": "string"
  }}
}}

IMPORTANT:
- Use ONLY residues from the candidate list above
- Use ONLY standard amino acids: A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
- Cite ONLY PMIDs from the literature list above
- Provide mechanism-based rationale, not just correlation
"""

    return prompt


def extract_json_from_response(response: str) -> Dict:
    """
    从LLM响应中提取JSON

    处理可能的格式问题：
    - 响应可能包含markdown代码块
    - 响应可能包含额外的文本
    - 响应可能被截断
    """
    import re

    # 尝试直接解析
    try:
        return json.loads(response)
    except json.JSONDecodeError:
        pass

    # 尝试提取markdown代码块中的JSON
    if "```json" in response:
        start = response.find("```json") + 7
        end = response.find("```", start)
        if end == -1:  # 没有找到结束标记，可能被截断
            # 尝试找到最后一个完整的 }
            json_str = response[start:].strip()
            # 找到最后一个 } 的位置
            last_brace = json_str.rfind("}")
            if last_brace != -1:
                json_str = json_str[:last_brace + 1]
        else:
            json_str = response[start:end].strip()

        try:
            return json.loads(json_str)
        except json.JSONDecodeError as e:
            # 尝试修复常见的 JSON 错误（尾部逗号）
            json_str = re.sub(r',(\s*[}\]])', r'\1', json_str)
            try:
                return json.loads(json_str)
            except json.JSONDecodeError:
                pass

    # 尝试提取{}之间的内容
    if "{" in response and "}" in response:
        start = response.find("{")
        end = response.rfind("}") + 1
        json_str = response[start:end]
        try:
            return json.loads(json_str)
        except json.JSONDecodeError:
            # 尝试修复
            json_str = re.sub(r',(\s*[}\]])', r'\1', json_str)
            try:
                return json.loads(json_str)
            except json.JSONDecodeError:
                pass

    raise ValueError(f"无法从响应中提取有效的JSON: {response[:200]}...")


# 使用示例
if __name__ == "__main__":
    # 示例数据
    example_analysis = {
        "system_id": "5c0a_run2",
        "metadata": {
            "peptide": "MVWGPDPLYV",
            "hla": "HLA-A",
            "hla_allele": "HLA-A*02:665"
        },
        "analysis_summary": {
            "quality": {"is_converged": True, "rmsd_mean": 0.676, "rmsd_std": 0.168},
            "interface": {"bsa_mean": 1500.0, "bsa_std": 50.0}
        },
        "candidate_residues": [
            {
                "residue": "ASP92",
                "chain": "alpha",
                "region": "CDR3",
                "rrcs_mean": 4.39,
                "occupancy": 0.90,
                "partner": "LYS66@HLA_alpha"
            }
        ]
    }

    example_literature = [
        {
            "pmid": "12345678",
            "title": "TCR engineering for enhanced affinity",
            "abstract": "We demonstrate that mutations in CDR3...",
            "year": 2023,
            "relevance_score": 0.85
        }
    ]

    prompt = build_user_prompt(
        analysis_data=example_analysis,
        literature=example_literature,
        design_goals=["affinity_enhancement"]
    )

    print("=" * 60)
    print("System Prompt:")
    print("=" * 60)
    print(SYSTEM_PROMPT)
    print("\n" + "=" * 60)
    print("User Prompt (example):")
    print("=" * 60)
    print(prompt[:1000] + "...")
