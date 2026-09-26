"""
突变推荐引擎 (DEPRECATED, 2026-05-27)

历史定位：把 Data Layer (`AnalysisDataFormatter`)、Knowledge Layer
(`LiteratureRetriever`)、Reasoning Layer (LLM) 拼成一次性 batch 推荐
流水线，给原先的 `ims recommend` CLI 子命令用。

为什么不再用：
    - 推荐任务已经全部迁到对话式 Design Copilot
      (`immunoscope/web/routers/design.py` +
      `immunoscope/agent/design_context.py`)。web agent 用工具调用 +
      持久化 session 把"挑候选 / 调权 / 写报告"这几步拆开，
      不像 batch engine 一次跑完就走。
    - CLI 子命令 `ims recommend` 已于 2026-05-27 从
      `immunoscope/cli/main.py` 撤掉。
    - 任务分支（TCR affinity vs peptide presentation）走 web agent 的
      `task` 字段；这里的 `MutationRecommendationEngine.recommend()`
      没有任务路由逻辑，继续用会把所有任务塞回 TCR-affinity 模板。

公共 API 状态：
    本模块**不再**通过 `immunoscope.recommendation` 顶层包 re-export。
    顶层 `__init__.py` 现在只 export data formatter / retriever / schema /
    TaskSpec。如果新代码需要走推荐，请走 web agent 路径。

如果你 *真的* 需要在脚本里复用 batch 流程，仍然可以
``from immunoscope.recommendation.engine import MutationRecommendationEngine``
显式导入；但这条路径不会再被维护，下一轮清理会整个删掉文件。
"""
import os
import json
import logging
from typing import List, Dict, Optional
from pathlib import Path
from datetime import datetime

from .data_formatter import AnalysisDataFormatter
from .retriever import LiteratureRetriever
from .schemas import RecommendationReport, MutationRecommendation
from .prompts import (
    SYSTEM_PROMPT,
    build_user_prompt,
    extract_json_from_response,
    system_prompt_for_task,
)
from .explainability_report import ExplainabilityReportGenerator
from .task_spec import TaskSpec
from immunoscope.llm_providers import (
    default_model as _provider_default_model,
    resolve_api_key as _provider_resolve_api_key,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class LLMAdapter:
    """LLM适配器 - 支持DeepSeek, Claude, OpenAI"""

    def __init__(
        self,
        provider: str = "deepseek",
        api_key: Optional[str] = None,
        model: Optional[str] = None
    ):
        """
        初始化LLM适配器

        Args:
            provider: LLM提供商 (deepseek/anthropic/openai)
            api_key: API密钥（如果为None，从环境变量读取）
            model: 模型名称（如果为None，使用默认）
        """
        self.provider = provider.lower()

        # 委托给共享的 provider 元数据表（immunoscope/llm_providers.py）：
        # 1) provider 校验：未知 provider 立即 ValueError，避免到 generate() 才报错；
        # 2) API key 解析：优先入参，回退到该 provider 对应的环境变量；
        # 3) 默认模型：调用方未指定时回退到表里的默认值。
        self.api_key = _provider_resolve_api_key(self.provider, api_key)
        self.model = model if model is not None else _provider_default_model(self.provider)

        logger.info(f"初始化LLM适配器: provider={self.provider}, model={self.model}")

    def generate(
        self,
        system_prompt: str,
        user_prompt: str,
        temperature: float = 0.3,
        max_tokens: int = 4000
    ) -> str:
        """
        调用LLM生成响应

        Args:
            system_prompt: 系统提示
            user_prompt: 用户提示
            temperature: 温度参数
            max_tokens: 最大token数

        Returns:
            LLM响应文本
        """
        if self.provider == "deepseek":
            return self._generate_deepseek(system_prompt, user_prompt, temperature, max_tokens)
        elif self.provider == "anthropic":
            return self._generate_anthropic(system_prompt, user_prompt, temperature, max_tokens)
        elif self.provider == "openai":
            return self._generate_openai(system_prompt, user_prompt, temperature, max_tokens)
        else:
            raise ValueError(f"Unsupported provider: {self.provider}")

    def _generate_deepseek(
        self,
        system_prompt: str,
        user_prompt: str,
        temperature: float,
        max_tokens: int
    ) -> str:
        """DeepSeek API调用（使用OpenAI兼容接口）"""
        import requests

        url = "https://api.deepseek.com/v1/chat/completions"
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json"
        }
        data = {
            "model": self.model,
            "messages": [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt}
            ],
            "temperature": temperature,
            "max_tokens": max_tokens
        }

        logger.info(f"调用DeepSeek API: model={self.model}, temp={temperature}")

        try:
            response = requests.post(url, headers=headers, json=data, timeout=120)
            response.raise_for_status()
        except requests.Timeout:
            logger.exception("DeepSeek 请求超时（120s）")
            raise
        except requests.HTTPError as exc:
            status = exc.response.status_code if exc.response is not None else "?"
            body = exc.response.text[:300] if exc.response is not None else ""
            logger.error("DeepSeek API HTTP %s: %s", status, body)
            raise
        except requests.RequestException:
            logger.exception("DeepSeek 网络错误")
            raise

        try:
            result = response.json()["choices"][0]["message"]["content"]
        except (ValueError, KeyError, IndexError, TypeError):
            logger.exception("DeepSeek 响应结构异常: %s", response.text[:300])
            raise
        logger.info(f"DeepSeek响应长度: {len(result)} 字符")
        return result

    def _generate_anthropic(
        self,
        system_prompt: str,
        user_prompt: str,
        temperature: float,
        max_tokens: int
    ) -> str:
        """Anthropic Claude API调用"""
        try:
            from anthropic import Anthropic, APIError
        except ImportError as exc:
            raise ImportError("请安装anthropic包: pip install anthropic") from exc

        client = Anthropic(api_key=self.api_key)
        logger.info(f"调用Anthropic API: model={self.model}, temp={temperature}")

        try:
            response = client.messages.create(
                model=self.model,
                max_tokens=max_tokens,
                temperature=temperature,
                system=system_prompt,
                messages=[{"role": "user", "content": user_prompt}],
            )
        except APIError:
            logger.exception("Anthropic API 错误")
            raise

        try:
            result = response.content[0].text
        except (IndexError, AttributeError):
            logger.exception("Anthropic 响应结构异常: %r", response)
            raise
        logger.info(f"Anthropic响应长度: {len(result)} 字符")
        return result

    def _generate_openai(
        self,
        system_prompt: str,
        user_prompt: str,
        temperature: float,
        max_tokens: int
    ) -> str:
        """OpenAI API调用"""
        try:
            from openai import OpenAI, OpenAIError
        except ImportError as exc:
            raise ImportError("请安装openai包: pip install openai") from exc

        client = OpenAI(api_key=self.api_key)
        logger.info(f"调用OpenAI API: model={self.model}, temp={temperature}")

        try:
            response = client.chat.completions.create(
                model=self.model,
                max_tokens=max_tokens,
                temperature=temperature,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt},
                ],
            )
        except OpenAIError:
            logger.exception("OpenAI API 错误")
            raise

        try:
            result = response.choices[0].message.content
        except (IndexError, AttributeError):
            logger.exception("OpenAI 响应结构异常: %r", response)
            raise
        logger.info(f"OpenAI响应长度: {len(result)} 字符")
        return result


class MutationRecommendationEngine:
    """突变推荐引擎"""

    def __init__(
        self,
        literature_db_path: str,
        llm_provider: str = "deepseek",
        llm_api_key: Optional[str] = None,
        llm_model: Optional[str] = None
    ):
        """
        初始化推荐引擎

        Args:
            literature_db_path: 文献向量数据库路径
            llm_provider: LLM提供商
            llm_api_key: LLM API密钥
            llm_model: LLM模型名称
        """
        self.retriever = LiteratureRetriever(literature_db_path)
        self.llm = LLMAdapter(
            provider=llm_provider,
            api_key=llm_api_key,
            model=llm_model
        )

        logger.info("突变推荐引擎初始化完成")

    def recommend(
        self,
        analysis_dir: str,
        design_goals: List[str],
        n_candidates: int = 15,
        n_papers: int = 10,
        n_recommendations: int = 12,
        output_file: Optional[str] = None,
        task: Optional[TaskSpec] = None,
    ) -> RecommendationReport:
        """
        生成突变推荐

        Args:
            analysis_dir: 分析结果目录
            design_goals: 设计目标列表
            n_candidates: 候选残基数量
            n_papers: 检索文献数量
            output_file: 输出文件路径（可选）
            task: TaskSpec declaring (subject, objective, target_subinterface).
                Defaults to TCR-affinity for backwards compatibility. Threads
                through the formatter (candidate filter) and prompt builder
                (system + user prompt wording) so the LLM sees the right task
                framing for non-default combinations.

        Returns:
            RecommendationReport对象
        """
        if task is None:
            task = TaskSpec.tcr_affinity()

        logger.info("=" * 60)
        logger.info("开始生成突变推荐")
        logger.info(f"任务: {task.slug}")
        logger.info("=" * 60)

        # 1. Data Layer - 格式化分析数据
        logger.info("[1/5] 格式化分析数据...")
        formatter = AnalysisDataFormatter(analysis_dir)
        analysis_summary = formatter.format(max_candidates=n_candidates, task=task)
        logger.info(f"✓ 提取了 {len(analysis_summary.candidate_residues)} 个候选残基")

        # 2. Knowledge Layer - 检索相关文献
        logger.info("[2/5] 检索相关文献...")
        papers = self.retriever.retrieve_for_mutation_design(
            analysis_data=analysis_summary.model_dump(),
            design_goal=design_goals[0] if design_goals else "affinity_enhancement",
            n_results=n_papers
        )
        logger.info(f"✓ 检索到 {len(papers)} 篇相关文献")

        # 3. Reasoning Layer - 构建Prompt
        logger.info("[3/5] 构建LLM Prompt...")
        user_prompt = build_user_prompt(
            analysis_data=analysis_summary.model_dump(),
            literature=papers,
            design_goals=design_goals,
            n_recommendations=n_recommendations,
            task=task,
        )
        logger.info(f"✓ Prompt长度: {len(user_prompt)} 字符")

        # 4. 调用LLM
        logger.info("[4/5] 调用LLM生成推荐...")
        response = self.llm.generate(
            system_prompt=system_prompt_for_task(task),
            user_prompt=user_prompt,
            temperature=0.3,
            max_tokens=8000  # 增加到 8000 以容纳可解释性字段
        )
        logger.info(f"✓ LLM响应长度: {len(response)} 字符")

        # 5. 解析和验证输出
        logger.info("[5/5] 解析和验证输出...")
        from pydantic import ValidationError

        try:
            result_dict = extract_json_from_response(response)
        except (ValueError, json.JSONDecodeError):
            logger.exception("LLM 输出 JSON 解析失败；原始响应: %s...", response[:500])
            raise

        # 添加元数据
        if "metadata" not in result_dict:
            result_dict["metadata"] = {}
        result_dict["metadata"].update({
            "generated_at": datetime.now().isoformat(),
            "model": self.llm.model,
            "provider": self.llm.provider,
            "version": "1.0",
        })

        try:
            report = RecommendationReport(**result_dict)
        except ValidationError:
            logger.exception("LLM 输出 schema 校验失败；原始响应: %s...", response[:500])
            raise

        # 验证推荐的残基
        candidate_residues = {r["residue"] for r in analysis_summary.candidate_residues}
        for rec in report.recommendations:
            if rec.residue not in candidate_residues:
                logger.warning(f"推荐的残基 {rec.residue} 不在候选列表中")
        logger.info(f"✓ 生成了 {len(report.recommendations)} 个突变推荐")

        # 保存输出
        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # 保存 JSON 报告
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(report.model_dump_json(indent=2, exclude_none=True))

            logger.info(f"✓ JSON报告已保存到: {output_file}")

            # 生成可解释性 Markdown 报告（次要输出：失败时跳过，但保留 JSON 报告）
            md_output = output_path.with_suffix('.explainability.md')
            try:
                explainer = ExplainabilityReportGenerator()
                explainer.generate_markdown_report(report, str(md_output))
                logger.info(f"✓ 可解释性报告已保存到: {md_output}")
            except Exception:
                # 广义捕获是有意的：可解释性报告是次要产物，不应阻塞主流程；
                # 用 logger.exception 保留 traceback 便于排查模板/字段缺失等问题。
                logger.exception("生成可解释性报告失败（次要输出，已跳过）")

        logger.info("=" * 60)
        logger.info("突变推荐生成完成！")
        logger.info("=" * 60)

        return report


# 使用示例
if __name__ == "__main__":
    # 检查API密钥
    if not os.environ.get("DEEPSEEK_API_KEY"):
        print("错误: 请设置环境变量 DEEPSEEK_API_KEY")
        print("export DEEPSEEK_API_KEY='your-api-key'")
        exit(1)

    # 初始化引擎
    engine = MutationRecommendationEngine(
        literature_db_path="development/literature_collection/knowledge_base/vector_db",
        llm_provider="deepseek"
    )

    # 生成推荐
    report = engine.recommend(
        analysis_dir="output/5c0a_run2_full_analysis",
        design_goals=["affinity_enhancement"],
        n_candidates=15,
        n_papers=5,
        output_file="output/5c0a_run2_full_analysis/mutation_recommendations.json"
    )

    # 打印结果
    print("\n" + "=" * 60)
    print("突变推荐结果:")
    print("=" * 60)
    print(f"系统ID: {report.system_id}")
    print(f"设计目标: {', '.join(report.design_goals)}")
    print(f"推荐数量: {len(report.recommendations)}")
    print(f"整体置信度: {report.confidence}")
    print("\n推荐列表:")
    for i, rec in enumerate(report.recommendations, 1):
        print(f"\n{i}. {rec.residue} ({rec.region}, {rec.chain})")
        print(f"   当前: {rec.current_aa} → 建议: {', '.join(rec.suggested_mutations)}")
        print(f"   优先级: {rec.priority}, 置信度: {rec.confidence}")
        print(f"   理由: {rec.rationale[:100]}...")
