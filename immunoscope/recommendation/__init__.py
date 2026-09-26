"""
ImmunoScope Mutation Recommendation System

基于MD分析数据和文献知识的TCR突变推荐系统。

历史/当前状态 (2026-05-27):
    一次性 batch 推荐入口（`MutationRecommendationEngine` +
    `LLMAdapter`，定义在 ``engine.py``）已经被对话式的 Design Copilot
    取代——见 ``immunoscope/web/routers/design.py`` 顶部那行
    "Replaces the old batch recommendation engine"。CLI 子命令
    `ims recommend` 已于 2026-05-27 撤掉，公共 API 不再 re-export
    那两个类，避免新代码继续接 batch 路径。

仍然 export 的部件：
    - ``AnalysisDataFormatter`` / ``LiteratureRetriever``：data + knowledge
      layer，web agent 路径还在用。
    - ``MutationRecommendation`` / ``RecommendationReport`` /
      ``AnalysisDataSummary``：序列化 schema，web agent
      ``save_recommendation`` 工具也在用同一份。
    - ``TaskSpec``：评测层 (``evaluation.benchmark.task_routing``) 和
      web agent 任务分支都依赖。

如果确实需要 batch engine（脚本里），仍然可以 ``from
immunoscope.recommendation.engine import MutationRecommendationEngine``
显式导入；这只是不再放进默认命名空间。
"""

from .data_formatter import AnalysisDataFormatter
from .retriever import LiteratureRetriever
from .schemas import MutationRecommendation, RecommendationReport, AnalysisDataSummary
from .task_spec import TaskSpec

__all__ = [
    'AnalysisDataFormatter',
    'LiteratureRetriever',
    'MutationRecommendation',
    'RecommendationReport',
    'AnalysisDataSummary',
    'TaskSpec',
]
