"""Reporter analysis support modules.

The reporter package provides the query layer for structured data retrieval
from analysis outputs. It does NOT include LLM reasoning — that's handled by
the Agent's run_query() entry point.

Public API:
  - QueryRouter: Routes natural language queries to appropriate data handlers
  - QueryAnswerBuilder: Builds structured JSON answers from analysis files
  - QueryAnswer: Data model for query responses
"""

from .query_router import QueryRoute, QueryRouter
from .query_answer_builder import QueryAnswerBuilder, QueryAnswer

__all__ = [
    'QueryRoute',
    'QueryRouter',
    'QueryAnswerBuilder',
    'QueryAnswer',
]
