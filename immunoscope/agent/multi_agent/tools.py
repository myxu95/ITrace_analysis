"""Multi-agent-only tools (v1, 2026-05-28).

These tools are exposed to specific multi-agent stages via
`BaseAgent(extra_tools=...)` — they are intentionally NOT registered in
the global `agent.tools.TOOL_REGISTRY` so they don't leak into the
generic conversational agent surface.

Currently exposed:

- `SearchLiteratureTool`  — Bio Agent's only external skill in v1.
                            Wraps `recommendation.retriever.LiteratureRetriever`
                            with a simple free-text query interface.

Future tools (post-v1) can be added here without touching the rest of
the pipeline; the orchestrator stays oblivious to which agent owns
which tool.
"""

from __future__ import annotations

import asyncio
import logging
import os
import threading
from pathlib import Path
from typing import Any, List, Optional

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult

log = logging.getLogger("immunoscope.agent.multi_agent.tools")

# Constructing a ChromaDB client is NOT safe to run concurrently from multiple
# threads on the same persist dir: the Rust bindings race during tenant
# validation ("'RustBindingsAPI' object has no attribute 'bindings'"). Under
# batch recommendation (multiple cases in asyncio.gather, each building its own
# SearchLiteratureTool / SearchFactBlocksTool), several `asyncio.to_thread`
# workers hit retriever construction at once. We serialize construction and
# share one retriever per resolved persist path across all tool instances —
# reads on the constructed client are then concurrency-safe.
_RETRIEVER_LOCK = threading.Lock()
_RETRIEVER_CACHE: dict = {}


_DEFAULT_LIT_DB_PATH = os.environ.get(
    "IMMUNOSCOPE_LITERATURE_DB",
    "development/literature_collection/knowledge_base/vector_db",
)

# The extracted-fact-block collection lives in the SAME ChromaDB persist
# directory as the abstract collection (`tcr_literature`), under a different
# collection name (`tcr_facts`). Built by
# development/literature_collection/build_fact_blocks.py.
_DEFAULT_FACTS_DB_PATH = os.environ.get(
    "IMMUNOSCOPE_FACTS_DB", _DEFAULT_LIT_DB_PATH
)
_FACTS_COLLECTION = "tcr_facts"
_FACTS_EMBED_MODEL = "all-MiniLM-L6-v2"
_FACT_BLOCK_TYPES = (
    "mutation_effect",
    "binding_measurement",
    "structural_observation",
    "design_heuristic",
)


def _resolve_repo_path(rel_or_abs: str) -> str:
    """Resolve a repo-relative path by walking up from this file (mirrors the
    convention used for the literature DB)."""
    if Path(rel_or_abs).is_absolute():
        return rel_or_abs
    here = Path(__file__).resolve()
    for ancestor in here.parents:
        candidate = ancestor / rel_or_abs
        if candidate.exists():
            return str(candidate)
    return rel_or_abs


class _SearchLiteratureInput(BaseModel):
    """Single free-text query against the literature RAG."""

    query: str = Field(
        ...,
        description=(
            "Free-text search query. Include peptide / HLA / TCR keywords "
            "plus the biological angle you care about (e.g. "
            "'NY-ESO-1 HLA-A*02:01 affinity-matured TCR engineering')."
        ),
    )
    n_results: int = Field(
        5,
        ge=1,
        le=20,
        description="Number of top papers to return (max 20).",
    )
    year_min: Optional[int] = Field(
        None,
        description="Optional minimum publication year filter.",
    )


class SearchLiteratureTool(Tool):
    """Thin wrapper around `LiteratureRetriever.db.search()` for Bio Agent.

    The retriever's high-level `retrieve_for_mutation_design` builds a
    query from an `analysis_data` dict — that's the batch-recommendation
    use case. Bio Agent issues free-text queries directly, so we skip
    the convenience wrapper and hit the vector DB's `search()`.

    The tool returns a compact markdown table of hits so the LLM can
    pick PMIDs to cite without parsing JSON.
    """

    name = "search_literature"
    description = (
        "Search the curated TCR / pMHC / HLA literature collection. "
        "Returns paper title, PMID, year, and abstract excerpt for each "
        "hit. Use this to gather biology background before emitting "
        "your final message. NEVER cite a PMID that did not come back "
        "from a call to this tool."
    )
    Input = _SearchLiteratureInput
    is_read_only = True
    is_concurrency_safe = True

    def __init__(self, db_path: Optional[str] = None) -> None:
        self._db_path = db_path or _DEFAULT_LIT_DB_PATH
        self._retriever: Any = None  # lazily initialized on first call

    def _ensure_retriever(self) -> Any:
        if self._retriever is not None:
            return self._retriever
        from immunoscope.recommendation.retriever import LiteratureRetriever

        db_path = self._db_path
        # Resolve relative path against the repo root so the same env
        # works from any cwd.
        if not Path(db_path).is_absolute():
            # Walk up until we find a directory containing
            # `development/literature_collection` (the convention
            # used by `recommendation/retriever.py`).
            here = Path(__file__).resolve()
            for ancestor in here.parents:
                candidate = ancestor / db_path
                if candidate.exists():
                    db_path = str(candidate)
                    break
        # Serialize ChromaDB client construction and share one retriever per
        # resolved path across all tool instances (see _RETRIEVER_LOCK note).
        with _RETRIEVER_LOCK:
            retriever = _RETRIEVER_CACHE.get(db_path)
            if retriever is None:
                retriever = LiteratureRetriever(db_path)
                _RETRIEVER_CACHE[db_path] = retriever
        self._retriever = retriever
        return self._retriever

    async def call(
        self, args: _SearchLiteratureInput, ctx: ToolContext
    ) -> ToolResult:
        try:
            retriever = await asyncio.to_thread(self._ensure_retriever)
        except Exception as e:
            log.exception("literature retriever init failed")
            return ToolResult(
                content=(
                    f"Literature retriever unavailable: {e}. "
                    "Skip the call and emit your message without "
                    "literature refs; flag biology claims with "
                    "[NEEDS_VERIFICATION]."
                ),
                is_error=True,
            )

        kwargs: dict = {"query": args.query, "n_results": args.n_results}
        if args.year_min is not None:
            kwargs["year_min"] = args.year_min
        try:
            raw = await asyncio.to_thread(retriever.db.search, **kwargs)
        except Exception as e:
            log.exception("literature search failed")
            return ToolResult(
                content=f"Literature search failed: {e}",
                is_error=True,
            )

        papers = retriever._format_results(raw)
        if not papers:
            return ToolResult(
                content=f"No literature hits for query: {args.query!r}",
                is_error=False,
            )
        return ToolResult(content=self._format_papers(papers), is_error=False)

    @staticmethod
    def _format_papers(papers: List[dict]) -> str:
        lines = [
            "| # | PMID | Year | Title | Abstract excerpt |",
            "|---|------|------|-------|------------------|",
        ]
        for i, p in enumerate(papers, 1):
            title = (p.get("title") or "").replace("|", "/")
            abstract = (p.get("abstract") or "").replace("|", "/").replace("\n", " ")
            if len(abstract) > 220:
                abstract = abstract[:220] + "…"
            lines.append(
                f"| {i} | {p.get('pmid', '?')} | {p.get('year', '?')} "
                f"| {title} | {abstract} |"
            )
        return "\n".join(lines)


class _TcrConservationInput(BaseModel):
    """Read the pre-computed TCR germline-conservation profile for a case."""

    case_dir: str = Field(
        ...,
        description=(
            "The case directory (same value passed to other analysis "
            "queries). The tool reads "
            "`analysis/conservation/conservation_tcr.csv` from it."
        ),
    )
    top_n: int = Field(
        15,
        ge=1,
        le=50,
        description="How many residues to list in the conserved / variable tables.",
    )


class TcrConservationTool(Tool):
    """Surface per-residue TCR germline conservation to the Bio Agent.

    Conservation is sequence-derived (IMGT germline V-gene reference via
    ANARCI), not an MD trajectory metric — it tells the Bio Agent which
    framework positions are evolutionarily load-bearing (do-not-touch) and
    which CDR1/CDR2 positions are germline-variable (tolerant). The CDR3
    junction is not covered by germline conservation and is reported as
    such.

    Implemented as a thin wrapper over the registered `conservation`
    presenter so the markdown the agent sees matches the standard view.
    """

    name = "query_tcr_conservation"
    description = (
        "Read the TCR germline-conservation profile for this case "
        "(per-residue conservation vs. the IMGT germline V-gene reference). "
        "Returns a per-region summary plus the most-conserved (avoid) and "
        "most-variable (tolerant) residues. Use it to ground "
        "`conserved_framework_positions` and anchor-avoidance claims in "
        "data rather than memory. Conservation does NOT cover the CDR3 "
        "junction."
    )
    Input = _TcrConservationInput
    is_read_only = True
    is_concurrency_safe = True

    async def call(
        self, args: _TcrConservationInput, ctx: ToolContext
    ) -> ToolResult:
        def _render() -> str:
            from immunoscope.analysis.features.core.locator import CaseLocator
            from immunoscope.agent.presenters import render

            locator = CaseLocator(args.case_dir)
            return render(
                "conservation", locator, filters={"top_n": args.top_n}
            )

        try:
            markdown = await asyncio.to_thread(_render)
        except Exception as e:  # noqa: BLE001
            log.exception("tcr conservation render failed")
            return ToolResult(
                content=(
                    f"TCR conservation unavailable: {e}. Proceed without "
                    "conservation evidence; do not fabricate conserved "
                    "positions."
                ),
                is_error=True,
            )
        return ToolResult(content=markdown, is_error=False)


class _SearchFactsInput(BaseModel):
    """Semantic search over the extracted biochemical fact-block store."""

    query: str = Field(
        ...,
        description=(
            "Free-text query describing the fact you want — e.g. "
            "'CDR3 mutations that increased affinity for NY-ESO-1' or "
            "'HLA-A*02:01 anchor residue substitution outcomes'."
        ),
    )
    n_results: int = Field(
        8, ge=1, le=20, description="Number of fact blocks to return (max 20)."
    )
    block_type: Optional[str] = Field(
        None,
        description=(
            "Optional filter to one fact type: 'mutation_effect', "
            "'binding_measurement', 'structural_observation', or "
            "'design_heuristic'."
        ),
    )


class SearchFactBlocksTool(Tool):
    """Search the structured biochemical fact blocks mined from full-text papers.

    Unlike `search_literature` (which returns whole abstracts), this returns
    atomic, citable facts extracted from the Results/Discussion/Methods/Tables
    of open-access papers: specific mutation outcomes, binding measurements,
    structural observations, and design heuristics. Every block carries a PMID
    and a VERBATIM evidence quote that was verified against the source — so the
    Bio Agent can populate `prior_engineering` with precisely-cited precedents
    rather than fuzzy abstract-level recall.
    """

    name = "search_fact_blocks"
    description = (
        "Search extracted biochemical fact blocks (specific mutation→effect "
        "outcomes, binding measurements, structural observations, design "
        "heuristics) mined from full-text TCR-pMHC papers. Each hit carries a "
        "PMID and a verbatim evidence quote. Use this to ground "
        "`prior_engineering` and quantitative biology claims in real, citable "
        "results. NEVER cite a PMID or number that did not come back from this "
        "tool. Optionally filter by block_type."
    )
    Input = _SearchFactsInput
    is_read_only = True
    is_concurrency_safe = True

    def __init__(self, db_path: Optional[str] = None) -> None:
        self._db_path = db_path or _DEFAULT_FACTS_DB_PATH
        self._collection: Any = None
        self._embedder: Any = None

    def _ensure_store(self) -> None:
        if self._collection is not None:
            return
        import chromadb
        from sentence_transformers import SentenceTransformer

        persist = _resolve_repo_path(self._db_path)
        # Serialize ChromaDB client construction + embedder load and share the
        # (collection, embedder) pair across tool instances (see
        # _RETRIEVER_LOCK note). get_collection raises if never built — that
        # propagates to call() which degrades gracefully and is not cached.
        cache_key = f"facts::{persist}"
        with _RETRIEVER_LOCK:
            cached = _RETRIEVER_CACHE.get(cache_key)
            if cached is None:
                client = chromadb.PersistentClient(path=persist)
                collection = client.get_collection(_FACTS_COLLECTION)
                embedder = SentenceTransformer(_FACTS_EMBED_MODEL)
                cached = (collection, embedder)
                _RETRIEVER_CACHE[cache_key] = cached
        self._collection, self._embedder = cached

    async def call(self, args: _SearchFactsInput, ctx: ToolContext) -> ToolResult:
        try:
            await asyncio.to_thread(self._ensure_store)
        except Exception as e:  # noqa: BLE001
            log.info("fact-block store unavailable: %s", e)
            return ToolResult(
                content=(
                    f"Fact-block store unavailable ({e}). Fall back to "
                    "`search_literature`; do not fabricate fact blocks."
                ),
                is_error=True,
            )

        where = None
        if args.block_type and args.block_type in _FACT_BLOCK_TYPES:
            where = {"block_type": args.block_type}

        def _query() -> dict:
            emb = self._embedder.encode(args.query).tolist()
            return self._collection.query(
                query_embeddings=[emb], n_results=args.n_results, where=where
            )

        try:
            res = await asyncio.to_thread(_query)
        except Exception as e:  # noqa: BLE001
            log.exception("fact-block query failed")
            return ToolResult(content=f"Fact-block search failed: {e}", is_error=True)

        return ToolResult(content=self._format(res), is_error=False)

    @staticmethod
    def _format(res: dict) -> str:
        docs = (res.get("documents") or [[]])[0]
        metas = (res.get("metadatas") or [[]])[0]
        dists = (res.get("distances") or [[]])[0]
        if not docs:
            return "No fact blocks matched the query."
        lines = ["| # | type | PMID | system | fact (with verbatim quote) |",
                 "|---|------|------|--------|----------------------------|"]
        for i, (doc, meta) in enumerate(zip(docs, metas), 1):
            meta = meta or {}
            bt = meta.get("block_type", "?")
            pmid = meta.get("pmid", "?")
            system = " ".join(
                x for x in (meta.get("tcr", ""), meta.get("hla", ""),
                            meta.get("region", "")) if x
            ) or "-"
            text = (doc or "").replace("|", "/").replace("\n", " ")
            if len(text) > 280:
                text = text[:280] + "…"
            lines.append(f"| {i} | {bt} | {pmid} | {system} | {text} |")
        return "\n".join(lines)


__all__ = ["SearchLiteratureTool", "TcrConservationTool", "SearchFactBlocksTool"]
