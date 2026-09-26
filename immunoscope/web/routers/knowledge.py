"""Knowledge Center — semantic search over the extracted TCR literature facts.

Surfaces the 4739 biochemical fact blocks mined from full-text papers
(`development/literature_collection/`, ChromaDB collection ``tcr_facts``) to web
users. The blocks are already consumed by the BioAgent's `search_fact_blocks`
tool; this exposes the same store as a browsable knowledge module.

Each ChromaDB id is ``{pmid}_{global_index}`` where the index is the line number
in ``fact_blocks.jsonl`` (verified 4739/4739). We use the vector store only for
ranking and recover the *structured* block (separate statement / verbatim quote)
from the jsonl by that index, so cards can render the claim and its evidence
quote distinctly.

Construction of the chromadb client + embedder is serialized + cached (the Rust
bindings race on concurrent first-init, same lesson as multi_agent/tools.py).
"""
from __future__ import annotations

import asyncio
import json
import re
import threading
from collections import Counter
from pathlib import Path
from typing import Any, Optional

from fastapi import APIRouter, Query

from immunoscope.agent.multi_agent.tools import (
    _DEFAULT_FACTS_DB_PATH,
    _FACT_BLOCK_TYPES,
    _FACTS_COLLECTION,
    _FACTS_EMBED_MODEL,
    _resolve_repo_path,
)

router = APIRouter(prefix="/api/knowledge", tags=["knowledge"])

_FACTS_JSONL = "development/literature_collection/knowledge_base/fact_blocks.jsonl"

_LOCK = threading.Lock()
_STATE: dict[str, Any] = {}
_ID_RE = re.compile(r"_(\d+)$")


def _ensure_state() -> dict[str, Any]:
    """Lazy, serialized init of (collection, embedder, blocks). Cached."""
    if _STATE:
        return _STATE
    with _LOCK:
        if _STATE:
            return _STATE
        import chromadb
        from sentence_transformers import SentenceTransformer

        client = chromadb.PersistentClient(path=_resolve_repo_path(_DEFAULT_FACTS_DB_PATH))
        collection = client.get_collection(_FACTS_COLLECTION)
        embedder = SentenceTransformer(_FACTS_EMBED_MODEL)

        blocks: list[dict] = []
        p = Path(_resolve_repo_path(_FACTS_JSONL))
        if p.exists():
            for line in p.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if not line:
                    continue
                try:
                    blocks.append(json.loads(line))
                except json.JSONDecodeError:
                    blocks.append({})
        _STATE.update(collection=collection, embedder=embedder, blocks=blocks)
    return _STATE


def _card(id_: str, meta: Optional[dict], doc: Optional[str],
          dist: Optional[float], blocks: list[dict]) -> dict:
    meta = meta or {}
    blk: dict = {}
    m = _ID_RE.search(id_ or "")
    if m:
        i = int(m.group(1))
        if 0 <= i < len(blocks):
            blk = blocks[i] or {}
    statement = blk.get("statement") or ""
    if not statement and doc:
        # Fallback: the embedded document is "[block_type] statement quote".
        statement = re.sub(r"^\[[^\]]+\]\s*", "", doc)
    system = " ".join(
        x for x in (meta.get("tcr", ""), meta.get("peptide", ""), meta.get("hla", ""))
        if x
    )
    return {
        "id": id_,
        "block_type": meta.get("block_type") or blk.get("block_type") or "",
        "statement": statement,
        "evidence_quote": blk.get("evidence_quote") or "",
        "pmid": str(meta.get("pmid") or blk.get("pmid") or ""),
        "confidence": meta.get("confidence") or blk.get("confidence") or "",
        "source_section": meta.get("source_section") or blk.get("source_section") or "",
        "region": meta.get("region") or "",
        "system": system,
        # cosine distance -> similarity-ish score in [0, 1] for display/sort
        "score": round(max(0.0, 1.0 - float(dist)), 3) if dist is not None else None,
    }


@router.get("/stats")
async def knowledge_stats() -> dict[str, Any]:
    """Total fact blocks + per-type counts (drives the landing + filter chips)."""
    st = await asyncio.to_thread(_ensure_state)
    by_type = Counter(b.get("block_type", "") for b in st["blocks"])
    return {
        "total": len(st["blocks"]),
        "by_type": {t: by_type.get(t, 0) for t in _FACT_BLOCK_TYPES},
        "types": list(_FACT_BLOCK_TYPES),
    }


@router.get("/search")
async def knowledge_search(
    q: str = Query("", description="free-text query"),
    type: str = Query("", description="optional block_type filter"),
    limit: int = Query(20, ge=1, le=50),
) -> dict[str, Any]:
    """Semantic search over the fact blocks; optional block_type filter."""
    st = await asyncio.to_thread(_ensure_state)
    if not q.strip():
        return {"results": [], "query": q, "type": type}

    where = {"block_type": type} if type in _FACT_BLOCK_TYPES else None

    def _run() -> dict:
        emb = st["embedder"].encode(q).tolist()
        return st["collection"].query(
            query_embeddings=[emb],
            n_results=limit,
            where=where,
            include=["metadatas", "documents", "distances"],
        )

    try:
        res = await asyncio.to_thread(_run)
    except Exception as exc:  # noqa: BLE001
        return {"results": [], "query": q, "type": type, "error": str(exc)}

    ids = (res.get("ids") or [[]])[0]
    metas = (res.get("metadatas") or [[]])[0]
    docs = (res.get("documents") or [[]])[0]
    dists = (res.get("distances") or [[]])[0]
    cards = [
        _card(i, m, d, dist, st["blocks"])
        for i, m, d, dist in zip(ids, metas, docs, dists)
    ]
    return {"results": cards, "query": q, "type": type, "count": len(cards)}
