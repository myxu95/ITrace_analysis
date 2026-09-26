"""
文献检索器

根据分析数据智能检索相关文献
"""
import importlib.util
import logging
from pathlib import Path
from typing import List, Dict, Optional

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Sister scripts directory inside this repo. Not on sys.path by default; we
# load vector_db.py by explicit file path from inside the constructor so
# that importing immunoscope.recommendation does not crash when chromadb /
# sentence-transformers are absent.
_DEFAULT_LITERATURE_COLLECTION_DIR = (
    Path(__file__).resolve().parent.parent.parent
    / "development"
    / "literature_collection"
)


def _load_vector_db_module(literature_collection_dir: Path):
    """Load the LiteratureVectorDB class without mutating sys.path."""
    module_path = literature_collection_dir / "vector_db.py"
    if not module_path.exists():
        raise FileNotFoundError(
            f"literature_collection/vector_db.py not found at {module_path}. "
            "Pass literature_collection_dir explicitly or install the "
            "literature collection scripts alongside this repo."
        )
    spec = importlib.util.spec_from_file_location(
        "immunoscope_literature_vector_db", module_path
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load vector_db spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.LiteratureVectorDB


class LiteratureRetriever:
    """文献检索器"""

    def __init__(
        self,
        db_path: str,
        literature_collection_dir: Optional[Path] = None,
    ):
        """
        初始化检索器

        Args:
            db_path: 向量数据库路径
            literature_collection_dir: optional override for the directory
                hosting vector_db.py. Defaults to the sibling
                ``development/literature_collection`` directory in this repo.
        """
        lit_dir = Path(literature_collection_dir or _DEFAULT_LITERATURE_COLLECTION_DIR)
        vector_db_cls = _load_vector_db_module(lit_dir)
        self.db = vector_db_cls(db_path)
        logger.info(f"初始化文献检索器: {db_path}")

    def retrieve_for_mutation_design(
        self,
        analysis_data: Dict,
        design_goal: str,
        n_results: int = 10
    ) -> List[Dict]:
        """
        根据分析数据和设计目标检索相关文献

        Args:
            analysis_data: 分析数据摘要（来自AnalysisDataFormatter）
            design_goal: 设计目标（affinity_enhancement/specificity/stability等）
            n_results: 返回结果数量

        Returns:
            文献列表，每篇文献包含title, abstract, pmid, year, relevance_score
        """
        logger.info(f"检索文献: 设计目标={design_goal}, 返回数量={n_results}")

        # 1. 从分析数据提取关键信息
        key_info = self._extract_key_info(analysis_data)

        # 2. 构建查询
        query = self._build_query(key_info, design_goal)
        logger.info(f"查询语句: {query}")

        # 3. 检索
        results = self.db.search(
            query=query,
            n_results=n_results,
            year_min=2015  # 只检索2015年后的文献
        )

        # 4. 格式化结果
        formatted_results = self._format_results(results)

        logger.info(f"检索到 {len(formatted_results)} 篇相关文献")
        return formatted_results

    def _extract_key_info(self, analysis_data: Dict) -> Dict:
        """从分析数据提取关键信息"""
        metadata = analysis_data.get("metadata", {})
        candidate_residues = analysis_data.get("candidate_residues", [])

        # 提取top残基的区域分布
        regions = [r.get("region", "") for r in candidate_residues[:5]]
        region_counts = {}
        for region in regions:
            if region:
                region_counts[region] = region_counts.get(region, 0) + 1

        # 找出主要区域
        dominant_region = max(region_counts.items(), key=lambda x: x[1])[0] if region_counts else "CDR3"

        return {
            "peptide": metadata.get("peptide", ""),
            "hla": metadata.get("hla", ""),
            "tcr_alpha_v": metadata.get("tcr_alpha_v", ""),
            "tcr_beta_v": metadata.get("tcr_beta_v", ""),
            "dominant_region": dominant_region,
            "top_residues": [r.get("residue", "") for r in candidate_residues[:3]],
            "n_candidates": len(candidate_residues)
        }

    def _build_query(self, key_info: Dict, design_goal: str) -> str:
        """构建智能查询"""
        # 设计目标映射
        goal_keywords = {
            "affinity_enhancement": "affinity enhancement binding improvement",
            "specificity": "specificity selectivity cross-reactivity",
            "stability": "stability thermostability expression",
            "immunogenicity": "immunogenicity immunogenic safety"
        }

        goal_text = goal_keywords.get(design_goal, design_goal)

        # 构建查询（优先考虑区域和设计目标）
        dominant_region = key_info.get("dominant_region", "CDR3")

        query_parts = [
            f"TCR {dominant_region} mutation",
            goal_text,
            "rational design",
            "structure-based"
        ]

        # 如果有HLA信息，添加到查询
        hla = key_info.get("hla", "")
        if hla:
            query_parts.append(hla)

        query = " ".join(query_parts)
        return query

    def _format_results(self, results: Dict) -> List[Dict]:
        """Normalize a ChromaDB query() result into a flat list of papers.

        ChromaDB returns one outer list per query string (we only ever send
        a single query, so we want index 0). Each field may also come back
        as a flat list when callers use lower-level APIs, so we tolerate
        both shapes rather than crashing with IndexError.
        """

        def _first_batch(payload, key: str) -> list:
            raw = payload.get(key)
            if raw is None:
                return []
            if not isinstance(raw, list) or not raw:
                return []
            # query()-style: list of lists, take batch 0.
            if isinstance(raw[0], list):
                return raw[0]
            # get()-style: already a flat list of items.
            return raw

        ids = _first_batch(results, "ids")
        distances = _first_batch(results, "distances")
        metadatas = _first_batch(results, "metadatas")
        documents = _first_batch(results, "documents")

        formatted: List[Dict] = []
        for i, doc_id in enumerate(ids):
            metadata = metadatas[i] if i < len(metadatas) else {}
            document = documents[i] if i < len(documents) else ""
            distance = distances[i] if i < len(distances) else None

            abstract = document.split("\n\n", 1)[1] if "\n\n" in document else document
            authors = metadata.get("authors", "") or ""

            formatted.append({
                "pmid": metadata.get("pmid", doc_id),
                "title": metadata.get("title", ""),
                "abstract": abstract[:500] + "..." if len(abstract) > 500 else abstract,
                "year": metadata.get("year", 0),
                "relevance_score": round(distance, 3) if distance is not None else 0.0,
                "journal": metadata.get("journal", ""),
                "authors": authors[:100] + "..." if len(authors) > 100 else authors,
            })

        return formatted


# 使用示例
if __name__ == "__main__":
    import json

    # 加载分析数据
    with open("output/5c0a_run2_full_analysis/analysis_summary_for_llm.json", "r") as f:
        analysis_data = json.load(f)

    # 初始化检索器
    retriever = LiteratureRetriever("development/literature_collection/knowledge_base/vector_db")

    # 检索文献
    papers = retriever.retrieve_for_mutation_design(
        analysis_data=analysis_data,
        design_goal="affinity_enhancement",
        n_results=5
    )

    # 打印结果
    print(f"\n检索到 {len(papers)} 篇相关文献:\n")
    for i, paper in enumerate(papers, 1):
        print(f"{i}. [{paper['year']}] {paper['title']}")
        print(f"   PMID: {paper['pmid']}, Relevance: {paper['relevance_score']}")
        print(f"   Abstract: {paper['abstract'][:200]}...")
        print()
