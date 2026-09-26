"""
Biological identity annotation module.

Current goals:
1. Extract basic identity information from complex structures
2. Provide local sequence-library alignment for HLA heavy chains
3. Output a unified JSON structure suitable for HTML reports
"""

from __future__ import annotations

import csv
import json
import logging
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio import Align

from ..structure.pdb_sequence_extractor import PDBSequenceExtractor
from .intelligent_chain_identifier import IntelligentChainIdentifier

logger = logging.getLogger(__name__)


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
HLA_LIBRARY_ROOT = PACKAGE_ROOT / "data" / "reference" / "hla"
HLA_RAW_DIR = HLA_LIBRARY_ROOT / "raw" / "class_i"
HLA_DERIVED_DIR = HLA_LIBRARY_ROOT / "derived"
HLA_DERIVED_FASTA = HLA_DERIVED_DIR / "class_i_extracellular.fasta"
HLA_DERIVED_METADATA = HLA_DERIVED_DIR / "class_i_metadata.csv"

_THREE_TO_ONE = PDBSequenceExtractor.THREE_TO_ONE
_HYDROPHOBIC = set("AILMFWVYCG")


@dataclass
class HLAReferenceEntry:
    allele: str
    locus: str
    source_file: str
    full_length: int
    extracellular_start_1based: int
    extracellular_end_1based: int
    extracellular_sequence: str


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    header: str | None = None
    seq_chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq_chunks)))
            header = line[1:]
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if header is not None:
        entries.append((header, "".join(seq_chunks)))
    return entries


def _extract_allele_from_header(header: str, fallback_locus: str) -> tuple[str, str]:
    tokens = header.split()
    allele_token = next((tok for tok in tokens if "*" in tok), f"{fallback_locus}*unknown")
    locus = allele_token.split("*", 1)[0]
    if not locus.startswith("HLA-"):
        locus = f"HLA-{locus}"
    if not allele_token.startswith("HLA-"):
        allele_token = f"HLA-{allele_token}"
    return allele_token, locus


def _find_mature_start(sequence: str) -> int:
    for motif in ("GSHS", "CSHS"):
        idx = sequence.find(motif)
        if 10 <= idx <= 45:
            return idx
    return 24 if len(sequence) > 300 else 0


def _find_transmembrane_start(sequence: str, search_start: int) -> int:
    best: int | None = None
    for start in range(max(search_start, 220), max(search_start, len(sequence) - 12)):
        window = sequence[start : start + 14]
        if len(window) < 14:
            break
        hydrophobic_fraction = sum(1 for aa in window if aa in _HYDROPHOBIC) / len(window)
        if hydrophobic_fraction >= 0.78:
            best = start
            break
    if best is not None:
        return best
    return min(len(sequence), search_start + 276)


def derive_hla_reference_library(
    raw_dir: Path = HLA_RAW_DIR,
    derived_dir: Path = HLA_DERIVED_DIR,
) -> tuple[Path, Path]:
    derived_dir.mkdir(parents=True, exist_ok=True)

    entries: list[HLAReferenceEntry] = []
    for fasta_path in sorted(raw_dir.glob("*_prot.fasta")):
        default_locus = f"HLA-{fasta_path.stem.split('_')[0]}"
        for header, sequence in _parse_fasta(fasta_path):
            allele, locus = _extract_allele_from_header(header, default_locus)
            mature_start = _find_mature_start(sequence)
            tm_start = _find_transmembrane_start(sequence, mature_start + 220)
            extracellular = sequence[mature_start:tm_start]
            entries.append(
                HLAReferenceEntry(
                    allele=allele,
                    locus=locus,
                    source_file=fasta_path.name,
                    full_length=len(sequence),
                    extracellular_start_1based=mature_start + 1,
                    extracellular_end_1based=tm_start,
                    extracellular_sequence=extracellular,
                )
            )

    with HLA_DERIVED_FASTA.open("w", encoding="utf-8") as handle:
        for entry in entries:
            handle.write(
                f">{entry.allele}|locus={entry.locus}|region=extracellular|start={entry.extracellular_start_1based}|end={entry.extracellular_end_1based}\n"
            )
            for start in range(0, len(entry.extracellular_sequence), 70):
                handle.write(entry.extracellular_sequence[start : start + 70] + "\n")

    with HLA_DERIVED_METADATA.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "allele",
                "locus",
                "source_file",
                "full_length",
                "extracellular_start_1based",
                "extracellular_end_1based",
                "extracellular_length",
                "sequence",
            ],
        )
        writer.writeheader()
        for entry in entries:
            writer.writerow(
                {
                    "allele": entry.allele,
                    "locus": entry.locus,
                    "source_file": entry.source_file,
                    "full_length": entry.full_length,
                    "extracellular_start_1based": entry.extracellular_start_1based,
                    "extracellular_end_1based": entry.extracellular_end_1based,
                    "extracellular_length": len(entry.extracellular_sequence),
                    "sequence": entry.extracellular_sequence,
                }
            )

    logger.info("Derived HLA library generated: %s", HLA_DERIVED_METADATA)
    return HLA_DERIVED_FASTA, HLA_DERIVED_METADATA


def ensure_hla_reference_library() -> tuple[Path, Path]:
    if HLA_DERIVED_FASTA.exists() and HLA_DERIVED_METADATA.exists():
        return HLA_DERIVED_FASTA, HLA_DERIVED_METADATA
    return derive_hla_reference_library()


def _build_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -1.0
    return aligner


def _alignment_metrics(query: str, target: str, alignment: Align.Alignment) -> tuple[float, float]:
    q_aln, t_aln = alignment[0], alignment[1]
    compared = 0
    matches = 0
    aligned_query = 0
    for qa, ta in zip(q_aln, t_aln):
        if qa != "-":
            aligned_query += 1
        if qa == "-" or ta == "-":
            continue
        compared += 1
        if qa == ta:
            matches += 1
    identity = matches / compared if compared else 0.0
    coverage = aligned_query / len(query) if query else 0.0
    return identity, coverage


class HLAIdentityAnnotator:
    """Local alignment annotator for HLA locus and candidate alleles."""

    def __init__(self, metadata_csv: Path | None = None):
        _, default_metadata = ensure_hla_reference_library()
        self.metadata_csv = metadata_csv or default_metadata
        self.aligner = _build_aligner()
        self.entries = self._load_entries()

    def _load_entries(self) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        with self.metadata_csv.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rows.append(row)
        return rows

    @staticmethod
    def _confidence(identity: float, coverage: float) -> str:
        if identity >= 0.98 and coverage >= 0.95:
            return "high"
        if identity >= 0.94 and coverage >= 0.90:
            return "medium"
        return "low"

    def annotate(self, query_sequence: str, top_n: int = 3) -> dict[str, Any]:
        scored: list[dict[str, Any]] = []
        for row in self.entries:
            target = row["sequence"]
            alignment = self.aligner.align(query_sequence, target)[0]
            identity, coverage = _alignment_metrics(query_sequence, target, alignment)
            scored.append(
                {
                    "allele": row["allele"],
                    "locus": row["locus"],
                    "identity": identity,
                    "coverage": coverage,
                    "score": float(alignment.score),
                }
            )
        scored.sort(key=lambda item: (item["identity"], item["coverage"], item["score"]), reverse=True)
        best = scored[0] if scored else {}
        return {
            "query_length": len(query_sequence),
            "best_locus": best.get("locus", ""),
            "best_candidate_allele": best.get("allele", ""),
            "identity": round(float(best.get("identity", 0.0)), 4),
            "coverage": round(float(best.get("coverage", 0.0)), 4),
            "confidence": self._confidence(float(best.get("identity", 0.0)), float(best.get("coverage", 0.0))),
            "top_candidates": [
                {
                    "allele": item["allele"],
                    "locus": item["locus"],
                    "identity": round(float(item["identity"]), 4),
                    "coverage": round(float(item["coverage"]), 4),
                    "score": round(float(item["score"]), 2),
                }
                for item in scored[:top_n]
            ],
        }


class TCRIdentityAnnotator:
    """TCR V/J annotator based on ANARCI germline assignment."""

    _GERMLINE_HEADER = "#|species|v_gene|v_identity|j_gene|j_identity|"

    def __init__(self, numbering_scheme: str = "imgt"):
        self.numbering_scheme = numbering_scheme
        self.anarci_executable = shutil.which("ANARCI")

    def _run_assign_germline(self, sequence: str, label: str) -> str:
        if not self.anarci_executable or not sequence:
            return ""

        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as handle:
            handle.write(f">{label}\n{sequence}\n")
            fasta_path = Path(handle.name)

        try:
            result = subprocess.run(
                [
                    self.anarci_executable,
                    "-i",
                    str(fasta_path),
                    "--scheme",
                    self.numbering_scheme,
                    "--assign_germline",
                ],
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            return result.stdout
        except (subprocess.SubprocessError, FileNotFoundError) as exc:
            logger.warning("ANARCI germline assignment failed (%s): %s", label, exc)
            return ""
        finally:
            fasta_path.unlink(missing_ok=True)

    def _parse_assign_germline(self, output: str) -> dict[str, Any]:
        if not output:
            return {}

        lines = [line.strip() for line in output.splitlines() if line.strip()]
        try:
            header_index = lines.index(self._GERMLINE_HEADER)
        except ValueError:
            return {}

        if header_index + 1 >= len(lines):
            return {}

        values = lines[header_index + 1].strip("#|").split("|")
        if len(values) < 5:
            return {}

        species, v_gene, v_identity, j_gene, j_identity = values[:5]
        parsed = {
            "species": species,
            "v_gene": v_gene,
            "j_gene": j_gene,
            "v_identity": float(v_identity) if v_identity else 0.0,
            "j_identity": float(j_identity) if j_identity else 0.0,
            "confidence": self._confidence(v_identity, j_identity),
            "source": "anarci_assign_germline",
        }
        return parsed

    @staticmethod
    def _confidence(v_identity: str, j_identity: str) -> str:
        try:
            v_val = float(v_identity)
            j_val = float(j_identity)
        except ValueError:
            return "low"
        if v_val >= 0.98 and j_val >= 0.98:
            return "high"
        if v_val >= 0.94 and j_val >= 0.94:
            return "medium"
        return "low"

    def annotate_pair(self, alpha_sequence: str, beta_sequence: str) -> dict[str, Any]:
        alpha_output = self._run_assign_germline(alpha_sequence, "TCR_alpha")
        beta_output = self._run_assign_germline(beta_sequence, "TCR_beta")
        return {
            "alpha": self._parse_assign_germline(alpha_output),
            "beta": self._parse_assign_germline(beta_output),
        }


class BiologicalIdentityAnnotator:
    """Integrated biological identity annotator for complexes."""

    def __init__(self):
        self.sequence_extractor = PDBSequenceExtractor()
        self.chain_identifier = IntelligentChainIdentifier(use_anarci=True)
        self.hla_annotator = HLAIdentityAnnotator()
        self.tcr_annotator = TCRIdentityAnnotator()

    @staticmethod
    def _build_chain_mapping(chain_results: dict[str, Any]) -> dict[str, str]:
        chain_mapping: dict[str, str] = {}
        for chain_id, chain_info in chain_results.items():
            chain_type = chain_info.chain_type
            if chain_type == "HLA_alpha":
                chain_mapping["mhc_alpha"] = chain_id
            elif chain_type == "beta2m":
                chain_mapping["b2m"] = chain_id
            elif chain_type == "peptide":
                chain_mapping["peptide"] = chain_id
            elif chain_type == "TCR_alpha":
                chain_mapping["tcr_alpha"] = chain_id
            elif chain_type == "TCR_beta":
                chain_mapping["tcr_beta"] = chain_id
        return chain_mapping

    @staticmethod
    def _load_cdr_metadata(cdr_metadata_path: Path | None) -> dict[str, Any]:
        if cdr_metadata_path and cdr_metadata_path.exists():
            return json.loads(cdr_metadata_path.read_text(encoding="utf-8"))
        return {}

    @staticmethod
    def _normalize_chain_mapping(chain_mapping: dict[str, Any] | None) -> dict[str, str]:
        if not chain_mapping:
            return {}
        return {str(key): str(value) for key, value in chain_mapping.items() if value}

    @staticmethod
    def _extract_cdr3_sequences(
        cdr_detection: dict[str, Any] | None,
        cdr_metadata: dict[str, Any],
    ) -> tuple[str, str, str]:
        if cdr_detection:
            chains = cdr_detection.get("chains", {}) or {}
            alpha_seq = str((chains.get("TCR_alpha", {}).get("cdrs", {}) or {}).get(3, {}).get("sequence", "") or "")
            beta_seq = str((chains.get("TCR_beta", {}).get("cdrs", {}) or {}).get(3, {}).get("sequence", "") or "")
            if alpha_seq or beta_seq:
                return alpha_seq, beta_seq, "pipeline_cdr_detection"

        cdr_chains = cdr_metadata.get("chains", {}) or {}
        alpha_meta = cdr_chains.get("TCR_alpha", {}) or {}
        beta_meta = cdr_chains.get("TCR_beta", {}) or {}
        return (
            str((alpha_meta.get("cdr3", {}) or {}).get("sequence", "") or ""),
            str((beta_meta.get("cdr3", {}) or {}).get("sequence", "") or ""),
            "anarci_cdr_metadata",
        )

    def annotate(
        self,
        structure_pdb: Path,
        cdr_metadata_path: Path | None = None,
        chain_mapping: dict[str, Any] | None = None,
        cdr_detection: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        if chain_mapping:
            normalized_chain_mapping = self._normalize_chain_mapping(chain_mapping)
        else:
            chain_results = self.chain_identifier.identify_chains(str(structure_pdb))
            normalized_chain_mapping = self._build_chain_mapping(chain_results)
        chains_data = self.sequence_extractor.extract_sequences_from_pdb(str(structure_pdb))
        cdr_metadata = self._load_cdr_metadata(cdr_metadata_path)

        peptide_chain = normalized_chain_mapping.get("peptide", "")
        mhc_chain = normalized_chain_mapping.get("mhc_alpha", "")
        tcr_alpha_chain = normalized_chain_mapping.get("tcr_alpha", "")
        tcr_beta_chain = normalized_chain_mapping.get("tcr_beta", "")

        peptide_sequence = chains_data.get(peptide_chain, {}).get("sequence", "")
        mhc_sequence = chains_data.get(mhc_chain, {}).get("sequence", "")
        tcr_alpha_sequence = chains_data.get(tcr_alpha_chain, {}).get("sequence", "")
        tcr_beta_sequence = chains_data.get(tcr_beta_chain, {}).get("sequence", "")
        cdr3_alpha, cdr3_beta, tcr_source = self._extract_cdr3_sequences(cdr_detection, cdr_metadata)

        hla_identity = self.hla_annotator.annotate(mhc_sequence) if mhc_sequence else {}
        tcr_germline = self.tcr_annotator.annotate_pair(tcr_alpha_sequence, tcr_beta_sequence)
        alpha_germline = tcr_germline.get("alpha", {}) or {}
        beta_germline = tcr_germline.get("beta", {}) or {}

        return {
            "complex_identity": {
                "mhc_alpha_chain": mhc_chain,
                "b2m_chain": normalized_chain_mapping.get("b2m", ""),
                "peptide_chain": peptide_chain,
                "tcr_alpha_chain": tcr_alpha_chain,
                "tcr_beta_chain": tcr_beta_chain,
            },
            "peptide_identity": {
                "sequence": peptide_sequence,
                "length": len(peptide_sequence),
            },
            "tcr_identity": {
                "alpha_length": chains_data.get(tcr_alpha_chain, {}).get("length", 0),
                "beta_length": chains_data.get(tcr_beta_chain, {}).get("length", 0),
                "alpha_v_gene": alpha_germline.get("v_gene", ""),
                "alpha_j_gene": alpha_germline.get("j_gene", ""),
                "alpha_v_identity": alpha_germline.get("v_identity", 0.0),
                "alpha_j_identity": alpha_germline.get("j_identity", 0.0),
                "alpha_genotype_confidence": alpha_germline.get("confidence", "low"),
                "beta_v_gene": beta_germline.get("v_gene", ""),
                "beta_j_gene": beta_germline.get("j_gene", ""),
                "beta_v_identity": beta_germline.get("v_identity", 0.0),
                "beta_j_identity": beta_germline.get("j_identity", 0.0),
                "beta_genotype_confidence": beta_germline.get("confidence", "low"),
                "cdr3_alpha_sequence": cdr3_alpha,
                "cdr3_beta_sequence": cdr3_beta,
                "source": tcr_source,
                "genotype_source": "anarci_assign_germline" if (alpha_germline or beta_germline) else "",
            },
            "hla_identity": hla_identity,
            "sequence_sources": {
                "structure_pdb": str(structure_pdb),
                "cdr_metadata": str(cdr_metadata_path) if cdr_metadata_path else "",
                "hla_reference_metadata": str(self.hla_annotator.metadata_csv),
            },
        }
