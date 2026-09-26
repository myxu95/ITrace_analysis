#!/usr/bin/env python3
"""
PDB Chain Standardization Module for Protein Complexes

This module provides functionality to standardize chain identifiers in PDB files
based on residue count ordering. Primarily designed for pHLA-TCR complexes but
can be extended to other multi-chain protein complexes.

Standard ordering for pHLA-TCR complexes (by residue count):
  Chain C: Peptide (shortest, 5-22 residues)
  Chain B: HLA-β/β2-microglobulin (~99-100 residues)
  Chain D: TCR-α (shorter TCR chain, 115-215 residues)
  Chain E: TCR-β (longer TCR chain, 209-252 residues)
  Chain A: HLA-α (longest, ~270-280 residues)

Author: ImmunoScope Team
License: MIT
"""

from pathlib import Path
from collections import defaultdict
from datetime import datetime
from multiprocessing import Pool
from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass
import csv
import logging
import tempfile
import shutil

try:
    from ..topology import IntelligentChainIdentifier
    INTELLIGENT_IDENTIFIER_AVAILABLE = True
except ImportError:
    INTELLIGENT_IDENTIFIER_AVAILABLE = False
    IntelligentChainIdentifier = None

__all__ = ['PDBChainStandardizer', 'ChainInfo', 'StandardizationResult']


@dataclass
class ChainInfo:
    """
    Information about a protein chain in a PDB file.

    Attributes:
        chain_id: Chain identifier (single character)
        residue_count: Number of residues in the chain
        atom_count: Number of atoms in the chain
    """
    chain_id: str
    residue_count: int
    atom_count: int

    def __repr__(self) -> str:
        return f"ChainInfo({self.chain_id}, {self.residue_count} res, {self.atom_count} atoms)"


@dataclass
class StandardizationResult:
    """
    Result of PDB chain standardization for a single task.

    Attributes:
        task_name: Name of the task/sample
        status: Processing status (OK, ALREADY_STANDARD, MULTICHAIN, INSUFFICIENT, ERROR, NO_PDB)
        num_chains: Number of protein chains detected
        chain_mapping: Dictionary mapping old chain IDs to new chain IDs
        input_file: Path to input PDB file
        output_file: Path to output PDB file (if created)
        processing_time: Time taken to process in seconds
        error_message: Error message if processing failed
    """
    task_name: str
    status: str
    num_chains: int
    chain_mapping: Dict[str, str]
    input_file: Optional[Path] = None
    output_file: Optional[Path] = None
    processing_time: float = 0.0
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Union[str, int, float]]:
        """Convert result to dictionary format for CSV export."""
        mapping_str = " / ".join([f"{old}→{new}" for old, new in self.chain_mapping.items()])
        return {
            "TaskName": self.task_name,
            "Status": self.status,
            "NumChains": self.num_chains,
            "ChainMapping": mapping_str,
            "ProcessingTime_sec": self.processing_time,
            "Error": self.error_message or ""
        }


class PDBChainStandardizer:
    """
    Standardize chain identifiers in PDB files based on residue count ordering.

    This class provides methods to analyze protein chains, create chain mappings,
    and standardize PDB files by renumbering chains according to a defined order.

    Attributes:
        standard_order: List of chain IDs in the desired standard order
        expected_chain_count: Expected number of chains in the complex
        max_residue_threshold: Maximum residues to consider as protein (filters out solvent)
        logger: Logger instance for this class

    Example:
        >>> standardizer = PDBChainStandardizer()
        >>> result = standardizer.standardize_file(
        ...     input_pdb="input.pdb",
        ...     output_pdb="output.pdb"
        ... )
        >>> print(result.status)  # 'OK' if successful
    """

    def __init__(
        self,
        standard_order: Optional[List[str]] = None,
        expected_chain_count: int = 5,
        max_residue_threshold: int = 1000,
        use_intelligent_identification: bool = False,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize PDB Chain Standardizer.

        Args:
            standard_order: List of chain IDs in desired order (default: ['C', 'B', 'D', 'E', 'A'])
            expected_chain_count: Expected number of chains (default: 5 for pHLA-TCR)
            max_residue_threshold: Max residues to consider as protein, filters solvent (default: 1000)
            use_intelligent_identification: Use ANARCI-based intelligent chain identification (default: False)
            logger: Custom logger instance (optional)
        """
        self.standard_order = standard_order or ['C', 'B', 'D', 'E', 'A']
        self.expected_chain_count = expected_chain_count
        self.max_residue_threshold = max_residue_threshold
        self.use_intelligent_identification = use_intelligent_identification
        self.logger = logger or self._setup_default_logger()

        # Initialize intelligent identifier if requested
        if use_intelligent_identification:
            if not INTELLIGENT_IDENTIFIER_AVAILABLE:
                self.logger.warning(
                    "Intelligent chain identification requested but module not available. "
                    "Falling back to length-based identification."
                )
                self.use_intelligent_identification = False
                self.intelligent_identifier = None
            else:
                try:
                    self.intelligent_identifier = IntelligentChainIdentifier(use_anarci=True)
                    self.logger.info("Intelligent chain identification enabled (ANARCI-based)")
                except Exception as e:
                    self.logger.warning(
                        f"Failed to initialize intelligent identifier: {e}. "
                        "Falling back to length-based identification."
                    )
                    self.use_intelligent_identification = False
                    self.intelligent_identifier = None
        else:
            self.intelligent_identifier = None

    def _setup_default_logger(self) -> logging.Logger:
        """Set up default logger for the class."""
        logger = logging.getLogger(__name__)
        if not logger.handlers:
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        return logger

    def analyze_chains(self, pdb_file: Union[str, Path]) -> List[ChainInfo]:
        """
        Analyze chain composition in a PDB file.

        Reads through the PDB file and extracts information about each protein chain,
        including residue count and atom count. Filters out solvent chains based on
        the max_residue_threshold.

        Args:
            pdb_file: Path to PDB file to analyze

        Returns:
            List of ChainInfo objects, sorted by residue count (ascending)

        Raises:
            FileNotFoundError: If PDB file does not exist
            ValueError: If PDB file is empty or malformed

        Example:
            >>> chains = standardizer.analyze_chains("protein.pdb")
            >>> for chain in chains:
            ...     print(f"Chain {chain.chain_id}: {chain.residue_count} residues")
        """
        pdb_file = Path(pdb_file)

        if not pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")

        chains = defaultdict(lambda: {'atoms': 0, 'residues': set()})

        with open(pdb_file) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain_id = line[21]
                    if chain_id.strip():  # Skip empty chain IDs (solvent)
                        try:
                            res_num = int(line[22:26].strip())
                            chains[chain_id]['atoms'] += 1
                            chains[chain_id]['residues'].add(res_num)
                        except ValueError:
                            # Skip malformed residue numbers
                            continue

        if not chains:
            raise ValueError(f"No protein chains found in {pdb_file}")

        # Convert to ChainInfo objects and filter solvent
        chain_list = []
        for chain_id, info in chains.items():
            residue_count = len(info['residues'])
            atom_count = info['atoms']

            # Filter out solvent chains
            if residue_count < self.max_residue_threshold:
                chain_list.append(ChainInfo(
                    chain_id=chain_id,
                    residue_count=residue_count,
                    atom_count=atom_count
                ))

        # Sort by residue count (ascending)
        chain_list.sort(key=lambda x: x.residue_count)

        return chain_list

    def create_mapping(
        self,
        input_pdb: Union[str, Path],
        chain_list: Optional[List[ChainInfo]] = None
    ) -> Tuple[Dict[str, str], str]:
        """
        Create chain ID mapping using intelligent or length-based method.

        If intelligent identification is enabled, uses ANARCI to identify TCR chains.
        Otherwise, falls back to simple residue count ordering.

        Args:
            input_pdb: Path to input PDB file (required for intelligent mode)
            chain_list: Optional pre-analyzed chain list (for length-based mode)

        Returns:
            Tuple of:
                - Dictionary mapping old chain IDs to new chain IDs
                - Status string: "OK", "MULTICHAIN", "INSUFFICIENT", or "ERROR"
        """
        # Backward compatibility: older callers passed only a pre-analyzed
        # chain list as the first positional argument.
        if chain_list is None and isinstance(input_pdb, list):
            chain_list = input_pdb
            input_pdb = ""

        # Use intelligent identification if enabled
        if self.use_intelligent_identification and self.intelligent_identifier:
            try:
                return self._create_mapping_intelligent(input_pdb)
            except Exception as e:
                self.logger.warning(
                    f"Intelligent identification failed: {e}. "
                    "Falling back to length-based method."
                )
                # Fall through to length-based method

        # Fall back to length-based method
        return self._create_mapping_by_length(chain_list)

    def _create_mapping_intelligent(
        self,
        input_pdb: Union[str, Path]
    ) -> Tuple[Dict[str, str], str]:
        """
        Create mapping using ANARCI-based intelligent identification.

        Args:
            input_pdb: Path to input PDB file

        Returns:
            Tuple of (mapping dict, status string)
        """
        # Use intelligent identifier
        identifications = self.intelligent_identifier.identify_chains(str(input_pdb))

        if len(identifications) != self.expected_chain_count:
            self.logger.warning(
                f"Expected {self.expected_chain_count} chains, identified {len(identifications)}"
            )
            if len(identifications) > self.expected_chain_count:
                return {}, "MULTICHAIN"
            else:
                return {}, "INSUFFICIENT"

        # Create mapping from identifications
        mapping = self.intelligent_identifier.create_standardization_mapping(identifications)

        # Log confidence scores
        for chain_id, info in identifications.items():
            self.logger.info(
                f"  {chain_id} → {mapping.get(chain_id, '?')}: "
                f"{info.chain_type} ({info.length} aa, confidence={info.confidence:.2f})"
            )

        return mapping, "OK"

    def _create_mapping_by_length(
        self,
        chain_list: List[ChainInfo]
    ) -> Tuple[Dict[str, str], str]:
        """
        Create mapping based on residue count ordering (legacy method).

        Maps original chain IDs to standard chain IDs according to the predefined
        standard_order, based on residue count (smallest to largest).

        Args:
            chain_list: List of ChainInfo objects, sorted by residue count

        Returns:
            Tuple of:
                - Dictionary mapping old chain IDs to new chain IDs
                - Status string: "OK", "MULTICHAIN", or "INSUFFICIENT"
        """
        num_chains = len(chain_list)

        if num_chains > self.expected_chain_count:
            self.logger.warning(
                f"Found {num_chains} chains, expected {self.expected_chain_count}. "
                "Possible multichain complex."
            )
            return {}, "MULTICHAIN"

        if num_chains < self.expected_chain_count:
            self.logger.warning(
                f"Found {num_chains} chains, expected {self.expected_chain_count}. "
                "Insufficient chains detected."
            )
            return {}, "INSUFFICIENT"

        # Create mapping: old_chain_id -> new_chain_id
        mapping = {}
        for i, chain_info in enumerate(chain_list):
            old_chain = chain_info.chain_id
            new_chain = self.standard_order[i]
            mapping[old_chain] = new_chain

        return mapping, "OK"

    def is_already_standard(
        self,
        chain_list: List[ChainInfo],
        chain_mapping: Dict[str, str]
    ) -> bool:
        """
        Check if PDB file already has standard chain ordering.

        Args:
            chain_list: List of ChainInfo objects (sorted by residue count)
            chain_mapping: Chain ID mapping dictionary

        Returns:
            True if chains are already in standard order, False otherwise
        """
        current_order = [chain.chain_id for chain in chain_list]
        return current_order == self.standard_order

    def standardize_file(
        self,
        input_pdb: Union[str, Path],
        output_pdb: Union[str, Path],
        chain_mapping: Dict[str, str]
    ) -> None:
        """
        Rewrite PDB file with standardized chain identifiers.

        Reads the input PDB file line by line and rewrites it with new chain IDs
        according to the provided mapping. Non-protein records (TER, CONECT, etc.)
        are preserved unchanged.

        IMPORTANT: Safely handles case where input and output are the same file.
        Uses temporary file to avoid data loss when input_pdb == output_pdb.

        Args:
            input_pdb: Path to input PDB file
            output_pdb: Path to output PDB file
            chain_mapping: Dictionary mapping old chain IDs to new chain IDs

        Raises:
            FileNotFoundError: If input PDB file does not exist
            PermissionError: If cannot write to output location
        """
        input_pdb = Path(input_pdb).resolve()  # Get absolute path
        output_pdb = Path(output_pdb).resolve()  # Get absolute path

        if not input_pdb.exists():
            raise FileNotFoundError(f"Input PDB file not found: {input_pdb}")

        # Ensure output directory exists
        output_pdb.parent.mkdir(parents=True, exist_ok=True)

        # Check if input and output are the same file
        is_same_file = input_pdb == output_pdb

        if is_same_file:
            # Use temporary file to avoid data loss
            temp_fd, temp_path = tempfile.mkstemp(
                suffix='.pdb',
                dir=output_pdb.parent,
                text=True
            )
            temp_file = Path(temp_path)

            try:
                # Write to temporary file
                with open(input_pdb) as f_in, open(temp_file, 'w') as f_out:
                    for line in f_in:
                        if line.startswith('ATOM') or line.startswith('HETATM'):
                            old_chain = line[21]
                            if old_chain in chain_mapping:
                                new_chain = chain_mapping[old_chain]
                                # Rewrite chain identifier (column 22, index 21)
                                new_line = line[:21] + new_chain + line[22:]
                                f_out.write(new_line)
                            else:
                                # Keep unchanged (solvent, etc.)
                                f_out.write(line)
                        else:
                            # Non-atom records preserved unchanged
                            f_out.write(line)

                # Verify temp file is not empty
                if temp_file.stat().st_size == 0:
                    raise ValueError("Standardized file is empty - possible processing error")

                # Replace original file with temp file
                shutil.move(str(temp_file), str(output_pdb))

            except Exception as e:
                # Clean up temp file on error
                if temp_file.exists():
                    temp_file.unlink()
                raise e
            finally:
                # Close the file descriptor
                import os
                try:
                    os.close(temp_fd)
                except OSError:
                    pass
        else:
            # Input and output are different - direct write is safe
            with open(input_pdb) as f_in, open(output_pdb, 'w') as f_out:
                for line in f_in:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        old_chain = line[21]
                        if old_chain in chain_mapping:
                            new_chain = chain_mapping[old_chain]
                            # Rewrite chain identifier (column 22, index 21)
                            new_line = line[:21] + new_chain + line[22:]
                            f_out.write(new_line)
                        else:
                            # Keep unchanged (solvent, etc.)
                            f_out.write(line)
                    else:
                        # Non-atom records preserved unchanged
                        f_out.write(line)

    def process_single(
        self,
        input_pdb: Union[str, Path],
        output_pdb: Union[str, Path],
        task_name: Optional[str] = None,
        skip_if_standard: bool = True
    ) -> StandardizationResult:
        """
        Process a single PDB file for chain standardization.

        Complete workflow: analyze chains -> create mapping -> check if needed -> standardize.

        Args:
            input_pdb: Path to input PDB file
            output_pdb: Path to output PDB file
            task_name: Optional task/sample name for result tracking
            skip_if_standard: If True, skip processing if already in standard order

        Returns:
            StandardizationResult object with processing details

        Example:
            >>> result = standardizer.process_single(
            ...     "1ao7.pdb",
            ...     "1ao7_standardized.pdb",
            ...     task_name="1ao7_run1"
            ... )
            >>> print(f"Status: {result.status}")
        """
        input_pdb = Path(input_pdb)
        output_pdb = Path(output_pdb)
        task_name = task_name or input_pdb.stem

        start_time = datetime.now()

        result = StandardizationResult(
            task_name=task_name,
            status="UNKNOWN",
            num_chains=0,
            chain_mapping={},
            input_file=input_pdb,
            output_file=output_pdb
        )

        try:
            # Check if input exists
            if not input_pdb.exists():
                result.status = "NO_PDB"
                result.error_message = f"PDB file not found: {input_pdb}"
                self.logger.error(result.error_message)
                return result

            # Analyze chains
            chain_list = self.analyze_chains(input_pdb)
            result.num_chains = len(chain_list)

            self.logger.info(f"\nTask: {task_name}")
            self.logger.info(f"  Found {len(chain_list)} protein chains:")
            for i, chain_info in enumerate(chain_list, 1):
                self.logger.info(
                    f"    {i}. Chain {chain_info.chain_id}: "
                    f"{chain_info.residue_count:3d} residues, "
                    f"{chain_info.atom_count:4d} atoms"
                )

            # Create mapping (intelligent or length-based)
            chain_mapping, status = self.create_mapping(input_pdb, chain_list)
            result.chain_mapping = chain_mapping
            result.status = status

            if status == "OK":
                # Display mapping
                mapping_str = " / ".join([f"{old}→{new}" for old, new in chain_mapping.items()])
                self.logger.info(f"  Chain mapping: {mapping_str}")

                # Check if already standard
                if skip_if_standard and self.is_already_standard(chain_list, chain_mapping):
                    self.logger.info("  ✓ Already in standard order, skipping")
                    result.status = "ALREADY_STANDARD"
                else:
                    # Standardize
                    self.standardize_file(input_pdb, output_pdb, chain_mapping)
                    self.logger.info(f"  ✓ Standardization complete: {output_pdb}")

            elif status == "MULTICHAIN":
                result.error_message = f"{len(chain_list)} chains, exceeds standard {self.expected_chain_count}"
                self.logger.warning(f"  ⚠️  {result.error_message}")

            elif status == "INSUFFICIENT":
                result.error_message = f"{len(chain_list)} chains, less than standard {self.expected_chain_count}"
                self.logger.warning(f"  ⚠️  {result.error_message}")

        except Exception as e:
            result.status = "ERROR"
            result.error_message = str(e)
            self.logger.error(f"  ✗ Processing failed: {e}")

        finally:
            end_time = datetime.now()
            result.processing_time = (end_time - start_time).total_seconds()

        return result

    def batch_process(
        self,
        input_output_pairs: List[Tuple[Path, Path, str]],
        n_processes: int = 4,
        skip_if_standard: bool = True
    ) -> List[StandardizationResult]:
        """
        Batch process multiple PDB files in parallel.

        Args:
            input_output_pairs: List of (input_pdb, output_pdb, task_name) tuples
            n_processes: Number of parallel processes (default: 4)
            skip_if_standard: Skip processing if already in standard order

        Returns:
            List of StandardizationResult objects

        Example:
            >>> pairs = [
            ...     (Path("1ao7.pdb"), Path("1ao7_std.pdb"), "1ao7"),
            ...     (Path("2vlk.pdb"), Path("2vlk_std.pdb"), "2vlk"),
            ... ]
            >>> results = standardizer.batch_process(pairs, n_processes=2)
            >>> successful = [r for r in results if r.status == "OK"]
        """
        self.logger.info(f"{'='*80}")
        self.logger.info("PDB Chain Standardization - Batch Processing")
        self.logger.info(f"{'='*80}")
        self.logger.info(f"Total tasks: {len(input_output_pairs)}")
        self.logger.info(f"Parallel processes: {n_processes}")
        self.logger.info(f"Standard order: {' -> '.join(self.standard_order)}")
        self.logger.info(f"{'='*80}\n")

        start_time = datetime.now()

        # Parallel processing
        with Pool(processes=n_processes) as pool:
            results = pool.starmap(
                self._batch_worker,
                [(inp, out, name, skip_if_standard)
                 for inp, out, name in input_output_pairs]
            )

        end_time = datetime.now()
        total_time = (end_time - start_time).total_seconds()

        # Statistics
        self._print_batch_summary(results, total_time)

        return results

    def _batch_worker(
        self,
        input_pdb: Path,
        output_pdb: Path,
        task_name: str,
        skip_if_standard: bool
    ) -> StandardizationResult:
        """Worker function for batch processing (used by multiprocessing)."""
        return self.process_single(
            input_pdb=input_pdb,
            output_pdb=output_pdb,
            task_name=task_name,
            skip_if_standard=skip_if_standard
        )

    def _print_batch_summary(
        self,
        results: List[StandardizationResult],
        total_time: float
    ) -> None:
        """Print summary statistics for batch processing."""
        self.logger.info(f"\n{'='*80}")
        self.logger.info("Processing Complete - Summary")
        self.logger.info(f"{'='*80}")

        ok = [r for r in results if r.status == "OK"]
        already_standard = [r for r in results if r.status == "ALREADY_STANDARD"]
        multichain = [r for r in results if r.status == "MULTICHAIN"]
        insufficient = [r for r in results if r.status == "INSUFFICIENT"]
        error = [r for r in results if r.status in ["ERROR", "NO_PDB"]]

        self.logger.info(f"\nTotal tasks: {len(results)}")
        self.logger.info(f"✓ Standardized: {len(ok)}")
        self.logger.info(f"✓ Already standard: {len(already_standard)}")
        self.logger.info(f"⚠️  Multichain: {len(multichain)}")
        self.logger.info(f"⚠️  Insufficient chains: {len(insufficient)}")
        self.logger.info(f"✗ Failed/Missing: {len(error)}")
        self.logger.info(f"Total time: {total_time:.1f}s")

        if multichain:
            self.logger.info(f"\n⚠️  Multichain tasks ({len(multichain)}):")
            for r in multichain[:5]:
                self.logger.info(f"  - {r.task_name}: {r.num_chains} chains")
            if len(multichain) > 5:
                self.logger.info(f"  ... and {len(multichain)-5} more")

        if error:
            self.logger.info(f"\n✗ Failed tasks ({len(error)}):")
            for r in error[:5]:
                self.logger.info(f"  - {r.task_name}: {r.error_message}")
            if len(error) > 5:
                self.logger.info(f"  ... and {len(error)-5} more")

    def save_report(
        self,
        results: List[StandardizationResult],
        output_csv: Union[str, Path]
    ) -> None:
        """
        Save processing results to CSV report.

        Args:
            results: List of StandardizationResult objects
            output_csv: Path to output CSV file

        Example:
            >>> standardizer.save_report(results, "standardization_report.csv")
        """
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)

        fieldnames = ["TaskName", "Status", "NumChains", "ChainMapping",
                      "ProcessingTime_sec", "Error"]

        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

            for result in results:
                writer.writerow(result.to_dict())

        self.logger.info(f"\n✓ Report saved: {output_csv}")

    def validate_standardized_pdb(
        self,
        pdb_file: Union[str, Path]
    ) -> Tuple[bool, str]:
        """
        Validate that a PDB file follows the standard chain ordering.

        Args:
            pdb_file: Path to PDB file to validate

        Returns:
            Tuple of (is_valid, message)

        Example:
            >>> is_valid, msg = standardizer.validate_standardized_pdb("output.pdb")
            >>> print(f"Valid: {is_valid}, Message: {msg}")
        """
        try:
            chain_list = self.analyze_chains(pdb_file)

            if len(chain_list) != self.expected_chain_count:
                return False, f"Expected {self.expected_chain_count} chains, found {len(chain_list)}"

            # Check if chain IDs match standard order
            actual_chains = [chain.chain_id for chain in chain_list]
            if actual_chains != self.standard_order:
                return False, f"Chain order {actual_chains} does not match standard {self.standard_order}"

            return True, "Valid standard PDB"

        except Exception as e:
            return False, f"Validation error: {str(e)}"

    def _standardize_protein_only(
        self,
        input_pdb: Union[str, Path],
        output_pdb: Union[str, Path],
        chain_mapping: Dict[str, str]
    ) -> None:
        """
        Standardize PDB keeping only protein chains.

        Args:
            input_pdb: Input PDB file
            output_pdb: Output PDB file (protein only)
            chain_mapping: Dictionary mapping old chain IDs to new chain IDs
        """
        input_pdb = Path(input_pdb)
        output_pdb = Path(output_pdb)

        output_pdb.parent.mkdir(parents=True, exist_ok=True)

        with open(input_pdb) as f_in, open(output_pdb, 'w') as f_out:
            for line in f_in:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain_id = line[21]

                    # Only process protein chains
                    if chain_id in chain_mapping:
                        new_chain = chain_mapping[chain_id]
                        new_line = line[:21] + new_chain + line[22:]
                        f_out.write(new_line)
                    # Skip water, ions, and other non-protein atoms

                elif line.startswith(('CRYST1', 'MODEL', 'ENDMDL', 'END')):
                    # Preserve structural information
                    f_out.write(line)

    def _standardize_with_solvent(
        self,
        input_pdb: Union[str, Path],
        output_pdb: Union[str, Path],
        chain_mapping: Dict[str, str]
    ) -> None:
        """
        Standardize PDB keeping all atoms (protein + solvent + ions).

        Args:
            input_pdb: Input PDB file
            output_pdb: Output PDB file (full system)
            chain_mapping: Dictionary mapping old chain IDs to new chain IDs
        """
        input_pdb = Path(input_pdb)
        output_pdb = Path(output_pdb)

        output_pdb.parent.mkdir(parents=True, exist_ok=True)

        with open(input_pdb) as f_in, open(output_pdb, 'w') as f_out:
            for line in f_in:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain_id = line[21]

                    if chain_id in chain_mapping:
                        # Standardize protein chain IDs
                        new_chain = chain_mapping[chain_id]
                        new_line = line[:21] + new_chain + line[22:]
                        f_out.write(new_line)
                    else:
                        # Keep water, ions unchanged
                        f_out.write(line)
                else:
                    # Preserve all other records
                    f_out.write(line)

    def process_single_dual_output(
        self,
        input_pdb: Union[str, Path],
        output_prefix: Union[str, Path],
        task_name: Optional[str] = None,
        skip_if_standard: bool = True
    ) -> Dict[str, StandardizationResult]:
        """
        Process single PDB file and generate two standardized versions.

        Generates:
        1. {output_prefix}_protein.pdb - Protein only (standardized chains)
        2. {output_prefix}_full.pdb - Full system (standardized protein chains, preserved solvent)

        Args:
            input_pdb: Input PDB file path
            output_prefix: Output file prefix (e.g., "1ao7_std")
            task_name: Optional task name for tracking
            skip_if_standard: Skip processing if already standardized

        Returns:
            Dictionary with keys 'protein' and 'full', each containing StandardizationResult

        Example:
            >>> standardizer = PDBChainStandardizer()
            >>> results = standardizer.process_single_dual_output(
            ...     input_pdb="complex.pdb",
            ...     output_prefix="complex_std"
            ... )
            >>> print(f"Protein: {results['protein'].status}")
            >>> print(f"Full: {results['full'].status}")
            >>> print(f"Chain mapping: {results['protein'].chain_mapping}")
        """
        input_pdb = Path(input_pdb)
        output_prefix = Path(output_prefix)
        task_name = task_name or input_pdb.stem

        start_time = datetime.now()

        # Output files
        output_protein = output_prefix.parent / f"{output_prefix.name}_protein.pdb"
        output_full = output_prefix.parent / f"{output_prefix.name}_full.pdb"

        # Initialize results
        result_protein = StandardizationResult(
            task_name=f"{task_name}_protein",
            status="UNKNOWN",
            num_chains=0,
            chain_mapping={},
            input_file=input_pdb,
            output_file=output_protein
        )

        result_full = StandardizationResult(
            task_name=f"{task_name}_full",
            status="UNKNOWN",
            num_chains=0,
            chain_mapping={},
            input_file=input_pdb,
            output_file=output_full
        )

        try:
            # Check if input exists
            if not input_pdb.exists():
                error_msg = f"PDB file not found: {input_pdb}"
                result_protein.status = result_full.status = "NO_PDB"
                result_protein.error_message = result_full.error_message = error_msg
                self.logger.error(error_msg)
                return {'protein': result_protein, 'full': result_full}

            # Analyze chains
            chain_list = self.analyze_chains(input_pdb)
            num_chains = len(chain_list)

            result_protein.num_chains = result_full.num_chains = num_chains

            # Check chain count
            if num_chains != self.expected_chain_count:
                self.logger.warning(
                    f"Expected {self.expected_chain_count} chains, found {num_chains}"
                )

            # Create chain mapping
            chain_mapping, method = self.create_mapping(input_pdb, chain_list)

            result_protein.chain_mapping = result_full.chain_mapping = chain_mapping

            # Check if already standardized
            if skip_if_standard and self._is_already_standard(chain_mapping):
                result_protein.status = result_full.status = "ALREADY_STANDARD"
                self.logger.info(f"PDB already in standard format: {input_pdb.name}")
                return {'protein': result_protein, 'full': result_full}

            # Generate protein-only version
            self._standardize_protein_only(input_pdb, output_protein, chain_mapping)
            result_protein.status = "OK"
            self.logger.info(f"Generated protein-only PDB: {output_protein.name}")

            # Generate full version
            self._standardize_with_solvent(input_pdb, output_full, chain_mapping)
            result_full.status = "OK"
            self.logger.info(f"Generated full system PDB: {output_full.name}")

            # Record processing time
            processing_time = (datetime.now() - start_time).total_seconds()
            result_protein.processing_time = result_full.processing_time = processing_time

            self.logger.info(
                f"\n✓ Dual output completed for {task_name}:\n"
                f"  - Protein-only: {output_protein}\n"
                f"  - Full system: {output_full}\n"
                f"  - Chain mapping: {chain_mapping}\n"
                f"  - Processing time: {processing_time:.2f}s"
            )

        except Exception as e:
            error_msg = f"Failed to process {input_pdb.name}: {str(e)}"
            result_protein.status = result_full.status = "ERROR"
            result_protein.error_message = result_full.error_message = error_msg
            self.logger.error(error_msg)

        return {'protein': result_protein, 'full': result_full}

    def _is_already_standard(self, chain_mapping: Dict[str, str]) -> bool:
        """Check if chain mapping indicates already standardized chains."""
        return all(old == new for old, new in chain_mapping.items())
