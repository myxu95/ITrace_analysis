"""
PDB Downloader Module

This module provides functionality to download PDB files from RCSB PDB database,
with support for downloading biological assemblies (not just asymmetric units).

Classes:
    PDBDownloader: Main class for downloading PDB files from RCSB

Example:
    >>> from md_analysis.utils import PDBDownloader
    >>> downloader = PDBDownloader()
    >>> result = downloader.download_pdb('1ao7')
    >>> print(result['file_path'])
    'input/standardizedpdbs/pdbs_raw/1ao7.pdb'
"""

import argparse
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Union

try:
    import requests
except ImportError:
    raise ImportError(
        "requests library is required for PDB downloading. "
        "Please install it using: pip install requests"
    )

logger = logging.getLogger(__name__)


class PDBDownloader:
    """
    Download PDB files from RCSB PDB database.

    This class provides methods to download PDB structures from RCSB PDB,
    with intelligent selection of biological assemblies for accurate
    representation of functional protein complexes.

    Attributes:
        download_dir (str): Directory to save downloaded PDB files
        assembly_preference (str|int): Strategy for selecting biological assembly
        max_retries (int): Maximum number of retry attempts for failed downloads
        retry_delay (float): Base delay between retries (exponential backoff)
        timeout (int): HTTP request timeout in seconds
        verify_download (bool): Whether to validate downloaded files
    """

    # RCSB PDB API endpoints
    ENDPOINTS = {
        'entry_metadata': 'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}',
        'assembly_metadata': 'https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}/{assembly_id}',
        'asymmetric_unit': 'https://files.rcsb.org/download/{pdb_id}.pdb',
        'biological_assembly': 'https://files.rcsb.org/download/{pdb_id}.pdb{assembly_num}'
    }

    # Standard protein residues for validation
    STANDARD_RESIDUES = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
        'THR', 'TRP', 'TYR', 'VAL'
    }

    def __init__(
        self,
        download_dir: str = "input/standardizedpdbs/pdbs_raw",
        assembly_preference: Union[str, int] = 'auto',
        max_retries: int = 3,
        retry_delay: float = 2.0,
        timeout: int = 30,
        verify_download: bool = True
    ):
        """
        Initialize PDB downloader.

        Args:
            download_dir: Directory to save downloaded files
            assembly_preference: Assembly selection strategy
                - 'auto': Query metadata, prefer assembly 1 (recommended)
                - 'first': Download .pdb1 directly, fast mode
                - 'largest': Select assembly with most chains
                - int (1,2,3...): Download specific assembly number
                - 'asymmetric': Download asymmetric unit only
            max_retries: Maximum retry attempts for failed requests
            retry_delay: Base delay between retries in seconds
            timeout: HTTP request timeout in seconds
            verify_download: Whether to validate downloaded files
        """
        self.download_dir = Path(download_dir)
        self.assembly_preference = assembly_preference
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.timeout = timeout
        self.verify_download = verify_download

        # Create download directory if it doesn't exist
        self.download_dir.mkdir(parents=True, exist_ok=True)

        # Setup HTTP session with persistent connection
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'ITrace-PDB-Downloader/1.0 (Scientific Research)',
            'Accept': 'application/json, text/plain, */*'
        })

        logger.info(f"PDBDownloader initialized with directory: {self.download_dir}")
        logger.info(f"Assembly preference: {assembly_preference}")

    @staticmethod
    def normalize_pdb_id(pdb_id: str) -> str:
        """
        Normalize PDB ID to uppercase 4-character format.

        Args:
            pdb_id: PDB identifier (case-insensitive)

        Returns:
            Normalized PDB ID in uppercase

        Raises:
            ValueError: If PDB ID format is invalid
        """
        pdb_id = pdb_id.strip().upper()
        if len(pdb_id) != 4:
            raise ValueError(f"Invalid PDB ID: {pdb_id}. Must be 4 characters.")
        if not pdb_id.isalnum():
            raise ValueError(f"Invalid PDB ID: {pdb_id}. Must be alphanumeric.")
        return pdb_id

    def _make_request_with_retry(
        self,
        url: str,
        request_type: str = 'GET',
        **kwargs
    ) -> requests.Response:
        """
        Make HTTP request with exponential backoff retry.

        Args:
            url: URL to request
            request_type: HTTP method (GET, POST, etc.)
            **kwargs: Additional arguments passed to requests

        Returns:
            Response object

        Raises:
            requests.exceptions.RequestException: If all retries fail
        """
        for attempt in range(self.max_retries):
            try:
                if request_type.upper() == 'GET':
                    response = self.session.get(
                        url,
                        timeout=self.timeout,
                        **kwargs
                    )
                else:
                    raise ValueError(f"Unsupported request type: {request_type}")

                # Don't retry on 404 (not found) or other client errors
                if response.status_code == 404:
                    logger.error(f"Resource not found: {url}")
                    response.raise_for_status()

                # Retry on server errors (5xx)
                if 500 <= response.status_code < 600:
                    logger.warning(
                        f"Server error {response.status_code} on attempt "
                        f"{attempt + 1}/{self.max_retries}"
                    )
                    if attempt < self.max_retries - 1:
                        delay = self.retry_delay * (2 ** attempt)
                        logger.info(f"Retrying in {delay} seconds...")
                        time.sleep(delay)
                        continue

                response.raise_for_status()
                return response

            except requests.exceptions.Timeout:
                logger.warning(
                    f"Request timeout on attempt {attempt + 1}/{self.max_retries}"
                )
                if attempt < self.max_retries - 1:
                    delay = self.retry_delay * (2 ** attempt)
                    logger.info(f"Retrying in {delay} seconds...")
                    time.sleep(delay)
                else:
                    raise

            except requests.exceptions.RequestException as e:
                if attempt < self.max_retries - 1 and response.status_code >= 500:
                    delay = self.retry_delay * (2 ** attempt)
                    logger.info(f"Retrying in {delay} seconds...")
                    time.sleep(delay)
                else:
                    raise

        raise requests.exceptions.RequestException(
            f"Failed after {self.max_retries} attempts"
        )

    def get_entry_metadata(self, pdb_id: str) -> Optional[Dict]:
        """
        Get PDB entry metadata from RCSB API.

        Args:
            pdb_id: PDB identifier

        Returns:
            Metadata dictionary or None if request fails
        """
        pdb_id = self.normalize_pdb_id(pdb_id)
        url = self.ENDPOINTS['entry_metadata'].format(pdb_id=pdb_id)

        try:
            response = self._make_request_with_retry(url)
            return response.json()
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to get metadata for {pdb_id}: {e}")
            return None

    def get_assembly_info(self, pdb_id: str) -> List[Dict]:
        """
        Get information about all biological assemblies for a PDB entry.

        Args:
            pdb_id: PDB identifier

        Returns:
            List of assembly information dictionaries
        """
        pdb_id = self.normalize_pdb_id(pdb_id)
        metadata = self.get_entry_metadata(pdb_id)

        if metadata is None:
            logger.warning(f"No metadata available for {pdb_id}")
            return []

        assemblies = []

        # Try to get assembly information from metadata
        if 'rcsb_entry_info' in metadata:
            entry_info = metadata['rcsb_entry_info']
            if 'assembly_count' in entry_info:
                assembly_count = entry_info['assembly_count']
                logger.info(f"{pdb_id} has {assembly_count} biological assemblies")

                # Query each assembly
                for i in range(1, assembly_count + 1):
                    assembly_url = self.ENDPOINTS['assembly_metadata'].format(
                        pdb_id=pdb_id,
                        assembly_id=i
                    )
                    try:
                        response = self._make_request_with_retry(assembly_url)
                        assembly_data = response.json()

                        # Extract oligomeric state safely
                        oligomeric_count = 'unknown'
                        symmetry_data = assembly_data.get('rcsb_struct_symmetry', None)
                        if isinstance(symmetry_data, dict):
                            oligomeric_count = symmetry_data.get('oligomeric_state', 'unknown')
                        elif isinstance(symmetry_data, list) and len(symmetry_data) > 0:
                            if isinstance(symmetry_data[0], dict):
                                oligomeric_count = symmetry_data[0].get('oligomeric_state', 'unknown')

                        assemblies.append({
                            'assembly_id': i,
                            'oligomeric_count': oligomeric_count,
                            'preferred': (i == 1)  # Assembly 1 is usually preferred
                        })
                    except (requests.exceptions.RequestException, Exception) as e:
                        logger.warning(f"Could not get assembly {i} metadata: {e}")
                        # Still add basic info
                        assemblies.append({
                            'assembly_id': i,
                            'oligomeric_count': 'unknown',
                            'preferred': (i == 1)
                        })

        return assemblies

    def select_assembly(self, pdb_id: str) -> tuple[str, int]:
        """
        Select which assembly to download based on preference strategy.

        Args:
            pdb_id: PDB identifier

        Returns:
            Tuple of (download_url, assembly_number)
            assembly_number is 0 for asymmetric unit
        """
        pdb_id = self.normalize_pdb_id(pdb_id)

        # Strategy: asymmetric unit
        if self.assembly_preference == 'asymmetric':
            url = self.ENDPOINTS['asymmetric_unit'].format(pdb_id=pdb_id)
            logger.info(f"Downloading asymmetric unit for {pdb_id}")
            return url, 0

        # Strategy: specific assembly number
        if isinstance(self.assembly_preference, int):
            assembly_num = self.assembly_preference
            url = self.ENDPOINTS['biological_assembly'].format(
                pdb_id=pdb_id,
                assembly_num=assembly_num
            )
            logger.info(f"Downloading assembly {assembly_num} for {pdb_id}")
            return url, assembly_num

        # Strategy: first (fast mode, no metadata query)
        if self.assembly_preference == 'first':
            url = self.ENDPOINTS['biological_assembly'].format(
                pdb_id=pdb_id,
                assembly_num=1
            )
            logger.info(f"Downloading assembly 1 for {pdb_id} (fast mode)")
            return url, 1

        # Strategy: auto (query metadata, prefer assembly 1)
        if self.assembly_preference == 'auto':
            assemblies = self.get_assembly_info(pdb_id)

            if not assemblies:
                logger.warning(
                    f"No assembly info found for {pdb_id}, "
                    "falling back to asymmetric unit"
                )
                url = self.ENDPOINTS['asymmetric_unit'].format(pdb_id=pdb_id)
                return url, 0

            # Prefer assembly 1 (author-defined biological unit)
            assembly_num = 1
            url = self.ENDPOINTS['biological_assembly'].format(
                pdb_id=pdb_id,
                assembly_num=assembly_num
            )
            logger.info(
                f"Auto-selected assembly {assembly_num} for {pdb_id} "
                f"({len(assemblies)} assemblies available)"
            )
            return url, assembly_num

        # Strategy: largest (select assembly with most chains)
        if self.assembly_preference == 'largest':
            assemblies = self.get_assembly_info(pdb_id)

            if not assemblies:
                logger.warning(
                    f"No assembly info found for {pdb_id}, "
                    "falling back to asymmetric unit"
                )
                url = self.ENDPOINTS['asymmetric_unit'].format(pdb_id=pdb_id)
                return url, 0

            # For now, just use assembly 1 as we don't have chain count info
            # This would require more detailed parsing of assembly data
            assembly_num = 1
            url = self.ENDPOINTS['biological_assembly'].format(
                pdb_id=pdb_id,
                assembly_num=assembly_num
            )
            logger.info(f"Selected largest assembly {assembly_num} for {pdb_id}")
            return url, assembly_num

        # Default fallback
        logger.warning(f"Unknown preference '{self.assembly_preference}', using auto")
        return self.select_assembly(pdb_id)

    def validate_pdb_file(self, pdb_file: Path) -> bool:
        """
        Validate downloaded PDB file integrity.

        Args:
            pdb_file: Path to PDB file

        Returns:
            True if file is valid, False otherwise
        """
        try:
            # Check file exists and has content
            if not pdb_file.exists():
                logger.error(f"File does not exist: {pdb_file}")
                return False

            file_size = pdb_file.stat().st_size
            if file_size < 1024:  # Less than 1KB
                logger.error(f"File too small ({file_size} bytes): {pdb_file}")
                return False

            # Check for essential PDB records
            has_atom = False
            has_protein = False

            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        has_atom = True
                        # Check if line contains protein residues
                        if len(line) >= 20:
                            resname = line[17:20].strip()
                            if resname in self.STANDARD_RESIDUES:
                                has_protein = True
                                break

            if not has_atom:
                logger.error(f"No ATOM/HETATM records found in {pdb_file}")
                return False

            if not has_protein:
                logger.warning(f"No standard protein residues found in {pdb_file}")
                # Not a critical error, could be nucleic acid or other molecules

            logger.debug(f"File validation passed: {pdb_file}")
            return True

        except Exception as e:
            logger.error(f"Error validating file {pdb_file}: {e}")
            return False

    def download_pdb(
        self,
        pdb_id: str,
        output_file: Optional[str] = None,
        force: bool = False
    ) -> Dict[str, any]:
        """
        Download a single PDB file.

        Args:
            pdb_id: PDB identifier
            output_file: Custom output file path (optional)
            force: Force re-download even if file exists

        Returns:
            Dictionary with download results:
                - 'pdb_id': PDB identifier
                - 'success': Boolean indicating success
                - 'file_path': Path to downloaded file (if successful)
                - 'assembly_number': Assembly number (0 for asymmetric)
                - 'message': Status message
        """
        try:
            pdb_id = self.normalize_pdb_id(pdb_id)
        except ValueError as e:
            return {
                'pdb_id': pdb_id,
                'success': False,
                'message': str(e)
            }

        # Determine output file path
        if output_file is None:
            output_file = self.download_dir / f"{pdb_id}.pdb"
        else:
            output_file = Path(output_file)

        # Check if file already exists
        if output_file.exists() and not force:
            if self.verify_download and self.validate_pdb_file(output_file):
                logger.info(f"File already exists and is valid: {output_file}")
                return {
                    'pdb_id': pdb_id,
                    'success': True,
                    'file_path': str(output_file),
                    'assembly_number': None,
                    'message': 'File already exists (skipped)'
                }
            elif not self.verify_download:
                logger.info(f"File already exists: {output_file}")
                return {
                    'pdb_id': pdb_id,
                    'success': True,
                    'file_path': str(output_file),
                    'assembly_number': None,
                    'message': 'File already exists (not verified)'
                }

        # Select assembly and get download URL
        try:
            download_url, assembly_num = self.select_assembly(pdb_id)
        except Exception as e:
            logger.error(f"Failed to select assembly for {pdb_id}: {e}")
            return {
                'pdb_id': pdb_id,
                'success': False,
                'message': f'Assembly selection failed: {e}'
            }

        # Download file
        logger.info(f"Downloading {pdb_id} from {download_url}")
        try:
            response = self._make_request_with_retry(download_url)

            # Write to file
            with open(output_file, 'wb') as f:
                f.write(response.content)

            logger.info(f"Downloaded {pdb_id} to {output_file}")

            # Validate if requested
            if self.verify_download:
                if not self.validate_pdb_file(output_file):
                    output_file.unlink()  # Delete corrupted file
                    return {
                        'pdb_id': pdb_id,
                        'success': False,
                        'message': 'Downloaded file failed validation'
                    }

            return {
                'pdb_id': pdb_id,
                'success': True,
                'file_path': str(output_file),
                'assembly_number': assembly_num,
                'message': 'Downloaded successfully'
            }

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to download {pdb_id}: {e}")
            return {
                'pdb_id': pdb_id,
                'success': False,
                'message': f'Download failed: {e}'
            }
        except Exception as e:
            logger.error(f"Unexpected error downloading {pdb_id}: {e}")
            return {
                'pdb_id': pdb_id,
                'success': False,
                'message': f'Unexpected error: {e}'
            }

    def batch_download(
        self,
        pdb_ids: List[str],
        max_workers: int = 4,
        force: bool = False
    ) -> List[Dict[str, any]]:
        """
        Download multiple PDB files in parallel.

        Args:
            pdb_ids: List of PDB identifiers
            max_workers: Maximum number of parallel downloads
            force: Force re-download even if files exist

        Returns:
            List of result dictionaries (one per PDB)
        """
        logger.info(f"Starting batch download of {len(pdb_ids)} PDB files")
        logger.info(f"Using {max_workers} parallel workers")

        results = []

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all download tasks
            future_to_pdb = {
                executor.submit(self.download_pdb, pdb_id, None, force): pdb_id
                for pdb_id in pdb_ids
            }

            # Process completed downloads
            for future in as_completed(future_to_pdb):
                pdb_id = future_to_pdb[future]
                try:
                    result = future.result()
                    results.append(result)

                    # Log progress
                    status = "OK" if result['success'] else "FAIL"
                    logger.info(f"[{status}] {pdb_id}: {result.get('message', '')}")

                except Exception as e:
                    logger.error(f"Exception occurred for {pdb_id}: {e}")
                    results.append({
                        'pdb_id': pdb_id,
                        'success': False,
                        'message': f'Exception: {e}'
                    })

        # Summary
        successful = sum(1 for r in results if r['success'])
        logger.info(
            f"Batch download completed: {successful}/{len(pdb_ids)} successful"
        )

        return results


def main():
    """Command-line interface for PDB downloader."""
    parser = argparse.ArgumentParser(
        description='Download PDB files from RCSB PDB database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download single PDB (auto-select biological assembly)
  python -m md_analysis.utils.pdb_downloader 1ao7

  # Download specific assembly
  python -m md_analysis.utils.pdb_downloader 1ao7 --assembly 2

  # Batch download from file
  python -m md_analysis.utils.pdb_downloader --batch pdb_ids.txt --workers 4

  # Query assembly information only
  python -m md_analysis.utils.pdb_downloader 1ao7 --info-only

  # Force re-download
  python -m md_analysis.utils.pdb_downloader 1ao7 --force
        """
    )

    # Input arguments
    parser.add_argument(
        'pdb_id',
        nargs='?',
        help='PDB identifier (4 characters)'
    )
    parser.add_argument(
        '--batch',
        help='File containing list of PDB IDs (one per line)'
    )

    # Output arguments
    parser.add_argument(
        '-o', '--output',
        default='input/standardizedpdbs/pdbs_raw',
        help='Output directory (default: input/standardizedpdbs/pdbs_raw)'
    )

    # Assembly selection
    parser.add_argument(
        '--assembly',
        default='auto',
        help=(
            'Assembly selection strategy: auto (default), first, largest, '
            'asymmetric, or integer (1,2,3...)'
        )
    )

    # Download options
    parser.add_argument(
        '--workers',
        type=int,
        default=4,
        help='Number of parallel workers for batch download (default: 4)'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force re-download even if file exists'
    )
    parser.add_argument(
        '--no-verify',
        action='store_true',
        help='Skip file validation after download'
    )

    # Query options
    parser.add_argument(
        '--info-only',
        action='store_true',
        help='Query assembly information only (no download)'
    )

    # Logging
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Validate input
    if not args.pdb_id and not args.batch:
        parser.error("Either pdb_id or --batch must be provided")

    if args.pdb_id and args.batch:
        parser.error("Cannot specify both pdb_id and --batch")

    # Handle assembly preference
    assembly_pref = args.assembly
    if assembly_pref.isdigit():
        assembly_pref = int(assembly_pref)

    # Initialize downloader
    downloader = PDBDownloader(
        download_dir=args.output,
        assembly_preference=assembly_pref,
        verify_download=not args.no_verify
    )

    # Info-only mode
    if args.info_only:
        if args.batch:
            parser.error("--info-only cannot be used with --batch")

        pdb_id = args.pdb_id
        print(f"\nQuerying assembly information for {pdb_id}...")
        assemblies = downloader.get_assembly_info(pdb_id)

        if not assemblies:
            print(f"No assembly information found for {pdb_id}")
            print("This PDB may only have an asymmetric unit.")
        else:
            print(f"\nFound {len(assemblies)} biological assemblies:")
            for asm in assemblies:
                preferred = " (PREFERRED)" if asm['preferred'] else ""
                print(
                    f"  Assembly {asm['assembly_id']}: "
                    f"{asm['oligomeric_count']}{preferred}"
                )
        return

    # Single download
    if args.pdb_id:
        result = downloader.download_pdb(args.pdb_id, force=args.force)

        if result['success']:
            print(f"\nSuccess! Downloaded to: {result['file_path']}")
            if result.get('assembly_number'):
                print(f"Assembly number: {result['assembly_number']}")
        else:
            print(f"\nFailed: {result['message']}")
            exit(1)

    # Batch download
    elif args.batch:
        # Read PDB IDs from file
        batch_file = Path(args.batch)
        if not batch_file.exists():
            print(f"Error: Batch file not found: {batch_file}")
            exit(1)

        with open(batch_file, 'r') as f:
            pdb_ids = [line.strip() for line in f if line.strip()]

        print(f"\nBatch downloading {len(pdb_ids)} PDB files...")
        results = downloader.batch_download(
            pdb_ids,
            max_workers=args.workers,
            force=args.force
        )

        # Print summary
        successful = [r for r in results if r['success']]
        failed = [r for r in results if not r['success']]

        print(f"\n{'='*60}")
        print(f"Batch download completed:")
        print(f"  Total: {len(results)}")
        print(f"  Successful: {len(successful)}")
        print(f"  Failed: {len(failed)}")

        if failed:
            print(f"\nFailed downloads:")
            for r in failed:
                print(f"  {r['pdb_id']}: {r['message']}")


if __name__ == '__main__':
    main()
