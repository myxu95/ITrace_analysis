"""
Structure Validator

Validates PDB structure quality and analyzes chain composition.
"""

import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import logging

logger = logging.getLogger(__name__)

try:
    import MDAnalysis as mda
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False
    logger.warning("MDAnalysis not available. Structure validation will be limited.")


class StructureValidator:
    """
    Validates PDB structure quality and composition.
    
    Note: This validator is specifically designed for PDB files which contain 
    chain information. GRO files do not contain chain data and should be 
    validated using different criteria.
    
    Checks:
    - Number of protein chains (PDB only)
    - Atom coordinate ranges and outliers
    - Missing residues and atoms
    - Structure completeness
    - Complex composition analysis
    """
    
    def __init__(self, 
                 expected_chain_count: int = 5,
                 max_coordinate_range: float = 1000.0,
                 max_missing_residue_ratio: float = 0.05):
        """
        Initialize structure validator.
        
        Args:
            expected_chain_count: Expected number of protein chains
            max_coordinate_range: Maximum reasonable coordinate range (Angstroms)
            max_missing_residue_ratio: Maximum allowed missing residue ratio
        """
        self.expected_chain_count = expected_chain_count
        self.max_coordinate_range = max_coordinate_range
        self.max_missing_residue_ratio = max_missing_residue_ratio
        
    def validate_structure(self, structure_file: Union[str, Path]) -> Dict:
        """
        Validate a single structure file.
        
        Args:
            structure_file: Path to PDB structure file
            
        Returns:
            Dictionary with validation results
        """
        structure_path = Path(structure_file)
        
        if not structure_path.exists():
            return {
                'status': 'error',
                'message': f"Structure file does not exist: {structure_path}",
                'file_path': str(structure_path),
                'issues': ['file_not_found']
            }
        
        # Check file format - only process PDB files for chain analysis
        file_extension = structure_path.suffix.lower()
        if file_extension not in ['.pdb', '.ent']:
            return {
                'status': 'error',
                'message': f"Unsupported file format: {file_extension}. StructureValidator only supports PDB files (.pdb, .ent) for chain analysis.",
                'file_path': str(structure_path),
                'file_format': file_extension,
                'issues': ['unsupported_format'],
                'recommendation': 'Use PDB files for chain analysis. GRO files do not contain chain information.'
            }
            
        logger.info(f"Validating PDB structure: {structure_path}")
        
        try:
            if HAS_MDANALYSIS:
                return self._validate_with_mdanalysis(structure_path)
            else:
                return self._validate_without_mdanalysis(structure_path)
                
        except Exception as e:
            logger.error(f"Error validating structure {structure_path}: {e}")
            return {
                'status': 'error',
                'message': f"Validation error: {e}",
                'file_path': str(structure_path),
                'issues': ['validation_error']
            }
    
    def _validate_with_mdanalysis(self, structure_path: Path) -> Dict:
        """Validate structure using MDAnalysis."""
        try:
            # Load structure
            universe = mda.Universe(str(structure_path))
            
            # Basic structure info
            structure_info = {
                'total_atoms': len(universe.atoms),
                'total_residues': len(universe.residues),
                'n_frames': len(universe.trajectory)
            }
            
            # Chain analysis
            chain_analysis = self._analyze_chains_mda(universe)
            
            # Coordinate analysis
            coord_analysis = self._analyze_coordinates_mda(universe)
            
            # Protein analysis
            protein_analysis = self._analyze_protein_mda(universe)
            
            # Overall assessment
            issues = []
            status = self._assess_structure_status(chain_analysis, coord_analysis, protein_analysis, issues)
            
            return {
                'status': status,
                'file_path': str(structure_path),
                'structure_info': structure_info,
                'chain_analysis': chain_analysis,
                'coordinate_analysis': coord_analysis,
                'protein_analysis': protein_analysis,
                'issues': issues,
                'validation_score': self._calculate_validation_score(chain_analysis, coord_analysis, protein_analysis)
            }
            
        except Exception as e:
            raise Exception(f"MDAnalysis validation failed: {e}")
    
    def _validate_without_mdanalysis(self, structure_path: Path) -> Dict:
        """Basic validation without MDAnalysis."""
        logger.warning("Using basic validation without MDAnalysis")
        
        # Basic file analysis
        file_stat = structure_path.stat()
        
        # Try to parse basic info from file
        basic_info = self._parse_structure_file_basic(structure_path)
        
        issues = []
        if basic_info['estimated_chains'] != self.expected_chain_count:
            issues.append('unexpected_chain_count')
            
        if file_stat.st_size < 1000:  # Very small file
            issues.append('file_too_small')
            
        status = 'partial' if not issues else 'incomplete'
        
        return {
            'status': status,
            'file_path': str(structure_path),
            'structure_info': basic_info,
            'chain_analysis': {'method': 'basic_parsing'},
            'coordinate_analysis': {'method': 'not_available'},
            'protein_analysis': {'method': 'not_available'},
            'issues': issues,
            'validation_score': 50.0 if not issues else 25.0
        }
    
    def _analyze_chains_mda(self, universe) -> Dict:
        """Analyze protein chains using MDAnalysis."""
        try:
            # Get protein selection
            protein = universe.select_atoms("protein")
            
            if len(protein) == 0:
                return {
                    'protein_chains': 0,
                    'chain_details': [],
                    'has_protein': False,
                    'chain_count_normal': False
                }
            
            # Analyze chains - handle potential None/empty chain IDs
            try:
                chain_ids = set(protein.chainIDs)
                # Filter out None, empty strings, and whitespace-only chain IDs
                valid_chain_ids = {cid for cid in chain_ids if cid is not None and str(cid).strip()}
            except (AttributeError, ValueError):
                # Fallback: if chainIDs attribute is not available or invalid
                logger.warning("ChainIDs not available or invalid, attempting alternative chain detection")
                valid_chain_ids = set()
                try:
                    # Try to get unique segment IDs as alternative
                    seg_ids = set(protein.segids)
                    valid_chain_ids = {sid for sid in seg_ids if sid is not None and str(sid).strip()}
                except:
                    # If all else fails, assume single chain
                    valid_chain_ids = {'A'}
                    logger.warning("Using default chain ID 'A' as fallback")
            
            chain_details = []
            
            for chain_id in sorted(valid_chain_ids):
                try:
                    # Use proper chain selection syntax
                    if str(chain_id).strip():
                        chain_atoms = protein.select_atoms(f"chainid {chain_id}")
                    else:
                        # Handle empty/whitespace chain IDs
                        chain_atoms = protein
                        
                    chain_residues = set(chain_atoms.resids)
                    
                    chain_details.append({
                        'chain_id': str(chain_id),
                        'atom_count': len(chain_atoms),
                        'residue_count': len(chain_residues),
                        'residue_range': (min(chain_residues), max(chain_residues)) if chain_residues else (0, 0)
                    })
                except Exception as chain_e:
                    logger.warning(f"Failed to analyze chain {chain_id}: {chain_e}")
                    continue
            
            # If no valid chains found but protein exists, create a default entry
            if not chain_details and len(protein) > 0:
                protein_residues = set(protein.resids)
                chain_details.append({
                    'chain_id': 'unknown',
                    'atom_count': len(protein),
                    'residue_count': len(protein_residues),
                    'residue_range': (min(protein_residues), max(protein_residues)) if protein_residues else (0, 0)
                })
                valid_chain_ids = {'unknown'}
            
            return {
                'protein_chains': len(valid_chain_ids),
                'chain_details': chain_details,
                'has_protein': True,
                'chain_count_normal': len(valid_chain_ids) == self.expected_chain_count,
                'chain_ids': sorted(list(valid_chain_ids))
            }
            
        except Exception as e:
            logger.warning(f"Chain analysis failed: {e}")
            return {
                'protein_chains': 0,
                'chain_details': [],
                'has_protein': False,
                'chain_count_normal': False,
                'error': str(e)
            }
    
    def _analyze_coordinates_mda(self, universe) -> Dict:
        """Analyze coordinate quality using MDAnalysis."""
        try:
            atoms = universe.atoms
            positions = atoms.positions
            
            # Calculate coordinate ranges
            coord_ranges = {
                'x_range': (float(positions[:, 0].min()), float(positions[:, 0].max())),
                'y_range': (float(positions[:, 1].min()), float(positions[:, 1].max())),
                'z_range': (float(positions[:, 2].min()), float(positions[:, 2].max()))
            }
            
            # Calculate overall coordinate span
            coord_span = max(
                coord_ranges['x_range'][1] - coord_ranges['x_range'][0],
                coord_ranges['y_range'][1] - coord_ranges['y_range'][0],
                coord_ranges['z_range'][1] - coord_ranges['z_range'][0]
            )
            
            # Check for coordinate outliers
            coord_reasonable = coord_span <= self.max_coordinate_range
            
            # Center of mass
            protein = universe.select_atoms("protein")
            if len(protein) > 0:
                center_of_mass = protein.center_of_mass()
                com = [float(x) for x in center_of_mass]
            else:
                com = [0.0, 0.0, 0.0]
            
            return {
                'coordinate_ranges': coord_ranges,
                'coordinate_span': coord_span,
                'coordinates_reasonable': coord_reasonable,
                'center_of_mass': com,
                'total_atoms_with_coords': len(atoms)
            }
            
        except Exception as e:
            logger.warning(f"Coordinate analysis failed: {e}")
            return {
                'coordinate_ranges': {},
                'coordinate_span': 0.0,
                'coordinates_reasonable': False,
                'error': str(e)
            }
    
    def _analyze_protein_mda(self, universe) -> Dict:
        """Analyze protein-specific properties using MDAnalysis."""
        try:
            protein = universe.select_atoms("protein")
            
            if len(protein) == 0:
                return {
                    'has_protein': False,
                    'protein_atoms': 0,
                    'protein_residues': 0
                }
            
            # Basic protein stats
            protein_residues = len(set(protein.resids))
            protein_atoms = len(protein)
            
            # Check for common amino acids
            amino_acids = set(protein.resnames)
            standard_aa = {
                'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
                'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
            }
            
            standard_aa_found = amino_acids.intersection(standard_aa)
            non_standard = amino_acids - standard_aa
            
            # Secondary structure elements (if available)
            has_backbone = len(universe.select_atoms("protein and backbone")) > 0
            
            return {
                'has_protein': True,
                'protein_atoms': protein_atoms,
                'protein_residues': protein_residues,
                'amino_acids_found': len(amino_acids),
                'standard_amino_acids': len(standard_aa_found),
                'non_standard_residues': list(non_standard),
                'has_backbone': has_backbone,
                'protein_completeness_ratio': len(standard_aa_found) / len(standard_aa) if standard_aa else 0.0
            }
            
        except Exception as e:
            logger.warning(f"Protein analysis failed: {e}")
            return {
                'has_protein': False,
                'protein_atoms': 0,
                'protein_residues': 0,
                'error': str(e)
            }
    
    def _parse_structure_file_basic(self, structure_path: Path) -> Dict:
        """Basic PDB file parsing without MDAnalysis."""
        info = {
            'file_size_bytes': structure_path.stat().st_size,
            'estimated_atoms': 0,
            'estimated_chains': 0,
            'file_format': structure_path.suffix.lower()
        }
        
        try:
            with open(structure_path, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()
                
            # PDB format parsing
            atom_lines = [l for l in lines if l.startswith('ATOM') or l.startswith('HETATM')]
            chain_ids = set()
            protein_atoms = 0
            
            for line in atom_lines:
                if len(line) > 21:
                    chain_id = line[21].strip()
                    if chain_id and chain_id != ' ':
                        chain_ids.add(chain_id)
                    
                    # Check if it's a protein atom (simple heuristic)
                    if len(line) > 17:
                        res_name = line[17:20].strip()
                        standard_aa = {
                            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
                            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
                        }
                        if res_name in standard_aa:
                            protein_atoms += 1
                            
            info['estimated_atoms'] = len(atom_lines)
            info['estimated_chains'] = len(chain_ids)
            info['protein_atoms'] = protein_atoms
            info['chain_ids'] = sorted(list(chain_ids))
            
            # Additional PDB-specific info
            header_lines = [l for l in lines if l.startswith('HEADER')]
            if header_lines:
                header = header_lines[0]
                if len(header) > 62:
                    info['pdb_id'] = header[62:66].strip()
                    
        except Exception as e:
            logger.warning(f"Basic PDB parsing failed for {structure_path}: {e}")
            
        return info
    
    def _assess_structure_status(self, chain_analysis: Dict, coord_analysis: Dict, 
                               protein_analysis: Dict, issues: List[str]) -> str:
        """Assess overall structure validation status."""
        # Check chain count
        if not chain_analysis.get('chain_count_normal', False):
            issues.append('unusual_chain_count')
            
        # Check coordinates
        if not coord_analysis.get('coordinates_reasonable', True):
            issues.append('unreasonable_coordinates')
            
        # Check protein presence
        if not protein_analysis.get('has_protein', True):
            issues.append('no_protein_found')
            
        # Check protein completeness
        completeness = protein_analysis.get('protein_completeness_ratio', 1.0)
        if completeness < (1.0 - self.max_missing_residue_ratio):
            issues.append('significant_missing_residues')
            
        # Determine status
        if not issues:
            return 'valid'
        elif len(issues) <= 1:
            return 'acceptable'
        else:
            return 'problematic'
    
    def _calculate_validation_score(self, chain_analysis: Dict, coord_analysis: Dict, 
                                  protein_analysis: Dict) -> float:
        """Calculate validation score (0-100)."""
        score = 0.0
        
        # Chain count score (30 points)
        if chain_analysis.get('chain_count_normal', False):
            score += 30
        elif chain_analysis.get('protein_chains', 0) > 0:
            score += 15  # Partial credit for having some chains
            
        # Coordinate quality score (25 points)
        if coord_analysis.get('coordinates_reasonable', False):
            score += 25
            
        # Protein presence score (25 points)
        if protein_analysis.get('has_protein', False):
            score += 25
            
        # Protein completeness score (20 points)
        completeness = protein_analysis.get('protein_completeness_ratio', 0.0)
        score += completeness * 20
        
        return min(100.0, score)
    
    def batch_validate(self, structure_files: List[Union[str, Path]], 
                      pdb_id_only: bool = True) -> List[Dict]:
        """
        Validate multiple PDB structure files in batch.
        
        Args:
            structure_files: List of PDB structure file paths (.pdb, .ent)
            pdb_id_only: If True, only validate files with PDB ID naming convention
            
        Returns:
            List of validation results
            
        Note:
            Only PDB files are processed for chain analysis. GRO files will be 
            skipped with an appropriate error message since they don't contain 
            chain information.
        """
        results = []
        pdb_files = []
        skipped_files = []
        
        # Pre-filter to identify PDB files and apply naming convention filter
        for structure_file in structure_files:
            file_path = Path(structure_file)
            
            # Check file extension
            if file_path.suffix.lower() not in ['.pdb', '.ent']:
                skipped_files.append({
                    'file': structure_file,
                    'reason': 'non_pdb_format',
                    'suffix': file_path.suffix.lower()
                })
                continue
            
            # Check PDB ID naming convention if required
            if pdb_id_only and not self._is_pdb_id_filename(file_path.name):
                skipped_files.append({
                    'file': structure_file,
                    'reason': 'non_pdb_id_naming',
                    'filename': file_path.name
                })
                continue
                
            pdb_files.append(structure_file)
        
        # Log skipped files
        if skipped_files:
            non_pdb_count = sum(1 for f in skipped_files if f['reason'] == 'non_pdb_format')
            non_id_count = sum(1 for f in skipped_files if f['reason'] == 'non_pdb_id_naming')
            
            if non_pdb_count > 0:
                logger.warning(f"Skipping {non_pdb_count} non-PDB files")
                
            if non_id_count > 0:
                logger.info(f"Skipping {non_id_count} PDB files with non-PDB ID naming convention")
                if pdb_id_only:
                    logger.info("Note: Only processing files with PDB ID naming (e.g., 1abc.pdb, 2xyz_clean.pdb)")
        
        logger.info(f"Validating {len(pdb_files)} PDB structure files")
        
        for structure_file in pdb_files:
            try:
                result = self.validate_structure(structure_file)
                results.append(result)
            except Exception as e:
                logger.error(f"Error validating structure {structure_file}: {e}")
                results.append({
                    'status': 'error',
                    'file_path': str(structure_file),
                    'message': str(e),
                    'issues': ['validation_error']
                })
        
        # Add results for skipped files
        for skipped_info in skipped_files:
            skipped_file = skipped_info['file']
            reason = skipped_info['reason']
            
            if reason == 'non_pdb_format':
                message = 'File skipped - StructureValidator only processes PDB files for chain analysis'
                issues = ['unsupported_format']
            else:  # non_pdb_id_naming
                message = 'File skipped - Does not follow PDB ID naming convention'
                issues = ['non_pdb_id_naming']
            
            results.append({
                'status': 'skipped',
                'file_path': str(skipped_file),
                'file_format': Path(skipped_file).suffix.lower(),
                'message': message,
                'issues': issues,
                'skip_reason': reason
            })
                
        return results
    
    def _is_pdb_id_filename(self, filename: str) -> bool:
        """
        Check if filename follows PDB ID naming convention.
        
        PDB ID format: 4-character alphanumeric code starting with digit
        Valid patterns:
        - 1abc.pdb, 2XYZ.pdb (4 chars + .pdb)
        - 1abc_clean.pdb, 2XYZ_processed.pdb (4 chars + suffix + .pdb)
        """
        import re
        
        # Remove extension
        name_without_ext = filename.lower()
        if name_without_ext.endswith('.pdb'):
            name_without_ext = name_without_ext[:-4]
        elif name_without_ext.endswith('.ent'):
            name_without_ext = name_without_ext[:-4]
        else:
            return False
        
        # Check for PDB ID pattern: digit + 3 alphanumeric (standard PDB format)
        pdb_id_pattern = r'^[0-9][a-z0-9]{3}(_.*)?$'
        
        # Additional validation: ensure it's exactly 4 characters before any underscore
        base_part = name_without_ext.split('_')[0]
        if len(base_part) != 4:
            return False
            
        return bool(re.match(pdb_id_pattern, name_without_ext))
    
    def get_validation_summary(self, results: List[Dict]) -> Dict:
        """Generate summary statistics from validation results."""
        if not results:
            return {}
            
        total = len(results)
        valid = sum(1 for r in results if r.get('status') == 'valid')
        acceptable = sum(1 for r in results if r.get('status') == 'acceptable')
        problematic = sum(1 for r in results if r.get('status') == 'problematic')
        errors = sum(1 for r in results if r.get('status') == 'error')
        
        avg_score = sum(r.get('validation_score', 0) for r in results) / total
        
        # Chain count statistics
        chain_counts = []
        for result in results:
            chain_analysis = result.get('chain_analysis', {})
            if 'protein_chains' in chain_analysis:
                chain_counts.append(chain_analysis['protein_chains'])
                
        # Issue statistics
        all_issues = []
        for result in results:
            all_issues.extend(result.get('issues', []))
            
        issue_counts = {}
        for issue in all_issues:
            issue_counts[issue] = issue_counts.get(issue, 0) + 1
            
        return {
            'total_structures': total,
            'valid': valid,
            'acceptable': acceptable,
            'problematic': problematic,
            'errors': errors,
            'validation_rate': ((valid + acceptable) / total) * 100 if total > 0 else 0,
            'average_validation_score': avg_score,
            'chain_count_distribution': {
                'mean': sum(chain_counts) / len(chain_counts) if chain_counts else 0,
                'expected': self.expected_chain_count,
                'counts': chain_counts
            },
            'common_issues': sorted(issue_counts.items(), key=lambda x: x[1], reverse=True)
        }