"""
MD Completeness Checker

Validates MD simulation completeness and output file integrity.
"""

import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import logging

logger = logging.getLogger(__name__)


class MDCompletenessChecker:
    """
    Checks MD simulation completeness and validates output files.
    
    Validates:
    - Presence of required output files (md.gro, trajectory, log files)
    - File sizes and basic integrity
    - Simulation completion status from log files
    - Simulation time vs expected duration
    """
    
    def __init__(self, 
                 min_trajectory_size_mb: float = 1.0,
                 min_simulation_time_ps: float = 1000.0,
                 required_files: Optional[List[str]] = None):
        """
        Initialize MD completeness checker.
        
        Args:
            min_trajectory_size_mb: Minimum expected trajectory file size in MB
            min_simulation_time_ps: Minimum expected simulation time in ps
            required_files: List of required files (default: md.gro, md.xtc, md.log)
        """
        self.min_trajectory_size = min_trajectory_size_mb * 1024 * 1024  # Convert to bytes
        self.min_simulation_time = min_simulation_time_ps
        
        if required_files is None:
            self.required_files = ['md.gro', 'md.xtc', 'md.log']
        else:
            self.required_files = required_files
            
        # Alternative file patterns
        self.trajectory_patterns = ['*.xtc', '*.trr']
        self.log_patterns = ['*.log', '*.out']
        self.structure_patterns = ['md.gro', 'prod.gro', 'production.gro']
        
    def check_directory(self, md_dir: Union[str, Path]) -> Dict:
        """
        Check MD simulation completeness for a single directory.
        
        Args:
            md_dir: MD simulation directory path
            
        Returns:
            Dictionary with completeness analysis results
        """
        md_path = Path(md_dir)
        
        if not md_path.exists() or not md_path.is_dir():
            return {
                'status': 'error',
                'message': f"Directory does not exist: {md_path}",
                'directory': str(md_path),
                'files_found': {},
                'issues': ['directory_not_found']
            }
        
        logger.info(f"Checking MD completeness for: {md_path}")
        
        # Check for required files
        files_found = self._find_required_files(md_path)
        
        # Analyze file integrity
        file_analysis = self._analyze_files(md_path, files_found)
        
        # Check simulation completion from log
        simulation_status = self._check_simulation_completion(md_path, files_found.get('log'))
        
        # Overall assessment
        issues = []
        status = self._assess_overall_status(files_found, file_analysis, simulation_status, issues)
        
        return {
            'status': status,
            'directory': str(md_path),
            'files_found': files_found,
            'file_analysis': file_analysis,
            'simulation_status': simulation_status,
            'issues': issues,
            'completeness_score': self._calculate_completeness_score(files_found, file_analysis, simulation_status)
        }
    
    def check_single_md(self, md_dir: Union[str, Path]) -> Dict:
        """
        Alias for check_directory method for backward compatibility.
        
        Args:
            md_dir: MD simulation directory path
            
        Returns:
            Dictionary with completeness analysis results
        """
        return self.check_directory(md_dir)
    
    def _find_required_files(self, md_path: Path) -> Dict[str, Optional[str]]:
        """Find required MD output files in the directory."""
        files_found = {}
        
        # Look for structure file (md.gro)
        structure_file = self._find_file_by_patterns(md_path, self.structure_patterns)
        files_found['structure'] = structure_file
        
        # Look for trajectory file
        trajectory_files = []
        for pattern in self.trajectory_patterns:
            trajectory_files.extend(list(md_path.glob(pattern)))
        
        if trajectory_files:
            # Choose the largest trajectory file (likely the main one)
            largest_traj = max(trajectory_files, key=lambda f: f.stat().st_size if f.exists() else 0)
            files_found['trajectory'] = str(largest_traj)
        else:
            files_found['trajectory'] = None
            
        # Look for log file
        log_files = []
        for pattern in self.log_patterns:
            log_files.extend(list(md_path.glob(pattern)))
            
        if log_files:
            # Choose the most recent log file
            latest_log = max(log_files, key=lambda f: f.stat().st_mtime if f.exists() else 0)
            files_found['log'] = str(latest_log)
        else:
            files_found['log'] = None
            
        # Additional files
        files_found['tpr'] = self._find_file_by_patterns(md_path, ['*.tpr', 'md.tpr'])
        files_found['edr'] = self._find_file_by_patterns(md_path, ['*.edr', 'md.edr'])
        
        return files_found
    
    def _find_file_by_patterns(self, directory: Path, patterns: List[str]) -> Optional[str]:
        """Find first matching file by patterns."""
        for pattern in patterns:
            matches = list(directory.glob(pattern))
            if matches:
                return str(matches[0])
        return None
    
    def _analyze_files(self, md_path: Path, files_found: Dict) -> Dict:
        """Analyze file sizes and basic integrity."""
        analysis = {}
        
        for file_type, file_path in files_found.items():
            if file_path and Path(file_path).exists():
                file_stat = Path(file_path).stat()
                analysis[file_type] = {
                    'path': file_path,
                    'size_bytes': file_stat.st_size,
                    'size_mb': file_stat.st_size / (1024 * 1024),
                    'exists': True,
                    'readable': os.access(file_path, os.R_OK)
                }
                
                # Specific checks for trajectory files
                if file_type == 'trajectory':
                    analysis[file_type]['size_adequate'] = file_stat.st_size >= self.min_trajectory_size
                    
            else:
                analysis[file_type] = {
                    'path': None,
                    'exists': False,
                    'size_bytes': 0,
                    'size_mb': 0.0,
                    'readable': False
                }
                
        return analysis
    
    def _check_simulation_completion(self, md_path: Path, log_file: Optional[str]) -> Dict:
        """Check simulation completion status from log file."""
        if not log_file or not Path(log_file).exists():
            return {
                'completed': False,
                'simulation_time_ps': 0.0,
                'time_adequate': False,
                'completion_message': 'No log file found'
            }
        
        try:
            with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
                log_content = f.read()
                
            # Look for completion indicators
            completion_patterns = [
                r'Finished mdrun',
                r'simulation completed',
                r'Thanks for using GROMACS',
                r'Performance:'
            ]
            
            completed = any(re.search(pattern, log_content, re.IGNORECASE) for pattern in completion_patterns)
            
            # Extract simulation time
            simulation_time = self._extract_simulation_time(log_content)
            
            return {
                'completed': completed,
                'simulation_time_ps': simulation_time,
                'time_adequate': simulation_time >= self.min_simulation_time,
                'completion_message': 'Simulation completed successfully' if completed else 'Simulation may be incomplete',
                'log_file': log_file
            }
            
        except Exception as e:
            logger.warning(f"Error reading log file {log_file}: {e}")
            return {
                'completed': False,
                'simulation_time_ps': 0.0,
                'time_adequate': False,
                'completion_message': f'Error reading log file: {e}'
            }
    
    def _extract_simulation_time(self, log_content: str) -> float:
        """Extract total simulation time from log content."""
        time_patterns = [
            r'Statistics over (\d+\.?\d*) ps',
            r'Total simulation time:\s*(\d+\.?\d*)\s*ps',
            r'Time:\s*(\d+\.?\d*)\s*ps'
        ]
        
        for pattern in time_patterns:
            match = re.search(pattern, log_content)
            if match:
                try:
                    return float(match.group(1))
                except ValueError:
                    continue
                    
        return 0.0
    
    def _assess_overall_status(self, files_found: Dict, file_analysis: Dict, 
                             simulation_status: Dict, issues: List[str]) -> str:
        """Assess overall completion status using simplified logic."""
        # Simplified logic: only check for md.gro and md.xtc existence
        has_gro = files_found.get('structure') is not None
        has_xtc = files_found.get('trajectory') is not None
        
        if has_gro and has_xtc:
            # Both files exist - simulation is complete
            return 'complete'
        elif has_xtc and not has_gro:
            # Has trajectory but no final structure - simulation didn't finish properly
            issues.append('missing_final_structure')
            return 'incomplete'
        else:
            # Missing trajectory file - simulation incomplete or failed
            if not has_xtc:
                issues.append('missing_trajectory_file')
            if not has_gro:
                issues.append('missing_structure_file')
            return 'incomplete'
    
    def _calculate_completeness_score(self, files_found: Dict, file_analysis: Dict, 
                                    simulation_status: Dict) -> float:
        """Calculate a completeness score (0-100) using simplified logic."""
        # Simplified scoring: only based on essential files
        has_gro = files_found.get('structure') is not None
        has_xtc = files_found.get('trajectory') is not None
        
        if has_gro and has_xtc:
            # Both essential files present - 100% complete
            return 100.0
        elif has_xtc and not has_gro:
            # Has trajectory but missing final structure - 50% 
            return 50.0
        else:
            # Missing trajectory or both files - 0%
            return 0.0
    
    def batch_check(self, root_directory: Union[str, Path], 
                   pattern: str = "*") -> List[Dict]:
        """
        Check multiple MD directories in batch.
        
        Args:
            root_directory: Root directory containing MD simulation directories
            pattern: Pattern to match MD directories (default: all subdirectories)
            
        Returns:
            List of completeness analysis results
        """
        root_path = Path(root_directory)
        
        if not root_path.exists():
            logger.error(f"Root directory does not exist: {root_path}")
            return []
            
        # Find MD directories
        md_directories = [d for d in root_path.glob(pattern) if d.is_dir()]
        
        logger.info(f"Found {len(md_directories)} directories to check")
        
        results = []
        for md_dir in md_directories:
            try:
                result = self.check_directory(md_dir)
                results.append(result)
            except Exception as e:
                logger.error(f"Error checking directory {md_dir}: {e}")
                results.append({
                    'status': 'error',
                    'directory': str(md_dir),
                    'message': str(e),
                    'issues': ['check_error']
                })
                
        return results
    
    def get_summary_statistics(self, results: List[Dict]) -> Dict:
        """Generate summary statistics from batch check results."""
        if not results:
            return {}
            
        total = len(results)
        complete = sum(1 for r in results if r.get('status') == 'complete')
        partial = sum(1 for r in results if r.get('status') == 'partial')
        incomplete = sum(1 for r in results if r.get('status') == 'incomplete')
        errors = sum(1 for r in results if r.get('status') == 'error')
        
        avg_score = sum(r.get('completeness_score', 0) for r in results) / total
        
        # Collect all issues
        all_issues = []
        for result in results:
            all_issues.extend(result.get('issues', []))
            
        issue_counts = {}
        for issue in all_issues:
            issue_counts[issue] = issue_counts.get(issue, 0) + 1
            
        return {
            'total_directories': total,
            'complete': complete,
            'partial': partial,
            'incomplete': incomplete,
            'errors': errors,
            'success_rate': (complete / total) * 100 if total > 0 else 0,
            'average_completeness_score': avg_score,
            'common_issues': sorted(issue_counts.items(), key=lambda x: x[1], reverse=True)
        }