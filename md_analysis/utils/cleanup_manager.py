#!/usr/bin/env python3
"""
Cleanup Manager for ITrace

Comprehensive cleanup utilities for temporary files and intermediate results
generated during MD trajectory processing.
"""

import os
import glob
import logging
from pathlib import Path
from typing import List, Optional, Set
import fnmatch

logger = logging.getLogger(__name__)


class CleanupManager:
    """Manages cleanup of temporary and intermediate files."""
    
    def __init__(self, keep_temp_files: bool = False, verbose: bool = False):
        """
        Initialize cleanup manager.
        
        Args:
            keep_temp_files: If True, skip cleanup for debugging
            verbose: If True, provide detailed cleanup information
        """
        self.keep_temp_files = keep_temp_files
        self.verbose = verbose
        self.registered_files = []
        self.cleanup_patterns = {
            # GROMACS temporary and backup files
            'gromacs_temp': [
                '*.temp_*.xtc',
                '*.temp_*.trr', 
                '*.temp_*.gro',
                '#*.xtc.*#',
                '#*.trr.*#',
                '#*.gro.*#',
                '*.backup.*'
            ],
            # Index files
            'index_files': [
                'temp_*.ndx',
                'shortest_chain_*.ndx',
                'all_chains_*.ndx'
            ],
            # Python temporary files
            'python_temp': [
                '*.tmp',
                '.temp_*',
                '__pycache__',
                '*.pyc'
            ],
            # Other MD processing temp files
            'md_temp': [
                '*.log.backup',
                'step*.pdb',
                '*.debug'
            ]
        }
    
    def register_file(self, filepath: str) -> str:
        """Register a file for cleanup."""
        if filepath and filepath not in self.registered_files:
            self.registered_files.append(filepath)
            if self.verbose:
                logger.debug(f"Registered file for cleanup: {filepath}")
        return filepath
    
    def register_pattern(self, pattern: str, directory: str = "."):
        """Register all files matching a pattern for cleanup."""
        matches = glob.glob(os.path.join(directory, pattern))
        for match in matches:
            self.register_file(match)
    
    def cleanup_registered_files(self) -> dict:
        """Clean up all registered files."""
        if self.keep_temp_files:
            logger.info(f"Keeping {len(self.registered_files)} registered files for debugging")
            if self.verbose:
                for file_path in self.registered_files:
                    logger.info(f"  Keeping: {file_path}")
            return {"kept": len(self.registered_files), "removed": 0}
        
        removed_count = 0
        failed_count = 0
        
        for file_path in self.registered_files:
            try:
                if os.path.exists(file_path):
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                        removed_count += 1
                        if self.verbose:
                            logger.debug(f"Removed: {file_path}")
                    elif os.path.isdir(file_path):
                        import shutil
                        shutil.rmtree(file_path)
                        removed_count += 1
                        if self.verbose:
                            logger.debug(f"Removed directory: {file_path}")
            except Exception as e:
                failed_count += 1
                logger.warning(f"Failed to remove {file_path}: {e}")
        
        self.registered_files.clear()
        
        if removed_count > 0:
            logger.info(f"Cleaned up {removed_count} registered files")
        if failed_count > 0:
            logger.warning(f"Failed to clean up {failed_count} files")
        
        return {"removed": removed_count, "failed": failed_count}
    
    def cleanup_by_patterns(self, directory: str = ".", 
                          categories: Optional[List[str]] = None) -> dict:
        """
        Clean up files matching predefined patterns.
        
        Args:
            directory: Directory to search for files
            categories: List of pattern categories to clean (default: all)
        
        Returns:
            Dictionary with cleanup statistics
        """
        if self.keep_temp_files:
            logger.info("Skipping pattern-based cleanup (keeping temp files)")
            return {"kept": 0, "removed": 0}
        
        if categories is None:
            categories = list(self.cleanup_patterns.keys())
        
        directory_path = Path(directory)
        removed_count = 0
        failed_count = 0
        
        for category in categories:
            if category not in self.cleanup_patterns:
                logger.warning(f"Unknown cleanup category: {category}")
                continue
            
            patterns = self.cleanup_patterns[category]
            logger.debug(f"Cleaning category '{category}' with {len(patterns)} patterns")
            
            for pattern in patterns:
                try:
                    # Handle special case for __pycache__ directories
                    if pattern == '__pycache__':
                        for pycache_dir in directory_path.rglob('__pycache__'):
                            if pycache_dir.is_dir():
                                import shutil
                                shutil.rmtree(pycache_dir)
                                removed_count += 1
                                if self.verbose:
                                    logger.debug(f"Removed __pycache__: {pycache_dir}")
                        continue
                    
                    # Regular pattern matching
                    matches = list(directory_path.glob(pattern))
                    for match in matches:
                        try:
                            if match.is_file():
                                match.unlink()
                                removed_count += 1
                                if self.verbose:
                                    logger.debug(f"Removed: {match}")
                            elif match.is_dir():
                                import shutil
                                shutil.rmtree(match)
                                removed_count += 1
                                if self.verbose:
                                    logger.debug(f"Removed directory: {match}")
                        except Exception as e:
                            failed_count += 1
                            logger.warning(f"Failed to remove {match}: {e}")
                            
                except Exception as e:
                    logger.warning(f"Error processing pattern '{pattern}': {e}")
        
        if removed_count > 0:
            logger.info(f"Pattern cleanup removed {removed_count} files/directories")
        if failed_count > 0:
            logger.warning(f"Pattern cleanup failed on {failed_count} items")
        
        return {"removed": removed_count, "failed": failed_count}
    
    def comprehensive_cleanup(self, directory: str = ".") -> dict:
        """
        Perform comprehensive cleanup of both registered files and pattern matches.
        
        Args:
            directory: Directory to clean
            
        Returns:
            Dictionary with detailed cleanup statistics
        """
        logger.info("Starting comprehensive cleanup...")
        
        # Clean registered files first
        registered_stats = self.cleanup_registered_files()
        
        # Then clean by patterns
        pattern_stats = self.cleanup_by_patterns(directory)
        
        total_stats = {
            "registered_removed": registered_stats.get("removed", 0),
            "registered_failed": registered_stats.get("failed", 0),
            "pattern_removed": pattern_stats.get("removed", 0),
            "pattern_failed": pattern_stats.get("failed", 0),
            "total_removed": registered_stats.get("removed", 0) + pattern_stats.get("removed", 0),
            "total_failed": registered_stats.get("failed", 0) + pattern_stats.get("failed", 0)
        }
        
        logger.info(f"Comprehensive cleanup completed: "
                   f"{total_stats['total_removed']} removed, "
                   f"{total_stats['total_failed']} failed")
        
        return total_stats
    
    def get_cleanup_preview(self, directory: str = ".") -> dict:
        """
        Preview what files would be cleaned without actually removing them.
        
        Args:
            directory: Directory to preview
            
        Returns:
            Dictionary with lists of files that would be cleaned
        """
        directory_path = Path(directory)
        preview = {
            "registered_files": list(self.registered_files),
            "pattern_matches": {}
        }
        
        for category, patterns in self.cleanup_patterns.items():
            category_matches = []
            for pattern in patterns:
                if pattern == '__pycache__':
                    category_matches.extend([str(p) for p in directory_path.rglob('__pycache__')])
                else:
                    category_matches.extend([str(p) for p in directory_path.glob(pattern)])
            preview["pattern_matches"][category] = category_matches
        
        return preview


def create_cleanup_manager(keep_temp_files: bool = False, 
                         verbose: bool = False) -> CleanupManager:
    """Factory function to create a cleanup manager with standard configuration."""
    return CleanupManager(keep_temp_files=keep_temp_files, verbose=verbose)