"""
Batch Tracker

Tracks MD simulation batches and detects duplicates/missing entries.
"""

import json
import os
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union
import logging

logger = logging.getLogger(__name__)


class BatchTracker:
    """
    Tracks MD simulation batches and analyzes PDB processing statistics.
    
    Features:
    - PDB ID extraction from directory names
    - Duplicate detection across batches
    - Missing PDB identification
    - Batch statistics and reporting
    - Processing history tracking
    """
    
    def __init__(self, 
                 pdb_id_pattern: str = r'^([A-Za-z0-9]{4})',
                 tracking_db_path: Optional[Union[str, Path]] = None):
        """
        Initialize batch tracker.
        
        Args:
            pdb_id_pattern: Regex pattern to extract PDB ID from directory names
            tracking_db_path: Path to tracking database file (JSON)
        """
        self.pdb_id_pattern = re.compile(pdb_id_pattern)
        
        if tracking_db_path is None:
            self.tracking_db_path = Path.home() / ".immunoscope_batch_tracking.json"
        else:
            self.tracking_db_path = Path(tracking_db_path)
            
        self.tracking_data = self._load_tracking_data()
        
    def _load_tracking_data(self) -> Dict:
        """Load existing tracking data from database."""
        if self.tracking_db_path.exists():
            try:
                with open(self.tracking_db_path, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Failed to load tracking data: {e}")
                
        return {
            'batches': {},
            'pdb_entries': {},
            'last_updated': None
        }
    
    def _save_tracking_data(self):
        """Save tracking data to database."""
        self.tracking_data['last_updated'] = datetime.now().isoformat()
        
        try:
            with open(self.tracking_db_path, 'w') as f:
                json.dump(self.tracking_data, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save tracking data: {e}")
    
    def extract_pdb_id(self, directory_name: str) -> Optional[str]:
        """
        Extract PDB ID from directory name.
        
        Args:
            directory_name: Name of MD simulation directory
            
        Returns:
            PDB ID if found, None otherwise
        """
        match = self.pdb_id_pattern.match(directory_name)
        if match:
            return match.group(1).upper()
        return None
    
    def analyze_batch(self, batch_directory: Union[str, Path], 
                     batch_name: Optional[str] = None) -> Dict:
        """
        Analyze a single batch directory.
        
        Args:
            batch_directory: Path to batch directory containing MD simulations
            batch_name: Optional name for the batch (default: directory name)
            
        Returns:
            Dictionary with batch analysis results
        """
        batch_path = Path(batch_directory)
        
        if not batch_path.exists() or not batch_path.is_dir():
            raise ValueError(f"Batch directory does not exist: {batch_path}")
            
        if batch_name is None:
            batch_name = batch_path.name
            
        logger.info(f"Analyzing batch: {batch_name}")
        
        # Find all MD simulation directories
        md_directories = [d for d in batch_path.iterdir() if d.is_dir()]
        
        # Extract PDB information
        batch_entries = []
        pdb_counts = {}
        
        for md_dir in md_directories:
            pdb_id = self.extract_pdb_id(md_dir.name)
            
            entry = {
                'directory_name': md_dir.name,
                'directory_path': str(md_dir),
                'pdb_id': pdb_id,
                'timestamp': datetime.fromtimestamp(md_dir.stat().st_mtime).isoformat()
            }
            
            batch_entries.append(entry)
            
            if pdb_id:
                if pdb_id not in pdb_counts:
                    pdb_counts[pdb_id] = []
                pdb_counts[pdb_id].append(entry)
        
        # Analyze duplicates within batch
        duplicates = {pdb_id: entries for pdb_id, entries in pdb_counts.items() 
                     if len(entries) > 1}
        
        # Batch statistics
        batch_stats = {
            'batch_name': batch_name,
            'batch_path': str(batch_path),
            'total_directories': len(md_directories),
            'identified_pdbs': len([e for e in batch_entries if e['pdb_id']]),
            'unidentified_directories': len([e for e in batch_entries if not e['pdb_id']]),
            'unique_pdbs': len(pdb_counts),
            'duplicate_pdbs': len(duplicates),
            'analysis_timestamp': datetime.now().isoformat()
        }
        
        batch_data = {
            'statistics': batch_stats,
            'entries': batch_entries,
            'pdb_counts': pdb_counts,
            'duplicates': duplicates
        }
        
        # Update tracking database
        self.tracking_data['batches'][batch_name] = batch_data
        self._update_pdb_entries(batch_name, batch_entries)
        self._save_tracking_data()
        
        return batch_data
    
    def _update_pdb_entries(self, batch_name: str, batch_entries: List[Dict]):
        """Update global PDB entry tracking."""
        for entry in batch_entries:
            pdb_id = entry['pdb_id']
            if pdb_id:
                if pdb_id not in self.tracking_data['pdb_entries']:
                    self.tracking_data['pdb_entries'][pdb_id] = []
                    
                self.tracking_data['pdb_entries'][pdb_id].append({
                    'batch_name': batch_name,
                    'directory_name': entry['directory_name'],
                    'directory_path': entry['directory_path'],
                    'timestamp': entry['timestamp']
                })
    
    def find_duplicates_across_batches(self) -> Dict[str, List[Dict]]:
        """
        Find PDB IDs that appear in multiple batches.
        
        Returns:
            Dictionary mapping PDB IDs to their batch occurrences
        """
        cross_batch_duplicates = {}
        
        for pdb_id, entries in self.tracking_data['pdb_entries'].items():
            if len(entries) > 1:
                # Group by batch
                batch_groups = {}
                for entry in entries:
                    batch_name = entry['batch_name']
                    if batch_name not in batch_groups:
                        batch_groups[batch_name] = []
                    batch_groups[batch_name].append(entry)
                
                if len(batch_groups) > 1:  # Appears in multiple batches
                    cross_batch_duplicates[pdb_id] = {
                        'total_occurrences': len(entries),
                        'batch_count': len(batch_groups),
                        'batches': batch_groups
                    }
        
        return cross_batch_duplicates
    
    def find_missing_pdbs(self, expected_pdbs: Set[str]) -> Dict:
        """
        Find PDB IDs that are expected but missing from all batches.
        
        Args:
            expected_pdbs: Set of expected PDB IDs
            
        Returns:
            Dictionary with missing PDB analysis
        """
        processed_pdbs = set(self.tracking_data['pdb_entries'].keys())
        missing_pdbs = expected_pdbs - processed_pdbs
        extra_pdbs = processed_pdbs - expected_pdbs
        
        return {
            'expected_total': len(expected_pdbs),
            'processed_total': len(processed_pdbs),
            'missing_count': len(missing_pdbs),
            'extra_count': len(extra_pdbs),
            'missing_pdbs': sorted(list(missing_pdbs)),
            'extra_pdbs': sorted(list(extra_pdbs)),
            'completion_rate': (len(processed_pdbs.intersection(expected_pdbs)) / len(expected_pdbs)) * 100 if expected_pdbs else 0
        }
    
    def get_pdb_processing_history(self, pdb_id: str) -> Optional[List[Dict]]:
        """
        Get complete processing history for a specific PDB ID.
        
        Args:
            pdb_id: PDB identifier
            
        Returns:
            List of processing entries for the PDB ID
        """
        return self.tracking_data['pdb_entries'].get(pdb_id.upper())
    
    def get_batch_summary(self, batch_name: str) -> Optional[Dict]:
        """
        Get summary for a specific batch.
        
        Args:
            batch_name: Name of the batch
            
        Returns:
            Batch summary data
        """
        batch_data = self.tracking_data['batches'].get(batch_name)
        if batch_data:
            return batch_data['statistics']
        return None
    
    def get_global_statistics(self) -> Dict:
        """
        Get global statistics across all tracked batches.
        
        Returns:
            Global processing statistics
        """
        all_batches = self.tracking_data['batches']
        
        if not all_batches:
            return {
                'total_batches': 0,
                'total_directories': 0,
                'total_unique_pdbs': 0,
                'total_duplicates': 0
            }
        
        total_directories = sum(batch['statistics']['total_directories'] 
                              for batch in all_batches.values())
        total_identified = sum(batch['statistics']['identified_pdbs'] 
                             for batch in all_batches.values())
        total_unidentified = sum(batch['statistics']['unidentified_directories'] 
                               for batch in all_batches.values())
        
        # Cross-batch duplicate analysis
        cross_duplicates = self.find_duplicates_across_batches()
        
        # PDB occurrence distribution
        occurrence_dist = {}
        for pdb_id, entries in self.tracking_data['pdb_entries'].items():
            count = len(entries)
            if count not in occurrence_dist:
                occurrence_dist[count] = 0
            occurrence_dist[count] += 1
        
        return {
            'total_batches': len(all_batches),
            'total_directories': total_directories,
            'total_identified_pdbs': total_identified,
            'total_unidentified': total_unidentified,
            'unique_pdbs': len(self.tracking_data['pdb_entries']),
            'cross_batch_duplicates': len(cross_duplicates),
            'pdb_occurrence_distribution': occurrence_dist,
            'identification_rate': (total_identified / total_directories) * 100 if total_directories > 0 else 0,
            'last_updated': self.tracking_data.get('last_updated')
        }
    
    def generate_batch_report(self, batch_name: str, output_file: Optional[Union[str, Path]] = None) -> Dict:
        """
        Generate detailed report for a specific batch.
        
        Args:
            batch_name: Name of the batch to report on
            output_file: Optional file path to save the report
            
        Returns:
            Detailed batch report
        """
        batch_data = self.tracking_data['batches'].get(batch_name)
        if not batch_data:
            raise ValueError(f"Batch '{batch_name}' not found in tracking data")
        
        # Enhanced analysis
        entries = batch_data['entries']
        pdb_counts = batch_data['pdb_counts']
        
        # Directory naming analysis
        naming_patterns = {}
        for entry in entries:
            dir_name = entry['directory_name']
            # Analyze naming patterns
            if '_' in dir_name:
                pattern = 'contains_underscore'
            elif '-' in dir_name:
                pattern = 'contains_dash'
            elif dir_name.isalnum():
                pattern = 'alphanumeric_only'
            else:
                pattern = 'other'
                
            naming_patterns[pattern] = naming_patterns.get(pattern, 0) + 1
        
        # PDB ID analysis
        pdb_analysis = {
            'single_occurrence': len([pdb for pdb, entries in pdb_counts.items() if len(entries) == 1]),
            'multiple_occurrence': len([pdb for pdb, entries in pdb_counts.items() if len(entries) > 1]),
            'max_occurrences': max([len(entries) for entries in pdb_counts.values()]) if pdb_counts else 0
        }
        
        report = {
            'batch_name': batch_name,
            'report_timestamp': datetime.now().isoformat(),
            'batch_statistics': batch_data['statistics'],
            'pdb_analysis': pdb_analysis,
            'naming_patterns': naming_patterns,
            'duplicates': batch_data['duplicates'],
            'unidentified_directories': [e for e in entries if not e['pdb_id']],
            'recommendations': self._generate_recommendations(batch_data)
        }
        
        if output_file:
            output_path = Path(output_file)
            try:
                with open(output_path, 'w') as f:
                    json.dump(report, f, indent=2)
                logger.info(f"Batch report saved to: {output_path}")
            except Exception as e:
                logger.error(f"Failed to save batch report: {e}")
        
        return report
    
    def _generate_recommendations(self, batch_data: Dict) -> List[str]:
        """Generate recommendations based on batch analysis."""
        recommendations = []
        stats = batch_data['statistics']
        
        if stats['unidentified_directories'] > 0:
            recommendations.append(
                f"Consider reviewing {stats['unidentified_directories']} directories with unidentifiable PDB IDs"
            )
        
        if stats['duplicate_pdbs'] > 0:
            recommendations.append(
                f"Found {stats['duplicate_pdbs']} PDB IDs with multiple runs - review for redundancy"
            )
        
        if stats['identified_pdbs'] / stats['total_directories'] < 0.9:
            recommendations.append(
                "Low PDB identification rate - consider updating PDB ID extraction pattern"
            )
        
        return recommendations
    
    def merge_batch_reports(self, batch_names: List[str], 
                           output_file: Optional[Union[str, Path]] = None) -> Dict:
        """
        Merge multiple batch reports into a consolidated analysis.
        
        Args:
            batch_names: List of batch names to merge
            output_file: Optional file path to save the merged report
            
        Returns:
            Merged analysis report
        """
        merged_data = {
            'merged_batches': batch_names,
            'merge_timestamp': datetime.now().isoformat(),
            'combined_statistics': {},
            'cross_batch_analysis': {},
            'recommendations': []
        }
        
        # Combine statistics
        total_dirs = sum(self.tracking_data['batches'][name]['statistics']['total_directories'] 
                        for name in batch_names if name in self.tracking_data['batches'])
        total_pdbs = sum(self.tracking_data['batches'][name]['statistics']['identified_pdbs'] 
                        for name in batch_names if name in self.tracking_data['batches'])
        
        merged_data['combined_statistics'] = {
            'total_batches': len(batch_names),
            'total_directories': total_dirs,
            'total_identified_pdbs': total_pdbs,
            'identification_rate': (total_pdbs / total_dirs) * 100 if total_dirs > 0 else 0
        }
        
        # Cross-batch analysis
        cross_duplicates = self.find_duplicates_across_batches()
        relevant_duplicates = {pdb_id: data for pdb_id, data in cross_duplicates.items() 
                             if any(batch in data['batches'] for batch in batch_names)}
        
        merged_data['cross_batch_analysis'] = {
            'cross_batch_duplicates': len(relevant_duplicates),
            'duplicate_details': relevant_duplicates
        }
        
        if output_file:
            output_path = Path(output_file)
            try:
                with open(output_path, 'w') as f:
                    json.dump(merged_data, f, indent=2)
                logger.info(f"Merged report saved to: {output_path}")
            except Exception as e:
                logger.error(f"Failed to save merged report: {e}")
        
        return merged_data
    
    def cleanup_tracking_data(self, keep_recent_days: int = 30):
        """
        Clean up old tracking data to prevent database bloat.
        
        Args:
            keep_recent_days: Number of recent days to keep in tracking data
        """
        from datetime import timedelta
        
        cutoff_date = datetime.now() - timedelta(days=keep_recent_days)
        
        # Clean old batch data
        batches_to_remove = []
        for batch_name, batch_data in self.tracking_data['batches'].items():
            batch_timestamp = datetime.fromisoformat(batch_data['statistics']['analysis_timestamp'])
            if batch_timestamp < cutoff_date:
                batches_to_remove.append(batch_name)
        
        for batch_name in batches_to_remove:
            del self.tracking_data['batches'][batch_name]
            logger.info(f"Removed old batch data: {batch_name}")
        
        # Clean PDB entries that no longer have associated batches
        active_batches = set(self.tracking_data['batches'].keys())
        for pdb_id, entries in list(self.tracking_data['pdb_entries'].items()):
            filtered_entries = [e for e in entries if e['batch_name'] in active_batches]
            if filtered_entries:
                self.tracking_data['pdb_entries'][pdb_id] = filtered_entries
            else:
                del self.tracking_data['pdb_entries'][pdb_id]
        
        self._save_tracking_data()
        logger.info(f"Cleaned up tracking data, keeping {keep_recent_days} recent days")