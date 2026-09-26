import os
import json
from pathlib import Path
from typing import Dict, Optional, Union
from datetime import datetime


class PathManager:
    """Path management utility for organizing MD analysis file paths."""
    
    def __init__(self, project_root: str, config_file: Optional[str] = None):
        """
        Initialize PathManager.
        
        Args:
            project_root: Root directory for the project
            config_file: Optional config file for custom path settings
        """
        self.project_root = Path(project_root)
        self.config_file = config_file
        self.paths = self._load_default_paths()
        
        if config_file and os.path.exists(config_file):
            self._load_config(config_file)
            
    def _load_default_paths(self) -> Dict[str, str]:
        """Load default directory structure."""
        return {
            "raw_data": "data/raw",
            "processed_data": "data/processed", 
            "trajectories": "data/trajectories",
            "structures": "data/structures",
            "results": "results",
            "plots": "results/plots",
            "analysis": "results/analysis",
            "logs": "logs",
            "temp": "temp",
            "scripts": "scripts"
        }
    
    def _load_config(self, config_file: str):
        """Load configuration from JSON file."""
        with open(config_file, 'r') as f:
            config = json.load(f)
            self.paths.update(config.get('paths', {}))
    
    def get_path(self, path_type: str, create: bool = True) -> Path:
        """
        Get path for specified type.
        
        Args:
            path_type: Type of path (e.g., 'results', 'plots')
            create: Whether to create directory if it doesn't exist
            
        Returns:
            Path object for the requested directory
        """
        if path_type not in self.paths:
            raise ValueError(f"Unknown path type: {path_type}")
            
        path = self.project_root / self.paths[path_type]
        
        if create:
            path.mkdir(parents=True, exist_ok=True)
            
        return path
    
    def get_timestamped_path(self, path_type: str, prefix: str = "") -> Path:
        """
        Get timestamped path for organizing results by time.
        
        Args:
            path_type: Base path type
            prefix: Optional prefix for the timestamped directory
            
        Returns:
            Timestamped path
        """
        base_path = self.get_path(path_type)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        if prefix:
            timestamped_dir = f"{prefix}_{timestamp}"
        else:
            timestamped_dir = timestamp
            
        timestamped_path = base_path / timestamped_dir
        timestamped_path.mkdir(parents=True, exist_ok=True)
        
        return timestamped_path
    
    def organize_file(self, 
                     file_path: Union[str, Path], 
                     destination_type: str,
                     subfolder: Optional[str] = None) -> Path:
        """
        Move or copy file to organized location.
        
        Args:
            file_path: Source file path
            destination_type: Destination path type
            subfolder: Optional subfolder within destination
            
        Returns:
            Final destination path
        """
        source = Path(file_path)
        dest_dir = self.get_path(destination_type)
        
        if subfolder:
            dest_dir = dest_dir / subfolder
            dest_dir.mkdir(parents=True, exist_ok=True)
            
        dest_path = dest_dir / source.name
        
        if source.exists():
            source.rename(dest_path)
            
        return dest_path
    
    def get_analysis_output_path(self, analysis_name: str, file_format: str = "png") -> Path:
        """
        Generate standardized output path for analysis results.
        
        Args:
            analysis_name: Name of the analysis
            file_format: Output file format
            
        Returns:
            Output file path
        """
        if file_format in ['png', 'jpg', 'svg', 'pdf']:
            output_dir = self.get_path('plots')
        else:
            output_dir = self.get_path('analysis')
            
        filename = f"{analysis_name}.{file_format}"
        return output_dir / filename
    
    def create_project_structure(self):
        """Create complete project directory structure."""
        for path_type in self.paths:
            self.get_path(path_type, create=True)
    
    def __str__(self) -> str:
        """String representation of path structure."""
        lines = [f"Project Root: {self.project_root}"]
        lines.append("Path Structure:")
        
        for path_type, relative_path in self.paths.items():
            full_path = self.project_root / relative_path
            exists = "✓" if full_path.exists() else "✗"
            lines.append(f"  {exists} {path_type}: {relative_path}")
            
        return "\n".join(lines)