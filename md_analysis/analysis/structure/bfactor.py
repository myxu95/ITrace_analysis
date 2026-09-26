import numpy as np
import pandas as pd
import MDAnalysis as mda
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class BFactorAnalyzer:
    """B-factor analysis module for PDB structures."""
    
    def __init__(self, structure_file: str):
        """
        Initialize BFactorAnalyzer.
        
        Args:
            structure_file: Structure file (.pdb, .gro, etc.)
        """
        self.structure_file = structure_file
        
        try:
            self.universe = mda.Universe(structure_file)
            logger.info(f"B-factor Analyzer initialized: {self.universe.atoms.n_atoms} atoms")
        except Exception as e:
            logger.error(f"Failed to load structure for B-factor analysis: {e}")
            raise
    
    def extract_bfactors(self, 
                        selection: str = "all",
                        output_file: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract B-factors from PDB structure.
        
        Args:
            selection: Atom selection
            output_file: Optional output file for B-factor data
            
        Returns:
            Tuple of (atom_indices, bfactors) arrays
        """
        atoms = self.universe.select_atoms(selection)
        
        try:
            bfactors = atoms.tempfactors
            atom_indices = atoms.indices
            logger.info(f"Extracted B-factors for {len(atoms)} atoms")
        except AttributeError:
            logger.warning("B-factors not available in structure file, using zeros")
            bfactors = np.zeros(len(atoms))
            atom_indices = atoms.indices
        
        if output_file:
            self._save_bfactor_data(atoms, bfactors, output_file)
        
        return atom_indices, bfactors
    
    def analyze_by_residue(self,
                          selection: str = "protein",
                          output_file: Optional[str] = None) -> pd.DataFrame:
        """
        Analyze B-factors by residue.
        
        Args:
            selection: Atom selection
            output_file: Optional output file
            
        Returns:
            DataFrame with residue-level B-factor statistics
        """
        atoms = self.universe.select_atoms(selection)
        
        try:
            bfactors = atoms.tempfactors
        except AttributeError:
            logger.warning("B-factors not available, using zeros")
            bfactors = np.zeros(len(atoms))
        
        df = pd.DataFrame({
            'Residue_ID': atoms.resids,
            'Residue_Name': atoms.resnames,
            'Chain': atoms.segids,
            'B_Factor': bfactors
        })
        
        residue_stats = df.groupby(['Chain', 'Residue_ID', 'Residue_Name']).agg({
            'B_Factor': ['mean', 'std', 'min', 'max', 'count']
        }).round(3)
        
        residue_stats.columns = ['Mean_BFactor', 'Std_BFactor', 'Min_BFactor', 
                               'Max_BFactor', 'Atom_Count']
        residue_stats.reset_index(inplace=True)
        
        if output_file:
            residue_stats.to_csv(output_file, index=False)
        
        logger.info(f"B-factor analysis completed for {len(residue_stats)} residues")
        return residue_stats
    
    def analyze_by_atom_type(self,
                           selection: str = "protein",
                           output_file: Optional[str] = None) -> pd.DataFrame:
        """
        Analyze B-factors by atom type.
        
        Args:
            selection: Atom selection
            output_file: Optional output file
            
        Returns:
            DataFrame with atom type B-factor statistics
        """
        atoms = self.universe.select_atoms(selection)
        
        try:
            bfactors = atoms.tempfactors
        except AttributeError:
            bfactors = np.zeros(len(atoms))
        
        df = pd.DataFrame({
            'Atom_Name': atoms.names,
            'Element': atoms.elements,
            'B_Factor': bfactors
        })
        
        atom_stats = df.groupby(['Atom_Name', 'Element']).agg({
            'B_Factor': ['mean', 'std', 'min', 'max', 'count']
        }).round(3)
        
        atom_stats.columns = ['Mean_BFactor', 'Std_BFactor', 'Min_BFactor', 
                             'Max_BFactor', 'Count']
        atom_stats.reset_index(inplace=True)
        
        if output_file:
            atom_stats.to_csv(output_file, index=False)
        
        logger.info(f"B-factor analysis by atom type completed")
        return atom_stats
    
    def find_high_bfactor_residues(self,
                                  selection: str = "protein",
                                  threshold: float = None,
                                  percentile: float = 90.0) -> pd.DataFrame:
        """
        Find residues with high B-factors.
        
        Args:
            selection: Atom selection
            threshold: B-factor threshold (if None, use percentile)
            percentile: Percentile threshold if threshold is None
            
        Returns:
            DataFrame with high B-factor residues
        """
        residue_stats = self.analyze_by_residue(selection)
        
        if threshold is None:
            threshold = np.percentile(residue_stats['Mean_BFactor'], percentile)
            logger.info(f"Using {percentile}th percentile threshold: {threshold:.2f}")
        
        high_bfactor = residue_stats[residue_stats['Mean_BFactor'] > threshold].copy()
        high_bfactor = high_bfactor.sort_values('Mean_BFactor', ascending=False)
        
        logger.info(f"Found {len(high_bfactor)} residues with B-factor > {threshold:.2f}")
        return high_bfactor
    
    def _save_bfactor_data(self, atoms, bfactors: np.ndarray, filename: str):
        """Save detailed B-factor data."""
        df = pd.DataFrame({
            'Atom_Index': atoms.indices,
            'Residue_ID': atoms.resids,
            'Residue_Name': atoms.resnames,
            'Atom_Name': atoms.names,
            'Chain': atoms.segids,
            'B_Factor': bfactors
        })
        df.to_csv(filename, index=False)
        logger.info(f"B-factor data saved to: {filename}")