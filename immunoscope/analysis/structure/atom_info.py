import numpy as np
import pandas as pd
import MDAnalysis as mda
from typing import Optional, Dict
import logging

logger = logging.getLogger(__name__)


class AtomInfoExtractor:
    """Atom information extraction module for structures."""
    
    def __init__(self, structure_file: str):
        """
        Initialize AtomInfoExtractor.
        
        Args:
            structure_file: Structure file (.pdb, .gro, etc.)
        """
        self.structure_file = structure_file
        
        try:
            self.universe = mda.Universe(structure_file)
            logger.info(f"Atom Info Extractor initialized: {self.universe.atoms.n_atoms} atoms")
        except Exception as e:
            logger.error(f"Failed to load structure for atom info extraction: {e}")
            raise
    
    def extract_comprehensive_info(self,
                                 selection: str = "all",
                                 output_file: Optional[str] = None) -> pd.DataFrame:
        """
        Extract comprehensive atom information.
        
        Args:
            selection: Atom selection
            output_file: Optional output file
            
        Returns:
            DataFrame with comprehensive atom information
        """
        atoms = self.universe.select_atoms(selection)
        
        # Try to get B-factors
        try:
            bfactors = atoms.tempfactors
        except AttributeError:
            bfactors = np.zeros(len(atoms))
            logger.warning("B-factors not available, using zeros")
        
        # Try to get occupancy
        try:
            occupancy = atoms.occupancies
        except AttributeError:
            occupancy = np.ones(len(atoms))
            logger.warning("Occupancy not available, using ones")
        
        # Try to get charges
        try:
            charges = atoms.charges
        except AttributeError:
            charges = np.zeros(len(atoms))
            logger.warning("Charges not available, using zeros")
        
        df = pd.DataFrame({
            'Atom_Index': atoms.indices,
            'Atom_Name': atoms.names,
            'Element': atoms.elements,
            'Residue_ID': atoms.resids,
            'Residue_Name': atoms.resnames,
            'Chain': atoms.segids,
            'X': atoms.positions[:, 0],
            'Y': atoms.positions[:, 1],
            'Z': atoms.positions[:, 2],
            'B_Factor': bfactors,
            'Occupancy': occupancy,
            'Charge': charges
        })
        
        if output_file:
            df.to_csv(output_file, index=False)
            logger.info(f"Comprehensive atom info saved to: {output_file}")
        
        logger.info(f"Comprehensive atom information extracted for {len(atoms)} atoms")
        return df
    
    def get_structure_summary(self) -> Dict:
        """Get summary statistics of the structure."""
        atoms = self.universe.atoms
        
        # Basic counts
        n_atoms = atoms.n_atoms
        n_residues = atoms.n_residues
        n_segments = atoms.n_segments
        
        # Atom types
        unique_elements = list(set(atoms.elements))
        unique_atom_names = list(set(atoms.names))
        unique_residues = list(set(atoms.resnames))
        
        # Element counts
        element_counts = {}
        for element in unique_elements:
            count = np.sum(atoms.elements == element)
            element_counts[element] = int(count)
        
        # Residue counts
        residue_counts = {}
        for residue in unique_residues:
            count = np.sum(atoms.resnames == residue)
            residue_counts[residue] = int(count)
        
        # Coordinate ranges
        positions = atoms.positions
        coord_ranges = {
            'x_range': [float(np.min(positions[:, 0])), float(np.max(positions[:, 0]))],
            'y_range': [float(np.min(positions[:, 1])), float(np.max(positions[:, 1]))],
            'z_range': [float(np.min(positions[:, 2])), float(np.max(positions[:, 2]))]
        }
        
        # Box dimensions if available
        dimensions = None
        if self.universe.dimensions is not None:
            dimensions = self.universe.dimensions.tolist()
        
        summary = {
            "file_path": self.structure_file,
            "n_atoms": n_atoms,
            "n_residues": n_residues,
            "n_segments": n_segments,
            "unique_elements": unique_elements,
            "unique_atom_names": unique_atom_names,
            "unique_residues": unique_residues,
            "element_counts": element_counts,
            "residue_counts": residue_counts,
            "coordinate_ranges": coord_ranges,
            "box_dimensions": dimensions
        }
        
        logger.info(f"Structure summary generated for {n_atoms} atoms")
        return summary
    
    def extract_residue_info(self,
                           selection: str = "protein",
                           output_file: Optional[str] = None) -> pd.DataFrame:
        """
        Extract residue-level information.
        
        Args:
            selection: Atom selection
            output_file: Optional output file
            
        Returns:
            DataFrame with residue information
        """
        atoms = self.universe.select_atoms(selection)
        
        # Group by residue
        residue_data = []
        
        for residue in atoms.residues:
            res_atoms = residue.atoms
            
            # Calculate center of mass
            com = res_atoms.center_of_mass()
            
            # Get atom counts by element
            elements = res_atoms.elements
            element_counts = {elem: int(np.sum(elements == elem)) for elem in set(elements)}
            
            # Try to get average B-factor
            try:
                avg_bfactor = float(np.mean(res_atoms.tempfactors))
            except AttributeError:
                avg_bfactor = 0.0
            
            residue_data.append({
                'Residue_ID': residue.resid,
                'Residue_Name': residue.resname,
                'Chain': residue.segment.segid,
                'N_Atoms': len(res_atoms),
                'COM_X': com[0],
                'COM_Y': com[1],
                'COM_Z': com[2],
                'Avg_BFactor': avg_bfactor,
                'Element_Counts': str(element_counts)
            })
        
        df = pd.DataFrame(residue_data)
        
        if output_file:
            df.to_csv(output_file, index=False)
            logger.info(f"Residue info saved to: {output_file}")
        
        logger.info(f"Residue information extracted for {len(df)} residues")
        return df
    
    def find_atoms_by_criteria(self,
                             element: Optional[str] = None,
                             residue_name: Optional[str] = None,
                             atom_name: Optional[str] = None,
                             bfactor_threshold: Optional[float] = None,
                             output_file: Optional[str] = None) -> pd.DataFrame:
        """
        Find atoms matching specific criteria.
        
        Args:
            element: Element type filter
            residue_name: Residue name filter
            atom_name: Atom name filter
            bfactor_threshold: B-factor threshold filter
            output_file: Optional output file
            
        Returns:
            DataFrame with matching atoms
        """
        atoms = self.universe.atoms
        
        # Start with all atoms
        mask = np.ones(len(atoms), dtype=bool)
        
        # Apply filters
        if element:
            mask &= (atoms.elements == element)
        
        if residue_name:
            mask &= (atoms.resnames == residue_name)
        
        if atom_name:
            mask &= (atoms.names == atom_name)
        
        if bfactor_threshold is not None:
            try:
                mask &= (atoms.tempfactors > bfactor_threshold)
            except AttributeError:
                logger.warning("B-factors not available for threshold filtering")
        
        # Extract matching atoms
        filtered_atoms = atoms[mask]
        
        if len(filtered_atoms) == 0:
            logger.warning("No atoms found matching the specified criteria")
            return pd.DataFrame()
        
        # Create DataFrame
        try:
            bfactors = filtered_atoms.tempfactors
        except AttributeError:
            bfactors = np.zeros(len(filtered_atoms))
        
        df = pd.DataFrame({
            'Atom_Index': filtered_atoms.indices,
            'Atom_Name': filtered_atoms.names,
            'Element': filtered_atoms.elements,
            'Residue_ID': filtered_atoms.resids,
            'Residue_Name': filtered_atoms.resnames,
            'Chain': filtered_atoms.segids,
            'X': filtered_atoms.positions[:, 0],
            'Y': filtered_atoms.positions[:, 1],
            'Z': filtered_atoms.positions[:, 2],
            'B_Factor': bfactors
        })
        
        if output_file:
            df.to_csv(output_file, index=False)
            logger.info(f"Filtered atom info saved to: {output_file}")
        
        logger.info(f"Found {len(df)} atoms matching criteria")
        return df