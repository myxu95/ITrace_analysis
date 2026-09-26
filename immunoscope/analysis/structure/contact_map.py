import numpy as np
import pandas as pd
import MDAnalysis as mda
from typing import Optional
import logging

logger = logging.getLogger(__name__)


class ContactMapCalculator:
    """Contact map calculation module for protein structures."""
    
    def __init__(self, structure_file: str):
        """
        Initialize ContactMapCalculator.
        
        Args:
            structure_file: Structure file (.pdb, .gro, etc.)
        """
        self.structure_file = structure_file
        
        try:
            self.universe = mda.Universe(structure_file)
            logger.info(f"Contact Map Calculator initialized: {self.universe.atoms.n_atoms} atoms")
        except Exception as e:
            logger.error(f"Failed to load structure for contact map: {e}")
            raise
    
    def calculate(self,
                 selection: str = "name CA",
                 cutoff: float = 8.0,
                 output_file: Optional[str] = None) -> np.ndarray:
        """
        Calculate contact map for the structure.
        
        Args:
            selection: Atom selection for contact calculation
            cutoff: Distance cutoff for contacts (Angstroms)
            output_file: Optional output file for contact map
            
        Returns:
            Contact map matrix
        """
        atoms = self.universe.select_atoms(selection)
        positions = atoms.positions
        
        logger.info(f"Calculating contact map for {len(atoms)} atoms with cutoff {cutoff} Å")
        
        # Calculate distance matrix
        diff = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
        distances = np.sqrt(np.sum(diff**2, axis=2))
        
        # Create contact map based on cutoff
        contact_map = (distances < cutoff).astype(int)
        np.fill_diagonal(contact_map, 0)  # Remove self-contacts
        
        if output_file:
            self._save_contact_map(contact_map, atoms, output_file)
        
        n_contacts = np.sum(contact_map) // 2  # Divide by 2 for symmetric matrix
        logger.info(f"Contact map calculated: {n_contacts} contacts within {cutoff} Å")
        
        return contact_map
    
    def calculate_residue_contacts(self,
                                 selection: str = "protein",
                                 cutoff: float = 8.0,
                                 output_file: Optional[str] = None) -> np.ndarray:
        """
        Calculate residue-residue contact map.
        
        Args:
            selection: Atom selection
            cutoff: Distance cutoff for contacts (Angstroms)
            output_file: Optional output file
            
        Returns:
            Residue contact map matrix
        """
        atoms = self.universe.select_atoms(selection)
        residue_ids = np.unique(atoms.resids)
        n_residues = len(residue_ids)
        
        contact_map = np.zeros((n_residues, n_residues), dtype=int)
        
        logger.info(f"Calculating residue contact map for {n_residues} residues")
        
        for i, res_i in enumerate(residue_ids):
            atoms_i = atoms.select_atoms(f"resid {res_i}")
            for j, res_j in enumerate(residue_ids):
                if i >= j:  # Skip diagonal and lower triangle
                    continue
                    
                atoms_j = atoms.select_atoms(f"resid {res_j}")
                
                # Calculate minimum distance between residues
                pos_i = atoms_i.positions
                pos_j = atoms_j.positions
                
                diff = pos_i[:, np.newaxis, :] - pos_j[np.newaxis, :, :]
                distances = np.sqrt(np.sum(diff**2, axis=2))
                min_dist = np.min(distances)
                
                if min_dist < cutoff:
                    contact_map[i, j] = 1
                    contact_map[j, i] = 1  # Symmetric
        
        if output_file:
            self._save_residue_contact_map(contact_map, residue_ids, output_file)
        
        n_contacts = np.sum(contact_map) // 2
        logger.info(f"Residue contact map calculated: {n_contacts} residue-residue contacts")
        
        return contact_map
    
    def analyze_contact_distribution(self,
                                   contact_map: np.ndarray,
                                   residue_ids: Optional[np.ndarray] = None) -> pd.DataFrame:
        """
        Analyze contact distribution from contact map.
        
        Args:
            contact_map: Contact map matrix
            residue_ids: Optional residue IDs for labeling
            
        Returns:
            DataFrame with contact statistics per residue
        """
        n_residues = contact_map.shape[0]
        
        if residue_ids is None:
            residue_ids = np.arange(1, n_residues + 1)
        
        contact_counts = np.sum(contact_map, axis=1)
        
        df = pd.DataFrame({
            'Residue_ID': residue_ids,
            'Contact_Count': contact_counts,
            'Contact_Fraction': contact_counts / (n_residues - 1)  # Exclude self
        })
        
        df = df.sort_values('Contact_Count', ascending=False)
        
        logger.info(f"Contact distribution analysis completed")
        return df
    
    def find_long_range_contacts(self,
                               contact_map: np.ndarray,
                               residue_ids: Optional[np.ndarray] = None,
                               sequence_separation: int = 5) -> pd.DataFrame:
        """
        Find long-range contacts (contacts between residues far in sequence).
        
        Args:
            contact_map: Contact map matrix
            residue_ids: Optional residue IDs
            sequence_separation: Minimum sequence separation for long-range contacts
            
        Returns:
            DataFrame with long-range contacts
        """
        n_residues = contact_map.shape[0]
        
        if residue_ids is None:
            residue_ids = np.arange(1, n_residues + 1)
        
        long_range_contacts = []
        
        for i in range(n_residues):
            for j in range(i + sequence_separation, n_residues):
                if contact_map[i, j] == 1:
                    long_range_contacts.append({
                        'Residue_1': residue_ids[i],
                        'Residue_2': residue_ids[j],
                        'Sequence_Separation': j - i
                    })
        
        df = pd.DataFrame(long_range_contacts)
        
        if len(df) > 0:
            df = df.sort_values('Sequence_Separation', ascending=False)
        
        logger.info(f"Found {len(df)} long-range contacts (sep >= {sequence_separation})")
        return df
    
    def _save_contact_map(self, contact_map: np.ndarray, atoms, filename: str):
        """Save contact map to file."""
        if filename.endswith('.csv'):
            # Save as matrix with residue labels
            residue_ids = [f"{atoms.resnames[i]}{atoms.resids[i]}" for i in range(len(atoms))]
            df = pd.DataFrame(contact_map, index=residue_ids, columns=residue_ids)
            df.to_csv(filename)
        else:
            np.savetxt(filename, contact_map, fmt='%d')
        
        logger.info(f"Contact map saved to: {filename}")
    
    def _save_residue_contact_map(self, contact_map: np.ndarray, residue_ids: np.ndarray, filename: str):
        """Save residue contact map to file."""
        if filename.endswith('.csv'):
            df = pd.DataFrame(contact_map, 
                            index=[f"Res{rid}" for rid in residue_ids],
                            columns=[f"Res{rid}" for rid in residue_ids])
            df.to_csv(filename)
        else:
            np.savetxt(filename, contact_map, fmt='%d')
        
        logger.info(f"Residue contact map saved to: {filename}")