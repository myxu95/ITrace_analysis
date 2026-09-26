import numpy as np
import pandas as pd
import MDAnalysis as mda
from typing import Dict, Optional
import logging

logger = logging.getLogger(__name__)


class GeometryAnalyzer:
    """Geometric analysis module for protein structures."""
    
    def __init__(self, structure_file: str):
        """
        Initialize GeometryAnalyzer.
        
        Args:
            structure_file: Structure file (.pdb, .gro, etc.)
        """
        self.structure_file = structure_file
        
        try:
            self.universe = mda.Universe(structure_file)
            logger.info(f"Geometry Analyzer initialized: {self.universe.atoms.n_atoms} atoms")
        except Exception as e:
            logger.error(f"Failed to load structure for geometry analysis: {e}")
            raise
    
    def analyze_basic_properties(self,
                                selection: str = "protein") -> Dict:
        """
        Analyze basic geometric properties of the structure.
        
        Args:
            selection: Atom selection
            
        Returns:
            Dictionary with geometric properties
        """
        atoms = self.universe.select_atoms(selection)
        
        # Center of mass
        com = atoms.center_of_mass()
        
        # Radius of gyration
        rg = atoms.radius_of_gyration()
        
        # Positions relative to center of mass
        positions = atoms.positions
        centered_pos = positions - com
        
        # Moments of inertia (simplified)
        xx = np.sum(centered_pos[:, 0]**2)
        yy = np.sum(centered_pos[:, 1]**2)
        zz = np.sum(centered_pos[:, 2]**2)
        xy = np.sum(centered_pos[:, 0] * centered_pos[:, 1])
        xz = np.sum(centered_pos[:, 0] * centered_pos[:, 2])
        yz = np.sum(centered_pos[:, 1] * centered_pos[:, 2])
        
        # Inertia tensor
        inertia_tensor = np.array([
            [xx, -xy, -xz],
            [-xy, yy, -yz],
            [-xz, -yz, zz]
        ])
        
        # Principal moments and axes
        eigenvals, eigenvecs = np.linalg.eigh(inertia_tensor)
        
        # Bounding box
        min_coords = np.min(positions, axis=0)
        max_coords = np.max(positions, axis=0)
        box_dimensions = max_coords - min_coords
        
        # Asphericity and other shape parameters
        rg_components = self._calculate_rg_components(atoms)
        asphericity = self._calculate_asphericity(rg_components)
        
        results = {
            "center_of_mass": com.tolist(),
            "radius_of_gyration": float(rg),
            "rg_components": rg_components,
            "asphericity": asphericity,
            "n_atoms": len(atoms),
            "bounding_box": {
                "min": min_coords.tolist(),
                "max": max_coords.tolist(),
                "dimensions": box_dimensions.tolist(),
                "volume": float(np.prod(box_dimensions))
            },
            "inertia_tensor": inertia_tensor.tolist(),
            "principal_moments": eigenvals.tolist(),
            "principal_axes": eigenvecs.tolist()
        }
        
        logger.info(f"Basic geometric analysis completed for {len(atoms)} atoms")
        return results
    
    def analyze_secondary_structure_geometry(self,
                                           selection: str = "protein") -> Dict:
        """
        Analyze geometric properties related to secondary structure.
        
        Args:
            selection: Protein selection
            
        Returns:
            Dictionary with secondary structure geometry
        """
        # This is a simplified analysis - for more sophisticated analysis,
        # consider integrating with DSSP or similar tools
        
        atoms = self.universe.select_atoms(selection)
        ca_atoms = self.universe.select_atoms(f"({selection}) and name CA")
        
        if len(ca_atoms) < 3:
            logger.warning("Not enough CA atoms for secondary structure analysis")
            return {"error": "Insufficient CA atoms"}
        
        # Calculate backbone angles (simplified)
        ca_positions = ca_atoms.positions
        n_residues = len(ca_positions)
        
        # Calculate consecutive CA-CA distances
        ca_distances = []
        for i in range(n_residues - 1):
            dist = np.linalg.norm(ca_positions[i+1] - ca_positions[i])
            ca_distances.append(dist)
        
        # Calculate backbone angles
        angles = []
        for i in range(n_residues - 2):
            v1 = ca_positions[i+1] - ca_positions[i]
            v2 = ca_positions[i+2] - ca_positions[i+1]
            
            cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Handle numerical errors
            angle = np.arccos(cos_angle) * 180 / np.pi
            angles.append(angle)
        
        results = {
            "n_residues": n_residues,
            "ca_distances": {
                "mean": float(np.mean(ca_distances)),
                "std": float(np.std(ca_distances)),
                "min": float(np.min(ca_distances)),
                "max": float(np.max(ca_distances))
            },
            "backbone_angles": {
                "mean": float(np.mean(angles)),
                "std": float(np.std(angles)),
                "min": float(np.min(angles)),
                "max": float(np.max(angles))
            }
        }
        
        logger.info(f"Secondary structure geometry analysis completed")
        return results
    
    def calculate_surface_area_approximation(self,
                                           selection: str = "protein",
                                           probe_radius: float = 1.4) -> Dict:
        """
        Calculate approximate surface area using a simple geometric approach.
        
        Args:
            selection: Atom selection
            probe_radius: Probe radius for surface calculation
            
        Returns:
            Dictionary with surface area estimates
        """
        atoms = self.universe.select_atoms(selection)
        positions = atoms.positions
        
        # Use van der Waals radii (simplified)
        vdw_radii = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "H": 1.2}
        default_radius = 1.5
        
        atom_radii = []
        for element in atoms.elements:
            radius = vdw_radii.get(element, default_radius) + probe_radius
            atom_radii.append(radius)
        
        atom_radii = np.array(atom_radii)
        
        # Simple surface area calculation (sum of spheres minus overlaps)
        total_area = 4 * np.pi * np.sum(atom_radii**2)
        
        # Estimate buried area (very simplified)
        n_atoms = len(atoms)
        buried_area = 0
        
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                dist = np.linalg.norm(positions[i] - positions[j])
                sum_radii = atom_radii[i] + atom_radii[j]
                
                if dist < sum_radii:
                    # Simplified overlap calculation
                    overlap_factor = (sum_radii - dist) / sum_radii
                    area_overlap = 4 * np.pi * min(atom_radii[i], atom_radii[j])**2 * overlap_factor
                    buried_area += area_overlap
        
        accessible_area = max(0, total_area - buried_area)
        
        results = {
            "total_sphere_area": float(total_area),
            "estimated_buried_area": float(buried_area),
            "estimated_accessible_area": float(accessible_area),
            "burial_fraction": float(buried_area / total_area) if total_area > 0 else 0,
            "probe_radius": probe_radius,
            "n_atoms": n_atoms
        }
        
        logger.info(f"Surface area approximation completed")
        return results
    
    def _calculate_rg_components(self, atoms) -> Dict:
        """Calculate radius of gyration components."""
        positions = atoms.positions
        com = atoms.center_of_mass()
        rel_pos = positions - com
        
        rg_total = atoms.radius_of_gyration()
        rgx = np.sqrt(np.mean(rel_pos[:, 0]**2))
        rgy = np.sqrt(np.mean(rel_pos[:, 1]**2))
        rgz = np.sqrt(np.mean(rel_pos[:, 2]**2))
        
        return {
            "total": float(rg_total),
            "x": float(rgx),
            "y": float(rgy),
            "z": float(rgz)
        }
    
    def _calculate_asphericity(self, rg_components: Dict) -> float:
        """Calculate asphericity parameter."""
        rg_vals = [rg_components["x"], rg_components["y"], rg_components["z"]]
        rg_vals.sort(reverse=True)  # Sort in descending order
        
        # Asphericity: measure of deviation from spherical shape
        rg_mean = np.mean(rg_vals)
        if rg_mean > 0:
            asphericity = sum((rg - rg_mean)**2 for rg in rg_vals) / (3 * rg_mean**2)
        else:
            asphericity = 0.0
        
        return float(asphericity)