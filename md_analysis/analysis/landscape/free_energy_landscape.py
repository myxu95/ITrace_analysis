"""
Free Energy Landscape (FEL) Builder

Constructs 2D free energy landscapes from MD trajectories using collective variables.
Supports RMSD, radius of gyration, contact number, and custom CVs.
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, contacts
from pathlib import Path
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass
import warnings


@dataclass
class FELResult:
    """Container for FEL computation results"""
    cv1_values: np.ndarray
    cv2_values: np.ndarray
    cv1_edges: np.ndarray
    cv2_edges: np.ndarray
    free_energy: np.ndarray
    cv1_name: str
    cv2_name: str
    temperature: float
    min_energy: float
    histogram: np.ndarray
    n_frames: int
    basins: List[Dict]


class FELBuilder:
    """
    Free Energy Landscape Builder

    Constructs 2D free energy landscapes from MD trajectories using:
    F(CV1, CV2) = -kT ln(P(CV1, CV2))

    where P is the probability distribution estimated from trajectory sampling.
    """

    def __init__(
        self,
        cv1: str = "rmsd",
        cv2: str = "contact_number",
        temperature: float = 300.0,
        bins: int = 50,
        kb: float = 0.001987204  # kcal/(mol·K)
    ):
        """
        Initialize FEL builder

        Args:
            cv1: First collective variable (rmsd, rg, contact_number)
            cv2: Second collective variable
            temperature: Temperature in Kelvin
            bins: Number of bins for 2D histogram
            kb: Boltzmann constant in kcal/(mol·K)
        """
        self.cv1 = cv1
        self.cv2 = cv2
        self.temperature = temperature
        self.bins = bins
        self.kb = kb
        self.kT = kb * temperature  # kcal/mol

    def compute_rmsd(
        self,
        universe: mda.Universe,
        reference: mda.Universe,
        selection: str = "backbone"
    ) -> np.ndarray:
        """
        Compute RMSD time series

        Args:
            universe: Trajectory universe
            reference: Reference structure
            selection: Atom selection for RMSD calculation

        Returns:
            RMSD values for each frame (Å)
        """
        rmsd_values = []
        ref_atoms = reference.select_atoms(selection)

        for ts in universe.trajectory:
            mobile_atoms = universe.select_atoms(selection)
            rmsd_val = rms.rmsd(
                mobile_atoms.positions,
                ref_atoms.positions,
                superposition=True
            )
            rmsd_values.append(rmsd_val)

        return np.array(rmsd_values)

    def compute_radius_of_gyration(
        self,
        universe: mda.Universe,
        selection: str = "protein"
    ) -> np.ndarray:
        """
        Compute radius of gyration time series

        Args:
            universe: Trajectory universe
            selection: Atom selection for Rg calculation

        Returns:
            Rg values for each frame (Å)
        """
        rg_values = []

        for ts in universe.trajectory:
            atoms = universe.select_atoms(selection)
            rg_values.append(atoms.radius_of_gyration())

        return np.array(rg_values)

    def compute_contact_number(
        self,
        universe: mda.Universe,
        selection_a: str = "segid A",
        selection_b: str = "segid B",
        cutoff: float = 4.5
    ) -> np.ndarray:
        """
        Compute number of contacts between two selections

        Args:
            universe: Trajectory universe
            selection_a: First selection (e.g., receptor)
            selection_b: Second selection (e.g., ligand)
            cutoff: Distance cutoff for contacts (Å)

        Returns:
            Contact number for each frame
        """
        contact_values = []
        group_a = universe.select_atoms(selection_a)
        group_b = universe.select_atoms(selection_b)

        for ts in universe.trajectory:
            # Compute distance matrix
            dist_matrix = contacts.distance_array(
                group_a.positions,
                group_b.positions
            )
            # Count contacts below cutoff
            n_contacts = np.sum(dist_matrix < cutoff)
            contact_values.append(n_contacts)

        return np.array(contact_values)

    def compute_cv(
        self,
        universe: mda.Universe,
        cv_name: str,
        reference: Optional[mda.Universe] = None,
        **kwargs
    ) -> np.ndarray:
        """
        Compute collective variable time series

        Args:
            universe: Trajectory universe
            cv_name: Name of CV (rmsd, rg, contact_number)
            reference: Reference structure (required for RMSD)
            **kwargs: Additional arguments for CV computation

        Returns:
            CV values for each frame
        """
        if cv_name == "rmsd":
            if reference is None:
                raise ValueError("Reference structure required for RMSD calculation")
            selection = kwargs.get("selection", "backbone")
            return self.compute_rmsd(universe, reference, selection)

        elif cv_name == "rg":
            selection = kwargs.get("selection", "protein")
            return self.compute_radius_of_gyration(universe, selection)

        elif cv_name == "contact_number":
            selection_a = kwargs.get("selection_a", "segid A")
            selection_b = kwargs.get("selection_b", "segid B")
            cutoff = kwargs.get("cutoff", 4.5)
            return self.compute_contact_number(
                universe, selection_a, selection_b, cutoff
            )

        else:
            raise ValueError(f"Unknown CV: {cv_name}")

    def build_2d_histogram(
        self,
        cv1_data: np.ndarray,
        cv2_data: np.ndarray,
        bins: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Build 2D histogram from CV data

        Args:
            cv1_data: First CV time series
            cv2_data: Second CV time series
            bins: Number of bins (uses self.bins if None)

        Returns:
            (histogram, cv1_edges, cv2_edges)
        """
        if bins is None:
            bins = self.bins

        # Remove NaN values
        valid_mask = ~(np.isnan(cv1_data) | np.isnan(cv2_data))
        cv1_clean = cv1_data[valid_mask]
        cv2_clean = cv2_data[valid_mask]

        if len(cv1_clean) == 0:
            raise ValueError("No valid data points after removing NaNs")

        # Build 2D histogram
        histogram, cv1_edges, cv2_edges = np.histogram2d(
            cv1_clean,
            cv2_clean,
            bins=bins
        )

        return histogram, cv1_edges, cv2_edges

    def histogram_to_free_energy(
        self,
        histogram: np.ndarray,
        min_count: int = 1
    ) -> np.ndarray:
        """
        Convert histogram to free energy landscape

        F = -kT ln(P) = -kT ln(N/N_total)

        Args:
            histogram: 2D histogram of CV sampling
            min_count: Minimum count to avoid log(0)

        Returns:
            Free energy landscape (kcal/mol)
        """
        # Avoid log(0) by setting minimum count
        histogram_safe = np.where(histogram < min_count, np.nan, histogram)

        # Normalize to probability
        total_counts = np.nansum(histogram_safe)
        probability = histogram_safe / total_counts

        # Compute free energy
        free_energy = -self.kT * np.log(probability)

        # Shift minimum to zero
        min_energy = np.nanmin(free_energy)
        free_energy = free_energy - min_energy

        return free_energy

    def identify_basins(
        self,
        free_energy: np.ndarray,
        cv1_edges: np.ndarray,
        cv2_edges: np.ndarray,
        energy_threshold: float = 2.0,
        min_size: int = 5
    ) -> List[Dict]:
        """
        Identify energy basins (local minima)

        Args:
            free_energy: Free energy landscape
            cv1_edges: CV1 bin edges
            cv2_edges: CV2 bin edges
            energy_threshold: Energy cutoff above minimum (kcal/mol)
            min_size: Minimum basin size (number of bins)

        Returns:
            List of basin dictionaries with location and properties
        """
        from scipy.ndimage import label, center_of_mass

        # Find regions below threshold
        valid_mask = ~np.isnan(free_energy)
        low_energy_mask = (free_energy <= energy_threshold) & valid_mask

        # Label connected regions
        labeled_array, num_features = label(low_energy_mask)

        basins = []
        for basin_id in range(1, num_features + 1):
            basin_mask = labeled_array == basin_id
            basin_size = np.sum(basin_mask)

            if basin_size < min_size:
                continue

            # Find minimum energy in basin
            basin_energies = free_energy[basin_mask]
            min_energy = np.min(basin_energies)

            # Find center of mass
            com = center_of_mass(basin_mask)
            cv1_center = cv1_edges[int(com[0])]
            cv2_center = cv2_edges[int(com[1])]

            basins.append({
                'id': basin_id,
                'size': basin_size,
                'min_energy': min_energy,
                'cv1_center': cv1_center,
                'cv2_center': cv2_center,
                'mask': basin_mask
            })

        # Sort by energy
        basins.sort(key=lambda x: x['min_energy'])

        return basins

    def build_fel(
        self,
        trajectory_path: str,
        topology_path: str,
        reference_pdb: Optional[str] = None,
        cv1_kwargs: Optional[Dict] = None,
        cv2_kwargs: Optional[Dict] = None
    ) -> FELResult:
        """
        Build complete free energy landscape

        Args:
            trajectory_path: Path to trajectory file
            topology_path: Path to topology file
            reference_pdb: Path to reference structure (for RMSD)
            cv1_kwargs: Additional arguments for CV1 computation
            cv2_kwargs: Additional arguments for CV2 computation

        Returns:
            FELResult object with all FEL data
        """
        if cv1_kwargs is None:
            cv1_kwargs = {}
        if cv2_kwargs is None:
            cv2_kwargs = {}

        # Load trajectory
        universe = mda.Universe(topology_path, trajectory_path)
        n_frames = len(universe.trajectory)

        # Load reference if needed
        reference = None
        if self.cv1 == "rmsd" or self.cv2 == "rmsd":
            if reference_pdb is None:
                raise ValueError("Reference PDB required for RMSD calculation")
            reference = mda.Universe(reference_pdb)

        # Compute CV1
        print(f"Computing {self.cv1}...")
        cv1_data = self.compute_cv(universe, self.cv1, reference, **cv1_kwargs)

        # Compute CV2
        print(f"Computing {self.cv2}...")
        cv2_data = self.compute_cv(universe, self.cv2, reference, **cv2_kwargs)

        # Build histogram
        print("Building 2D histogram...")
        histogram, cv1_edges, cv2_edges = self.build_2d_histogram(cv1_data, cv2_data)

        # Convert to free energy
        print("Computing free energy...")
        free_energy = self.histogram_to_free_energy(histogram)
        min_energy = np.nanmin(free_energy)

        # Identify basins
        print("Identifying energy basins...")
        basins = self.identify_basins(free_energy, cv1_edges, cv2_edges)

        # Compute bin centers
        cv1_centers = (cv1_edges[:-1] + cv1_edges[1:]) / 2
        cv2_centers = (cv2_edges[:-1] + cv2_edges[1:]) / 2

        return FELResult(
            cv1_values=cv1_centers,
            cv2_values=cv2_centers,
            cv1_edges=cv1_edges,
            cv2_edges=cv2_edges,
            free_energy=free_energy,
            cv1_name=self.cv1,
            cv2_name=self.cv2,
            temperature=self.temperature,
            min_energy=min_energy,
            histogram=histogram,
            n_frames=n_frames,
            basins=basins
        )
