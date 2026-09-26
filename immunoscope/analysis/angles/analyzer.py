"""
Docking Angle Analyzer - Complete Implementation

This module contains all logic for TCR-pMHC docking angle calculation.
Internal helper classes are prefixed with underscore and not exported.

Public API:
    - DockingAngleAnalyzer: Main analyzer class conforming to ImmunoScope 6 design principles

Internal classes (not exported):
    - _GeometryUtils: Geometry utilities (plane fitting, vector angles)
    - _SequenceAligner: MHC sequence alignment
    - _MHCGrooveDetector: MHC groove geometry detection
    - _TCRAxisCalculator: TCR axis calculation using conserved cysteines

Author: ImmunoScope Development Team
Date: 2026-03-18 (Refactored from 9 files to 1 file)
"""

import logging
import csv
from typing import Tuple, Optional, Callable, Dict, List, Any
from datetime import datetime
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, Superimposer

from .angle_data_structures import DockingAngleInput, DockingAngleResult

logger = logging.getLogger(__name__)


# ============================================================================
# Internal Helper Classes (Not Exported)
# ============================================================================

class _GeometryUtils:
    """
    Internal geometry utilities (plane fitting, vector angles, etc.).

    This class combines functionality from:
    - plane_fitting.py (PlaneFitter)
    - vector_angles.py (utility functions)
    """

    @staticmethod
    def fit_plane_svd(atoms: mda.AtomGroup) -> Tuple[np.ndarray, np.ndarray]:
        """
        Fit plane to atom positions using SVD.

        Algorithm:
        1. Calculate centroid (center of mass) of atoms
        2. Center coordinates by subtracting centroid
        3. Apply SVD to centered coordinates
        4. Normal = smallest singular vector (3rd column of V matrix)

        Parameters
        ----------
        atoms : MDAnalysis.AtomGroup
            Atoms to fit plane through

        Returns
        -------
        centroid : np.ndarray, shape (3,)
            Center point of the plane (center of mass of atoms)
        normal : np.ndarray, shape (3,)
            Unit normal vector to the plane

        Raises
        ------
        ValueError
            If atoms group is empty or has insufficient atoms
        """
        if len(atoms) == 0:
            raise ValueError("Cannot fit plane to empty atom group")

        if len(atoms) < 3:
            raise ValueError(
                f"Need at least 3 atoms to fit a plane, got {len(atoms)}"
            )

        # Calculate centroid
        centroid = atoms.center_of_mass()

        # Center coordinates
        positions = atoms.positions
        centered = positions - centroid

        # Apply SVD
        U, S, Vt = np.linalg.svd(centered, full_matrices=False)

        # Extract normal (3rd row of Vt = 3rd column of V)
        normal = Vt[2, :]

        # Ensure unit vector
        normal = normal / np.linalg.norm(normal)

        return centroid, normal

    @staticmethod
    def validate_plane_fit(
        atoms: mda.AtomGroup,
        plane_centroid: np.ndarray,
        plane_normal: np.ndarray,
        threshold: float = 1.0
    ) -> Tuple[bool, float]:
        """
        Check if atoms are roughly coplanar.

        Parameters
        ----------
        atoms : MDAnalysis.AtomGroup
            Atoms to validate
        plane_centroid : np.ndarray, shape (3,)
            Center point of plane
        plane_normal : np.ndarray, shape (3,)
            Unit normal vector to plane
        threshold : float, default=1.0
            Maximum acceptable RMS deviation in Angstroms

        Returns
        -------
        is_valid : bool
            True if RMS deviation < threshold
        rms_deviation : float
            RMS deviation in Angstroms
        """
        if len(atoms) == 0:
            return False, np.inf

        positions = atoms.positions
        vectors = positions - plane_centroid
        distances = np.abs(np.dot(vectors, plane_normal))
        rms_deviation = np.sqrt(np.mean(distances ** 2))

        is_valid = rms_deviation < threshold

        return is_valid, rms_deviation

    @staticmethod
    def angle_between_vectors(
        v1: np.ndarray,
        v2: np.ndarray,
        degrees: bool = True
    ) -> float:
        """
        Calculate angle between two vectors.

        Parameters
        ----------
        v1 : np.ndarray
            First 3D vector
        v2 : np.ndarray
            Second 3D vector
        degrees : bool, default=True
            Return angle in degrees (True) or radians (False)

        Returns
        -------
        angle : float
            Angle between vectors (0-180 degrees or 0-pi radians)
        """
        v1_norm = v1 / np.linalg.norm(v1)
        v2_norm = v2 / np.linalg.norm(v2)

        cos_theta = np.dot(v1_norm, v2_norm)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)

        angle = np.arccos(cos_theta)

        if degrees:
            angle = np.degrees(angle)

        return float(angle)

    @staticmethod
    def canonicalize_axial_angle(angle: float) -> float:
        """
        Convert a 0-180 degree axial angle into acute-angle form.

        For geometric objects with equivalent opposite directions, such as the groove axis
        and groove normal, theta and 180-theta represent the same axial relationship.
        Normalize to 0-90 degrees to avoid sign-driven pseudo jumps.
        """
        angle = float(angle)
        if angle < 0.0:
            angle = -angle
        if angle > 180.0:
            angle = angle % 180.0
        return min(angle, 180.0 - angle)

    @staticmethod
    def project_vector_onto_plane(
        vector: np.ndarray,
        plane_normal: np.ndarray
    ) -> np.ndarray:
        """
        Project vector onto plane.

        The projection is calculated as: projection = v - (v·n)n
        where n is the normalized plane normal.

        Parameters
        ----------
        vector : np.ndarray
            Vector to project
        plane_normal : np.ndarray
            Plane normal vector

        Returns
        -------
        projected : np.ndarray
            Projected vector
        """
        normal_norm = plane_normal / np.linalg.norm(plane_normal)
        projection_length = np.dot(vector, normal_norm)
        projected = vector - projection_length * normal_norm

        return projected

    @staticmethod
    def principal_axis_from_coords(coords: np.ndarray) -> np.ndarray:
        """
        Extract the first principal-component direction from a coordinate point cloud.

        Parameters
        ----------
        coords : np.ndarray
            Point-cloud coordinates with shape=(N, 3)

        Returns
        -------
        np.ndarray
            Unit principal-axis vector
        """
        centered = coords - coords.mean(axis=0)
        _, _, vt = np.linalg.svd(centered, full_matrices=False)
        axis = vt[0]
        axis = axis / np.linalg.norm(axis)
        return axis

    @staticmethod
    def align_vector_orientation(
        vector: np.ndarray,
        reference: Optional[np.ndarray]
    ) -> np.ndarray:
        """
        Align vector orientation to a reference vector.

        For PCA axes, normal vectors, and other sign-equivalent geometric objects,
        flip the direction when it points opposite to the reference vector.
        """
        if reference is None:
            return vector
        if np.dot(vector, reference) < 0:
            return -vector
        return vector


class _StructureGrooveFallback:
    """
    Structure-based groove fallback detector.

    When reference sequence mapping fails, identify the two groove helices from local
    CA helix geometry within the first 180 residues of the MHC alpha chain.
    """

    HELIX_I_I3_MIN = 4.6
    HELIX_I_I3_MAX = 6.4
    MIN_SEGMENT_LENGTH = 8
    SEARCH_MAX_INDEX = 180
    HELIX1_INDEX_RANGE = (35, 105)
    HELIX2_INDEX_RANGE = (95, 185)
    MAX_AXIS_ANGLE_DEG = 35.0
    MIN_CENTROID_DISTANCE = 8.0
    MAX_CENTROID_DISTANCE = 35.0
    GOLDEN_REFERENCE_PATH = (
        Path('/home/xumy/work/development/ImmunoScope/input/database/parallel2/1OGA_sd_run2/1OGA_sd.pdb')
    )
    GOLDEN_ALIGN_MAX_INDEX = 180
    GOLDEN_HELIX1_RANGE = (50, 86)
    GOLDEN_HELIX2_RANGE = (140, 176)
    _golden_cache: Optional[Dict[str, Any]] = None

    def __init__(self, universe: mda.Universe, mhc_chain_selection: str):
        self.universe = universe
        self.mhc_chain_selection = mhc_chain_selection

    def get_alpha1_alpha2_residues(self) -> Tuple[List[int], List[int], Dict]:
        residues, coords = self._extract_ca_trace()
        helix_mask = self._build_helix_like_mask(coords)
        segments = self._collect_segments(residues, coords, helix_mask)
        helix1, helix2, pair_score = self._pick_best_helix_pair(segments)

        if helix1 is None or helix2 is None:
            raise ValueError(
                "Structure-based groove fallback failed to identify both helix segments."
            )

        alpha1_resids = [resid for _, resid in helix1['residues']]
        alpha2_resids = [resid for _, resid in helix2['residues']]

        quality = {
            'method': 'structure_fallback',
            'identity': None,
            'score': None,
            'alpha1_gaps': [],
            'alpha2_gaps': [],
            'alpha1_gap_fraction': 0.0,
            'alpha2_gap_fraction': 0.0,
            'alpha1_mapped_count': len(alpha1_resids),
            'alpha2_mapped_count': len(alpha2_resids),
            'helix1_index_range': [helix1['start_index'], helix1['end_index']],
            'helix2_index_range': [helix2['start_index'], helix2['end_index']],
            'helix1_score': helix1['score'],
            'helix2_score': helix2['score'],
            'pair_score': pair_score,
            'helix_axis_angle_deg': self._parallel_angle_deg(helix1['axis'], helix2['axis']),
            'helix_centroid_distance': float(np.linalg.norm(helix2['centroid'] - helix1['centroid'])),
            'golden_reference_used': self._golden_reference_available(),
        }

        return alpha1_resids, alpha2_resids, quality

    def _extract_ca_trace(self) -> Tuple[List[Tuple[int, int]], np.ndarray]:
        mhc_chain = self.universe.select_atoms(f"({self.mhc_chain_selection}) and protein and name CA")
        if len(mhc_chain) == 0:
            raise ValueError(f"No CA atoms found for MHC selection: {self.mhc_chain_selection}")

        residues = []
        coords = []
        seen_resids = set()
        for atom in mhc_chain:
            if atom.resid in seen_resids:
                continue
            seen_resids.add(atom.resid)
            residues.append((len(residues), int(atom.resid)))
            coords.append(atom.position.copy())
            if len(residues) >= self.SEARCH_MAX_INDEX:
                break

        if len(coords) < self.MIN_SEGMENT_LENGTH + 3:
            raise ValueError("Insufficient MHC CA trace for structure fallback")

        return residues, np.asarray(coords, dtype=float)

    def _build_helix_like_mask(self, coords: np.ndarray) -> np.ndarray:
        mask = np.zeros(len(coords), dtype=bool)
        for i in range(len(coords) - 3):
            dist_i3 = np.linalg.norm(coords[i + 3] - coords[i])
            if self.HELIX_I_I3_MIN <= dist_i3 <= self.HELIX_I_I3_MAX:
                mask[i:i + 4] = True
        return mask

    def _collect_segments(
        self,
        residues: List[Tuple[int, int]],
        coords: np.ndarray,
        helix_mask: np.ndarray
    ) -> List[Dict[str, Any]]:
        segments = []
        start = None
        for i, is_helix in enumerate(helix_mask):
            if is_helix and start is None:
                start = i
            elif not is_helix and start is not None:
                self._append_segment(segments, residues, coords, start, i - 1)
                start = None
        if start is not None:
            self._append_segment(segments, residues, coords, start, len(helix_mask) - 1)
        return segments

    def _append_segment(
        self,
        segments: List[Dict[str, Any]],
        residues: List[Tuple[int, int]],
        coords: np.ndarray,
        start: int,
        end: int
    ) -> None:
        length = end - start + 1
        if length < self.MIN_SEGMENT_LENGTH:
            return
        segment_residues = residues[start:end + 1]
        segment_coords = coords[start:end + 1]
        axis = self._principal_axis(segment_coords)
        centroid = segment_coords.mean(axis=0)
        linearity = self._linearity_score(segment_coords, axis, centroid)
        score = float(length) * (1.0 + linearity)
        segments.append(
            {
                'start_index': segment_residues[0][0] + 1,
                'end_index': segment_residues[-1][0] + 1,
                'residues': segment_residues,
                'coords': segment_coords,
                'axis': axis,
                'centroid': centroid,
                'linearity': linearity,
                'score': score,
            }
        )

    def _pick_best_helix_pair(
        self,
        segments: List[Dict[str, Any]]
    ) -> Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]], Optional[float]]:
        helix1_candidates = self._filter_candidates(segments, self.HELIX1_INDEX_RANGE)
        helix2_candidates = self._filter_candidates(segments, self.HELIX2_INDEX_RANGE)

        best_pair = (None, None, None)
        best_score = -np.inf

        for helix1 in helix1_candidates:
            for helix2 in helix2_candidates:
                if helix2['start_index'] <= helix1['end_index']:
                    continue
                axis_angle = self._parallel_angle_deg(helix1['axis'], helix2['axis'])
                centroid_distance = float(np.linalg.norm(helix2['centroid'] - helix1['centroid']))
                pair_score = self._score_pair(helix1, helix2, axis_angle, centroid_distance)
                if pair_score > best_score:
                    best_score = pair_score
                    best_pair = (helix1, helix2, pair_score)

        return best_pair

    def _filter_candidates(
        self,
        segments: List[Dict[str, Any]],
        index_range: Tuple[int, int]
    ) -> List[Dict[str, Any]]:
        start_limit, end_limit = index_range
        return [
            seg for seg in segments
            if seg['start_index'] >= start_limit and seg['end_index'] <= end_limit
        ]

    def _score_pair(
        self,
        helix1: Dict[str, Any],
        helix2: Dict[str, Any],
        axis_angle: float,
        centroid_distance: float
    ) -> float:
        if axis_angle > self.MAX_AXIS_ANGLE_DEG:
            return -np.inf
        if not (self.MIN_CENTROID_DISTANCE <= centroid_distance <= self.MAX_CENTROID_DISTANCE):
            return -np.inf

        length_score = helix1['score'] + helix2['score']
        parallel_bonus = max(0.0, 1.0 - axis_angle / self.MAX_AXIS_ANGLE_DEG) * 20.0
        distance_bonus = max(
            0.0,
            1.0 - abs(centroid_distance - 18.0) / 18.0
        ) * 10.0
        golden_bonus = self._golden_reference_bonus(helix1, helix2)
        return length_score + parallel_bonus + distance_bonus + golden_bonus

    def _principal_axis(self, coords: np.ndarray) -> np.ndarray:
        centered = coords - coords.mean(axis=0)
        _, _, vt = np.linalg.svd(centered, full_matrices=False)
        axis = vt[0]
        start_to_end = coords[-1] - coords[0]
        if np.dot(axis, start_to_end) < 0:
            axis = -axis
        return axis / np.linalg.norm(axis)

    def _linearity_score(self, coords: np.ndarray, axis: np.ndarray, centroid: np.ndarray) -> float:
        projected = np.dot(coords - centroid, axis)
        reconstructed = np.outer(projected, axis) + centroid
        residual = coords - reconstructed
        rms = np.sqrt(np.mean(np.sum(residual ** 2, axis=1)))
        return max(0.0, 1.0 - rms / 2.5)

    def _parallel_angle_deg(self, axis1: np.ndarray, axis2: np.ndarray) -> float:
        cos_theta = np.clip(abs(np.dot(axis1, axis2)), 0.0, 1.0)
        return float(np.degrees(np.arccos(cos_theta)))

    @classmethod
    def _golden_reference_available(cls) -> bool:
        return cls.GOLDEN_REFERENCE_PATH.exists()

    @classmethod
    def _load_golden_reference(cls) -> Optional[Dict[str, Any]]:
        if cls._golden_cache is not None:
            return cls._golden_cache

        if not cls.GOLDEN_REFERENCE_PATH.exists():
            cls._golden_cache = None
            return None

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('golden_reference', str(cls.GOLDEN_REFERENCE_PATH))

        chain_a = None
        for chain in structure.get_chains():
            if chain.id == 'A':
                chain_a = chain
                break

        if chain_a is None:
            cls._golden_cache = None
            return None

        residues = [res for res in chain_a if res.id[0] == ' ' and 'CA' in res]
        if len(residues) < cls.GOLDEN_ALIGN_MAX_INDEX:
            cls._golden_cache = None
            return None

        ca_coords = np.asarray([res['CA'].coord for res in residues[:cls.GOLDEN_ALIGN_MAX_INDEX]], dtype=float)
        helix1_slice = slice(cls.GOLDEN_HELIX1_RANGE[0] - 1, cls.GOLDEN_HELIX1_RANGE[1])
        helix2_slice = slice(cls.GOLDEN_HELIX2_RANGE[0] - 1, cls.GOLDEN_HELIX2_RANGE[1])

        cls._golden_cache = {
            'align_coords': ca_coords,
            'helix1_center': ca_coords[helix1_slice].mean(axis=0),
            'helix2_center': ca_coords[helix2_slice].mean(axis=0),
        }
        return cls._golden_cache

    def _golden_reference_bonus(self, helix1: Dict[str, Any], helix2: Dict[str, Any]) -> float:
        reference = self._load_golden_reference()
        if reference is None:
            return 0.0

        current_trace = self._extract_current_trace_for_alignment()
        if current_trace is None:
            return 0.0

        reference_coords = reference['align_coords']
        common_len = min(len(current_trace), len(reference_coords))
        if common_len < 40:
            return 0.0

        moving = [self._pseudo_atom(coord) for coord in current_trace[:common_len]]
        fixed = [self._pseudo_atom(coord) for coord in reference_coords[:common_len]]
        sup = Superimposer()
        sup.set_atoms(fixed, moving)

        transformed_h1 = np.dot(helix1['centroid'], sup.rotran[0]) + sup.rotran[1]
        transformed_h2 = np.dot(helix2['centroid'], sup.rotran[0]) + sup.rotran[1]

        dist_h1 = float(np.linalg.norm(transformed_h1 - reference['helix1_center']))
        dist_h2 = float(np.linalg.norm(transformed_h2 - reference['helix2_center']))
        score_h1 = max(0.0, 1.0 - dist_h1 / 20.0)
        score_h2 = max(0.0, 1.0 - dist_h2 / 20.0)
        return (score_h1 + score_h2) * 10.0

    def _extract_current_trace_for_alignment(self) -> Optional[np.ndarray]:
        mhc_chain = self.universe.select_atoms(f"({self.mhc_chain_selection}) and protein and name CA")
        if len(mhc_chain) == 0:
            return None
        coords = []
        seen_resids = set()
        for atom in mhc_chain:
            if atom.resid in seen_resids:
                continue
            seen_resids.add(atom.resid)
            coords.append(atom.position.copy())
            if len(coords) >= self.GOLDEN_ALIGN_MAX_INDEX:
                break
        if len(coords) < 40:
            return None
        return np.asarray(coords, dtype=float)

    @staticmethod
    def _pseudo_atom(coord: np.ndarray):
        class _Atom:
            def __init__(self, vector):
                self.coord = np.asarray(vector, dtype=float)
            def get_coord(self):
                return self.coord
        return _Atom(coord)


class _SequenceAligner:
    """
    Internal MHC sequence alignment.

    Aligns PDB MHC sequences to HLA-A*02:01 reference to identify α1/α2 helix regions.
    Handles non-standard PDB numbering via BioPython global alignment.

    This class is a refactored version of mhc_sequence_aligner.py (MHCSequenceAligner).
    """

    # HLA-A*02:01 reference sequence (mature protein, signal peptide cleaved)
    REFERENCE_SEQUENCE = (
        "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAAS"
        "PRMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEAGSHTLQRMYGCDVGPDGR"
        "LLRGHDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAARVAEQLRAYLEGTCVEWLRRYLE"
        "NGKDTLERADPPKTHVTHHPVFDYEATLRCWALSFYPAEITLTWQRDGEDQTQDVELVETRPAGDR"
        "TFQKWAAVVVPSGEEQRYTCHVQHEGLPEPLMLRWKQSSLPTIPIVGIVAGLAVLAVVVIGAVVAA"
        "VMWRRKSSDRKGGSYTQAACSDSAQGSDVSLTACKV"
    )

    ALPHA1_POSITIONS = (50, 86)
    ALPHA2_POSITIONS = (140, 176)

    MIN_ALIGNMENT_SCORE = 150.0
    MIN_SEQUENCE_IDENTITY = 0.65
    MAX_GAP_FRACTION = 0.15

    AA_3TO1 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        'MSE': 'M', 'UNK': 'X',
    }

    def __init__(self, universe: mda.Universe, mhc_chain_selection: str):
        """
        Initialize MHC sequence aligner.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            MDAnalysis Universe containing the structure
        mhc_chain_selection : str
            Selection string for MHC chain (e.g., 'chainID A')
        """
        self.universe = universe
        self.mhc_chain = universe.select_atoms(mhc_chain_selection)

        if len(self.mhc_chain) == 0:
            raise ValueError(f"No atoms found for MHC selection: {mhc_chain_selection}")

    def extract_pdb_sequence(self) -> Tuple[str, Dict[int, int]]:
        """
        Extract amino acid sequence from PDB MHC chain.

        Returns
        -------
        sequence : str
            Single-letter amino acid codes
        residue_map : dict
            Maps sequence index (0-based) to PDB residue ID
        """
        residues = []
        seen_resids = set()

        for atom in self.mhc_chain:
            resid = atom.resid
            if resid not in seen_resids:
                residues.append((resid, atom.resname))
                seen_resids.add(resid)

        residues.sort(key=lambda x: x[0])

        sequence = []
        residue_map = {}

        for seq_idx, (resid, resname) in enumerate(residues):
            aa_code = self.AA_3TO1.get(resname.upper(), 'X')
            sequence.append(aa_code)
            residue_map[seq_idx] = resid

        sequence_str = ''.join(sequence)

        return sequence_str, residue_map

    def align_to_reference(self, pdb_sequence: str) -> Tuple:
        """
        Align PDB sequence to HLA-A*02:01 reference using BioPython.

        Parameters
        ----------
        pdb_sequence : str
            Amino acid sequence extracted from PDB

        Returns
        -------
        alignment : Bio.Align.PairwiseAlignment
            Alignment object (first result from aligner)
        score : float
            Alignment score (BLOSUM62-based)
        identity : float
            Sequence identity (0.0-1.0)
        """
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        alignments = aligner.align(self.REFERENCE_SEQUENCE, pdb_sequence)

        if len(alignments) == 0:
            raise ValueError("No alignment found between reference and PDB sequence")

        best_alignment = alignments[0]
        score = best_alignment.score

        ref_aligned = best_alignment[0]
        pdb_aligned = best_alignment[1]

        matches = sum(1 for r, p in zip(ref_aligned, pdb_aligned) if r == p and r != '-')
        aligned_length = len(ref_aligned)
        identity = matches / aligned_length if aligned_length > 0 else 0.0

        return best_alignment, score, identity

    def map_reference_to_pdb(
        self,
        alignment,
        residue_map: Dict[int, int]
    ) -> Tuple[Dict[int, int], List[int]]:
        """
        Map reference sequence positions to PDB residue IDs.

        Parameters
        ----------
        alignment : Bio.Align.PairwiseAlignment
            Alignment between reference and PDB sequence
        residue_map : dict
            Maps PDB sequence index to PDB residue ID

        Returns
        -------
        position_map : dict
            Maps reference position (1-based) to PDB residue ID
        gaps : list
            Reference positions with no PDB match (gaps)
        """
        ref_aligned = alignment[0]
        pdb_aligned = alignment[1]

        position_map = {}
        gaps = []

        ref_pos = 0
        pdb_pos = 0

        for ref_char, pdb_char in zip(ref_aligned, pdb_aligned):
            if ref_char != '-':
                ref_pos_1based = ref_pos + 1

                if pdb_char != '-':
                    pdb_resid = residue_map[pdb_pos]
                    position_map[ref_pos_1based] = pdb_resid
                    pdb_pos += 1
                else:
                    gaps.append(ref_pos_1based)

                ref_pos += 1
            else:
                if pdb_char != '-':
                    pdb_pos += 1

        return position_map, gaps

    def get_alpha1_alpha2_residues(self) -> Tuple[List[int], List[int], Dict]:
        """
        Get PDB residue IDs corresponding to α1 and α2 helices.

        Full workflow:
        1. Extract PDB sequence
        2. Align to HLA-A*02:01 reference
        3. Map α1 (50-86) and α2 (140-176) positions to PDB residue IDs
        4. Validate alignment quality

        Returns
        -------
        alpha1_resids : list of int
            PDB residue IDs for α1 helix region
        alpha2_resids : list of int
            PDB residue IDs for α2 helix region
        alignment_quality : dict
            Contains 'score', 'identity', 'alpha1_gaps', 'alpha2_gaps', etc.

        Raises
        ------
        ValueError
            If alignment quality is too low or too many gaps in critical regions
        """
        pdb_sequence, residue_map = self.extract_pdb_sequence()
        alignment, score, identity = self.align_to_reference(pdb_sequence)
        position_map, gaps = self.map_reference_to_pdb(alignment, residue_map)

        alpha1_start, alpha1_end = self.ALPHA1_POSITIONS
        alpha2_start, alpha2_end = self.ALPHA2_POSITIONS

        alpha1_resids = []
        alpha1_gaps = []
        for ref_pos in range(alpha1_start, alpha1_end + 1):
            if ref_pos in position_map:
                alpha1_resids.append(position_map[ref_pos])
            else:
                alpha1_gaps.append(ref_pos)

        alpha2_resids = []
        alpha2_gaps = []
        for ref_pos in range(alpha2_start, alpha2_end + 1):
            if ref_pos in position_map:
                alpha2_resids.append(position_map[ref_pos])
            else:
                alpha2_gaps.append(ref_pos)

        alpha1_length = alpha1_end - alpha1_start + 1
        alpha2_length = alpha2_end - alpha2_start + 1

        alpha1_gap_fraction = len(alpha1_gaps) / alpha1_length
        alpha2_gap_fraction = len(alpha2_gaps) / alpha2_length

        alignment_quality = {
            'score': score,
            'identity': identity,
            'alpha1_gaps': alpha1_gaps,
            'alpha2_gaps': alpha2_gaps,
            'alpha1_gap_fraction': alpha1_gap_fraction,
            'alpha2_gap_fraction': alpha2_gap_fraction,
            'alpha1_mapped_count': len(alpha1_resids),
            'alpha2_mapped_count': len(alpha2_resids)
        }

        # Quality checks
        if score < self.MIN_ALIGNMENT_SCORE:
            raise ValueError(
                f"Alignment score too low: {score:.1f} < {self.MIN_ALIGNMENT_SCORE}. "
                f"PDB sequence may not be HLA Class I MHC."
            )

        if identity < self.MIN_SEQUENCE_IDENTITY:
            raise ValueError(
                f"Sequence identity too low: {identity:.2%} < {self.MIN_SEQUENCE_IDENTITY:.0%}. "
                f"PDB sequence diverges significantly from HLA-A*02:01."
            )

        if alpha1_gap_fraction > self.MAX_GAP_FRACTION:
            raise ValueError(
                f"Too many gaps in α1 region: {alpha1_gap_fraction:.1%} > {self.MAX_GAP_FRACTION:.0%}. "
                f"Missing reference positions: {alpha1_gaps}"
            )

        if alpha2_gap_fraction > self.MAX_GAP_FRACTION:
            raise ValueError(
                f"Too many gaps in α2 region: {alpha2_gap_fraction:.1%} > {self.MAX_GAP_FRACTION:.0%}. "
                f"Missing reference positions: {alpha2_gaps}"
            )

        if len(alpha1_resids) < 20:
            raise ValueError(
                f"Too few α1 residues mapped: {len(alpha1_resids)} < 20. "
                f"Alignment may be unreliable."
            )

        if len(alpha2_resids) < 20:
            raise ValueError(
                f"Too few α2 residues mapped: {len(alpha2_resids)} < 20. "
                f"Alignment may be unreliable."
            )

        return alpha1_resids, alpha2_resids, alignment_quality


class _MHCGrooveDetector:
    """
    Internal MHC groove geometry detector.

    Calculates MHC groove axis and plane using sequence-aligned α1/α2 residues.
    Combines sequence alignment with geometric calculations.

    This class is a refactored version of mhc_groove_detector.py (MHCGrooveDetector).
    """

    def __init__(self, universe: mda.Universe, mhc_chain_selection: str):
        """
        Initialize MHC groove detector.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            MDAnalysis Universe containing the structure
        mhc_chain_selection : str
            Selection string for MHC chain (e.g., 'chainID A')
        """
        self.universe = universe
        self.mhc_chain_selection = mhc_chain_selection
        self.aligner = _SequenceAligner(universe, mhc_chain_selection)

        # Cache alignment results
        self._alpha1_resids = None
        self._alpha2_resids = None
        self._alignment_quality = None

    def _ensure_alignment(self):
        """Perform alignment if not already done."""
        if self._alpha1_resids is None:
            try:
                (
                    self._alpha1_resids,
                    self._alpha2_resids,
                    self._alignment_quality
                ) = self.aligner.get_alpha1_alpha2_residues()
            except Exception as exc:
                logger.warning(
                    "Sequence-based MHC groove mapping failed, switching to structure fallback: %s",
                    exc,
                )
                fallback = _StructureGrooveFallback(self.universe, self.mhc_chain_selection)
                (
                    self._alpha1_resids,
                    self._alpha2_resids,
                    self._alignment_quality
                ) = fallback.get_alpha1_alpha2_residues()
                self._alignment_quality['fallback_reason'] = str(exc)

    def get_alignment_quality(self) -> Dict:
        """
        Get alignment quality metrics.

        Returns
        -------
        alignment_quality : dict
            Contains alignment score, identity, gaps, etc.
        """
        self._ensure_alignment()
        return self._alignment_quality

    def get_reference_resids(self) -> Tuple[list, list]:
        """
        Get fixed α1/α2 residue IDs for reference structure.

        Should be called once to determine which residues to use for
        groove geometry calculation across all frames.

        Returns
        -------
        alpha1_resids : list
            Residue IDs for α1 helix
        alpha2_resids : list
            Residue IDs for α2 helix
        """
        self._ensure_alignment()
        return self._alpha1_resids, self._alpha2_resids

    def detect_groove_axis(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict]:
        """
        Calculate groove axis using sequence-aligned α1/α2 positions.

        Workflow:
        1. Get α1/α2 residue IDs from sequence alignment
        2. Select Cα atoms
        3. Calculate COM(α1) and COM(α2)
        4. Groove axis = COM(α2) - COM(α1)

        Returns
        -------
        groove_axis : np.ndarray, shape (3,)
            Unit vector from N-end (α1) to C-end (α2)
        com_n_end : np.ndarray, shape (3,)
            Center of mass of α1 region
        com_c_end : np.ndarray, shape (3,)
            Center of mass of α2 region
        alignment_info : dict
            Alignment quality metrics
        """
        self._ensure_alignment()

        alpha1_resid_str = " ".join(map(str, self._alpha1_resids))
        alpha2_resid_str = " ".join(map(str, self._alpha2_resids))

        alpha1_ca = self.universe.select_atoms(
            f"({self.mhc_chain_selection}) and (resid {alpha1_resid_str}) and name CA"
        )
        alpha2_ca = self.universe.select_atoms(
            f"({self.mhc_chain_selection}) and (resid {alpha2_resid_str}) and name CA"
        )

        if len(alpha1_ca) == 0:
            raise ValueError(
                f"No Cα atoms found for α1 region (resids: {self._alpha1_resids[:5]}...)"
            )

        if len(alpha2_ca) == 0:
            raise ValueError(
                f"No Cα atoms found for α2 region (resids: {self._alpha2_resids[:5]}...)"
            )

        com_n_end = alpha1_ca.center_of_mass()
        com_c_end = alpha2_ca.center_of_mass()

        groove_axis_raw = com_c_end - com_n_end
        groove_axis = groove_axis_raw / np.linalg.norm(groove_axis_raw)

        return groove_axis, com_n_end, com_c_end, self._alignment_quality

    def detect_groove_plane(self) -> Tuple[np.ndarray, np.ndarray, float, Dict]:
        """
        Fit groove plane using α1/α2 backbone atoms.

        Workflow:
        1. Get α1/α2 residue IDs from sequence alignment
        2. Select backbone atoms (CA, C, N, O)
        3. Use SVD fitting for plane
        4. Validate fit quality

        Returns
        -------
        plane_centroid : np.ndarray, shape (3,)
            Center point of plane
        plane_normal : np.ndarray, shape (3,)
            Unit normal vector to the plane
        fit_quality : float
            RMS deviation of atoms from plane (Angstroms)
        alignment_info : dict
            Alignment quality metrics
        """
        self._ensure_alignment()

        all_resids = sorted(self._alpha1_resids + self._alpha2_resids)
        resid_str = " ".join(map(str, all_resids))

        backbone = self.universe.select_atoms(
            f"({self.mhc_chain_selection}) and "
            f"(resid {resid_str}) and "
            f"(name CA or name C or name N or name O)"
        )

        if len(backbone) == 0:
            raise ValueError(
                f"No backbone atoms found for α1/α2 regions "
                f"(resids: {all_resids[:5]}...{all_resids[-5:]})"
            )

        if len(backbone) < 3:
            raise ValueError(
                f"Insufficient backbone atoms for plane fit: {len(backbone)} < 3"
            )

        plane_centroid, plane_normal = _GeometryUtils.fit_plane_svd(backbone)

        is_valid, rms_deviation = _GeometryUtils.validate_plane_fit(
            backbone, plane_centroid, plane_normal, threshold=10.0
        )

        if not is_valid:
            import warnings
            warnings.warn(
                f"MHC groove plane fit quality is very poor: RMS = {rms_deviation:.2f} Å > 10.0 Å. "
                f"This may indicate severe structural distortion."
            )

        return plane_centroid, plane_normal, rms_deviation, self._alignment_quality


class _TCRAxisCalculator:
    """
    Internal TCR axis calculator using conserved cysteines.

    Uses ANARCI IMGT numbering to identify conserved disulfide bond cysteines
    (positions 23 and 104) in TCR V domains, with intelligent neighbor search fallback.

    This class is a refactored version of conserved_cysteine_detector.py (ConservedCysteineDetector).
    """

    IMGT_CYS_N_TERMINAL = 23
    IMGT_CYS_C_TERMINAL = 104

    DISULFIDE_MAX_CA_DISTANCE = 10.0
    DISULFIDE_MAX_SG_DISTANCE = 3.0

    CYS_SEARCH_RANGE = 20

    def __init__(self, universe: mda.Universe):
        """
        Initialize detector.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            MDAnalysis universe with loaded structure
        """
        self.universe = universe

        from immunoscope.analysis.topology import ANARCIWrapper
        self.anarci_wrapper = ANARCIWrapper(
            allow_fallback=True,  # Use distance fallback when ANARCI is unavailable.
            numbering_scheme="imgt"
        )

        # Cache for reference residue IDs (set once for trajectory analysis)
        self._reference_cys_resids = {}

    def detect_conserved_cysteines(
        self,
        chain_selection: str,
        chain_type: str = "TCR"
    ) -> Dict[str, Any]:
        """
        Detect conserved cysteines for a TCR chain using ANARCI with neighbor search.

        Strategy:
        1. Run ANARCI to get approximate IMGT positions 23 and 104
        2. Check if located positions are CYS; if not, search for nearest CYS
        3. Validate discovered CYS pair by CA-CA and SG-SG distances

        Parameters
        ----------
        chain_selection : str
            MDAnalysis selection string for the chain (e.g., 'chainID D')
        chain_type : str, default='TCR'
            Chain type identifier

        Returns
        -------
        result : dict
            Detection results with keys:
            - 'cys_23_resid', 'cys_104_resid'
            - 'cys_23_ca_selection', 'cys_104_ca_selection', 'cys_pair_selection'
            - 'ca_distance', 'sg_distance'
            - 'disulfide_validated'
            - 'numbering_scheme', 'chain_type_detected'
            - 'search_offset_23', 'search_offset_104'

        Raises
        ------
        ValueError
            If conserved cysteines cannot be found or validated
        """
        logger.info(f"Detecting conserved cysteines for selection: {chain_selection}")

        chain_atoms = self.universe.select_atoms(chain_selection)
        if len(chain_atoms) == 0:
            raise ValueError(f"No atoms found for selection: {chain_selection}")

        sequence = self._extract_sequence(chain_atoms)
        logger.info(f"Chain sequence: {len(sequence)} residues")

        anarci_result = self.anarci_wrapper.run_anarci(sequence, chain_type)

        if anarci_result['method'] == 'regex_fallback':
            logger.warning("ANARCI not available, using distance-based fallback")
            return self._distance_based_fallback(chain_atoms, chain_selection)

        numbering = anarci_result.get('numbering')
        if not numbering:
            logger.warning("ANARCI returned no numbering, using distance-based fallback")
            return self._distance_based_fallback(chain_atoms, chain_selection)

        cys_23_info = self._find_imgt_position(numbering, self.IMGT_CYS_N_TERMINAL)
        cys_104_info = self._find_imgt_position(numbering, self.IMGT_CYS_C_TERMINAL)

        if not cys_23_info or not cys_104_info:
            logger.warning("IMGT positions not found, using distance-based fallback")
            return self._distance_based_fallback(chain_atoms, chain_selection)

        ca_atoms = chain_atoms.select_atoms("name CA")

        seq_idx_23, aa_23 = cys_23_info
        seq_idx_104, aa_104 = cys_104_info

        imgt_23_resid = int(ca_atoms[seq_idx_23].resid)
        imgt_104_resid = int(ca_atoms[seq_idx_104].resid)

        logger.info(f"ANARCI located IMGT-23 → residue {imgt_23_resid} ({aa_23}), "
                   f"IMGT-104 → residue {imgt_104_resid} ({aa_104})")

        cys_23_resid, search_offset_23 = self._find_cysteine_near_position(
            chain_atoms, imgt_23_resid, "IMGT-23"
        )
        cys_104_resid, search_offset_104 = self._find_cysteine_near_position(
            chain_atoms, imgt_104_resid, "IMGT-104"
        )

        if not cys_23_resid or not cys_104_resid:
            raise ValueError("Cannot find conserved cysteines near IMGT positions")

        cys_23_atoms = chain_atoms.select_atoms(f"resid {cys_23_resid}")
        cys_104_atoms = chain_atoms.select_atoms(f"resid {cys_104_resid}")

        cys_23_ca = cys_23_atoms.select_atoms("name CA")
        cys_104_ca = cys_104_atoms.select_atoms("name CA")

        if len(cys_23_ca) == 0 or len(cys_104_ca) == 0:
            raise ValueError("Cannot find CA atoms for discovered cysteines")

        ca_distance = np.linalg.norm(cys_23_ca[0].position - cys_104_ca[0].position)

        cys_23_sg = cys_23_atoms.select_atoms("name SG")
        cys_104_sg = cys_104_atoms.select_atoms("name SG")

        sg_distance = None
        disulfide_validated = False

        if len(cys_23_sg) > 0 and len(cys_104_sg) > 0:
            sg_distance = np.linalg.norm(cys_23_sg[0].position - cys_104_sg[0].position)
            disulfide_validated = (
                ca_distance <= self.DISULFIDE_MAX_CA_DISTANCE and
                sg_distance <= self.DISULFIDE_MAX_SG_DISTANCE
            )
            logger.info(f"Disulfide validation: CA distance={ca_distance:.2f} Å, "
                       f"SG distance={sg_distance:.2f} Å")

            if disulfide_validated:
                logger.info(f"✓ Disulfide bond validated: {cys_23_resid}-{cys_104_resid}")
            else:
                logger.warning(
                    f"Distance validation failed: CA={ca_distance:.2f} Å (max {self.DISULFIDE_MAX_CA_DISTANCE}), "
                    f"SG={sg_distance:.2f} Å (max {self.DISULFIDE_MAX_SG_DISTANCE})"
                )
        else:
            logger.warning("Missing SG atoms, using CA distance only")
            disulfide_validated = ca_distance <= self.DISULFIDE_MAX_CA_DISTANCE

            if not disulfide_validated:
                logger.warning(
                    f"CA distance {ca_distance:.2f} Å exceeds maximum {self.DISULFIDE_MAX_CA_DISTANCE} Å. "
                    "This may indicate incorrect detection or unusual structure."
                )

        chain_id = self._extract_chain_id(chain_selection, chain_atoms)

        cys_23_selection = f"{chain_id} and resid {cys_23_resid} and name CA"
        cys_104_selection = f"{chain_id} and resid {cys_104_resid} and name CA"

        result = {
            'cys_23_resid': cys_23_resid,
            'cys_104_resid': cys_104_resid,
            'cys_23_ca_selection': cys_23_selection,
            'cys_104_ca_selection': cys_104_selection,
            'cys_pair_selection': f"{chain_id} and resid {cys_23_resid} {cys_104_resid} and name CA",
            'ca_distance': float(ca_distance),
            'sg_distance': float(sg_distance) if sg_distance else None,
            'disulfide_validated': disulfide_validated,
            'numbering_scheme': 'imgt_with_neighbor_search',
            'chain_type_detected': anarci_result.get('chain_type', 'unknown'),
            'imgt_positions': (self.IMGT_CYS_N_TERMINAL, self.IMGT_CYS_C_TERMINAL),
            'search_offset_23': search_offset_23,
            'search_offset_104': search_offset_104
        }

        logger.info(
            f"✓ Detection complete: residue {cys_23_resid} (offset {search_offset_23:+d}) - "
            f"residue {cys_104_resid} (offset {search_offset_104:+d})"
        )

        return result

    def _extract_sequence(self, atoms: mda.AtomGroup) -> str:
        """Extract amino acid sequence from atom group."""
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
            'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
            'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
            'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }

        ca_atoms = atoms.select_atoms("name CA")
        sequence = ''.join([aa_map.get(res, 'X') for res in ca_atoms.resnames])

        return sequence

    def _find_imgt_position(
        self,
        numbering: list,
        target_position: int
    ) -> Optional[Tuple[int, str]]:
        """Find sequence index for a given IMGT position."""
        for seq_idx, item in enumerate(numbering):
            if isinstance(item, tuple) and len(item) == 2:
                pos_info, aa = item
                if isinstance(pos_info, tuple) and len(pos_info) == 2:
                    pos, insertion = pos_info
                    if pos == target_position and aa != '-':
                        return (seq_idx, aa)

        return None

    def _find_cysteine_near_position(
        self,
        chain_atoms: mda.AtomGroup,
        target_resid: int,
        position_name: str
    ) -> Tuple[Optional[int], int]:
        """
        Find CYS near target position (IMGT-defined).

        Strategy:
        1. Check if target position itself is CYS (return immediately)
        2. If not, search for nearest CYS within ±CYS_SEARCH_RANGE residues
        3. Return the discovered CYS and offset from target
        """
        target_res = chain_atoms.select_atoms(f"resid {target_resid}")
        if len(target_res.residues) > 0:
            target_resname = target_res.residues[0].resname
            if target_resname == 'CYS':
                logger.info(f"  {position_name} (residue {target_resid}): ✓ already CYS")
                return target_resid, 0

        logger.info(f"  {position_name} (residue {target_resid}): {target_resname} ≠ CYS, "
                   f"searching within ±{self.CYS_SEARCH_RANGE} residues...")

        all_cys = chain_atoms.select_atoms("resname CYS")
        if len(all_cys.residues) == 0:
            logger.error(f"  No CYS residues found in chain")
            return None, 0

        min_diff = float('inf')
        nearest_cys_resid = None

        for cys_res in all_cys.residues:
            diff = abs(cys_res.resid - target_resid)
            if diff < min_diff and diff <= self.CYS_SEARCH_RANGE:
                min_diff = diff
                nearest_cys_resid = cys_res.resid

        if nearest_cys_resid:
            offset = nearest_cys_resid - target_resid
            logger.info(f"  ✓ Found nearest CYS: residue {nearest_cys_resid} (offset {offset:+d})")
            return nearest_cys_resid, offset
        else:
            logger.error(f"  ✗ No CYS found within ±{self.CYS_SEARCH_RANGE} residues of {target_resid}")
            return None, 0

    def _distance_based_fallback(
        self,
        chain_atoms: mda.AtomGroup,
        chain_selection: str
    ) -> Dict:
        """Complete distance-based fallback when ANARCI fails."""
        logger.info("Using complete distance-based fallback (ANARCI unavailable)")

        all_cys = chain_atoms.select_atoms("resname CYS")
        if len(all_cys.residues) < 2:
            raise ValueError(f"Insufficient CYS residues for fallback (found {len(all_cys.residues)}, need ≥2)")

        min_ca_dist = float('inf')
        best_pair = None

        for i, res1 in enumerate(all_cys.residues):
            for j, res2 in enumerate(all_cys.residues):
                if i < j:
                    ca1 = res1.atoms.select_atoms("name CA")
                    ca2 = res2.atoms.select_atoms("name CA")

                    if len(ca1) > 0 and len(ca2) > 0:
                        dist = np.linalg.norm(ca1[0].position - ca2[0].position)
                        if dist < min_ca_dist:
                            min_ca_dist = dist
                            best_pair = (res1, res2)

        if not best_pair:
            raise ValueError("No valid CYS pair found in fallback")

        res1, res2 = best_pair
        cys_23_resid = min(res1.resid, res2.resid)
        cys_104_resid = max(res1.resid, res2.resid)

        logger.info(f"Fallback discovered CYS pair: {cys_23_resid}-{cys_104_resid}, "
                   f"CA distance={min_ca_dist:.2f} Å")

        sg1 = res1.atoms.select_atoms("name SG")
        sg2 = res2.atoms.select_atoms("name SG")

        sg_distance = None
        if len(sg1) > 0 and len(sg2) > 0:
            sg_distance = np.linalg.norm(sg1[0].position - sg2[0].position)
            logger.info(f"  SG distance: {sg_distance:.2f} Å")

        disulfide_validated = (
            min_ca_dist <= self.DISULFIDE_MAX_CA_DISTANCE and
            (sg_distance is None or sg_distance <= self.DISULFIDE_MAX_SG_DISTANCE)
        )

        chain_id = self._extract_chain_id(chain_selection, chain_atoms)

        result = {
            'cys_23_resid': cys_23_resid,
            'cys_104_resid': cys_104_resid,
            'cys_23_ca_selection': f"{chain_id} and resid {cys_23_resid} and name CA",
            'cys_104_ca_selection': f"{chain_id} and resid {cys_104_resid} and name CA",
            'cys_pair_selection': f"{chain_id} and resid {cys_23_resid} {cys_104_resid} and name CA",
            'ca_distance': float(min_ca_dist),
            'sg_distance': float(sg_distance) if sg_distance else None,
            'disulfide_validated': disulfide_validated,
            'numbering_scheme': 'distance_based_fallback',
            'chain_type_detected': 'unknown',
            'search_offset_23': 0,
            'search_offset_104': 0
        }

        logger.info(f"✓ Fallback complete: {cys_23_resid}-{cys_104_resid}")

        return result

    def _extract_chain_id(self, selection: str, atoms: mda.AtomGroup) -> str:
        """Extract chain identifier from selection string."""
        if 'chainID' in selection:
            import re
            match = re.search(r'chainID\s+([A-Z])', selection)
            if match:
                return f"chainID {match.group(1)}"
        elif 'segname' in selection:
            import re
            match = re.search(r'segname\s+(\w+)', selection)
            if match:
                return f"segname {match.group(1)}"

        first_atom = atoms[0]
        if hasattr(first_atom, 'chainID') and first_atom.chainID:
            return f"chainID {first_atom.chainID}"
        elif hasattr(first_atom, 'segid') and first_atom.segid:
            return f"segname {first_atom.segid}"

        logger.warning("Could not extract chain ID, using generic selection")
        return selection

    def detect_reference_cysteines(
        self,
        chain_selection: str,
        chain_type: str = "TCR"
    ) -> tuple:
        """
        Detect conserved cysteines ONCE for reference structure.

        This method should be called only once (e.g., on first frame) to determine
        which residue IDs correspond to conserved cysteines. These IDs are then
        cached and reused for all subsequent frames.

        Parameters
        ----------
        chain_selection : str
            MDAnalysis selection string for the chain
        chain_type : str
            Chain type identifier

        Returns
        -------
        cys_23_resid : int
            Residue ID for N-terminal conserved cysteine
        cys_104_resid : int
            Residue ID for C-terminal conserved cysteine
        chain_id_str : str
            Chain identifier string (e.g., 'chainID D')
        """
        result = self.detect_conserved_cysteines(chain_selection, chain_type)
        return (
            result['cys_23_resid'],
            result['cys_104_resid'],
            self._extract_chain_id(chain_selection, self.universe.select_atoms(chain_selection))
        )

    def get_cys_selection_from_resids(
        self,
        cys_23_resid: int,
        cys_104_resid: int,
        chain_id_str: str
    ) -> str:
        """
        Build CA selection string from fixed residue IDs.

        Parameters
        ----------
        cys_23_resid : int
            N-terminal cysteine residue ID
        cys_104_resid : int
            C-terminal cysteine residue ID
        chain_id_str : str
            Chain identifier (e.g., 'chainID D')

        Returns
        -------
        selection : str
            MDAnalysis selection string for CA atoms of both cysteines
        """
        return (
            f"({chain_id_str} and resid {cys_23_resid} and name CA) or "
            f"({chain_id_str} and resid {cys_104_resid} and name CA)"
        )

    def detect_vonly_structure(
        self,
        alpha_selection: str,
        beta_selection: str
    ) -> bool:
        """
        Detect if TCR structure contains only V-domain (no C-domain).

        V-only structures typically have ~110-120 residues per chain.
        Full TCR extracellular domains have ~190-250 residues per chain.

        Parameters
        ----------
        alpha_selection : str
            Selection string for TCR alpha chain
        beta_selection : str
            Selection string for TCR beta chain

        Returns
        -------
        is_vonly : bool
            True if structure appears to be V-domain only
        """
        alpha_atoms = self.universe.select_atoms(alpha_selection)
        beta_atoms = self.universe.select_atoms(beta_selection)

        alpha_length = len(set(alpha_atoms.resids))
        beta_length = len(set(beta_atoms.resids))

        # Threshold: <150 residues indicates V-only
        is_vonly = (alpha_length < 150) or (beta_length < 150)

        if is_vonly:
            logger.info(
                f"V-only TCR structure detected: "
                f"alpha={alpha_length} aa, beta={beta_length} aa"
            )

        return is_vonly

    def get_com_based_tcr_axis(
        self,
        alpha_selection: str,
        beta_selection: str,
        alpha_cys_23: int,
        alpha_cys_104: int,
        beta_cys_23: int,
        beta_cys_104: int
    ) -> np.ndarray:
        """
        Calculate TCR axis using center-of-mass vector between Vα and Vβ.

        This method is more stable for V-only structures where PCA is unreliable.
        Uses COM(Vα) - COM(Vβ) as the TCR direction vector (pointing from β to α).

        Parameters
        ----------
        alpha_selection : str
            Selection string for TCR alpha chain
        beta_selection : str
            Selection string for TCR beta chain
        alpha_cys_23, alpha_cys_104 : int
            Conserved cysteine residue IDs for alpha chain
        beta_cys_23, beta_cys_104 : int
            Conserved cysteine residue IDs for beta chain

        Returns
        -------
        axis : np.ndarray
            Normalized TCR axis vector
        """
        alpha_atoms = self.universe.select_atoms(alpha_selection)
        beta_atoms = self.universe.select_atoms(beta_selection)

        # Select V-domain regions
        alpha_vd = alpha_atoms.select_atoms(
            f"resid {min(alpha_cys_23, alpha_cys_104)}:{max(alpha_cys_23, alpha_cys_104)}"
        )
        beta_vd = beta_atoms.select_atoms(
            f"resid {min(beta_cys_23, beta_cys_104)}:{max(beta_cys_23, beta_cys_104)}"
        )

        if len(alpha_vd) == 0 or len(beta_vd) == 0:
            raise ValueError("Cannot find V-domain atoms for COM calculation")

        alpha_com = alpha_vd.center_of_mass()
        beta_com = beta_vd.center_of_mass()

        # Try alpha - beta direction (pointing from β to α)
        axis_raw = alpha_com - beta_com
        axis = axis_raw / np.linalg.norm(axis_raw)

        return axis

    def get_vdomain_pca_axis(
        self,
        chain_selection: str,
        cys_23_resid: int,
        cys_104_resid: int,
        ref_axis: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Define TCR orientation using the first PCA axis of all V-domain CA atoms.

        Compared with a two-point Cys-Cys vector, this is less sensitive to local
        thermal noise and better suited to per-frame trajectory analysis.
        """
        chain_atoms = self.universe.select_atoms(chain_selection)
        if len(chain_atoms) == 0:
            raise ValueError(f"No atoms found for selection: {chain_selection}")

        resid_lo = min(cys_23_resid, cys_104_resid)
        resid_hi = max(cys_23_resid, cys_104_resid)
        vdomain_ca = chain_atoms.select_atoms(f"name CA and resid {resid_lo}:{resid_hi}")

        if len(vdomain_ca) < 10:
            logger.warning(
                f"V-domain Calpha count too low ({len(vdomain_ca)}); "
                "falling back to two-point Cys vector."
            )
            cys_23_ca = chain_atoms.select_atoms(f"resid {cys_23_resid} and name CA")
            cys_104_ca = chain_atoms.select_atoms(f"resid {cys_104_resid} and name CA")
            if len(cys_23_ca) == 0 or len(cys_104_ca) == 0:
                raise ValueError("Cannot find Cys CA atoms for fallback axis calculation.")
            raw = cys_104_ca.positions[0] - cys_23_ca.positions[0]
            axis = raw / np.linalg.norm(raw)
            if ref_axis is not None and np.dot(axis, ref_axis) < 0:
                axis = -axis
            return axis

        coords = vdomain_ca.positions
        centroid = coords.mean(axis=0)
        centered = coords - centroid
        _, _, vt = np.linalg.svd(centered, full_matrices=False)
        axis = vt[0]

        if ref_axis is not None:
            if np.dot(axis, ref_axis) < 0:
                axis = -axis
        else:
            cys_23_ca = chain_atoms.select_atoms(f"resid {cys_23_resid} and name CA")
            cys_104_ca = chain_atoms.select_atoms(f"resid {cys_104_resid} and name CA")
            if len(cys_23_ca) > 0 and len(cys_104_ca) > 0:
                cys_vec = cys_104_ca.positions[0] - cys_23_ca.positions[0]
                if np.dot(axis, cys_vec) < 0:
                    axis = -axis

        return axis / np.linalg.norm(axis)

    def get_robust_tcr_axis(
        self,
        alpha_selection: str,
        beta_selection: str,
        groove_com: np.ndarray
    ) -> np.ndarray:
        """
        Calculate TCR axis using robust COM-based method.

        This method calculates the TCR axis as the vector from MHC groove COM
        to TCR overall COM. This approach:
        - Does not depend on conserved Cys residue detection
        - Is not affected by PCA direction ambiguity
        - Works for both V-only and full-length TCR structures
        - Has clear biological meaning: TCR orientation relative to MHC

        Parameters
        ----------
        alpha_selection : str
            Selection string for TCR alpha chain
        beta_selection : str
            Selection string for TCR beta chain
        groove_com : np.ndarray
            Center of mass of MHC groove (α1 + α2 helices)

        Returns
        -------
        axis : np.ndarray
            Normalized TCR axis vector (pointing from groove to TCR)
        """
        # Get TCR overall COM
        tcr_atoms = self.universe.select_atoms(f"({alpha_selection}) or ({beta_selection})")
        if len(tcr_atoms) == 0:
            raise ValueError("No TCR atoms found")

        tcr_com = tcr_atoms.center_of_mass()

        # Calculate axis: from groove to TCR
        axis_raw = tcr_com - groove_com
        axis = axis_raw / np.linalg.norm(axis_raw)

        return axis


# ============================================================================
# Public API
# ============================================================================

class DockingAngleAnalyzer:
    """
    TCR-pMHC Docking Angle Analyzer.

    Calculates Crossing and Incident angles for TCR-pMHC docking geometry.
    Conforms to ImmunoScope 6 design principles:
    1. Clear Inputs - DockingAngleInput dataclass
    2. Clear Outputs - DockingAngleResult dataclass
    3. Clear Side Effects - Track all output files
    4. Clear Errors - Structured error handling
    5. Testable - Input validation + mock-friendly
    6. Schedulable - Progress callbacks + cancellation

    Angle definitions:
    - Crossing angle: angle(TCR_axis, groove_axis)
    - Incident angle: angle(TCR_axis, groove_normal)

    TCR axis: COM(Vβ disulfide) - COM(Vα disulfide) via ANARCI detection
    Groove axis: COM(α2) - COM(α1) via sequence alignment
    Groove normal: SVD plane fit of α1/α2 backbone

    Examples
    --------
    Single frame analysis:
    >>> from immunoscope.analysis.angles import DockingAngleAnalyzer, DockingAngleInput
    >>> analyzer = DockingAngleAnalyzer()
    >>> result = analyzer.analyze(DockingAngleInput(
    ...     topology='structure.pdb'
    ... ))
    >>> print(result.get_summary())

    Trajectory analysis:
    >>> result = analyzer.analyze(DockingAngleInput(
    ...     topology='md.tpr',
    ...     trajectory='md_pbc.xtc',
    ...     stride=10,
    ...     output_dir='./results'
    ... ))
    >>> print(f"Crossing: {result.statistics['crossing_mean']:.2f}°")

    With progress tracking:
    >>> def progress_callback(progress, message):
    ...     print(f"[{progress*100:.0f}%] {message}")
    >>> analyzer.set_progress_callback(progress_callback)
    >>> result = analyzer.analyze(input_params)
    """

    def __init__(self):
        """Initialize analyzer."""
        self._cancel_flag = False
        self.progress_callback: Optional[Callable[[float, str], None]] = None

    def set_progress_callback(self, callback: Callable[[float, str], None]):
        """
        Set progress callback for schedulable execution.

        Parameters
        ----------
        callback : callable
            Function with signature: callback(progress: float, message: str)
            where progress is 0.0 to 1.0
        """
        self.progress_callback = callback

    def _report_progress(self, progress: float, message: str):
        """Report progress if callback is set."""
        if self.progress_callback:
            self.progress_callback(progress, message)

    def cancel(self):
        """Cancel ongoing analysis."""
        self._cancel_flag = True
        logger.info("Analysis cancellation requested")

    def is_cancelled(self) -> bool:
        """Check if cancellation was requested."""
        return self._cancel_flag

    def analyze(self, input_params: DockingAngleInput) -> DockingAngleResult:
        """
        Execute docking angle analysis.

        Parameters
        ----------
        input_params : DockingAngleInput
            Standardized input parameters

        Returns
        -------
        DockingAngleResult
            Standardized output result
        """
        try:
            # 0. Check cancellation
            if self.is_cancelled():
                return DockingAngleResult(
                    success=False,
                    error_message="Analysis cancelled by user"
                )

            # 1. Input validation
            self._report_progress(0.0, "Validating input parameters")
            logger.info("Starting docking angle analysis")
            logger.info(f"  Topology: {input_params.topology}")
            logger.info(f"  Trajectory: {input_params.trajectory or 'None (single frame)'}")
            logger.info(f"  Auto-identify chains: {input_params.auto_identify_chains}")

            input_params.validate()

            # 2. Initialize universe
            self._report_progress(0.1, "Loading structure and trajectory")
            if input_params.trajectory:
                universe = mda.Universe(input_params.topology, input_params.trajectory)
            else:
                universe = mda.Universe(input_params.topology)

            # 3. Auto chain identification (if enabled)
            auto_selections = None
            chain_identifications = None

            if input_params.auto_identify_chains:
                self._report_progress(0.15, "Identifying chains automatically")
                auto_selections, chain_identifications = self._auto_identify_chains(
                    input_params.topology,
                    input_params.use_anarci
                )

            # 4. Determine selections
            mhc_sel, tcr_alpha_sel, tcr_beta_sel = self._determine_selections(
                input_params, auto_selections
            )

            # 5. Execute analysis
            if input_params.trajectory:
                result = self._analyze_trajectory(
                    universe, input_params, mhc_sel, tcr_alpha_sel, tcr_beta_sel,
                    chain_identifications
                )
            else:
                result = self._analyze_single_frame(
                    universe, input_params, mhc_sel, tcr_alpha_sel, tcr_beta_sel,
                    chain_identifications
                )

            self._report_progress(1.0, "Analysis complete")
            logger.info("Docking angle analysis completed successfully")

            return result

        except Exception as e:
            logger.exception("Docking angle analysis failed")
            return DockingAngleResult(
                success=False,
                error_message=f"{type(e).__name__}: {str(e)}"
            )

    def _auto_identify_chains(
        self,
        topology_file: str,
        use_anarci: bool
    ) -> Tuple[Dict, Dict]:
        """
        Perform automatic chain identification.

        Returns
        -------
        auto_selections : dict
            Selection strings for each chain type
        identifications : dict
            Full chain identification results
        """
        from immunoscope.analysis import ChainIdentificationAdapter
        from immunoscope.utils import SelectionStringBuilder

        logger.info("=" * 60)
        logger.info("Initializing automatic chain identification...")
        logger.info("=" * 60)

        adapter = ChainIdentificationAdapter(use_anarci=use_anarci)

        try:
            identifications = adapter.identify_chains(
                topology_file=topology_file,
                structure_file=None
            )
        except Exception as e:
            logger.error(f"Chain identification failed: {e}")
            raise ValueError(f"Failed to identify chains: {e}")

        is_valid, messages = adapter.validate_identification(identifications, strict=False)

        if not is_valid:
            error_messages = [msg for msg in messages if msg.startswith("ERROR")]
            logger.error("Chain identification validation failed:")
            for msg in error_messages:
                logger.error(f"  {msg}")
            raise ValueError(
                f"Chain identification failed validation. "
                f"Errors: {'; '.join(error_messages)}"
            )

        warning_messages = [msg for msg in messages if msg.startswith("WARNING")]
        for msg in warning_messages:
            logger.warning(msg)

        builder = SelectionStringBuilder()
        auto_selections = builder.build_all_selections(identifications)

        logger.info("=" * 60)
        logger.info("✓ Automatic chain identification successful")
        logger.info("=" * 60)
        logger.info(adapter.get_chain_summary(identifications))
        logger.info("=" * 60)

        return auto_selections, identifications

    def _determine_selections(
        self,
        input_params: DockingAngleInput,
        auto_selections: Optional[Dict]
    ) -> Tuple[str, str, str]:
        """
        Determine which selections to use (manual vs auto).

        Returns
        -------
        mhc_selection : str
        tcr_alpha_selection : str
        tcr_beta_selection : str
        """
        if input_params.mhc_selection and input_params.tcr_alpha_selection and input_params.tcr_beta_selection:
            logger.info("Using manually provided chain selections")
            return (
                input_params.mhc_selection,
                input_params.tcr_alpha_selection,
                input_params.tcr_beta_selection
            )

        elif input_params.auto_identify_chains and auto_selections:
            logger.info("Using auto-identified chain selections")
            mhc_sel = auto_selections.get('mhc_selection')
            tcr_alpha_sel = auto_selections.get('tcr_alpha_selection')
            tcr_beta_sel = auto_selections.get('tcr_beta_selection')

            if not tcr_alpha_sel and not tcr_beta_sel:
                raise ValueError(
                    "No TCR chains were auto-identified. Cannot calculate docking angles."
                )

            if not mhc_sel:
                raise ValueError(
                    "No HLA-alpha chain was auto-identified. Cannot calculate docking angles."
                )

            return mhc_sel, tcr_alpha_sel, tcr_beta_sel

        else:
            raise ValueError(
                "Chain selections are required. Either:\n"
                "  1. Provide all three parameters (mhc_selection, tcr_alpha_selection, tcr_beta_selection), or\n"
                "  2. Enable auto_identify_chains=True"
            )

    def _analyze_trajectory(
        self,
        universe: mda.Universe,
        input_params: DockingAngleInput,
        mhc_sel: str,
        tcr_alpha_sel: str,
        tcr_beta_sel: str,
        chain_identifications: Optional[Dict]
    ) -> DockingAngleResult:
        """
        Perform trajectory analysis with reference structure strategy.

        Key improvement: Detect conserved cysteines and MHC regions ONCE on
        the first frame, then use fixed residue IDs for all subsequent frames.
        This eliminates instability from frame-by-frame detection.
        """
        self._report_progress(0.2, "Initializing geometry detectors")

        # Initialize detectors
        groove_detector = _MHCGrooveDetector(universe, mhc_sel)
        tcr_calculator = _TCRAxisCalculator(universe)

        # === REFERENCE STRUCTURE INITIALIZATION (FIRST FRAME) ===
        logger.info("=" * 60)
        logger.info("Initializing reference geometry (first frame)")
        logger.info("=" * 60)

        universe.trajectory[0]  # Move to first frame

        # Detect conserved cysteines on reference structure
        logger.info("Detecting conserved cysteines for TCR alpha chain...")
        alpha_cys_23, alpha_cys_104, alpha_chain_id = tcr_calculator.detect_reference_cysteines(
            tcr_alpha_sel, "TCR"
        )
        logger.info(f"  Alpha chain: Cys {alpha_cys_23}, {alpha_cys_104}")

        logger.info("Detecting conserved cysteines for TCR beta chain...")
        beta_cys_23, beta_cys_104, beta_chain_id = tcr_calculator.detect_reference_cysteines(
            tcr_beta_sel, "TCR"
        )
        logger.info(f"  Beta chain: Cys {beta_cys_23}, {beta_cys_104}")

        # Get MHC α1/α2 regions from sequence alignment
        logger.info("Detecting MHC groove regions...")
        alpha1_resids, alpha2_resids = groove_detector.get_reference_resids()
        logger.info(f"  α1 helix: {len(alpha1_resids)} residues")
        logger.info(f"  α2 helix: {len(alpha2_resids)} residues")

        logger.info("=" * 60)
        logger.info("Reference geometry initialized")
        logger.info("Will use FIXED residue IDs for all frames")
        logger.info("=" * 60)

        # Build fixed selection strings
        alpha_cys_selection = tcr_calculator.get_cys_selection_from_resids(
            alpha_cys_23, alpha_cys_104, alpha_chain_id
        )
        beta_cys_selection = tcr_calculator.get_cys_selection_from_resids(
            beta_cys_23, beta_cys_104, beta_chain_id
        )

        # Extract chain identifier from mhc_sel
        mhc_chain_atoms = universe.select_atoms(mhc_sel)
        first_mhc_atom = mhc_chain_atoms[0]
        if hasattr(first_mhc_atom, 'chainID') and first_mhc_atom.chainID:
            mhc_chain_id = f"chainID {first_mhc_atom.chainID}"
        elif hasattr(first_mhc_atom, 'segid') and first_mhc_atom.segid:
            mhc_chain_id = f"segname {first_mhc_atom.segid}"
        else:
            mhc_chain_id = mhc_sel

        # Reuse AtomGroups in the hot loop to avoid repeatedly parsing selection strings.
        alpha1_ca_ag = universe.select_atoms(
            f"{mhc_chain_id} and (resid {' '.join(map(str, alpha1_resids))}) and name CA"
        )
        alpha2_ca_ag = universe.select_atoms(
            f"{mhc_chain_id} and (resid {' '.join(map(str, alpha2_resids))}) and name CA"
        )

        ref_all_resids = sorted(alpha1_resids + alpha2_resids)
        ref_backbone = universe.select_atoms(
            f"{mhc_chain_id} and "
            f"(resid {' '.join(map(str, ref_all_resids))}) and "
            f"(name CA or name C or name N or name O)"
        )

        groove_frame = self._build_reference_groove_frame(
            ref_alpha1_ca=alpha1_ca_ag,
            ref_alpha2_ca=alpha2_ca_ag,
            ref_backbone=ref_backbone,
        )
        ref_groove_longitudinal = groove_frame['longitudinal_axis']
        ref_groove_lateral = groove_frame['lateral_axis']
        ref_groove_normal = groove_frame['normal_axis']
        logger.info("Groove local frame fixed at reference frame (PCA helix axes + reference plane)")

        # Calculate groove COM for robust TCR axis calculation
        groove_com = ref_backbone.center_of_mass()

        # Use robust TCR axis calculation (works for both V-only and full TCR)
        logger.info("Using robust TCR axis calculation (TCR COM - Groove COM)")
        ref_tcr_axis = tcr_calculator.get_robust_tcr_axis(
            tcr_alpha_sel, tcr_beta_sel, groove_com
        )

        # === TRAJECTORY ANALYSIS (ALL FRAMES) ===
        stride = input_params.stride
        n_frames = len(range(0, len(universe.trajectory), stride))
        times = np.zeros(n_frames)
        crossing_angles = np.zeros(n_frames)
        incident_angles = np.zeros(n_frames)
        previous_tcr_axis = ref_tcr_axis.copy()

        logger.info(f"\nCalculating docking angles for {n_frames} frames (stride={stride})")
        logger.info("Using fixed reference groove frame and frame-continuous TCR axis")
        logger.info("-" * 60)

        # Iterate through trajectory
        for i, ts in enumerate(universe.trajectory[::stride]):
            if self.is_cancelled():
                return DockingAngleResult(
                    success=False,
                    error_message="Analysis cancelled by user"
                )

            progress = 0.2 + (0.6 * (i + 1) / n_frames)
            self._report_progress(progress, f"Processing frame {i+1}/{n_frames}")

            # Calculate TCR axis using robust method (works for both V-only and full TCR)
            current_groove_com = ref_backbone.center_of_mass()
            tcr_axis = tcr_calculator.get_robust_tcr_axis(
                tcr_alpha_sel, tcr_beta_sel, current_groove_com
            )

            tcr_axis = _GeometryUtils.align_vector_orientation(tcr_axis, previous_tcr_axis)
            previous_tcr_axis = tcr_axis.copy()

            # Use a fixed groove reference frame; rebuilding the groove per frame amplifies local noise.
            groove_longitudinal = ref_groove_longitudinal
            groove_lateral = ref_groove_lateral
            groove_normal = ref_groove_normal

            # Crossing: project the TCR direction onto the groove plane, then compare it with the groove axis.
            tcr_projected = _GeometryUtils.project_vector_onto_plane(tcr_axis, groove_normal)
            projected_norm = np.linalg.norm(tcr_projected)
            if projected_norm < 1e-8:
                crossing = 90.0
            else:
                tcr_projected = tcr_projected / projected_norm
                tcr_projected = _GeometryUtils.align_vector_orientation(
                    tcr_projected, ref_groove_longitudinal
                )
                crossing = _GeometryUtils.angle_between_vectors(
                    tcr_projected, groove_longitudinal
                )

            # Incident: angle between TCR direction and groove normal, reported in acute-angle form.
            incident_raw = _GeometryUtils.angle_between_vectors(tcr_axis, groove_normal)
            incident = _GeometryUtils.canonicalize_axial_angle(incident_raw)

            times[i] = ts.time
            crossing_angles[i] = crossing
            incident_angles[i] = incident

            if input_params.print_each_frame:
                logger.info(
                    f"  Frame {i + 1:>4}/{n_frames:<4} | "
                    f"Time {ts.time:>10.2f} ps | "
                    f"Crossing {crossing:>8.2f}° | "
                    f"Incident {incident:>8.2f}°"
                )

            if (i + 1) % 100 == 0:
                logger.info(f"  Processed {i + 1}/{n_frames} frames")

        logger.info("-" * 60)

        self._report_progress(0.8, "Calculating statistics")

        # Calculate statistics
        statistics = {
            'crossing_mean': float(np.mean(crossing_angles)),
            'crossing_std': float(np.std(crossing_angles)),
            'crossing_min': float(np.min(crossing_angles)),
            'crossing_max': float(np.max(crossing_angles)),
            'incident_mean': float(np.mean(incident_angles)),
            'incident_std': float(np.std(incident_angles)),
            'incident_min': float(np.min(incident_angles)),
            'incident_max': float(np.max(incident_angles)),
            'n_frames': len(times)
        }

        logger.info(
            f"Trajectory statistics: Crossing={statistics['crossing_mean']:.2f}±{statistics['crossing_std']:.2f}°, "
            f"Incident={statistics['incident_mean']:.2f}±{statistics['incident_std']:.2f}°"
        )

        # Save output files
        output_files = []
        if input_params.output_dir:
            output_dir = Path(input_params.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            csv_file = output_dir / "docking_angles.csv"
            self._save_csv(csv_file, times, crossing_angles, incident_angles)
            output_files.append(str(csv_file))
            logger.info(f"Saved results to: {csv_file}")

        # Collect metadata
        metadata = {
            'module_version': '4.3.0',
            'algorithm': 'reference_groove_frame_projected_crossing',
            'angle_convention': 'crossing_projected_raw__incident_acute',
            'timestamp': datetime.now().isoformat(),
            'topology': input_params.topology,
            'trajectory': input_params.trajectory,
            'stride': input_params.stride,
            'print_each_frame': input_params.print_each_frame,
            'auto_identify_chains': input_params.auto_identify_chains,
            'reference_residues': {
                'tcr_alpha_cys': [alpha_cys_23, alpha_cys_104],
                'tcr_beta_cys': [beta_cys_23, beta_cys_104],
                'mhc_alpha1_count': len(alpha1_resids),
                'mhc_alpha2_count': len(alpha2_resids),
                'groove_frame': {
                    'longitudinal_axis': ref_groove_longitudinal.tolist(),
                    'lateral_axis': ref_groove_lateral.tolist(),
                    'normal_axis': ref_groove_normal.tolist(),
                }
            }
        }

        if chain_identifications:
            metadata['chain_identifications'] = chain_identifications

        try:
            alignment_quality = groove_detector.get_alignment_quality()
            metadata['alignment_quality'] = alignment_quality
        except Exception:
            pass

        return DockingAngleResult(
            success=True,
            times=times,
            crossing_angles=crossing_angles,
            incident_angles=incident_angles,
            statistics=statistics,
            output_files=output_files,
            metadata=metadata
        )

    def _build_reference_groove_frame(
        self,
        ref_alpha1_ca: mda.AtomGroup,
        ref_alpha2_ca: mda.AtomGroup,
        ref_backbone: mda.AtomGroup,
    ) -> Dict[str, np.ndarray]:
        """
        Build a local coordinate frame from the two groove walls in the reference frame.

        Design inspired by STCRpy:
        - use the PCA axes of the two helices to define the longitudinal groove axis
        - use the centroid connector between the two helices to define the lateral axis
        - use the reference-frame plane normal to orient the final normal vector
        """
        alpha1_coords = ref_alpha1_ca.positions.copy()
        alpha2_coords = ref_alpha2_ca.positions.copy()

        alpha1_axis = _GeometryUtils.principal_axis_from_coords(alpha1_coords)
        alpha2_axis = _GeometryUtils.principal_axis_from_coords(alpha2_coords)
        alpha2_axis = _GeometryUtils.align_vector_orientation(alpha2_axis, alpha1_axis)

        longitudinal = alpha1_axis + alpha2_axis
        longitudinal = longitudinal / np.linalg.norm(longitudinal)

        c1 = alpha1_coords.mean(axis=0)
        c2 = alpha2_coords.mean(axis=0)
        lateral = c2 - c1
        lateral = lateral - np.dot(lateral, longitudinal) * longitudinal
        lateral = lateral / np.linalg.norm(lateral)

        _, plane_normal = _GeometryUtils.fit_plane_svd(ref_backbone)
        normal = np.cross(longitudinal, lateral)
        normal = normal / np.linalg.norm(normal)
        normal = _GeometryUtils.align_vector_orientation(normal, plane_normal)

        # Re-orthogonalize to avoid numerical error propagation.
        lateral = np.cross(normal, longitudinal)
        lateral = lateral / np.linalg.norm(lateral)

        return {
            'alpha1_axis': alpha1_axis,
            'alpha2_axis': alpha2_axis,
            'longitudinal_axis': longitudinal,
            'lateral_axis': lateral,
            'normal_axis': normal,
            'alpha1_centroid': c1,
            'alpha2_centroid': c2,
        }

    def _analyze_single_frame(
        self,
        universe: mda.Universe,
        input_params: DockingAngleInput,
        mhc_sel: str,
        tcr_alpha_sel: str,
        tcr_beta_sel: str,
        chain_identifications: Optional[Dict]
    ) -> DockingAngleResult:
        """Perform single frame analysis."""
        self._report_progress(0.5, "Calculating docking angles for single frame")

        # Initialize detectors
        groove_detector = _MHCGrooveDetector(universe, mhc_sel)
        tcr_calculator = _TCRAxisCalculator(universe)

        # Detect conserved cysteines
        alpha_cys_23, alpha_cys_104, _ = tcr_calculator.detect_reference_cysteines(
            tcr_alpha_sel, "TCR"
        )
        beta_cys_23, beta_cys_104, _ = tcr_calculator.detect_reference_cysteines(
            tcr_beta_sel, "TCR"
        )

        # Calculate groove geometry first (needed for robust TCR axis)
        alpha1_resids, alpha2_resids = groove_detector.get_reference_resids()
        mhc_chain_atoms = universe.select_atoms(mhc_sel)
        first_mhc_atom = mhc_chain_atoms[0]
        if hasattr(first_mhc_atom, 'chainID') and first_mhc_atom.chainID:
            mhc_chain_id = f"chainID {first_mhc_atom.chainID}"
        elif hasattr(first_mhc_atom, 'segid') and first_mhc_atom.segid:
            mhc_chain_id = f"segname {first_mhc_atom.segid}"
        else:
            mhc_chain_id = mhc_sel
        alpha1_ca = universe.select_atoms(
            f"{mhc_chain_id} and (resid {' '.join(map(str, alpha1_resids))}) and name CA"
        )
        alpha2_ca = universe.select_atoms(
            f"{mhc_chain_id} and (resid {' '.join(map(str, alpha2_resids))}) and name CA"
        )
        all_resids = sorted(alpha1_resids + alpha2_resids)
        backbone = universe.select_atoms(
            f"{mhc_chain_id} and "
            f"(resid {' '.join(map(str, all_resids))}) and "
            f"(name CA or name C or name N or name O)"
        )
        groove_frame = self._build_reference_groove_frame(alpha1_ca, alpha2_ca, backbone)
        groove_longitudinal = groove_frame['longitudinal_axis']
        groove_normal = groove_frame['normal_axis']

        # Calculate groove COM for robust TCR axis calculation
        groove_com = backbone.center_of_mass()

        # Use robust TCR axis calculation (works for both V-only and full TCR)
        logger.info("Using robust TCR axis calculation (TCR COM - Groove COM)")
        tcr_axis = tcr_calculator.get_robust_tcr_axis(
            tcr_alpha_sel, tcr_beta_sel, groove_com
        )

        tcr_projected = _GeometryUtils.project_vector_onto_plane(tcr_axis, groove_normal)
        projected_norm = np.linalg.norm(tcr_projected)
        if projected_norm < 1e-8:
            crossing = 90.0
        else:
            tcr_projected = tcr_projected / projected_norm
            tcr_projected = _GeometryUtils.align_vector_orientation(
                tcr_projected, groove_longitudinal
            )
            crossing = _GeometryUtils.angle_between_vectors(tcr_projected, groove_longitudinal)
        incident_raw = _GeometryUtils.angle_between_vectors(tcr_axis, groove_normal)
        incident = _GeometryUtils.canonicalize_axial_angle(incident_raw)

        logger.info(f"Single frame result: Crossing={crossing:.2f}, Incident={incident:.2f}")

        # Collect metadata
        metadata = {
            'module_version': '4.3.0',
            'algorithm': 'reference_groove_frame_projected_crossing',
            'angle_convention': 'crossing_projected_raw__incident_acute',
            'timestamp': datetime.now().isoformat(),
            'topology': input_params.topology,
            'auto_identify_chains': input_params.auto_identify_chains
        }

        if chain_identifications:
            metadata['chain_identifications'] = chain_identifications

        return DockingAngleResult(
            success=True,
            crossing_angle=crossing,
            incident_angle=incident,
            metadata=metadata
        )

    def _save_csv(
        self,
        filepath: Path,
        times: np.ndarray,
        crossing: np.ndarray,
        incident: np.ndarray
    ):
        """Save results to CSV file."""
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Time(ps)', 'Crossing(deg)', 'Incident(deg)'])

            for t, c, i in zip(times, crossing, incident):
                writer.writerow([f'{t:.2f}', f'{c:.4f}', f'{i:.4f}'])
