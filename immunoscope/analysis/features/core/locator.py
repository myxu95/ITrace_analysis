"""
Case directory locator for resolving module output paths.

Originally lived under immunoscope/agent/presenters/. Promoted to analysis/features/core/
so that both the presenters layer (Markdown views for the Agent) and the features
layer (design feature engineering) can share a single implementation.

The agent/presenters/locator.py file re-exports this class to preserve backward
compatibility.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Optional


class CaseLocatorError(Exception):
    """Raised when a case directory cannot be resolved."""
    pass


class CaseLocator:
    """
    Resolves module output paths from a case directory.

    Reads run_summary.json to find module roots, with fallback to
    hardcoded glob patterns for legacy/incomplete runs.
    """

    # Fallback directory patterns for each analysis module, in priority order.
    _FALLBACK_PATTERNS: Dict[str, list[str]] = {
        "rrcs":         ["analysis/rrcs", "overview/rrcs", "rrcs"],
        "bsa":          ["analysis/bsa", "overview/bsa", "bsa"],
        "interface":    ["analysis/interface", "overview/interface", "interface"],
        "rmsf":         ["analysis/rmsf", "overview/rmsf", "rmsf"],
        "contact":      ["analysis/contact", "overview/contact", "contact"],
        "inter_cluster":["analysis/inter_cluster", "overview/inter_cluster", "inter_cluster"],
        "landscape":    ["analysis/landscape", "overview/landscape", "landscape"],
        "quality":      ["analysis/quality", "preparation/quality", "overview/quality", "quality"],
        "identity":     ["analysis/identity", "overview/identity", "identity"],
        "cdr3_geometry":["analysis/angles/cdr3_geometry", "analysis/geometry/cdr3_geometry", "analysis/cdr3_geometry"],
        "hbond":        ["analysis/interactions/hydrogen_bonds",         "analysis/hbond",       "hbond"],
        "saltbridge":   ["analysis/interactions/salt_bridges",           "analysis/saltbridge",  "saltbridge"],
        "hydrophobic":  ["analysis/interactions/hydrophobic_contacts",   "analysis/hydrophobic", "hydrophobic"],
        "pipi":         ["analysis/interactions/pi_interactions",        "analysis/interactions/pi_pi", "analysis/pipi", "pipi"],
        "cationpi":     ["analysis/interactions/cation_pi_interactions", "analysis/interactions/cation_pi", "analysis/cationpi", "cationpi"],
        "structure":    ["analysis/structure", "overview/structure", "structure"],
        "chi_dihedrals":["analysis/geometry/chi_dihedrals", "analysis/chi_dihedrals", "chi_dihedrals"],
        "dihedrals":    ["analysis/dihedrals", "overview/dihedrals", "dihedrals"],
        "conservation": ["analysis/conservation", "overview/conservation", "conservation"],
    }

    def __init__(self, case_dir: Path | str):
        """
        Initialize locator.

        Args:
            case_dir: Path to analysis case directory

        Raises:
            CaseLocatorError: If case_dir does not exist
        """
        self.case_dir = Path(case_dir).resolve()

        if not self.case_dir.exists():
            raise CaseLocatorError(f"Case directory not found: {self.case_dir}")

        self._module_roots: Dict[str, Optional[Path]] = {}
        self._run_summary: Optional[Dict] = None

        self._load_run_summary()

    # ------------------------------------------------------------------ #
    # Initialization helpers
    # ------------------------------------------------------------------ #

    def _load_run_summary(self) -> None:
        """Load run_summary.json if it exists (non-fatal on failure)."""
        summary_path = self.case_dir / "run_summary.json"
        if not summary_path.exists():
            return
        try:
            with open(summary_path) as f:
                self._run_summary = json.load(f)
        except Exception:
            # Fall back to glob patterns if parsing fails.
            self._run_summary = None

    # ------------------------------------------------------------------ #
    # Module resolution
    # ------------------------------------------------------------------ #

    def get_module_root(self, module: str) -> Optional[Path]:
        """
        Get the root directory for a module's output.

        Args:
            module: Module name (rrcs, bsa, rmsf, contact, etc.)

        Returns:
            Path to module root, or None if not found
        """
        if module in self._module_roots:
            return self._module_roots[module]

        # Try run_summary.json first
        if self._run_summary:
            module_results = self._run_summary.get("module_results", {})
            if module in module_results:
                root_str = module_results[module].get("root")
                if root_str:
                    root = Path(root_str)
                    if root.exists():
                        self._module_roots[module] = root
                        return root

        # Fallback: hardcoded patterns
        root = self._find_module_root_fallback(module)
        self._module_roots[module] = root
        return root

    def _find_module_root_fallback(self, module: str) -> Optional[Path]:
        for pattern in self._FALLBACK_PATTERNS.get(module, []):
            candidate = self.case_dir / pattern
            if candidate.exists():
                return candidate
        return None

    # ------------------------------------------------------------------ #
    # File lookup
    # ------------------------------------------------------------------ #

    def find_file(self, module: str, *path_parts: str) -> Optional[Path]:
        """
        Find a file within a module's output.

        Args:
            module: Module name
            *path_parts: Path components relative to module root

        Returns:
            Full path to file, or None if not found
        """
        root = self.get_module_root(module)
        if not root:
            return None

        file_path = root / Path(*path_parts)
        return file_path if file_path.exists() else None

    def first_existing(self, *paths: Path) -> Optional[Path]:
        """
        Return the first existing path from a list.

        Args:
            *paths: Candidate paths (absolute or relative to case_dir)

        Returns:
            First existing path, or None
        """
        for path in paths:
            if not path.is_absolute():
                path = self.case_dir / path

            if path.exists():
                return path

        return None

    # ------------------------------------------------------------------ #
    # Metadata
    # ------------------------------------------------------------------ #

    def get_case_id(self) -> str:
        """Get case/system ID from directory name or run_summary."""
        if self._run_summary:
            job_id = self._run_summary.get("job_id")
            if job_id:
                return job_id

        return self.case_dir.name

    def available_modules(self) -> list[str]:
        """
        List modules that have resolvable output roots in this case.

        Returns:
            Sorted list of module names whose output directories exist.
        """
        found: list[str] = []

        # From run_summary if available
        if self._run_summary:
            for name in self._run_summary.get("module_results", {}).keys():
                if self.get_module_root(name) is not None:
                    found.append(name)

        # From fallback patterns
        for name in self._FALLBACK_PATTERNS.keys():
            if name not in found and self.get_module_root(name) is not None:
                found.append(name)

        return sorted(set(found))
