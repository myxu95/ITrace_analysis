import subprocess
from pathlib import Path
from typing import List, Optional, Dict, Any
import logging

from ...utils.group_selector import GroupSelector
from ...utils.cleanup_manager import CleanupManager

logger = logging.getLogger(__name__)


class PBCProcessor:
    """PBC (Periodic Boundary Conditions) processing for MD trajectories."""

    def __init__(self, gmx_executable: str = "gmx", keep_temp_files: bool = False):
        """
        Initialize PBCProcessor.

        Args:
            gmx_executable: GROMACS executable command
            keep_temp_files: Whether to keep temporary files for debugging (default: False)
        """
        self.gmx = gmx_executable
        self.group_selector = None  # Will be initialized when topology is provided
        self.keep_temp_files = keep_temp_files
        self.cleanup_manager = CleanupManager(keep_temp_files=keep_temp_files, verbose=True)
        self._check_gromacs()

    def _check_gromacs(self):
        """Check if GROMACS is available."""
        try:
            subprocess.run([self.gmx, "--version"], capture_output=True, text=True, check=True)
            logger.info("GROMACS found and accessible")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("GROMACS not found or not accessible")

    def _register_temp_file(self, filepath: str) -> str:
        """Register a temporary file for later cleanup."""
        return self.cleanup_manager.register_file(filepath)

    def _cleanup_temp_files(self, output_dir: Optional[str] = None):
        """Clean up all registered temporary files and apply pattern-based cleanup."""
        self.cleanup_manager.cleanup_registered_files()
        if output_dir and not self.keep_temp_files:
            self.cleanup_manager.cleanup_by_patterns(
                directory=output_dir,
                categories=['gromacs_temp', 'index_files']
            )

    def _run_gmx_command(self,
                        command: List[str],
                        stdin_input: Optional[str] = None,
                        cwd: Optional[str] = None) -> subprocess.CompletedProcess:
        """
        Run GROMACS command with error handling.

        Args:
            command: GROMACS command list
            stdin_input: Input to pass to stdin
            cwd: Working directory

        Returns:
            CompletedProcess result
        """
        full_command = [self.gmx] + command
        logger.info(f"Running: {' '.join(full_command)}")

        try:
            return subprocess.run(
                full_command,
                input=stdin_input,
                text=True,
                capture_output=True,
                check=True,
                cwd=cwd
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"GROMACS command failed: {e}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            raise

    def remove_pbc_2step(self,
                         trajectory: str,
                         topology: str,
                         output: str,
                         fit_group: Optional[str] = None,
                         output_group: Optional[str] = None,
                         dt: Optional[float] = None) -> str:
        """
        Remove PBC using simplified 2-step process: nojump -> fit rot+trans.
        """
        if self.group_selector is None:
            self.group_selector = GroupSelector(topology, self.gmx)

        if fit_group is None:
            fit_group = self.group_selector.select_group("fit_backbone")
        if output_group is None:
            output_group = self.group_selector.select_group("output_system")

        temp_nojump = self._register_temp_file(str(Path(output).with_suffix('.temp_nojump.xtc')))

        try:
            logger.info("Step 1: Applying PBC nojump to prevent atom jumps")
            nojump_cmd = [
                "trjconv",
                "-f", trajectory,
                "-s", topology,
                "-o", temp_nojump,
                "-pbc", "nojump"
            ]

            if dt is not None:
                nojump_cmd.extend(["-dt", str(dt)])
                logger.info(f"Using dt={dt} ps for trajectory downsampling")

            self._run_gmx_command(nojump_cmd, stdin_input=f"{output_group}\n")
            logger.info("PBC nojump completed")

            logger.info("Step 2: Fitting rotational and translational motion")
            fit_cmd = [
                "trjconv",
                "-f", temp_nojump,
                "-s", topology,
                "-o", output,
                "-fit", "rot+trans"
            ]

            fit_input = f"{fit_group}\n{output_group}\n"
            self._run_gmx_command(fit_cmd, stdin_input=fit_input)
            logger.info(f"Applied rot+trans fitting using group: {fit_group}")

            logger.info(f"2-step PBC processing completed: {output}")
            logger.info(f"Used groups - Fit: {fit_group}, Output: {output_group}")
            return output
        finally:
            self._cleanup_temp_files(str(Path(output).parent))

    def remove_pbc(self,
                   trajectory: str,
                   topology: str,
                   output: str,
                   center_group: Optional[str] = None,
                   output_group: Optional[str] = None,
                   fit_group: Optional[str] = None,
                   dt: Optional[float] = None,
                   use_nojump: bool = False) -> str:
        """
        Remove PBC using three-step process: center -> pbc whole -> fit rot+trans.
        """
        if self.group_selector is None:
            self.group_selector = GroupSelector(topology, self.gmx)

        if center_group is None:
            center_group = self.group_selector.select_group("center_protein")
        if output_group is None:
            output_group = self.group_selector.select_group("output_system")
        if fit_group is None:
            fit_group = self.group_selector.select_group("fit_backbone")

        temp_centered = self._register_temp_file(str(Path(output).with_suffix('.temp_centered.xtc')))
        temp_whole = self._register_temp_file(str(Path(output).with_suffix('.temp_whole.xtc')))

        try:
            pbc_method = "nojump" if use_nojump else "atom"
            if dt is not None:
                logger.info(f"Step 1: Centering trajectory (center pbc {pbc_method}) with dt={dt} ps")
            else:
                logger.info(f"Step 1: Centering trajectory (center pbc {pbc_method})")

            center_cmd = [
                "trjconv",
                "-f", trajectory,
                "-s", topology,
                "-o", temp_centered,
                "-center",
                "-pbc", pbc_method
            ]

            if self.group_selector and self.group_selector.has_shortest_chain_group():
                index_file = self.group_selector.get_shortest_chain_index_file()
                if index_file:
                    center_cmd.extend(["-n", index_file])
                    logger.info(f"Using custom index file for centering: {index_file}")

            if dt is not None:
                center_cmd.extend(["-dt", str(dt)])

            center_input = f"{center_group}\n{output_group}\n"
            self._run_gmx_command(center_cmd, stdin_input=center_input)
            logger.info(f"Centered trajectory using group: {center_group}")

            logger.info("Step 2: Applying PBC whole")
            whole_cmd = [
                "trjconv",
                "-f", temp_centered,
                "-s", topology,
                "-o", temp_whole,
                "-pbc", "whole"
            ]
            self._run_gmx_command(whole_cmd, stdin_input=f"{output_group}\n")
            logger.info("Applied PBC whole")

            logger.info("Step 3: Fitting rotational and translational motion")
            fit_cmd = [
                "trjconv",
                "-f", temp_whole,
                "-s", topology,
                "-o", output,
                "-fit", "rot+trans"
            ]
            fit_input = f"{fit_group}\n{output_group}\n"
            self._run_gmx_command(fit_cmd, stdin_input=fit_input)
            logger.info(f"Applied rot+trans fitting using group: {fit_group}")

            logger.info(f"3-step PBC processing completed: {output}")
            logger.info(f"Used groups - Center: {center_group}, Fit: {fit_group}, Output: {output_group}")
            return output
        finally:
            self._cleanup_temp_files(str(Path(output).parent))

    def comprehensive_pbc_process(self,
                                 trajectory: str,
                                 topology: str,
                                 output_dir: str,
                                 method: str = "2step",
                                 dt: Optional[float] = None,
                                 fit_group: Optional[str] = None,
                                 output_group: Optional[str] = None,
                                 center_group: Optional[str] = None) -> Dict[str, Any]:
        """
        Comprehensive PBC processing with trajectory conversion, PBC removal, and PDB generation.
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        base_name = Path(trajectory).stem
        processed_file = output_path / f"{base_name}_processed.xtc"
        pdb_file = output_path / f"{processed_file.stem}_converted.pdb"

        result = {
            "input_trajectory": trajectory,
            "input_topology": topology,
            "output_dir": str(output_path),
            "processed": None,
            "pdb": None,
            "success": False,
            "method": method,
            "shortest_chain_id": None,
        }

        try:
            if self.group_selector is None:
                self.group_selector = GroupSelector(topology, self.gmx)

            if self.group_selector.has_shortest_chain_group():
                result["shortest_chain_id"] = self.group_selector.get_shortest_chain_id()

            if method == "2step":
                final_trajectory = self.remove_pbc_2step(
                    trajectory=trajectory,
                    topology=topology,
                    output=str(processed_file),
                    fit_group=fit_group,
                    output_group=output_group,
                    dt=dt,
                )
            elif method == "3step":
                final_trajectory = self.remove_pbc(
                    trajectory=trajectory,
                    topology=topology,
                    output=str(processed_file),
                    center_group=center_group,
                    output_group=output_group,
                    fit_group=fit_group,
                    dt=dt,
                )
            else:
                raise ValueError(f"Unsupported PBC method: {method}")

            try:
                convert_cmd = [
                    "trjconv",
                    "-f", final_trajectory,
                    "-s", topology,
                    "-o", str(pdb_file),
                    "-dump", "0"
                ]
                pdb_output_group = output_group or (
                    self.group_selector.select_group("output_system")
                    if self.group_selector else "System"
                )
                self._run_gmx_command(convert_cmd, stdin_input=f"{pdb_output_group}\n")
                result["pdb"] = str(pdb_file)
            except Exception as exc:
                logger.warning(f"Failed to generate converted PDB: {exc}")

            result["processed"] = final_trajectory
            result["success"] = True
            return result
        except Exception as exc:
            logger.error(f"Comprehensive PBC processing failed: {exc}")
            result["error"] = str(exc)
            return result
