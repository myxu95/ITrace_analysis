import subprocess
import re
import json
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class GroupSelector:
    """GROMACS group selection manager with predefined mappings and dynamic detection."""
    
    def __init__(self, topology_file: str, gmx_executable: str = "gmx"):
        """
        Initialize GroupSelector.

        Args:
            topology_file: Topology file (.tpr, .gro, .pdb)
            gmx_executable: GROMACS executable command
        """
        self.topology = topology_file
        self.gmx = gmx_executable
        self.available_groups = {}
        self.group_mappings = self._load_default_mappings()
        self.config_file = Path.home() / ".itrace_groups.json"
        self.shortest_chain_index_file = None
        self.shortest_chain_group_id = None

        # Load user configurations if available
        self._load_user_config()

        # Detect available groups
        self._detect_available_groups()

        # Try to create shortest chain index from md.gro
        self._try_create_shortest_chain_index()
    
    def _load_default_mappings(self) -> Dict[str, List[str]]:
        """Load default group mappings based on standard GROMACS group indices."""
        return {
            # RMSD analysis - fit groups (for structural alignment/superposition)
            "rmsd_fit_backbone": ["4", "5", "3"],  # Backbone → MainChain → C-alpha
            "rmsd_fit_protein": ["4", "5", "3"],   # Usually fit on backbone even for protein RMSD
            "rmsd_fit_calpha": ["3"],              # C-alpha
            "rmsd_fit_heavy": ["4", "5", "2"],     # Backbone → MainChain → Protein-H
            
            # RMSD analysis - calculation groups (for actual RMSD calculation)
            "rmsd_calc_backbone": ["4", "5"],      # Backbone → MainChain
            "rmsd_calc_protein": ["1", "2"],       # Protein → Protein-H
            "rmsd_calc_calpha": ["3"],             # C-alpha
            "rmsd_calc_heavy": ["2", "1"],         # Protein-H → Protein
            
            # Centering operations (will be updated if shortest chain is available)
            "center_protein": ["1", "2"],          # Protein → Protein-H (fallback)
            "center_backbone": ["4", "5"],         # Backbone → MainChain
            
            # Output selections
            "output_system": ["0"],               # System
            "output_protein": ["1", "2"],         # Protein → Protein-H
            "output_water": ["12", "13"],         # Water → SOL
            "output_non_water": ["14"],           # non-Water
            
            # RDF analysis
            "rdf_water": ["12", "13"],            # Water → SOL
            "rdf_protein": ["1", "2"],            # Protein → Protein-H
            "rdf_ions": ["16", "17"],             # SOD → CLA
            
            # Distance analysis
            "distance_com": ["1", "2"],           # Protein → Protein-H
            "distance_backbone": ["4", "5"],      # Backbone → MainChain
            "distance_calpha": ["3"],             # C-alpha
            
            # Energy groups (for trajectory conversion)
            "energy_protein": ["1", "2"],         # Protein → Protein-H
            "energy_water": ["12", "13"],         # Water → SOL
            "energy_ions": ["16", "17"],          # SOD → CLA
            
            # Special selections for fitting and PBC
            "fit_protein": ["1", "2"],            # Protein → Protein-H
            "fit_backbone": ["4", "5"],           # Backbone → MainChain
            "pbc_protein": ["1", "2"],            # Protein → Protein-H
            "pbc_system": ["0"]                   # System
        }
    
    def _detect_available_groups(self):
        """Detect available groups from topology file."""
        try:
            # Method 1: Try to get groups from gmx make_ndx
            self._detect_groups_from_make_ndx()
        except Exception as e:
            logger.warning(f"Could not detect groups from gmx make_ndx: {e}")
            # Method 2: Use default GROMACS groups
            self._use_default_gromacs_groups()
    
    def _detect_groups_from_make_ndx(self):
        """Detect groups by running gmx make_ndx."""
        cmd = [self.gmx, "make_ndx", "-f", self.topology, "-o", "/dev/null"]
        
        try:
            result = subprocess.run(
                cmd,
                input="q\n",  # Quit immediately
                text=True,
                capture_output=True,
                check=True,
                timeout=30
            )
            
            # Parse the output to extract group information
            self._parse_make_ndx_output(result.stdout)
            
        except subprocess.TimeoutExpired:
            logger.warning("gmx make_ndx timed out")
        except subprocess.CalledProcessError as e:
            logger.warning(f"gmx make_ndx failed: {e}")
            logger.debug(f"make_ndx stderr: {e.stderr}")
            logger.debug(f"make_ndx stdout: {e.stdout}")
            
    def _parse_make_ndx_output(self, output: str):
        """Parse gmx make_ndx output to extract group information."""
        lines = output.split('\n')
        
        for line in lines:
            # Look for group definitions like "  0 System              : 12345 atoms"
            match = re.match(r'\s*(\d+)\s+(\w+(?:[-_]\w+)*)\s*:\s*(\d+)\s+atoms', line)
            if match:
                group_id = int(match.group(1))
                group_name = match.group(2)
                atom_count = int(match.group(3))
                
                self.available_groups[group_id] = {
                    'name': group_name,
                    'atoms': atom_count
                }
        
        logger.info(f"Detected {len(self.available_groups)} groups")
    
    def _use_default_gromacs_groups(self):
        """Use standard GROMACS group layout when detection fails."""
        default_groups = {
            0: {'name': 'System', 'atoms': 'all'},
            1: {'name': 'Protein', 'atoms': 'all_protein'},
            2: {'name': 'Protein-H', 'atoms': 'protein_heavy'},
            3: {'name': 'C-alpha', 'atoms': 'ca_atoms'},
            4: {'name': 'Backbone', 'atoms': 'backbone'},
            5: {'name': 'MainChain', 'atoms': 'mainchain'},
            6: {'name': 'MainChain+Cb', 'atoms': 'mainchain_cb'},
            7: {'name': 'MainChain+H', 'atoms': 'mainchain_h'},
            8: {'name': 'SideChain', 'atoms': 'sidechain'},
            9: {'name': 'SideChain-H', 'atoms': 'sidechain_heavy'},
            10: {'name': 'Prot-Masses', 'atoms': 'protein_masses'},
            11: {'name': 'non-Protein', 'atoms': 'non_protein'},
            12: {'name': 'Water', 'atoms': 'water'},
            13: {'name': 'SOL', 'atoms': 'sol'},
            14: {'name': 'non-Water', 'atoms': 'non_water'},
            15: {'name': 'Other', 'atoms': 'other'},
            16: {'name': 'SOD', 'atoms': 'sodium'},
            17: {'name': 'CLA', 'atoms': 'chloride'}
        }
        
        self.available_groups = default_groups
        logger.info("Using standard GROMACS group layout (0-17)")
    
    def select_group(self, purpose: str, preference: Optional[str] = None) -> str:
        """
        Select appropriate group for given purpose.
        
        Args:
            purpose: Purpose of the selection (e.g., 'rmsd_protein', 'center_protein')
            preference: User preference if not using auto-selection
            
        Returns:
            Group name or number as string
        """
        # 1. Try user preference first
        if preference:
            if self._validate_group(preference):
                return preference
            else:
                logger.warning(f"Preferred group '{preference}' not found")
        
        # 2. Special handling for center_protein: try shortest chain first
        if purpose == "center_protein" and self.has_shortest_chain_group():
            # Use the precise group ID from the generated index file
            shortest_chain_id = self.get_shortest_chain_group_id()
            if shortest_chain_id:
                logger.info(f"Using shortest chain for centering: group {shortest_chain_id} (from index file)")
                return shortest_chain_id
        
        # 3. Try auto-selection based on purpose
        auto_selected = self._auto_select_group(purpose)
        if auto_selected:
            return auto_selected
        
        # 4. Fallback to first mapping option
        if purpose in self.group_mappings:
            return self.group_mappings[purpose][0]
        
        # 5. Ultimate fallback
        logger.warning(f"No suitable group found for purpose '{purpose}', using 'System'")
        return "System"
    
    def _auto_select_group(self, purpose: str) -> Optional[str]:
        """Auto-select group based on purpose and available groups."""
        if purpose not in self.group_mappings:
            return None
        
        candidates = self.group_mappings[purpose]
        
        # Try to find the first available candidate
        for candidate in candidates:
            if self._validate_group(candidate):
                logger.info(f"Auto-selected group '{candidate}' for purpose '{purpose}'")
                return candidate
        
        return None
    
    def _validate_group(self, group: str) -> bool:
        """Validate if a group exists in available groups."""
        # Check by group number
        if group.isdigit():
            group_id = int(group)
            return group_id in self.available_groups
        
        # Check by group name
        for group_info in self.available_groups.values():
            if group_info['name'] == group:
                return True
        
        return False
    
    def get_group_input_string(self, *purposes) -> str:
        """
        Generate input string for GROMACS commands requiring multiple group selections.
        
        Args:
            *purposes: Multiple purposes for group selection
            
        Returns:
            Input string with newlines for stdin
        """
        selections = []
        for purpose in purposes:
            group = self.select_group(purpose)
            selections.append(group)
        
        input_string = '\n'.join(selections) + '\n'
        logger.info(f"Generated input string for {purposes}: {repr(input_string)}")
        return input_string
    
    def get_rmsd_groups(self, rmsd_type: str = "backbone") -> Tuple[str, str]:
        """
        Get fit and calculation groups for RMSD analysis.
        
        Args:
            rmsd_type: Type of RMSD analysis ('backbone', 'protein', 'calpha', 'heavy')
            
        Returns:
            Tuple of (fit_group, calc_group)
        """
        fit_purpose = f"rmsd_fit_{rmsd_type}"
        calc_purpose = f"rmsd_calc_{rmsd_type}"
        
        fit_group = self.select_group(fit_purpose)
        calc_group = self.select_group(calc_purpose)
        
        logger.info(f"RMSD groups for {rmsd_type}: fit='{fit_group}', calc='{calc_group}'")
        return fit_group, calc_group
    
    def get_rmsd_input_string(self, rmsd_type: str = "backbone") -> str:
        """
        Generate input string specifically for gmx rms command.
        
        Args:
            rmsd_type: Type of RMSD analysis
            
        Returns:
            Input string for gmx rms (fit_group\ncalc_group\n)
        """
        fit_group, calc_group = self.get_rmsd_groups(rmsd_type)
        input_string = f"{fit_group}\n{calc_group}\n"
        
        logger.info(f"RMSD input string for {rmsd_type}: {repr(input_string)}")
        return input_string
    
    def list_available_groups(self) -> Dict[int, Dict[str, str]]:
        """List all available groups."""
        return self.available_groups
    
    def add_custom_mapping(self, purpose: str, groups: List[str]):
        """Add custom group mapping for specific purpose."""
        self.group_mappings[purpose] = groups
        self._save_user_config()
        logger.info(f"Added custom mapping for '{purpose}': {groups}")
    
    def _load_user_config(self):
        """Load user-specific group configurations."""
        if self.config_file.exists():
            try:
                with open(self.config_file, 'r') as f:
                    user_config = json.load(f)
                    self.group_mappings.update(user_config.get('mappings', {}))
                logger.info(f"Loaded user configuration from {self.config_file}")
            except Exception as e:
                logger.warning(f"Could not load user config: {e}")
    
    def _save_user_config(self):
        """Save user-specific group configurations."""
        try:
            config = {'mappings': self.group_mappings}
            with open(self.config_file, 'w') as f:
                json.dump(config, f, indent=2)
            logger.info(f"Saved user configuration to {self.config_file}")
        except Exception as e:
            logger.warning(f"Could not save user config: {e}")
    
    def create_selection_preset(self, name: str, selections: Dict[str, str]):
        """
        Create a preset for common analysis workflows.
        
        Args:
            name: Preset name
            selections: Dictionary mapping purposes to group selections
        """
        if not hasattr(self, 'presets'):
            self.presets = {}
        
        self.presets[name] = selections
        self._save_user_config()
        logger.info(f"Created preset '{name}' with {len(selections)} selections")
    
    def use_preset(self, name: str) -> Dict[str, str]:
        """Use a predefined selection preset."""
        if not hasattr(self, 'presets') or name not in self.presets:
            raise ValueError(f"Preset '{name}' not found")
        
        return self.presets[name]
    
    def suggest_groups_for_purpose(self, purpose: str) -> List[Tuple[str, str]]:
        """
        Suggest available groups for a given purpose.
        
        Returns:
            List of (group_name, reason) tuples
        """
        suggestions = []
        
        if purpose in self.group_mappings:
            candidates = self.group_mappings[purpose]
            
            for candidate in candidates:
                if self._validate_group(candidate):
                    suggestions.append((candidate, "Direct mapping"))
        
        # Add other potential matches based on purpose keywords
        purpose_keywords = {
            'protein': ['Protein', 'Protein_noH', 'C-alpha', 'Backbone'],
            'backbone': ['Backbone', 'MainChain', 'C-alpha'],
            'water': ['Water', 'SOL', 'HOH'],
            'system': ['System'],
            'center': ['Protein', 'Backbone'],
            'heavy': ['Protein_noH', 'Heavy']
        }
        
        for keyword, group_names in purpose_keywords.items():
            if keyword in purpose.lower():
                for group_name in group_names:
                    if self._validate_group(group_name) and group_name not in [s[0] for s in suggestions]:
                        suggestions.append((group_name, f"Keyword match: {keyword}"))
        
        return suggestions
    
    def _try_create_shortest_chain_index(self):
        """Strict shortest-chain detection required for safe complex processing."""
        try:
            from ..analysis.topology.strict_shortest_chain_detector import detect_shortest_chain_strict

            # Locate md.gro if available.
            md_gro_path = self._find_md_gro_file()

            logger.info("Starting strict shortest-chain detection")
            logger.info(f"   TPR file: {self.topology}")
            if md_gro_path:
                logger.info(f"   GRO file: {md_gro_path}")
            else:
                logger.info("   GRO file not found; using TPR only")

            # Run strict shortest-chain detection.
            group_id, index_file = detect_shortest_chain_strict(
                topology_file=self.topology,
                gro_file=md_gro_path,
                gmx_executable=self.gmx
            )

            self.shortest_chain_index_file = index_file
            self.shortest_chain_group_id = group_id

            # Force centering to use the detected shortest chain.
            self.group_mappings["center_protein"] = [group_id]

            logger.info("Shortest-chain detection succeeded")
            logger.info(f"   Group ID: {group_id}")
            logger.info(f"   Index file: {index_file}")
            logger.info("   Using the shortest chain for centering")

        except Exception as e:
            logger.error(f"Shortest-chain detection failed: {e}")
            logger.error("Complex systems require shortest-chain peptide detection for centering")
            logger.error("No fallback group is used because that would compromise scientific correctness")
            raise RuntimeError(f"Shortest-chain detection failed; complex processing is unsafe: {e}")

    def _find_md_gro_file(self) -> Optional[str]:
        """Find an md.gro file."""
        # Search the topology directory and the current working directory.
        topo_path = Path(self.topology)
        search_dirs = [topo_path.parent, Path.cwd()]

        for search_dir in search_dirs:
            md_gro = search_dir / "md.gro"
            if md_gro.exists() and md_gro.is_file():
                # Basic file validity check.
                try:
                    file_size = md_gro.stat().st_size
                    if file_size > 1000:
                        return str(md_gro)
                except Exception:
                    continue

        return None
    
    def _run_splitch_analysis(self, temp_ndx_file: str):
        """Run gmx make_ndx with splitch command to analyze chains."""
        cmd = [self.gmx, "make_ndx", "-f", self.topology, "-o", temp_ndx_file]
        
        # Use splitch command to split protein by chains
        # splitch will create separate groups for each chain
        input_commands = "splitch\nq\n"
        
        try:
            result = subprocess.run(
                cmd,
                input=input_commands,
                text=True,
                capture_output=True,
                check=True,
                timeout=60
            )
            
            logger.debug("Successfully ran splitch analysis")
            logger.debug(f"gmx make_ndx output: {result.stdout}")
            
            # Parse stdout to get chain information
            self._parse_splitch_output(result.stdout)
            
        except subprocess.TimeoutExpired:
            logger.warning("gmx make_ndx splitch timed out")
            raise
        except subprocess.CalledProcessError as e:
            logger.warning(f"gmx make_ndx splitch failed: {e}")
            logger.debug(f"stderr: {e.stderr}")
            raise
    
    def _parse_splitch_output(self, output: str):
        """Parse gmx make_ndx splitch output to extract chain information."""
        lines = output.split('\n')
        
        # Look for chain groups created by splitch
        # Format: "  18 ch0_Protein         :  1234 atoms"
        chain_pattern = re.compile(r'\s*(\d+)\s+(ch\d+_\w+)\s*:\s*(\d+)\s+atoms')
        
        for line in lines:
            match = chain_pattern.match(line)
            if match:
                group_id = int(match.group(1))
                chain_name = match.group(2)
                atom_count = int(match.group(3))
                
                self.peptide_chains[chain_name] = {
                    'group_id': group_id,
                    'atom_count': atom_count,
                    'name': chain_name
                }
                
                logger.info(f"Found chain: {chain_name} with {atom_count} atoms (group {group_id})")
        
        if self.peptide_chains:
            logger.info(f"Detected {len(self.peptide_chains)} peptide chains")
        else:
            logger.warning("No peptide chains found by splitch command")
    
    def _parse_chain_analysis(self, ndx_file: str):
        """Parse the generated index file to confirm chain groups."""
        try:
            with open(ndx_file, 'r') as f:
                content = f.read()
            
            # Look for chain group definitions in index file
            # Format: [ ch0_Protein ]
            chain_sections = re.findall(r'\[\s*(ch\d+_\w+)\s*\]', content)
            
            for chain_name in chain_sections:
                if chain_name not in self.peptide_chains:
                    # If we found a chain in index file but not in stdout parsing
                    logger.info(f"Found additional chain in index file: {chain_name}")
                    # Add with placeholder info, will be updated if needed
                    self.peptide_chains[chain_name] = {
                        'group_id': None,
                        'atom_count': 0,
                        'name': chain_name
                    }
                    
        except Exception as e:
            logger.debug(f"Could not parse index file: {e}")
    
    def _generate_shortest_chain_index(self):
        """Generate index file containing the shortest peptide chain."""
        if not self.peptide_chains:
            logger.warning("No peptide chains available for shortest chain selection")
            return
        
        # Find the chain with minimum atom count
        shortest_chain = min(
            self.peptide_chains.values(),
            key=lambda x: x['atom_count']
        )
        
        logger.info(f"Shortest chain: {shortest_chain['name']} with {shortest_chain['atom_count']} atoms")
        
        # Create index file for the shortest chain
        self._create_shortest_chain_index_file(shortest_chain)
    
    def _create_shortest_chain_index_file(self, shortest_chain: Dict):
        """Create index file specifically for the shortest chain."""
        try:
            # Create temporary index file
            temp_file = tempfile.NamedTemporaryFile(
                mode='w', suffix='.ndx', delete=False, 
                prefix='shortest_chain_'
            )
            
            # Generate index file using gmx make_ndx
            # We'll select the shortest chain group and save it as "Shortest_Chain"
            cmd = [self.gmx, "make_ndx", "-f", self.topology, "-o", temp_file.name]
            
            # Commands: select the shortest chain group, rename it, and quit
            chain_group_id = shortest_chain['group_id']
            new_group_id = len(self.available_groups) + 18
            
            if chain_group_id is not None:
                input_commands = f"splitch\n{chain_group_id}\nname {new_group_id} Shortest_Chain\nq\n"
            else:
                # Fallback: try splitch and select first chain
                input_commands = f"splitch\nch0_Protein\nname {new_group_id} Shortest_Chain\nq\n"
            
            result = subprocess.run(
                cmd,
                input=input_commands,
                text=True,
                capture_output=True,
                check=True,
                timeout=30
            )
            
            self.shortest_chain_index_file = temp_file.name
            temp_file.close()
            
            # Add to available groups
            new_group_id = len(self.available_groups) + 18
            self.available_groups[new_group_id] = {
                'name': 'Shortest_Chain',
                'atoms': shortest_chain['atom_count']
            }
            
            # Update group mappings to use shortest chain for centering
            self.group_mappings["center_shortest_chain"] = [str(new_group_id)]
            
            logger.info(f"Generated shortest chain index file: {self.shortest_chain_index_file}")
            logger.info(f"Shortest chain available as group {new_group_id}: Shortest_Chain")
            
        except Exception as e:
            logger.error(f"Failed to create shortest chain index file: {e}")
            self.shortest_chain_index_file = None
    
    def get_shortest_chain_index_file(self) -> Optional[str]:
        """Get path to the shortest chain index file."""
        return self.shortest_chain_index_file
    
    def get_shortest_chain_info(self) -> Optional[Dict]:
        """Get information about the shortest peptide chain."""
        if not self.peptide_chains:
            return None
        
        shortest_chain = min(
            self.peptide_chains.values(),
            key=lambda x: x['atom_count']
        )
        return shortest_chain
    
    def has_shortest_chain_group(self) -> bool:
        """Check if shortest chain group is available."""
        return self.shortest_chain_index_file is not None and Path(self.shortest_chain_index_file).exists()
    
    def _generate_all_chains_index(self) -> Optional[str]:
        """Generate index file containing all peptide chains using splitch."""
        try:
            # Create temporary file for all chains
            temp_file = tempfile.NamedTemporaryFile(
                mode='w', suffix='.ndx', delete=False, 
                prefix='all_chains_'
            )
            temp_file.close()
            
            # Run gmx make_ndx with splitch 1 to split protein into chains
            cmd = [self.gmx, "make_ndx", "-f", self.topology, "-o", temp_file.name]
            
            # Commands: splitch protein group (group 1), confirm with enter, quit
            input_commands = "splitch 1\n\nq\n"
            
            result = subprocess.run(
                cmd,
                input=input_commands,
                text=True,
                capture_output=True,
                check=True,
                timeout=60
            )
            
            logger.info(f"Generated all chains index file: {temp_file.name}")
            logger.debug(f"gmx make_ndx output: {result.stdout}")
            
            return temp_file.name
            
        except subprocess.TimeoutExpired:
            logger.warning("gmx make_ndx (all chains) timed out")
            return None
        except subprocess.CalledProcessError as e:
            logger.warning(f"gmx make_ndx (all chains) failed: {e}")
            logger.debug(f"stderr: {e.stderr}")
            return None
        except Exception as e:
            logger.error(f"Failed to generate all chains index: {e}")
            return None
    
    def _parse_chains_from_index_file(self, index_file: str) -> Dict[str, int]:
        """Parse index file to extract chain information and atom counts."""
        chain_info = {}
        
        try:
            with open(index_file, 'r') as f:
                content = f.read()
            
            # Find all chain group sections
            # Pattern: [ Protein_chain_A ] or similar
            sections = re.split(r'\[\s*([^\]]+)\s*\]', content)
            
            for i in range(1, len(sections), 2):  # Skip empty sections
                group_name = sections[i].strip()
                if i + 1 < len(sections):
                    atom_data = sections[i + 1].strip()
                    
                    # Log all non-empty groups for debugging
                    atom_count = self._count_atoms_in_group_data(atom_data)
                    if atom_count > 0:
                        logger.debug(f"Found group: {group_name} with {atom_count} atoms")
                        
                        # Check if this looks like a protein chain
                        if self._is_chain_group_name(group_name):
                            chain_info[group_name] = atom_count
                            logger.info(f"Identified as chain: {group_name} with {atom_count} atoms")
                        else:
                            logger.debug(f"Not identified as chain: {group_name}")
            
            logger.info(f"Analyzed {len(chain_info)} chains from index file")
            return chain_info
            
        except Exception as e:
            logger.error(f"Failed to parse chains from index file: {e}")
            return {}
    
    def _is_chain_group_name(self, group_name: str) -> bool:
        """Check if group name represents a protein chain from splitch command."""
        # splitch 1 can generate various naming patterns depending on GROMACS version and system
        group_lower = group_name.lower()
        
        # Common splitch-generated chain patterns
        chain_patterns = [
            r'^protein_chain_[a-z]$',       # Protein_chain_A, Protein_chain_B, ...
            r'^protein_chain_\d+$',         # Protein_chain_1, Protein_chain_2, ...
            r'^protein_chain[a-z]$',        # Protein_chainA, Protein_chainB, ...
            r'^protein_chain\d+$',          # Protein_chain1, Protein_chain2, ...
            r'^ch\d+_protein$',             # ch0_Protein, ch1_Protein, ...
            r'^chain_[a-z]$',               # Chain_A, Chain_B, ...
            r'^chain_\d+$',                 # Chain_1, Chain_2, ...
            r'^protein[a-z]$',              # ProteinA, ProteinB, ...
            r'^protein\d+$',                # Protein1, Protein2, ...
        ]
        
        # Check all patterns
        for pattern in chain_patterns:
            if re.match(pattern, group_lower):
                logger.debug(f"Chain group matched pattern '{pattern}': {group_name}")
                return True
        
        # Additional flexible matching: any group containing both "protein" and "chain"
        if 'protein' in group_lower and 'chain' in group_lower:
            logger.debug(f"Chain group matched flexible pattern (protein+chain): {group_name}")
            return True
        
        # Check if it's not an obvious non-chain group
        non_chain_keywords = ['system', 'water', 'sol', 'ion', 'sod', 'cla', 'backbone', 'mainchain']
        for keyword in non_chain_keywords:
            if keyword in group_lower:
                return False
        
        return False
    
    def _count_atoms_in_group_data(self, atom_data: str) -> int:
        """Count number of atoms from group data section."""
        try:
            # Split by whitespace and count valid atom indices
            atom_indices = atom_data.split()
            # Filter out non-numeric entries
            valid_atoms = [x for x in atom_indices if x.isdigit()]
            return len(valid_atoms)
        except Exception:
            return 0
    
    def _find_shortest_valid_chain(self, chain_info: Dict[str, int]) -> Optional[Tuple[str, int]]:
        """Find the shortest valid peptide chain, excluding ions and small molecules."""
        if not chain_info:
            return None
        
        valid_chains = []
        
        for chain_name, atom_count in chain_info.items():
            if self._is_valid_peptide_chain(chain_name, atom_count):
                valid_chains.append((chain_name, atom_count))
        
        if not valid_chains:
            logger.warning("No valid peptide chains found")
            return None
        
        # Find shortest among valid chains
        shortest = min(valid_chains, key=lambda x: x[1])
        logger.info(f"Selected shortest valid chain: {shortest[0]} with {shortest[1]} atoms")
        
        return shortest
    
    def _is_valid_peptide_chain(self, chain_name: str, atom_count: int) -> bool:
        """Check if this is a valid peptide chain from splitch command."""
        
        # Since splitch 1 only splits protein chains, groups matching our patterns are valid
        # Just need basic sanity checks
        
        # 1. Minimum atom count (avoid empty or tiny groups)
        if atom_count < 10:  # Very small threshold since splitch should give valid chains
            logger.debug(f"Chain {chain_name} has too few atoms ({atom_count})")
            return False
        
        # 2. Maximum atom count (avoid errors)
        if atom_count > 100000:  # Very large threshold
            logger.debug(f"Chain {chain_name} has suspiciously many atoms ({atom_count})")
            return False
        
        # 3. If it matches our chain pattern, it's valid (splitch already filtered for us)
        return self._is_chain_group_name(chain_name)
    
    def _find_group_id_in_index(self, index_file: str, target_chain_name: str) -> Optional[int]:
        """Find the group ID for a specific chain name in the index file."""
        try:
            with open(index_file, 'r') as f:
                content = f.read()
            
            # Parse the index file to find group IDs
            sections = re.split(r'\[\s*([^\]]+)\s*\]', content)
            group_id = 0
            
            for i in range(1, len(sections), 2):  # Skip empty sections
                if i + 1 < len(sections):
                    group_name = sections[i].strip()
                    atom_data = sections[i + 1].strip()
                    
                    if atom_data:  # Only count non-empty groups
                        if group_name == target_chain_name:
                            logger.info(f"Found target chain '{target_chain_name}' at group ID {group_id}")
                            return group_id
                        group_id += 1
            
            logger.warning(f"Target chain '{target_chain_name}' not found in index file")
            return None
            
        except Exception as e:
            logger.error(f"Failed to find group ID in index file: {e}")
            return None
    
    def _generate_shortest_chain_only_index(self, shortest_chain: Tuple[str, int], all_chains_index: str):
        """Generate final index file containing only the shortest chain."""
        try:
            chain_name, atom_count = shortest_chain
            
            # Create final index file
            final_file = tempfile.NamedTemporaryFile(
                mode='w', suffix='.ndx', delete=False, 
                prefix='shortest_chain_'
            )
            final_file.close()
            
            # Run gmx make_ndx to select only the shortest chain
            cmd = [self.gmx, "make_ndx", "-f", self.topology, "-n", all_chains_index, "-o", final_file.name]
            
            # First find the group number for our target chain by parsing the all_chains_index
            target_group_id = self._find_group_id_in_index(all_chains_index, chain_name)
            
            if target_group_id is None:
                logger.error(f"Could not find group ID for chain {chain_name}")
                return
            
            # Commands: select the target chain group, rename it, quit
            input_commands = f"{target_group_id}\nname {target_group_id} Shortest_Chain\nq\n"
            
            result = subprocess.run(
                cmd,
                input=input_commands,
                text=True,
                capture_output=True,
                check=True,
                timeout=30
            )
            
            self.shortest_chain_index_file = final_file.name
            
            # Parse the generated index file to get accurate group mappings
            self._parse_generated_index_mappings(final_file.name)
            
            logger.info(f"Generated shortest chain index file: {final_file.name}")
            logger.info(f"Shortest chain '{chain_name}' ready for use")
            
        except Exception as e:
            logger.error(f"Failed to generate shortest chain only index: {e}")
            self.shortest_chain_index_file = None
    
    def _parse_generated_index_mappings(self, index_file: str):
        """Parse the generated index file to create accurate group ID mappings."""
        try:
            with open(index_file, 'r') as f:
                content = f.read()
            
            # Find all group sections and their actual positions
            group_mappings = {}
            group_id = 0
            
            # Split content by group headers
            sections = re.split(r'\[\s*([^\]]+)\s*\]', content)
            
            for i in range(1, len(sections), 2):  # Skip empty sections
                if i + 1 < len(sections):
                    group_name = sections[i].strip()
                    atom_data = sections[i + 1].strip()
                    
                    if atom_data:  # Only count non-empty groups
                        group_mappings[group_name] = str(group_id)
                        
                        # Update available groups
                        if group_name == "Shortest_Chain":
                            atom_count = self._count_atoms_in_group_data(atom_data)
                            self.available_groups[group_id] = {
                                'name': 'Shortest_Chain',
                                'atoms': atom_count
                            }
                            
                            # Update group mappings for center selection
                            self.group_mappings["center_shortest_chain"] = [str(group_id)]
                            
                            logger.info(f"Shortest_Chain mapped to group {group_id} with {atom_count} atoms")
                        
                        group_id += 1
            
            # Store the complete mapping for reference
            self.index_file_mappings = group_mappings
            
            logger.debug(f"Generated index file mappings: {group_mappings}")
            
        except Exception as e:
            logger.error(f"Failed to parse generated index mappings: {e}")
    
    def get_group_id_from_index(self, group_name: str) -> Optional[str]:
        """Get the actual group ID for a group name from the generated index file."""
        if hasattr(self, 'index_file_mappings') and group_name in self.index_file_mappings:
            return self.index_file_mappings[group_name]
        return None
    
    def get_shortest_chain_group_id(self) -> Optional[str]:
        """Get the group ID for the shortest chain."""
        return self.shortest_chain_group_id

    def get_shortest_chain_id(self) -> Optional[str]:
        """Compatibility alias returning the shortest-chain group ID."""
        return self.get_shortest_chain_group_id()

    def get_direct_shortest_chain_group_id(self) -> Optional[str]:
        """Get the group ID for the shortest chain from the generated index file."""
        return self.get_group_id_from_index("Shortest_Chain")
    
    def _find_group_id_in_index(self, index_file: str, target_chain_name: str) -> Optional[int]:
        """Find the group ID number for a specific chain name in the index file."""
        try:
            with open(index_file, 'r') as f:
                content = f.read()
            
            # Look for the target chain section
            # Pattern: [ target_chain_name ]
            escaped_name = re.escape(target_chain_name)
            pattern = rf'\[\s*{escaped_name}\s*\]'
            
            if re.search(pattern, content, re.IGNORECASE):
                # If we find the section, we need to determine its group number
                # This is tricky because gmx make_ndx assigns numbers sequentially
                # Let's count how many groups appear before our target
                
                all_groups = re.findall(r'\[\s*([^\]]+)\s*\]', content)
                for i, group_name in enumerate(all_groups):
                    if group_name.strip().lower() == target_chain_name.lower():
                        # Group numbers start from 0, so the index is the group number
                        return i
            
            logger.warning(f"Could not find group ID for chain {target_chain_name}")
            return None
            
        except Exception as e:
            logger.error(f"Error finding group ID: {e}")
            return None
