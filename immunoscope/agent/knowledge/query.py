"""
Knowledge base query interface for TCR-pMHC mutation design.

Provides structured access to domain knowledge including:
- System knowledge (pMHC-TCR structure and function)
- Design principles (mutation strategies)
- Region definitions (CDR, framework, interface zones)
- Mutation rules (design patterns and success rates)
- Interaction types (H-bond, π-π, salt bridge, etc.)
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass


KNOWLEDGE_DIR = Path(__file__).parent


@dataclass
class RegionInfo:
    """Information about a TCR or peptide region."""
    name: str
    full_name: str
    function: str
    mutation_tolerance: str
    design_priority: str
    design_notes: str
    recommended_mutations: List[str]
    avoid_mutations: List[str]
    typical_rmsf: Optional[str] = None
    typical_length: Optional[str] = None


@dataclass
class MutationRule:
    """A mutation design rule with rationale and examples."""
    rule_id: str
    pattern: str
    region: str
    interaction_type: str
    rationale: str
    success_rate: str
    success_rate_numeric: Optional[float]
    affinity_impact: str
    specificity_impact: str
    examples: List[Dict[str, Any]]
    when_to_use: List[str]
    cautions: List[str]
    contraindications: List[str]


@dataclass
class InteractionType:
    """Information about an interaction type."""
    name: str
    abbreviation: str
    definition: str
    typical_distance: str
    strength: str
    common_residue_pairs: List[Dict[str, Any]]
    design_strategies: Dict[str, str]
    cautions: List[str]
    interpretation: Dict[str, str]


class KnowledgeBase:
    """Query interface for TCR-pMHC mutation design knowledge."""

    def __init__(self):
        self._regions = None
        self._mutation_rules = None
        self._interaction_types = None
        self._system_knowledge = None
        self._design_principles = None
        self._peptide_design_principles = None
        self._amino_acid_properties = None
        self._substitution_strategies = None
        self._case_studies = None

    @property
    def regions(self) -> Dict[str, Dict]:
        """Load region definitions."""
        if self._regions is None:
            with open(KNOWLEDGE_DIR / "regions.json") as f:
                data = json.load(f)
                self._regions = data["tcr_regions"]
        return self._regions

    @property
    def mutation_rules(self) -> List[Dict]:
        """Load mutation rules."""
        if self._mutation_rules is None:
            with open(KNOWLEDGE_DIR / "mutation_rules.json") as f:
                data = json.load(f)
                self._mutation_rules = data["mutation_rules"]
        return self._mutation_rules

    @property
    def interaction_types(self) -> Dict[str, Dict]:
        """Load interaction type definitions."""
        if self._interaction_types is None:
            with open(KNOWLEDGE_DIR / "interaction_types.json") as f:
                data = json.load(f)
                self._interaction_types = data["interaction_types"]
        return self._interaction_types

    @property
    def system_knowledge(self) -> str:
        """Load system knowledge markdown."""
        if self._system_knowledge is None:
            with open(KNOWLEDGE_DIR / "system_knowledge.md") as f:
                self._system_knowledge = f.read()
        return self._system_knowledge

    @property
    def design_principles(self) -> str:
        """Load design principles markdown."""
        if self._design_principles is None:
            with open(KNOWLEDGE_DIR / "design_principles.md") as f:
                self._design_principles = f.read()
        return self._design_principles

    @property
    def peptide_design_principles(self) -> str:
        """Load peptide design principles markdown."""
        if self._peptide_design_principles is None:
            with open(KNOWLEDGE_DIR / "peptide_design_principles.md") as f:
                self._peptide_design_principles = f.read()
        return self._peptide_design_principles

    @property
    def amino_acid_properties(self) -> Dict[str, Any]:
        """Load amino acid properties."""
        if self._amino_acid_properties is None:
            with open(KNOWLEDGE_DIR / "amino_acid_properties.json") as f:
                self._amino_acid_properties = json.load(f)
        return self._amino_acid_properties

    @property
    def substitution_strategies(self) -> Dict[str, Any]:
        """Load substitution strategies."""
        if self._substitution_strategies is None:
            with open(KNOWLEDGE_DIR / "substitution_strategies.json") as f:
                self._substitution_strategies = json.load(f)
        return self._substitution_strategies

    @property
    def case_studies(self) -> Dict[str, Any]:
        """Load case studies."""
        if self._case_studies is None:
            with open(KNOWLEDGE_DIR / "case_studies.json") as f:
                self._case_studies = json.load(f)
        return self._case_studies

    # Query methods

    def get_region_info(self, region_name: str) -> Optional[RegionInfo]:
        """Get information about a specific region."""
        region_data = self.regions.get(region_name)
        if not region_data:
            return None

        return RegionInfo(
            name=region_name,
            full_name=region_data.get("full_name", ""),
            function=region_data.get("function", ""),
            mutation_tolerance=region_data.get("mutation_tolerance", ""),
            design_priority=region_data.get("design_priority", ""),
            design_notes=region_data.get("design_notes", ""),
            recommended_mutations=region_data.get("recommended_mutations", []),
            avoid_mutations=region_data.get("avoid_mutations", []),
            typical_rmsf=region_data.get("typical_rmsf"),
            typical_length=region_data.get("typical_length")
        )

    def get_mutation_rules_for_region(self, region: str) -> List[MutationRule]:
        """Get all mutation rules applicable to a region."""
        rules = []
        for rule_data in self.mutation_rules:
            if region in rule_data.get("region", ""):
                rules.append(MutationRule(
                    rule_id=rule_data["rule_id"],
                    pattern=rule_data["pattern"],
                    region=rule_data["region"],
                    interaction_type=rule_data["interaction_type"],
                    rationale=rule_data["rationale"],
                    success_rate=rule_data["success_rate"],
                    success_rate_numeric=rule_data.get("success_rate_numeric"),
                    affinity_impact=rule_data["affinity_impact"],
                    specificity_impact=rule_data["specificity_impact"],
                    examples=rule_data["examples"],
                    when_to_use=rule_data["when_to_use"],
                    cautions=rule_data["cautions"],
                    contraindications=rule_data["contraindications"]
                ))
        return rules

    def get_mutation_rules_for_interaction(self, interaction_type: str) -> List[MutationRule]:
        """Get mutation rules for a specific interaction type."""
        rules = []
        for rule_data in self.mutation_rules:
            if interaction_type in rule_data.get("interaction_type", ""):
                rules.append(MutationRule(
                    rule_id=rule_data["rule_id"],
                    pattern=rule_data["pattern"],
                    region=rule_data["region"],
                    interaction_type=rule_data["interaction_type"],
                    rationale=rule_data["rationale"],
                    success_rate=rule_data["success_rate"],
                    success_rate_numeric=rule_data.get("success_rate_numeric"),
                    affinity_impact=rule_data["affinity_impact"],
                    specificity_impact=rule_data["specificity_impact"],
                    examples=rule_data["examples"],
                    when_to_use=rule_data["when_to_use"],
                    cautions=rule_data["cautions"],
                    contraindications=rule_data["contraindications"]
                ))
        return rules

    def get_interaction_info(self, interaction_type: str) -> Optional[InteractionType]:
        """Get information about an interaction type."""
        int_data = self.interaction_types.get(interaction_type)
        if not int_data:
            return None

        return InteractionType(
            name=int_data["name"],
            abbreviation=int_data["abbreviation"],
            definition=int_data["definition"],
            typical_distance=int_data["typical_distance"],
            strength=int_data["strength"],
            common_residue_pairs=int_data["common_residue_pairs"],
            design_strategies=int_data["design_strategies"],
            cautions=int_data["cautions"],
            interpretation=int_data["interpretation"]
        )

    def find_rules_by_criteria(
        self,
        region: Optional[str] = None,
        interaction_type: Optional[str] = None,
        min_success_rate: Optional[float] = None
    ) -> List[MutationRule]:
        """Find mutation rules matching criteria."""
        rules = []
        for rule_data in self.mutation_rules:
            # Filter by region
            if region and region not in rule_data.get("region", ""):
                continue

            # Filter by interaction type
            if interaction_type and interaction_type not in rule_data.get("interaction_type", ""):
                continue

            # Filter by success rate
            if min_success_rate is not None:
                numeric_rate = rule_data.get("success_rate_numeric")
                if numeric_rate is None or numeric_rate < min_success_rate:
                    continue

            rules.append(MutationRule(
                rule_id=rule_data["rule_id"],
                pattern=rule_data["pattern"],
                region=rule_data["region"],
                interaction_type=rule_data["interaction_type"],
                rationale=rule_data["rationale"],
                success_rate=rule_data["success_rate"],
                success_rate_numeric=rule_data.get("success_rate_numeric"),
                affinity_impact=rule_data["affinity_impact"],
                specificity_impact=rule_data["specificity_impact"],
                examples=rule_data["examples"],
                when_to_use=rule_data["when_to_use"],
                cautions=rule_data["cautions"],
                contraindications=rule_data["contraindications"]
            ))

        return rules

    def get_design_strategy(self, goal: str) -> Optional[Dict[str, Any]]:
        """Get design strategy for a specific goal."""
        with open(KNOWLEDGE_DIR / "mutation_rules.json") as f:
            data = json.load(f)
            strategies = data.get("mutation_strategies", {})
            return strategies.get(goal)

    def interpret_metric(self, metric_name: str, value: float) -> str:
        """Interpret a metric value (RRCS, occupancy, BSA, RMSF)."""
        with open(KNOWLEDGE_DIR / "interaction_types.json") as f:
            data = json.load(f)
            metric_info = data.get("metric_interpretation", {}).get(metric_name)

            if not metric_info:
                return f"Unknown metric: {metric_name}"

            interpretation = metric_info.get("interpretation", {})

            # Find matching interpretation
            for key, desc in interpretation.items():
                if ">" in desc:
                    threshold = float(desc.split(">")[1].split()[0])
                    if value > threshold:
                        return desc
                elif "<" in desc:
                    threshold = float(desc.split("<")[1].split()[0])
                    if value < threshold:
                        return desc
                elif "-" in desc and desc[0].isdigit():
                    parts = desc.split("-")
                    low = float(parts[0])
                    high = float(parts[1].split()[0])
                    if low <= value <= high:
                        return desc

            return f"{metric_name} = {value}"

    def get_system_knowledge_section(self, section: str) -> str:
        """Extract a specific section from system knowledge."""
        content = self.system_knowledge
        lines = content.split("\n")

        # Find section header
        section_start = None
        section_level = None
        for i, line in enumerate(lines):
            if line.startswith("#") and section.lower() in line.lower():
                section_start = i
                section_level = len(line.split()[0])  # Count # characters
                break

        if section_start is None:
            return f"Section '{section}' not found"

        # Extract until next section of same or higher level
        section_lines = [lines[section_start]]
        for i in range(section_start + 1, len(lines)):
            line = lines[i]
            if line.startswith("#"):
                current_level = len(line.split()[0])
                if current_level <= section_level:
                    break
            section_lines.append(line)

        return "\n".join(section_lines)

    def get_design_principles_section(self, section: str) -> str:
        """Extract a specific section from design principles."""
        content = self.design_principles
        lines = content.split("\n")

        # Find section header
        section_start = None
        section_level = None
        for i, line in enumerate(lines):
            if line.startswith("#") and section.lower() in line.lower():
                section_start = i
                section_level = len(line.split()[0])
                break

        if section_start is None:
            return f"Section '{section}' not found"

        # Extract until next section of same or higher level
        section_lines = [lines[section_start]]
        for i in range(section_start + 1, len(lines)):
            line = lines[i]
            if line.startswith("#"):
                current_level = len(line.split()[0])
                if current_level <= section_level:
                    break
            section_lines.append(line)

        return "\n".join(section_lines)

    def get_peptide_design_section(self, section: str) -> str:
        """Extract a specific section from peptide design principles."""
        content = self.peptide_design_principles
        lines = content.split("\n")

        # Find section header
        section_start = None
        section_level = None
        for i, line in enumerate(lines):
            if line.startswith("#") and section.lower() in line.lower():
                section_start = i
                section_level = len(line.split()[0])
                break

        if section_start is None:
            return f"Section '{section}' not found"

        # Extract until next section of same or higher level
        section_lines = [lines[section_start]]
        for i in range(section_start + 1, len(lines)):
            line = lines[i]
            if line.startswith("#"):
                current_level = len(line.split()[0])
                if current_level <= section_level:
                    break
            section_lines.append(line)

        return "\n".join(section_lines)

    # Amino acid knowledge query methods

    def _normalize_residue(self, residue: str) -> str:
        """Normalize residue code to 3-letter uppercase."""
        one_to_three = {
            "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
            "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
            "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
            "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
        }
        r = residue.upper().strip()
        return one_to_three.get(r, r)

    def get_amino_acid_properties(self, residue: str) -> Optional[Dict[str, Any]]:
        """Get complete properties for an amino acid.

        Args:
            residue: Amino acid code (1-letter or 3-letter, e.g., "Y" or "TYR")

        Returns:
            Dictionary with properties, classification, and substitution info,
            or None if residue not found.
        """
        code = self._normalize_residue(residue)
        return self.amino_acid_properties["amino_acids"].get(code)

    def get_conservative_substitutions(self, residue: str) -> List[str]:
        """Get conservative substitution options for an amino acid.

        Args:
            residue: Amino acid code (1-letter or 3-letter)

        Returns:
            List of 3-letter codes for conservative substitutions.
        """
        info = self.get_amino_acid_properties(residue)
        if info is None:
            return []
        return info.get("conservative_substitutions", [])

    def get_avoid_substitutions(self, residue: str) -> List[str]:
        """Get substitutions to avoid for an amino acid.

        Args:
            residue: Amino acid code (1-letter or 3-letter)

        Returns:
            List of 3-letter codes for substitutions to avoid.
        """
        info = self.get_amino_acid_properties(residue)
        if info is None:
            return []
        return info.get("avoid_substitutions", [])

    def get_substitution_score(self, res1: str, res2: str) -> Optional[float]:
        """Get BLOSUM-like substitution score between two amino acids.

        Higher scores indicate more conservative (favorable) substitutions.
        Scores typically range from -4 (very unfavorable) to +11 (identical).

        Args:
            res1: Original amino acid code
            res2: Substitute amino acid code

        Returns:
            Substitution score, or None if matrix unavailable.
        """
        code1 = self._normalize_residue(res1)
        code2 = self._normalize_residue(res2)
        matrix_data = self.amino_acid_properties.get("substitution_matrix")
        if matrix_data is None:
            return None
        # Matrix scores are nested under "scores" key
        matrix = matrix_data.get("scores", matrix_data)
        row = matrix.get(code1)
        if row is None:
            return None
        return row.get(code2)

    def find_substitution_by_property(
        self,
        original: str,
        maintain_properties: List[str],
    ) -> List[str]:
        """Find substitutions that maintain specific properties.

        Args:
            original: Original amino acid (1-letter or 3-letter)
            maintain_properties: Properties to maintain, e.g.,
                ["polarity", "charge", "aromaticity", "size"]

        Returns:
            List of 3-letter codes that maintain the specified properties.
        """
        orig_info = self.get_amino_acid_properties(original)
        if orig_info is None:
            return []

        orig_props = orig_info.get("properties", {})
        candidates = []

        for code, info in self.amino_acid_properties["amino_acids"].items():
            if code == orig_info["three_letter"]:
                continue
            cand_props = info.get("properties", {})
            matches = True
            for prop in maintain_properties:
                if prop not in orig_props or prop not in cand_props:
                    matches = False
                    break
                if orig_props[prop] != cand_props[prop]:
                    matches = False
                    break
            if matches:
                candidates.append(code)
        return candidates

    def get_substitution_strategy(self, rule_id: str) -> Optional[Dict[str, Any]]:
        """Get a specific substitution strategy by rule_id.

        Args:
            rule_id: e.g., "hydrophobic_conservative", "aromatic_tuning",
                "charged_conservative", "polar_conservative", "size_matching"

        Returns:
            Strategy dictionary with pattern, rationale, success rate,
            and examples, or None if not found.
        """
        for rule in self.substitution_strategies.get("property_based_rules", []):
            if rule.get("rule_id") == rule_id:
                return rule
        return None

    def get_forbidden_substitutions(self) -> List[Dict[str, Any]]:
        """Get list of forbidden substitution patterns.

        Returns:
            List of forbidden substitution dictionaries with pattern,
            reason, and failure rate.
        """
        return self.substitution_strategies.get("forbidden_substitutions", [])

    def get_case_studies(
        self,
        mutation_type: Optional[str] = None,
        region: Optional[str] = None,
        outcome: Optional[str] = None,
        design_principle: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Query case study database.

        Args:
            mutation_type: Filter by type, e.g., "TCR" or "peptide"
            region: Filter by region, e.g., "CDR3alpha", "CDR3beta"
            outcome: Filter by outcome, "success" or "failure"
            design_principle: Filter by design principle,
                e.g., "aromatic_tuning", "charged_conservative"

        Returns:
            List of matching case study dictionaries.
        """
        cases = []
        if outcome is None or outcome == "success":
            cases.extend(self.case_studies.get("success_cases", []))
        if outcome is None or outcome == "failure":
            cases.extend(self.case_studies.get("failure_cases", []))

        filtered = []
        for case in cases:
            if mutation_type and case.get("mutation_type") != mutation_type:
                continue
            if region and region not in case.get("region", ""):
                continue
            if design_principle and case.get("design_principle") != design_principle:
                continue
            filtered.append(case)
        return filtered

    def explain_substitution(self, original: str, mutant: str) -> Dict[str, Any]:
        """Explain why a substitution might work or fail.

        Combines amino acid properties, BLOSUM score, and case studies
        to provide a comprehensive assessment.

        Args:
            original: Original amino acid (1-letter or 3-letter)
            mutant: Substitute amino acid (1-letter or 3-letter)

        Returns:
            Dictionary with assessment fields:
                - original / mutant: 3-letter codes
                - blosum_score: Substitution score
                - is_conservative: Boolean
                - is_forbidden: Boolean
                - property_changes: Dict of changed properties
                - maintained_properties: List of preserved properties
                - relevant_cases: Matching case studies
                - verdict: Summary assessment string
        """
        orig_code = self._normalize_residue(original)
        mut_code = self._normalize_residue(mutant)

        orig_info = self.get_amino_acid_properties(orig_code)
        mut_info = self.get_amino_acid_properties(mut_code)

        if orig_info is None or mut_info is None:
            return {
                "original": orig_code,
                "mutant": mut_code,
                "error": "Unknown amino acid code",
            }

        score = self.get_substitution_score(orig_code, mut_code)
        conservative_subs = orig_info.get("conservative_substitutions", [])
        avoid_subs = orig_info.get("avoid_substitutions", [])
        is_conservative = mut_code in conservative_subs
        is_forbidden = mut_code in avoid_subs

        # Compare properties
        orig_props = orig_info.get("properties", {})
        mut_props = mut_info.get("properties", {})
        changes = {}
        maintained = []
        for key in orig_props:
            if key not in mut_props:
                continue
            if orig_props[key] != mut_props[key]:
                changes[key] = {"from": orig_props[key], "to": mut_props[key]}
            else:
                maintained.append(key)

        # Find relevant cases
        relevant_cases = []
        for case in self.case_studies.get("success_cases", []) + self.case_studies.get("failure_cases", []):
            orig_res = case.get("original_residue", "")
            mut_res = case.get("mutant_residue", "")
            if orig_res == orig_code and mut_res == mut_code:
                relevant_cases.append(case)

        # Build verdict
        if is_forbidden:
            verdict = f"AVOID: {orig_code}->{mut_code} is a high-risk substitution"
        elif is_conservative:
            verdict = f"RECOMMENDED: {orig_code}->{mut_code} is a conservative substitution"
        elif score is not None and score >= 1:
            verdict = f"ACCEPTABLE: {orig_code}->{mut_code} has positive BLOSUM score ({score})"
        elif score is not None and score < 0:
            verdict = f"CAUTION: {orig_code}->{mut_code} has negative BLOSUM score ({score})"
        else:
            verdict = f"NEUTRAL: {orig_code}->{mut_code} has no strong recommendation"

        return {
            "original": orig_code,
            "mutant": mut_code,
            "blosum_score": score,
            "is_conservative": is_conservative,
            "is_forbidden": is_forbidden,
            "property_changes": changes,
            "maintained_properties": maintained,
            "relevant_cases": relevant_cases,
            "verdict": verdict,
        }


# Global instance
_kb = None

def get_knowledge_base() -> KnowledgeBase:
    """Get global knowledge base instance."""
    global _kb
    if _kb is None:
        _kb = KnowledgeBase()
    return _kb
