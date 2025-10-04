"""
Filter classes for filtering molecules based on structural patterns and properties.
"""

from typing import List, Dict, Optional, Union, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors

from .properties import MolecularProperties
from .data.pains_patterns import PAINS_PATTERNS
from .data.reactive_patterns import REACTIVE_PATTERNS
from .data.heterocycle_patterns import HETEROCYCLE_PATTERNS


class BaseFilter:
    """Base class for molecular filters."""

    def __init__(self, patterns: Dict[str, str]):
        """
        Initialize the filter with SMARTS patterns.
        
        Args:
            patterns: Dictionary mapping pattern names to SMARTS strings
        """
        self.patterns = patterns
        self._compiled_patterns = {}
        self._compile_patterns()

    def _compile_patterns(self):
        """Compile SMARTS patterns for efficient matching."""
        for name, smarts in self.patterns.items():
            try:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern is not None:
                    self._compiled_patterns[name] = pattern
            except Exception as e:
                print(f"Warning: Could not compile pattern '{name}': {e}")

    def check_molecule(self, mol: Chem.Mol) -> Tuple[bool, List[str]]:
        """
        Check if a molecule matches any filter patterns.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            Tuple of (passes_filter, list_of_matched_patterns)
            passes_filter is True if molecule passes (no matches), False if it fails
        """
        if mol is None:
            return False, ["Invalid molecule"]

        matched_patterns = []
        for name, pattern in self._compiled_patterns.items():
            if mol.HasSubstructMatch(pattern):
                matched_patterns.append(name)

        passes = len(matched_patterns) == 0
        return passes, matched_patterns

    def filter_molecules(
        self, molecules: List[Chem.Mol]
    ) -> Tuple[List[Chem.Mol], List[Chem.Mol], Dict[int, List[str]]]:
        """
        Filter a list of molecules.
        
        Args:
            molecules: List of RDKit molecule objects
            
        Returns:
            Tuple of (passed_molecules, failed_molecules, failure_reasons)
            failure_reasons is a dict mapping failed molecule indices to lists of matched patterns
        """
        passed = []
        failed = []
        failure_reasons = {}

        for idx, mol in enumerate(molecules):
            passes, matched = self.check_molecule(mol)
            if passes:
                passed.append(mol)
            else:
                failed.append(mol)
                failure_reasons[idx] = matched

        return passed, failed, failure_reasons


class PAINSFilter(BaseFilter):
    """
    Filter for Pan-Assay Interference Compounds (PAINS).
    
    PAINS are compounds that frequently show up as hits in high-throughput screens
    but are likely to be false positives due to various interference mechanisms.
    """

    def __init__(self, custom_patterns: Optional[Dict[str, str]] = None):
        """
        Initialize PAINS filter.
        
        Args:
            custom_patterns: Optional custom SMARTS patterns to use instead of defaults
        """
        patterns = custom_patterns if custom_patterns is not None else PAINS_PATTERNS
        super().__init__(patterns)


class ReactiveFilter(BaseFilter):
    """
    Filter for reactive functional groups.
    
    Identifies chemically reactive groups that may cause issues in biological assays
    or drug development.
    """

    def __init__(self, custom_patterns: Optional[Dict[str, str]] = None):
        """
        Initialize reactive group filter.
        
        Args:
            custom_patterns: Optional custom SMARTS patterns to use instead of defaults
        """
        patterns = custom_patterns if custom_patterns is not None else REACTIVE_PATTERNS
        super().__init__(patterns)


class HeterocycleFilter(BaseFilter):
    """
    Filter for heterocyclic scaffolds.
    
    This filter can be used to identify or require the presence of specific
    heterocyclic patterns common in medicinal chemistry.
    """

    def __init__(
        self,
        custom_patterns: Optional[Dict[str, str]] = None,
        require_heterocycle: bool = False,
    ):
        """
        Initialize heterocycle filter.
        
        Args:
            custom_patterns: Optional custom SMARTS patterns to use instead of defaults
            require_heterocycle: If True, pass only molecules WITH heterocycles.
                               If False, pass only molecules WITHOUT heterocycles.
        """
        patterns = custom_patterns if custom_patterns is not None else HETEROCYCLE_PATTERNS
        super().__init__(patterns)
        self.require_heterocycle = require_heterocycle

    def check_molecule(self, mol: Chem.Mol) -> Tuple[bool, List[str]]:
        """
        Check if molecule passes heterocycle filter.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            Tuple of (passes_filter, list_of_matched_patterns)
        """
        if mol is None:
            return False, ["Invalid molecule"]

        matched_patterns = []
        for name, pattern in self._compiled_patterns.items():
            if mol.HasSubstructMatch(pattern):
                matched_patterns.append(name)

        if self.require_heterocycle:
            passes = len(matched_patterns) > 0
        else:
            passes = len(matched_patterns) == 0

        return passes, matched_patterns


class PropertyFilter:
    """
    Filter molecules based on calculated molecular properties.
    
    Supports filtering by logP, TPSA, HBD, HBA, molecular weight, etc.
    """

    def __init__(
        self,
        logp_range: Optional[Tuple[float, float]] = None,
        tpsa_range: Optional[Tuple[float, float]] = None,
        mw_range: Optional[Tuple[float, float]] = None,
        hbd_range: Optional[Tuple[int, int]] = None,
        hba_range: Optional[Tuple[int, int]] = None,
        rotatable_bonds_max: Optional[int] = None,
        aromatic_rings_range: Optional[Tuple[int, int]] = None,
    ):
        """
        Initialize property filter with ranges.
        
        Args:
            logp_range: (min, max) tuple for logP
            tpsa_range: (min, max) tuple for TPSA in Ų
            mw_range: (min, max) tuple for molecular weight in g/mol
            hbd_range: (min, max) tuple for number of H-bond donors
            hba_range: (min, max) tuple for number of H-bond acceptors
            rotatable_bonds_max: Maximum number of rotatable bonds
            aromatic_rings_range: (min, max) tuple for number of aromatic rings
        """
        self.logp_range = logp_range
        self.tpsa_range = tpsa_range
        self.mw_range = mw_range
        self.hbd_range = hbd_range
        self.hba_range = hba_range
        self.rotatable_bonds_max = rotatable_bonds_max
        self.aromatic_rings_range = aromatic_rings_range

    def check_molecule(self, mol: Chem.Mol) -> Tuple[bool, List[str]]:
        """
        Check if molecule passes property filters.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            Tuple of (passes_filter, list_of_failed_properties)
        """
        if mol is None:
            return False, ["Invalid molecule"]

        failed_properties = []
        props = MolecularProperties.calculate_all_properties(mol)

        if self.logp_range is not None:
            if not (self.logp_range[0] <= props['logP'] <= self.logp_range[1]):
                failed_properties.append(
                    f"logP={props['logP']:.2f} outside range {self.logp_range}"
                )

        if self.tpsa_range is not None:
            if not (self.tpsa_range[0] <= props['tpsa'] <= self.tpsa_range[1]):
                failed_properties.append(
                    f"TPSA={props['tpsa']:.2f} outside range {self.tpsa_range}"
                )

        if self.mw_range is not None:
            if not (self.mw_range[0] <= props['molecular_weight'] <= self.mw_range[1]):
                failed_properties.append(
                    f"MW={props['molecular_weight']:.2f} outside range {self.mw_range}"
                )

        if self.hbd_range is not None:
            if not (self.hbd_range[0] <= props['hbd'] <= self.hbd_range[1]):
                failed_properties.append(f"HBD={props['hbd']} outside range {self.hbd_range}")

        if self.hba_range is not None:
            if not (self.hba_range[0] <= props['hba'] <= self.hba_range[1]):
                failed_properties.append(f"HBA={props['hba']} outside range {self.hba_range}")

        if self.rotatable_bonds_max is not None:
            if props['rotatable_bonds'] > self.rotatable_bonds_max:
                failed_properties.append(
                    f"Rotatable bonds={props['rotatable_bonds']} > max {self.rotatable_bonds_max}"
                )

        if self.aromatic_rings_range is not None:
            if not (
                self.aromatic_rings_range[0]
                <= props['aromatic_rings']
                <= self.aromatic_rings_range[1]
            ):
                failed_properties.append(
                    f"Aromatic rings={props['aromatic_rings']} "
                    f"outside range {self.aromatic_rings_range}"
                )

        passes = len(failed_properties) == 0
        return passes, failed_properties

    def filter_molecules(
        self, molecules: List[Chem.Mol]
    ) -> Tuple[List[Chem.Mol], List[Chem.Mol], Dict[int, List[str]]]:
        """
        Filter a list of molecules based on properties.
        
        Args:
            molecules: List of RDKit molecule objects
            
        Returns:
            Tuple of (passed_molecules, failed_molecules, failure_reasons)
        """
        passed = []
        failed = []
        failure_reasons = {}

        for idx, mol in enumerate(molecules):
            passes, reasons = self.check_molecule(mol)
            if passes:
                passed.append(mol)
            else:
                failed.append(mol)
                failure_reasons[idx] = reasons

        return passed, failed, failure_reasons


class FilterGroup:
    """
    Combine multiple filters and apply them sequentially.
    
    This allows for complex filtering pipelines with multiple criteria.
    """

    def __init__(self, filters: List[Union[BaseFilter, PropertyFilter]]):
        """
        Initialize filter group.
        
        Args:
            filters: List of filter objects to apply in sequence
        """
        self.filters = filters

    def filter_molecules(
        self, molecules: List[Chem.Mol], verbose: bool = False
    ) -> Tuple[List[Chem.Mol], Dict[int, Dict[str, List[str]]]]:
        """
        Apply all filters in sequence to a list of molecules.
        
        Args:
            molecules: List of RDKit molecule objects
            verbose: If True, print progress information
            
        Returns:
            Tuple of (passed_molecules, failure_details)
            failure_details is a dict mapping molecule indices to filter names and reasons
        """
        current_molecules = list(molecules)
        all_failure_details = {}

        for filter_idx, filter_obj in enumerate(self.filters):
            filter_name = filter_obj.__class__.__name__
            if verbose:
                print(f"Applying filter {filter_idx + 1}/{len(self.filters)}: {filter_name}")
                print(f"  Input: {len(current_molecules)} molecules")

            passed, failed, reasons = filter_obj.filter_molecules(current_molecules)

            # Track failures
            for idx, mol in enumerate(current_molecules):
                if mol in failed:
                    original_idx = molecules.index(mol)
                    if original_idx not in all_failure_details:
                        all_failure_details[original_idx] = {}
                    
                    failed_idx = failed.index(mol)
                    # Find the corresponding reason from the reasons dict
                    reason_key = None
                    for k in reasons.keys():
                        if failed_idx <= k:
                            reason_key = k
                            break
                    if reason_key is not None:
                        all_failure_details[original_idx][filter_name] = reasons[reason_key]

            current_molecules = passed

            if verbose:
                print(f"  Passed: {len(passed)} molecules")
                print(f"  Failed: {len(failed)} molecules")

        return current_molecules, all_failure_details

    def check_molecule(self, mol: Chem.Mol) -> Tuple[bool, Dict[str, List[str]]]:
        """
        Check a single molecule against all filters.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            Tuple of (passes_all_filters, dict_of_filter_failures)
        """
        failure_details = {}

        for filter_obj in self.filters:
            filter_name = filter_obj.__class__.__name__
            passes, reasons = filter_obj.check_molecule(mol)

            if not passes:
                failure_details[filter_name] = reasons

        passes_all = len(failure_details) == 0
        return passes_all, failure_details
