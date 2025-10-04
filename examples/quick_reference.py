"""
Quick Reference Card for medchem-molecule-filter

Common usage patterns and examples.
"""

from rdkit import Chem
from medchem_filter import (
    PAINSFilter,
    ReactiveFilter,
    HeterocycleFilter,
    PropertyFilter,
    FilterGroup,
    MolecularProperties,
)

# =============================================================================
# 1. BASIC FILTERING
# =============================================================================

# Check single molecule against PAINS
pains_filter = PAINSFilter()
mol = Chem.MolFromSmiles("c1ccccc1")
passes, matched_patterns = pains_filter.check_molecule(mol)

# Filter list of molecules
molecules = [Chem.MolFromSmiles(s) for s in ["CCO", "c1ccccc1", "CC=O"]]
passed, failed, reasons = pains_filter.filter_molecules(molecules)

# =============================================================================
# 2. PROPERTY-BASED FILTERING (Lipinski's Rule of Five)
# =============================================================================

ro5_filter = PropertyFilter(
    mw_range=(0, 500),      # Molecular weight <= 500 Da
    logp_range=(-float('inf'), 5),  # LogP <= 5
    hbd_range=(0, 5),       # HBD <= 5
    hba_range=(0, 10),      # HBA <= 10
)

mol = Chem.MolFromSmiles("CCN(CC)CC")
passes, reasons = ro5_filter.check_molecule(mol)

# =============================================================================
# 3. LEAD-LIKE PROPERTIES
# =============================================================================

lead_like_filter = PropertyFilter(
    mw_range=(200, 350),
    logp_range=(1, 3),
    hbd_range=(0, 3),
    hba_range=(0, 6),
)

# =============================================================================
# 4. COMBINING FILTERS (Drug-like pipeline)
# =============================================================================

drug_like_pipeline = FilterGroup([
    PAINSFilter(),              # Remove PAINS
    ReactiveFilter(),           # Remove reactive groups
    PropertyFilter(             # Apply Lipinski-like rules
        mw_range=(200, 500),
        logp_range=(-1, 5),
        hbd_range=(0, 5),
        hba_range=(0, 10),
        rotatable_bonds_max=10,
    ),
])

mol = Chem.MolFromSmiles("c1ccccc1CCO")
passes, failure_details = drug_like_pipeline.check_molecule(mol)

# =============================================================================
# 5. HETEROCYCLE FILTERING
# =============================================================================

# Require at least one heterocycle
require_het = HeterocycleFilter(require_heterocycle=True)
mol = Chem.MolFromSmiles("c1ccncc1")  # Pyridine
passes, matched = require_het.check_molecule(mol)

# Exclude heterocycles
exclude_het = HeterocycleFilter(require_heterocycle=False)
passes, matched = exclude_het.check_molecule(mol)

# =============================================================================
# 6. CALCULATING PROPERTIES
# =============================================================================

mol = Chem.MolFromSmiles("c1ccccc1CCO")

# Calculate individual properties
logp = MolecularProperties.calculate_logp(mol)
tpsa = MolecularProperties.calculate_tpsa(mol)
mw = MolecularProperties.calculate_molecular_weight(mol)
hbd = MolecularProperties.calculate_hbd(mol)
hba = MolecularProperties.calculate_hba(mol)

# Calculate all properties at once
props = MolecularProperties.calculate_all_properties(mol)
# Returns: {'logP': 1.22, 'tpsa': 20.23, 'hbd': 1, 'hba': 1, 
#           'molecular_weight': 122.17, 'rotatable_bonds': 2, 
#           'aromatic_rings': 1}

# =============================================================================
# 7. CUSTOM SMARTS PATTERNS
# =============================================================================

# Define custom PAINS patterns
custom_pains = PAINSFilter(custom_patterns={
    "my_pattern_1": "c1ccccc1N=N",
    "my_pattern_2": "C(=O)O[N+](=O)[O-]",
})

# Define custom reactive patterns
custom_reactive = ReactiveFilter(custom_patterns={
    "my_reactive_1": "[Si]",  # Silicon-containing
    "my_reactive_2": "[B]",   # Boron-containing
})

# =============================================================================
# 8. BATCH FILTERING WITH VERBOSE OUTPUT
# =============================================================================

smiles_list = ["CCO", "c1ccccc1", "CC=O", "c1ccncc1"]
molecules = [Chem.MolFromSmiles(s) for s in smiles_list]

filter_group = FilterGroup([PAINSFilter(), ReactiveFilter()])
passed, failed_details = filter_group.filter_molecules(molecules, verbose=True)

print(f"Passed: {len(passed)}/{len(molecules)} molecules")
for idx, details in failed_details.items():
    print(f"Molecule {idx} failed:")
    for filter_name, reasons in details.items():
        print(f"  {filter_name}: {reasons}")

# =============================================================================
# 9. COMMON PROPERTY RANGES
# =============================================================================

# Lipinski's Rule of Five
lipinski = PropertyFilter(mw_range=(0, 500), logp_range=(-float('inf'), 5), 
                          hbd_range=(0, 5), hba_range=(0, 10))

# Veber's Rules
veber = PropertyFilter(rotatable_bonds_max=10, tpsa_range=(0, 140))

# Lead-like
lead_like = PropertyFilter(mw_range=(200, 350), logp_range=(1, 3), 
                           hbd_range=(0, 3), hba_range=(0, 6))

# Fragment-like
fragment_like = PropertyFilter(mw_range=(100, 250), logp_range=(0, 3))

# =============================================================================
# 10. FILTERING WORKFLOW EXAMPLE
# =============================================================================

def filter_compound_library(smiles_list):
    """Complete workflow for filtering a compound library."""
    
    # Parse SMILES
    molecules = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            molecules.append(mol)
    
    # Define filter pipeline
    pipeline = FilterGroup([
        PAINSFilter(),
        ReactiveFilter(),
        PropertyFilter(
            mw_range=(200, 500),
            logp_range=(-1, 5),
            hbd_range=(0, 5),
            hba_range=(0, 10),
            rotatable_bonds_max=10,
            tpsa_range=(20, 140),
        ),
    ])
    
    # Apply filters
    passed, failed_details = pipeline.filter_molecules(molecules, verbose=True)
    
    # Calculate properties for passed molecules
    results = []
    for mol in passed:
        props = MolecularProperties.calculate_all_properties(mol)
        results.append({
            'smiles': Chem.MolToSmiles(mol),
            'properties': props,
        })
    
    return results, failed_details


if __name__ == "__main__":
    print("medchem-molecule-filter Quick Reference")
    print("=" * 60)
    print("\nSee the code above for common usage patterns.")
    print("\nFor more examples, run: python examples/basic_usage.py")
    print("\nFor documentation, visit:")
    print("https://frequencykg.github.io/medchem-molecule-filter")
