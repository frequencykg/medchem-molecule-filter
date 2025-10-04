#!/usr/bin/env python3
"""
Example demonstrating the use of all four filter types with reference to 
the comprehensive filter documentation.

This example shows how to:
1. Use each filter type individually
2. Understand what patterns are being matched
3. Refer to the documentation for visual representations

For detailed information about each pattern including SMARTS strings and 
visual representations, see the documentation at:
docs/_build/html/filter_rules.html
"""

from rdkit import Chem
from medchem_filter import PAINSFilter, ReactiveFilter, HeterocycleFilter, PropertyFilter, FilterGroup

# Sample molecules for testing
test_molecules = {
    "ethanol": "CCO",
    "benzene": "c1ccccc1",
    "quinone": "O=C1C=CC(=O)C=C1",  # PAINS pattern
    "catechol": "c1cc(O)c(O)cc1",  # PAINS pattern
    "acyl_chloride": "CC(=O)Cl",  # Reactive pattern
    "epoxide": "C1OC1",  # Reactive pattern
    "pyridine": "c1ccncc1",  # Heterocycle pattern
    "indole": "c1ccc2[nH]ccc2c1",  # Heterocycle pattern
}

print("=" * 80)
print("Comprehensive Filter Example")
print("=" * 80)
print("\nFor detailed pattern information with visual representations, see:")
print("  docs/_build/html/filter_rules.html")
print("  or build the docs with: cd docs && make html")
print("\n")

# 1. PAINS Filter (19 patterns)
print("1. PAINS Filter")
print("-" * 80)
print("Identifies 19 Pan-Assay Interference patterns including:")
print("  - Quinones (2 patterns)")
print("  - Catechols (2 patterns)")
print("  - Michael acceptors (2 patterns)")
print("  - Rhodanines (3 patterns)")
print("  - And 10 other problematic scaffolds")
print()

pains_filter = PAINSFilter()
for name, smiles in test_molecules.items():
    mol = Chem.MolFromSmiles(smiles)
    passes, matched = pains_filter.check_molecule(mol)
    status = "PASS" if passes else f"FAIL (matched: {', '.join(matched)})"
    print(f"  {name:20s}: {status}")

# 2. Reactive Filter (22 patterns)
print("\n2. Reactive Filter")
print("-" * 80)
print("Identifies 22 reactive functional groups including:")
print("  - Electrophiles (7 patterns)")
print("  - Michael acceptors (3 patterns)")
print("  - Peroxides (2 patterns)")
print("  - Diazo compounds (2 patterns)")
print("  - And 8 other reactive groups")
print()

reactive_filter = ReactiveFilter()
for name, smiles in test_molecules.items():
    mol = Chem.MolFromSmiles(smiles)
    passes, matched = reactive_filter.check_molecule(mol)
    status = "PASS" if passes else f"FAIL (matched: {', '.join(matched)})"
    print(f"  {name:20s}: {status}")

# 3. Heterocycle Filter (31 patterns)
print("\n3. Heterocycle Filter (require heterocycle)")
print("-" * 80)
print("Identifies 31 heterocyclic scaffolds including:")
print("  - 5-membered rings with 1 heteroatom (3 patterns)")
print("  - 5-membered rings with 2 heteroatoms (5 patterns)")
print("  - 6-membered rings (1-3 heteroatoms, 9 patterns)")
print("  - Fused bicyclic systems (9 patterns)")
print("  - Saturated heterocycles (5 patterns)")
print()

het_filter = HeterocycleFilter(require_heterocycle=True)
for name, smiles in test_molecules.items():
    mol = Chem.MolFromSmiles(smiles)
    passes, matched = het_filter.check_molecule(mol)
    status = "PASS" if passes else f"FAIL (no heterocycle)"
    if matched:
        status = f"PASS (contains: {', '.join(matched)})"
    print(f"  {name:20s}: {status}")

# 4. Property Filter
print("\n4. Property Filter (Lipinski's Rule of Five)")
print("-" * 80)
print("Filters based on calculated molecular properties:")
print("  - Molecular Weight ≤ 500 Da")
print("  - LogP ≤ 5")
print("  - HBD (Hydrogen Bond Donors) ≤ 5")
print("  - HBA (Hydrogen Bond Acceptors) ≤ 10")
print()

prop_filter = PropertyFilter(
    mw_range=(0, 500),
    logp_range=(-5, 5),
    hbd_range=(0, 5),
    hba_range=(0, 10)
)
for name, smiles in test_molecules.items():
    mol = Chem.MolFromSmiles(smiles)
    passes, failed_props = prop_filter.check_molecule(mol)
    status = "PASS" if passes else f"FAIL ({', '.join(failed_props)})"
    print(f"  {name:20s}: {status}")

# 5. Combined Filter Pipeline
print("\n5. Combined Filter Pipeline")
print("-" * 80)
print("Applying all filters in sequence...")
print()

filter_group = FilterGroup([
    PAINSFilter(),
    ReactiveFilter(),
    HeterocycleFilter(require_heterocycle=True),
    PropertyFilter(mw_range=(0, 500), logp_range=(-5, 5))
])

molecules = [Chem.MolFromSmiles(s) for s in test_molecules.values()]
passed, failed_details = filter_group.filter_molecules(molecules, verbose=False)

print(f"Input molecules: {len(molecules)}")
print(f"Passed all filters: {len(passed)}")
print(f"Failed: {len(failed_details)}")

if failed_details:
    print("\nFailure details:")
    mol_names = list(test_molecules.keys())
    for idx, details in failed_details.items():
        print(f"  {mol_names[idx]:20s}:")
        for filter_name, reasons in details.items():
            print(f"    - {filter_name}: {reasons}")

print("\n" + "=" * 80)
print("Summary")
print("=" * 80)
print(f"Total SMARTS patterns available: 72")
print(f"  - PAINS patterns: 19")
print(f"  - Reactive patterns: 22")
print(f"  - Heterocycle patterns: 31")
print()
print("All patterns are documented with:")
print("  - Pattern name")
print("  - SMARTS string")
print("  - Description")
print("  - Visual representation (RDKit-generated)")
print()
print("Documentation location: docs/_build/html/filter_rules.html")
print("=" * 80)
