"""
Example usage of the medchem-molecule-filter library.

This script demonstrates how to use various filters to screen molecules.
"""

from rdkit import Chem

from medchem_filter import (
    FilterGroup,
    HeterocycleFilter,
    MolecularProperties,
    PAINSFilter,
    PropertyFilter,
    ReactiveFilter,
)


def example_basic_filtering():
    """Example of basic PAINS filtering."""
    print("=" * 60)
    print("Example 1: Basic PAINS Filtering")
    print("=" * 60)

    pains_filter = PAINSFilter()

    # Test molecules
    test_smiles = [
        ("CCO", "Ethanol"),
        ("c1ccccc1", "Benzene"),
        ("O=C1C=CC(=O)C=C1", "Quinone (PAINS)"),
    ]

    for smiles, name in test_smiles:
        mol = Chem.MolFromSmiles(smiles)
        passes, matched = pains_filter.check_molecule(mol)

        print(f"\n{name} ({smiles}):")
        if passes:
            print("  ✓ Passed PAINS filter")
        else:
            print(f"  ✗ Failed PAINS filter - matched: {matched}")


def example_property_filtering():
    """Example of property-based filtering."""
    print("\n" + "=" * 60)
    print("Example 2: Property-Based Filtering")
    print("=" * 60)

    # Create a Lipinski-like filter
    prop_filter = PropertyFilter(
        mw_range=(0, 500),
        logp_range=(-1, 5),
        hbd_range=(0, 5),
        hba_range=(0, 10),
    )

    test_smiles = [
        ("CCO", "Ethanol"),
        ("c1ccc(cc1)c2ccccc2", "Biphenyl"),
        ("CCCCCCCCCCCCCCCCCC", "Octadecane (high MW/logP)"),
    ]

    for smiles, name in test_smiles:
        mol = Chem.MolFromSmiles(smiles)
        passes, reasons = prop_filter.check_molecule(mol)

        print(f"\n{name} ({smiles}):")
        if passes:
            print("  ✓ Passed property filter")
            props = MolecularProperties.calculate_all_properties(mol)
            print(f"    MW: {props['molecular_weight']:.1f}, logP: {props['logP']:.2f}")
        else:
            print(f"  ✗ Failed property filter:")
            for reason in reasons:
                print(f"    - {reason}")


def example_combined_filtering():
    """Example of combining multiple filters."""
    print("\n" + "=" * 60)
    print("Example 3: Combined Filter Pipeline")
    print("=" * 60)

    # Create a comprehensive filter
    filter_group = FilterGroup(
        [
            PAINSFilter(),
            ReactiveFilter(),
            PropertyFilter(mw_range=(0, 500), logp_range=(-1, 5)),
        ]
    )

    test_smiles = [
        ("CCO", "Ethanol"),
        ("c1ccncc1", "Pyridine"),
        ("CC=O", "Acetaldehyde (reactive)"),
        ("O=C1C=CC(=O)C=C1", "Quinone (PAINS)"),
    ]

    molecules = [Chem.MolFromSmiles(smiles) for smiles, _ in test_smiles]

    print("\nFiltering molecules:")
    passed, failed_details = filter_group.filter_molecules(molecules, verbose=True)

    print(f"\n\nResults:")
    for idx, (smiles, name) in enumerate(test_smiles):
        if idx in failed_details:
            print(f"\n{name} ({smiles}): FAILED")
            for filter_name, reasons in failed_details[idx].items():
                print(f"  - {filter_name}: {reasons}")
        else:
            print(f"\n{name} ({smiles}): PASSED ✓")


def example_heterocycle_filtering():
    """Example of heterocycle filtering."""
    print("\n" + "=" * 60)
    print("Example 4: Heterocycle Filtering")
    print("=" * 60)

    # Require heterocycles
    filter_require = HeterocycleFilter(require_heterocycle=True)

    # Exclude heterocycles
    filter_exclude = HeterocycleFilter(require_heterocycle=False)

    test_smiles = [
        ("c1ccccc1", "Benzene"),
        ("c1ccncc1", "Pyridine"),
        ("c1ccc2[nH]ccc2c1", "Indole"),
    ]

    print("\nRequiring heterocycles:")
    for smiles, name in test_smiles:
        mol = Chem.MolFromSmiles(smiles)
        passes, matched = filter_require.check_molecule(mol)
        status = "✓" if passes else "✗"
        print(f"  {status} {name} - matched: {matched if matched else 'none'}")

    print("\nExcluding heterocycles:")
    for smiles, name in test_smiles:
        mol = Chem.MolFromSmiles(smiles)
        passes, matched = filter_exclude.check_molecule(mol)
        status = "✓" if passes else "✗"
        print(f"  {status} {name} - matched: {matched if matched else 'none'}")


def example_property_calculation():
    """Example of calculating molecular properties."""
    print("\n" + "=" * 60)
    print("Example 5: Property Calculation")
    print("=" * 60)

    test_smiles = [
        ("CCO", "Ethanol"),
        ("c1ccccc1CCO", "Phenethyl alcohol"),
        ("c1ccc(cc1)c2ccccc2", "Biphenyl"),
    ]

    for smiles, name in test_smiles:
        mol = Chem.MolFromSmiles(smiles)
        props = MolecularProperties.calculate_all_properties(mol)

        print(f"\n{name} ({smiles}):")
        print(f"  Molecular Weight: {props['molecular_weight']:.2f} g/mol")
        print(f"  LogP: {props['logP']:.2f}")
        print(f"  TPSA: {props['tpsa']:.2f} Ų")
        print(f"  HBD: {props['hbd']}")
        print(f"  HBA: {props['hba']}")
        print(f"  Rotatable Bonds: {props['rotatable_bonds']}")
        print(f"  Aromatic Rings: {props['aromatic_rings']}")


if __name__ == "__main__":
    print("\nmedchem-molecule-filter - Example Usage\n")

    example_basic_filtering()
    example_property_filtering()
    example_combined_filtering()
    example_heterocycle_filtering()
    example_property_calculation()

    print("\n" + "=" * 60)
    print("Examples completed!")
    print("=" * 60)
