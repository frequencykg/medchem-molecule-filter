"""
Pandas Integration Examples for medchem-molecule-filter

This script demonstrates how to use the library with pandas DataFrames
containing SMILES strings for batch processing of molecular data.
"""

import pandas as pd
from rdkit import Chem

# Import the main library and pandas utilities
from medchem_filter import (
    FilterGroup,
    HeterocycleFilter,
    PAINSFilter,
    PropertyFilter,
    ReactiveFilter,
    add_mol_column,
    apply_filters_to_dataframe,
    calculate_properties_for_dataframe,
    filter_by_properties,
    filter_dataframe,
)


def example_1_basic_filtering():
    """Example 1: Basic filtering of a DataFrame with SMILES."""
    print("=" * 70)
    print("Example 1: Basic Filtering with DataFrame")
    print("=" * 70)

    # Create a sample DataFrame
    df = pd.DataFrame(
        {
            "compound_id": ["COMP001", "COMP002", "COMP003", "COMP004", "COMP005"],
            "smiles": [
                "CCO",  # Ethanol - clean
                "c1ccccc1",  # Benzene - clean
                "O=C1C=CC(=O)C=C1",  # Quinone - PAINS
                "CC=O",  # Acetaldehyde - reactive
                "c1ccncc1",  # Pyridine - heterocycle
            ],
            "name": ["Ethanol", "Benzene", "Quinone", "Acetaldehyde", "Pyridine"],
        }
    )

    print("\nOriginal DataFrame:")
    print(df)

    # Apply PAINS filter
    pains_filter = PAINSFilter()
    filtered_df = filter_dataframe(df, pains_filter, keep_failures=True)

    print("\n\nAfter PAINS Filter:")
    print(filtered_df[["compound_id", "name", "smiles", "passes_filter", "failure_reasons"]])

    # Keep only passing molecules
    clean_df = filter_dataframe(df, pains_filter, keep_failures=False)
    print(f"\n\nClean molecules: {len(clean_df)}/{len(df)}")
    print(clean_df[["compound_id", "name", "smiles"]])


def example_2_multiple_filters():
    """Example 2: Apply multiple filters with detailed results."""
    print("\n\n" + "=" * 70)
    print("Example 2: Multiple Filters with Detailed Results")
    print("=" * 70)

    # Create DataFrame
    df = pd.DataFrame(
        {
            "id": range(1, 6),
            "smiles": [
                "CCO",  # Clean
                "O=C1C=CC(=O)C=C1",  # PAINS
                "CC=O",  # Reactive
                "c1ccccc1CCO",  # Clean
                "invalid_smiles",  # Invalid
            ],
        }
    )

    print("\nOriginal DataFrame:")
    print(df)

    # Apply multiple filters
    filters = [
        PAINSFilter(),
        ReactiveFilter(),
    ]

    result = apply_filters_to_dataframe(df, filters, detailed=True)

    print("\n\nFiltered with Details:")
    print(result)

    # Show only failures
    failures = result[~result["passes_all"]]
    print("\n\nFailed Molecules:")
    print(failures[["id", "smiles", "failed_PAINSFilter", "failed_ReactiveFilter"]])


def example_3_property_calculation():
    """Example 3: Calculate molecular properties for a DataFrame."""
    print("\n\n" + "=" * 70)
    print("Example 3: Calculate Molecular Properties")
    print("=" * 70)

    # Create DataFrame
    df = pd.DataFrame(
        {
            "compound": ["Drug A", "Drug B", "Drug C"],
            "smiles": [
                "CCO",  # Ethanol
                "c1ccccc1CCO",  # Phenethyl alcohol
                "c1ccc(cc1)c2ccccc2",  # Biphenyl
            ],
        }
    )

    print("\nOriginal DataFrame:")
    print(df)

    # Calculate all properties
    df_with_props = calculate_properties_for_dataframe(df)

    print("\n\nWith Molecular Properties:")
    print(
        df_with_props[
            [
                "compound",
                "smiles",
                "molecular_weight",
                "logP",
                "tpsa",
                "hbd",
                "hba",
                "rotatable_bonds",
                "aromatic_rings",
            ]
        ].to_string(index=False)
    )

    # Calculate only specific properties
    df_selected = calculate_properties_for_dataframe(
        df, properties=["molecular_weight", "logP", "hbd", "hba"]
    )
    print("\n\nSelected Properties Only:")
    print(
        df_selected[["compound", "molecular_weight", "logP", "hbd", "hba"]].to_string(index=False)
    )


def example_4_property_filtering():
    """Example 4: Filter by molecular properties (Lipinski's Rule of Five)."""
    print("\n\n" + "=" * 70)
    print("Example 4: Property-Based Filtering (Lipinski's Rule of Five)")
    print("=" * 70)

    # Create DataFrame with various molecules
    df = pd.DataFrame(
        {
            "compound": ["Ethanol", "Aspirin", "Caffeine", "Large Molecule", "Lipophilic"],
            "smiles": [
                "CCO",  # Small, passes
                "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin, passes
                "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine, passes
                "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" * 10,  # Too large
                "CCCCCCCCCCCCCCCCCC",  # Too lipophilic
            ],
        }
    )

    print("\nOriginal DataFrame:")
    print(df)

    # Apply Lipinski's Rule of Five
    lipinski_filtered = filter_by_properties(
        df,
        mw_range=(0, 500),
        logp_range=(-float("inf"), 5),
        hbd_range=(0, 5),
        hba_range=(0, 10),
        keep_failures=True,
    )

    print("\n\nLipinski Rule of Five Results:")
    print(lipinski_filtered[["compound", "smiles", "passes_filter"]])

    # Show only passing compounds
    passing = lipinski_filtered[lipinski_filtered["passes_filter"]]
    print(f"\n\nPassing compounds: {len(passing)}/{len(df)}")
    print(passing[["compound", "smiles"]].to_string(index=False))


def example_5_comprehensive_pipeline():
    """Example 5: Comprehensive drug-like filtering pipeline."""
    print("\n\n" + "=" * 70)
    print("Example 5: Comprehensive Drug-Like Filtering Pipeline")
    print("=" * 70)

    # Create DataFrame with diverse compounds
    df = pd.DataFrame(
        {
            "id": ["CMPD_" + str(i).zfill(3) for i in range(1, 11)],
            "smiles": [
                "CCO",  # Ethanol - passes
                "c1ccccc1CCO",  # Phenethyl alcohol - passes
                "O=C1C=CC(=O)C=C1",  # Quinone - PAINS
                "CC=O",  # Acetaldehyde - reactive
                "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin - passes
                "CCCCCCCCCCCCCCCC",  # Long chain - fails properties
                "c1ccncc1",  # Pyridine - passes
                "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine - passes
                "invalid_smiles",  # Invalid
                "c1ccc2ccccc2c1",  # Naphthalene - passes
            ],
        }
    )

    print("\nOriginal DataFrame (10 compounds):")
    print(df)

    # Step 1: Calculate properties first
    print("\n\nStep 1: Calculating molecular properties...")
    df_props = calculate_properties_for_dataframe(df)

    print("\nMolecular Properties:")
    print(df_props[["id", "molecular_weight", "logP", "hbd", "hba"]].head(5).to_string(index=False))

    # Step 2: Apply structural filters
    print("\n\nStep 2: Applying structural filters (PAINS, Reactive)...")
    filters = [PAINSFilter(), ReactiveFilter()]
    df_filtered = apply_filters_to_dataframe(df_props, filters, detailed=True)

    structural_pass = df_filtered[df_filtered["passes_all"]]
    print(f"Passed structural filters: {len(structural_pass)}/{len(df)}")

    # Step 3: Apply property filters
    print("\n\nStep 3: Applying property filters (Lipinski's Rule)...")
    df_final = filter_by_properties(
        structural_pass,
        mw_range=(0, 500),
        logp_range=(-1, 5),
        hbd_range=(0, 5),
        hba_range=(0, 10),
        keep_failures=False,
    )

    print(f"\nFinal drug-like compounds: {len(df_final)}/{len(df)}")
    print("\nFinal Results:")
    print(
        df_final[["id", "smiles", "molecular_weight", "logP", "hbd", "hba"]].to_string(index=False)
    )


def example_6_csv_workflow():
    """Example 6: Complete workflow with CSV file I/O."""
    print("\n\n" + "=" * 70)
    print("Example 6: CSV File Workflow")
    print("=" * 70)

    # Create sample data
    df = pd.DataFrame(
        {
            "compound_id": ["C001", "C002", "C003", "C004", "C005"],
            "smiles": [
                "CCO",
                "c1ccccc1",
                "CC=O",
                "c1ccncc1",
                "CC(=O)Oc1ccccc1C(=O)O",
            ],
            "source": ["ChEMBL", "PubChem", "ChEMBL", "ChEMBL", "DrugBank"],
        }
    )

    print("\nOriginal data:")
    print(df)

    # Save to CSV
    csv_file = "/tmp/compounds.csv"
    df.to_csv(csv_file, index=False)
    print(f"\nSaved to: {csv_file}")

    # Read from CSV
    df_loaded = pd.read_csv(csv_file)
    print("\nLoaded from CSV:")
    print(df_loaded)

    # Apply filters
    pipeline = FilterGroup([PAINSFilter(), ReactiveFilter()])
    df_filtered = filter_dataframe(df_loaded, pipeline, keep_failures=False)

    # Calculate properties for passing molecules
    df_filtered = calculate_properties_for_dataframe(df_filtered)

    # Save filtered results
    output_file = "/tmp/compounds_filtered.csv"
    df_filtered.to_csv(output_file, index=False)
    print(f"\n\nFiltered results saved to: {output_file}")
    print(f"Passing compounds: {len(df_filtered)}/{len(df)}")
    print("\nFiltered data:")
    print(df_filtered[["compound_id", "smiles", "molecular_weight", "logP"]].to_string(index=False))


def example_7_merge_with_activity_data():
    """Example 7: Merge filtered compounds with activity data."""
    print("\n\n" + "=" * 70)
    print("Example 7: Merge with Activity Data")
    print("=" * 70)

    # Compound library
    compounds = pd.DataFrame(
        {
            "compound_id": ["A", "B", "C", "D", "E"],
            "smiles": [
                "CCO",
                "c1ccccc1",
                "O=C1C=CC(=O)C=C1",  # PAINS
                "CC=O",  # Reactive
                "c1ccncc1",
            ],
        }
    )

    # Activity data
    activity = pd.DataFrame(
        {
            "compound_id": ["A", "B", "C", "D", "E"],
            "ic50_nm": [1000, 500, 10, 5, 2000],
            "assay": ["Kinase", "Kinase", "Kinase", "Kinase", "Kinase"],
        }
    )

    print("\nCompounds:")
    print(compounds)
    print("\nActivity Data:")
    print(activity)

    # Filter compounds
    filters = [PAINSFilter(), ReactiveFilter()]
    compounds_filtered = apply_filters_to_dataframe(compounds, filters, detailed=True)

    # Merge with activity data
    merged = compounds_filtered.merge(activity, on="compound_id")

    print("\n\nMerged Data:")
    print(merged[["compound_id", "smiles", "passes_all", "ic50_nm"]])

    # Flag problematic actives
    merged["reliable_hit"] = (merged["passes_all"]) & (merged["ic50_nm"] < 1000)

    print("\n\nReliable Hits (IC50 < 1000 nM, passes filters):")
    reliable = merged[merged["reliable_hit"]]
    print(reliable[["compound_id", "smiles", "ic50_nm"]])


if __name__ == "__main__":
    print("\n")
    print("=" * 70)
    print("MEDCHEM-MOLECULE-FILTER: PANDAS INTEGRATION EXAMPLES")
    print("=" * 70)

    # Run all examples
    example_1_basic_filtering()
    example_2_multiple_filters()
    example_3_property_calculation()
    example_4_property_filtering()
    example_5_comprehensive_pipeline()
    example_6_csv_workflow()
    example_7_merge_with_activity_data()

    print("\n" + "=" * 70)
    print("All examples completed!")
    print("=" * 70)
