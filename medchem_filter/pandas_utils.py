"""
Pandas helper utilities for working with DataFrames containing SMILES strings.

This module provides convenient functions for applying molecular filters to
pandas DataFrames containing SMILES strings.
"""

from typing import Any, Dict, List, Optional, Union

import pandas as pd
from rdkit import Chem

from .filters import BaseFilter, FilterGroup, PropertyFilter


def smiles_to_mol(smiles: str) -> Optional[Chem.Mol]:
    """
    Convert a SMILES string to an RDKit molecule object.

    Args:
        smiles: SMILES string representation of the molecule

    Returns:
        RDKit Mol object or None if SMILES is invalid
    """
    if pd.isna(smiles) or smiles == "":
        return None
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def add_mol_column(
    df: pd.DataFrame, smiles_column: str = "smiles", mol_column: str = "mol"
) -> pd.DataFrame:
    """
    Add a column with RDKit molecule objects from SMILES strings.

    Args:
        df: DataFrame containing SMILES strings
        smiles_column: Name of the column containing SMILES strings
        mol_column: Name for the new column containing molecule objects

    Returns:
        DataFrame with added molecule column

    Example:
        >>> df = pd.DataFrame({"smiles": ["CCO", "c1ccccc1"]})
        >>> df = add_mol_column(df)
        >>> df["mol"].iloc[0]
        <rdkit.Chem.rdchem.Mol object at ...>
    """
    df = df.copy()
    df[mol_column] = df[smiles_column].apply(smiles_to_mol)
    return df


def filter_dataframe(
    df: pd.DataFrame,
    filter_obj: Union[BaseFilter, PropertyFilter, FilterGroup],
    smiles_column: str = "smiles",
    keep_failures: bool = False,
    add_failure_column: bool = True,
) -> pd.DataFrame:
    """
    Filter a DataFrame based on molecular filters.

    Args:
        df: DataFrame containing SMILES strings
        filter_obj: Filter or FilterGroup object to apply
        smiles_column: Name of the column containing SMILES strings
        keep_failures: If True, return all rows with pass/fail info. If False, return only passing rows
        add_failure_column: If True, add columns with failure reasons

    Returns:
        Filtered DataFrame with optional failure information columns

    Example:
        >>> from medchem_filter import PAINSFilter
        >>> df = pd.DataFrame({"smiles": ["CCO", "O=C1C=CC(=O)C=C1"]})
        >>> filtered = filter_dataframe(df, PAINSFilter())
        >>> len(filtered)
        1
    """
    df = df.copy()

    # Add molecule objects
    df["_mol_temp"] = df[smiles_column].apply(smiles_to_mol)

    # Apply filter
    results = []
    for idx, row in df.iterrows():
        mol = row["_mol_temp"]
        if mol is None:
            passes = False
            details = ["Invalid SMILES"]
        else:
            passes, details = filter_obj.check_molecule(mol)

        results.append(
            {
                "index": idx,
                "passes": passes,
                "failure_reasons": str(details) if not passes else "",
            }
        )

    results_df = pd.DataFrame(results).set_index("index")

    # Add results to original dataframe
    df["passes_filter"] = results_df["passes"]
    if add_failure_column:
        df["failure_reasons"] = results_df["failure_reasons"]

    # Remove temporary molecule column
    df = df.drop(columns=["_mol_temp"])

    # Filter rows if requested
    if not keep_failures:
        df = df[df["passes_filter"]].drop(columns=["passes_filter", "failure_reasons"])

    return df


def apply_filters_to_dataframe(
    df: pd.DataFrame,
    filters: List[Union[BaseFilter, PropertyFilter]],
    smiles_column: str = "smiles",
    detailed: bool = False,
) -> pd.DataFrame:
    """
    Apply multiple filters sequentially to a DataFrame.

    Args:
        df: DataFrame containing SMILES strings
        filters: List of filter objects to apply sequentially
        smiles_column: Name of the column containing SMILES strings
        detailed: If True, add columns for each filter's results

    Returns:
        DataFrame with filter results

    Example:
        >>> from medchem_filter import PAINSFilter, ReactiveFilter
        >>> df = pd.DataFrame({"smiles": ["CCO", "CC=O", "c1ccccc1"]})
        >>> filters = [PAINSFilter(), ReactiveFilter()]
        >>> result = apply_filters_to_dataframe(df, filters, detailed=True)
    """
    df = df.copy()
    filter_group = FilterGroup(filters)

    # Add molecule objects
    df["_mol_temp"] = df[smiles_column].apply(smiles_to_mol)

    # Apply filters
    results = []
    for idx, row in df.iterrows():
        mol = row["_mol_temp"]
        if mol is None:
            passes = False
            details = {"Error": ["Invalid SMILES"]}
        else:
            passes, details = filter_group.check_molecule(mol)

        result = {
            "index": idx,
            "passes_all": passes,
        }

        if detailed:
            for filter_name, reasons in details.items():
                result[f"failed_{filter_name}"] = str(reasons) if reasons else ""

        results.append(result)

    results_df = pd.DataFrame(results).set_index("index")

    # Merge results with original dataframe
    df = df.join(results_df)
    df = df.drop(columns=["_mol_temp"])

    return df


def calculate_properties_for_dataframe(
    df: pd.DataFrame,
    smiles_column: str = "smiles",
    properties: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Calculate molecular properties for all molecules in a DataFrame.

    Args:
        df: DataFrame containing SMILES strings
        smiles_column: Name of the column containing SMILES strings
        properties: List of property names to calculate. If None, calculates all.
                   Available: 'logP', 'tpsa', 'hbd', 'hba', 'molecular_weight',
                   'rotatable_bonds', 'aromatic_rings'

    Returns:
        DataFrame with added property columns

    Example:
        >>> df = pd.DataFrame({"smiles": ["CCO", "c1ccccc1"]})
        >>> df = calculate_properties_for_dataframe(df)
        >>> df[["smiles", "molecular_weight", "logP"]]
    """
    from .properties import MolecularProperties

    df = df.copy()
    df["_mol_temp"] = df[smiles_column].apply(smiles_to_mol)

    # Calculate properties for each molecule
    prop_results = []
    for idx, row in df.iterrows():
        mol = row["_mol_temp"]
        if mol is None:
            props = {
                "logP": None,
                "tpsa": None,
                "hbd": None,
                "hba": None,
                "molecular_weight": None,
                "rotatable_bonds": None,
                "aromatic_rings": None,
            }
        else:
            props = MolecularProperties.calculate_all_properties(mol)

        prop_results.append(props)

    prop_df = pd.DataFrame(prop_results, index=df.index)

    # Filter properties if specified
    if properties is not None:
        prop_df = prop_df[properties]

    # Add to original dataframe
    df = pd.concat([df, prop_df], axis=1)
    df = df.drop(columns=["_mol_temp"])

    return df


def filter_by_properties(
    df: pd.DataFrame,
    smiles_column: str = "smiles",
    logp_range: Optional[tuple] = None,
    tpsa_range: Optional[tuple] = None,
    mw_range: Optional[tuple] = None,
    hbd_range: Optional[tuple] = None,
    hba_range: Optional[tuple] = None,
    keep_failures: bool = False,
) -> pd.DataFrame:
    """
    Filter DataFrame by molecular property ranges.

    Args:
        df: DataFrame containing SMILES strings
        smiles_column: Name of the column containing SMILES strings
        logp_range: (min, max) tuple for logP
        tpsa_range: (min, max) tuple for TPSA
        mw_range: (min, max) tuple for molecular weight
        hbd_range: (min, max) tuple for H-bond donors
        hba_range: (min, max) tuple for H-bond acceptors
        keep_failures: If True, keep all rows with pass/fail info

    Returns:
        Filtered DataFrame

    Example:
        >>> df = pd.DataFrame({"smiles": ["CCO", "c1ccccc1"]})
        >>> filtered = filter_by_properties(df, mw_range=(0, 100))
    """
    prop_filter = PropertyFilter(
        logp_range=logp_range,
        tpsa_range=tpsa_range,
        mw_range=mw_range,
        hbd_range=hbd_range,
        hba_range=hba_range,
    )

    return filter_dataframe(df, prop_filter, smiles_column, keep_failures=keep_failures)
