"""Tests for pandas utility functions."""

import pytest
import pandas as pd
from rdkit import Chem

# Check if pandas is available
pandas_available = True
try:
    from medchem_filter.pandas_utils import (
        add_mol_column,
        filter_dataframe,
        apply_filters_to_dataframe,
        calculate_properties_for_dataframe,
        filter_by_properties,
    )
    from medchem_filter import PAINSFilter, ReactiveFilter, PropertyFilter, FilterGroup
except ImportError:
    pandas_available = False


@pytest.mark.skipif(not pandas_available, reason="pandas not installed")
class TestPandasUtils:
    """Test pandas utility functions."""

    def test_add_mol_column(self):
        """Test adding molecule column to DataFrame."""
        df = pd.DataFrame({"smiles": ["CCO", "c1ccccc1", "invalid"]})
        df_with_mol = add_mol_column(df)

        assert "mol" in df_with_mol.columns
        assert df_with_mol["mol"].iloc[0] is not None
        assert df_with_mol["mol"].iloc[1] is not None
        assert df_with_mol["mol"].iloc[2] is None  # Invalid SMILES

    def test_filter_dataframe_basic(self):
        """Test basic DataFrame filtering."""
        df = pd.DataFrame({
            "id": [1, 2, 3],
            "smiles": ["CCO", "O=C1C=CC(=O)C=C1", "c1ccccc1"]
        })

        pains_filter = PAINSFilter()
        filtered = filter_dataframe(df, pains_filter, keep_failures=False)

        # Should filter out quinone (PAINS)
        assert len(filtered) == 2
        assert "CCO" in filtered["smiles"].values
        assert "c1ccccc1" in filtered["smiles"].values

    def test_filter_dataframe_with_failures(self):
        """Test DataFrame filtering keeping failure information."""
        df = pd.DataFrame({
            "smiles": ["CCO", "O=C1C=CC(=O)C=C1"]
        })

        pains_filter = PAINSFilter()
        filtered = filter_dataframe(df, pains_filter, keep_failures=True)

        assert len(filtered) == 2
        assert "passes_filter" in filtered.columns
        assert "failure_reasons" in filtered.columns
        assert filtered["passes_filter"].iloc[0] == True
        assert filtered["passes_filter"].iloc[1] == False

    def test_apply_filters_to_dataframe(self):
        """Test applying multiple filters to DataFrame."""
        df = pd.DataFrame({
            "smiles": ["CCO", "CC=O", "c1ccccc1"]
        })

        filters = [PAINSFilter(), ReactiveFilter()]
        result = apply_filters_to_dataframe(df, filters, detailed=True)

        assert "passes_all" in result.columns
        assert result["passes_all"].iloc[0] == True  # Ethanol passes
        assert result["passes_all"].iloc[1] == False  # Acetaldehyde fails (reactive)
        assert result["passes_all"].iloc[2] == True  # Benzene passes

    def test_calculate_properties_for_dataframe(self):
        """Test calculating properties for DataFrame."""
        df = pd.DataFrame({
            "smiles": ["CCO", "c1ccccc1"]
        })

        df_with_props = calculate_properties_for_dataframe(df)

        assert "molecular_weight" in df_with_props.columns
        assert "logP" in df_with_props.columns
        assert "hbd" in df_with_props.columns
        assert "hba" in df_with_props.columns

        # Check ethanol properties
        assert df_with_props["hbd"].iloc[0] == 1
        assert df_with_props["hba"].iloc[0] == 1

    def test_calculate_properties_selected(self):
        """Test calculating only selected properties."""
        df = pd.DataFrame({
            "smiles": ["CCO"]
        })

        df_with_props = calculate_properties_for_dataframe(
            df, properties=["molecular_weight", "logP"]
        )

        assert "molecular_weight" in df_with_props.columns
        assert "logP" in df_with_props.columns
        assert "hbd" not in df_with_props.columns

    def test_filter_by_properties(self):
        """Test property-based filtering."""
        df = pd.DataFrame({
            "smiles": ["CCO", "CCCCCCCCCCCCCCCC"]  # Small and large molecules
        })

        filtered = filter_by_properties(df, mw_range=(0, 100), keep_failures=False)

        # Only ethanol should pass (MW < 100)
        assert len(filtered) == 1
        assert filtered["smiles"].iloc[0] == "CCO"

    def test_invalid_smiles_handling(self):
        """Test handling of invalid SMILES."""
        df = pd.DataFrame({
            "smiles": ["CCO", "invalid_smiles", "c1ccccc1"]
        })

        pains_filter = PAINSFilter()
        filtered = filter_dataframe(df, pains_filter, keep_failures=True)

        # Invalid SMILES should fail
        invalid_row = filtered[filtered["smiles"] == "invalid_smiles"]
        assert invalid_row["passes_filter"].iloc[0] == False
        assert "Invalid" in invalid_row["failure_reasons"].iloc[0]

    def test_empty_dataframe(self):
        """Test handling of empty DataFrame."""
        df = pd.DataFrame({"smiles": []})

        pains_filter = PAINSFilter()
        filtered = filter_dataframe(df, pains_filter)

        assert len(filtered) == 0

    def test_custom_column_names(self):
        """Test using custom column names."""
        df = pd.DataFrame({
            "my_smiles": ["CCO", "c1ccccc1"]
        })

        pains_filter = PAINSFilter()
        filtered = filter_dataframe(df, pains_filter, smiles_column="my_smiles")

        assert len(filtered) == 2

    def test_filter_group_integration(self):
        """Test integration with FilterGroup."""
        df = pd.DataFrame({
            "smiles": ["CCO", "O=C1C=CC(=O)C=C1", "CC=O"]
        })

        filter_group = FilterGroup([PAINSFilter(), ReactiveFilter()])
        filtered = filter_dataframe(df, filter_group, keep_failures=False)

        # Only ethanol should pass
        assert len(filtered) == 1
        assert filtered["smiles"].iloc[0] == "CCO"

    def test_property_filter_integration(self):
        """Test integration with PropertyFilter."""
        df = pd.DataFrame({
            "smiles": ["CCO", "CCCCCCCCCCCCCCCC"]
        })

        prop_filter = PropertyFilter(mw_range=(0, 100))
        filtered = filter_dataframe(df, prop_filter, keep_failures=True)

        assert filtered["passes_filter"].iloc[0] == True  # Ethanol passes
        assert filtered["passes_filter"].iloc[1] == False  # Long chain fails

    def test_none_and_empty_smiles(self):
        """Test handling of None and empty SMILES."""
        df = pd.DataFrame({
            "smiles": ["CCO", None, "", "c1ccccc1"]
        })

        pains_filter = PAINSFilter()
        filtered = filter_dataframe(df, pains_filter, keep_failures=False)

        # Only valid SMILES should pass
        assert len(filtered) == 2

    def test_dataframe_preserves_other_columns(self):
        """Test that other columns are preserved during filtering."""
        df = pd.DataFrame({
            "id": ["A", "B", "C"],
            "smiles": ["CCO", "c1ccccc1", "CC=O"],
            "name": ["Ethanol", "Benzene", "Acetaldehyde"],
            "activity": [10.5, 20.3, 5.1]
        })

        pains_filter = PAINSFilter()
        filtered = filter_dataframe(df, pains_filter, keep_failures=False)

        # Check that all original columns are preserved
        assert "id" in filtered.columns
        assert "name" in filtered.columns
        assert "activity" in filtered.columns

        # Check that data is correct
        assert len(filtered) == 3  # All pass PAINS
        assert all(col in filtered.columns for col in df.columns)

    def test_lipinski_rule_of_five(self):
        """Test Lipinski's Rule of Five filtering."""
        df = pd.DataFrame({
            "smiles": [
                "CCO",  # Passes
                "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin - passes
                "CCCCCCCCCCCCCCCCCC",  # Too lipophilic - fails
            ]
        })

        filtered = filter_by_properties(
            df,
            mw_range=(0, 500),
            logp_range=(-float('inf'), 5),
            hbd_range=(0, 5),
            hba_range=(0, 10),
            keep_failures=False
        )

        # First two should pass
        assert len(filtered) >= 2
