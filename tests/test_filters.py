"""Tests for molecular filters."""

import pytest
from rdkit import Chem

from medchem_filter.filters import (
    FilterGroup,
    HeterocycleFilter,
    PAINSFilter,
    PropertyFilter,
    ReactiveFilter,
)


class TestPAINSFilter:
    """Test PAINS filter."""

    def test_init(self):
        """Test filter initialization."""
        filter_obj = PAINSFilter()
        assert filter_obj is not None
        assert len(filter_obj._compiled_patterns) > 0

    def test_quinone_detection(self):
        """Test detection of quinone PAINS pattern."""
        filter_obj = PAINSFilter()
        mol = Chem.MolFromSmiles("O=C1C=CC(=O)C=C1")  # p-Benzoquinone
        passes, matched = filter_obj.check_molecule(mol)
        assert not passes  # Should fail (contains PAINS)
        assert len(matched) > 0

    def test_clean_molecule(self):
        """Test that clean molecules pass."""
        filter_obj = PAINSFilter()
        mol = Chem.MolFromSmiles("CCO")  # Ethanol - clean
        passes, matched = filter_obj.check_molecule(mol)
        assert passes
        assert len(matched) == 0


class TestReactiveFilter:
    """Test reactive group filter."""

    def test_init(self):
        """Test filter initialization."""
        filter_obj = ReactiveFilter()
        assert filter_obj is not None
        assert len(filter_obj._compiled_patterns) > 0

    def test_aldehyde_detection(self):
        """Test detection of aldehyde."""
        filter_obj = ReactiveFilter()
        mol = Chem.MolFromSmiles("CC=O")  # Acetaldehyde
        passes, matched = filter_obj.check_molecule(mol)
        assert not passes  # Should fail (contains reactive group)
        assert len(matched) > 0

    def test_clean_molecule(self):
        """Test that non-reactive molecules pass."""
        filter_obj = ReactiveFilter()
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene - clean
        passes, matched = filter_obj.check_molecule(mol)
        assert passes
        assert len(matched) == 0


class TestHeterocycleFilter:
    """Test heterocycle filter."""

    def test_init_default(self):
        """Test filter initialization."""
        filter_obj = HeterocycleFilter()
        assert filter_obj is not None
        assert len(filter_obj._compiled_patterns) > 0

    def test_pyridine_detection_exclude(self):
        """Test detection of pyridine when excluding heterocycles."""
        filter_obj = HeterocycleFilter(require_heterocycle=False)
        mol = Chem.MolFromSmiles("c1ccncc1")  # Pyridine
        passes, matched = filter_obj.check_molecule(mol)
        assert not passes  # Should fail (contains heterocycle)
        assert len(matched) > 0

    def test_pyridine_detection_require(self):
        """Test detection of pyridine when requiring heterocycles."""
        filter_obj = HeterocycleFilter(require_heterocycle=True)
        mol = Chem.MolFromSmiles("c1ccncc1")  # Pyridine
        passes, matched = filter_obj.check_molecule(mol)
        assert passes  # Should pass (contains required heterocycle)
        assert len(matched) > 0

    def test_benzene_exclude(self):
        """Test that benzene passes when excluding heterocycles."""
        filter_obj = HeterocycleFilter(require_heterocycle=False)
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
        passes, matched = filter_obj.check_molecule(mol)
        assert passes  # Should pass (no heterocycle)
        assert len(matched) == 0

    def test_benzene_require(self):
        """Test that benzene fails when requiring heterocycles."""
        filter_obj = HeterocycleFilter(require_heterocycle=True)
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
        passes, matched = filter_obj.check_molecule(mol)
        assert not passes  # Should fail (no heterocycle)


class TestPropertyFilter:
    """Test property-based filter."""

    def test_init(self):
        """Test filter initialization."""
        filter_obj = PropertyFilter(logp_range=(-1, 5))
        assert filter_obj is not None

    def test_logp_filter_pass(self):
        """Test logP filter with passing molecule."""
        filter_obj = PropertyFilter(logp_range=(-1, 5))
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        passes, reasons = filter_obj.check_molecule(mol)
        assert passes

    def test_logp_filter_fail(self):
        """Test logP filter with failing molecule."""
        filter_obj = PropertyFilter(logp_range=(-10, -5))
        mol = Chem.MolFromSmiles("CCCCCCCCCC")  # Decane (high logP)
        passes, reasons = filter_obj.check_molecule(mol)
        assert not passes
        assert len(reasons) > 0

    def test_mw_filter(self):
        """Test molecular weight filter."""
        filter_obj = PropertyFilter(mw_range=(0, 50))
        mol = Chem.MolFromSmiles("CCO")  # Ethanol, MW ~46
        passes, reasons = filter_obj.check_molecule(mol)
        assert passes

        mol = Chem.MolFromSmiles("CCCCCCCCCC")  # Decane, MW ~142
        passes, reasons = filter_obj.check_molecule(mol)
        assert not passes

    def test_hbd_filter(self):
        """Test HBD filter."""
        filter_obj = PropertyFilter(hbd_range=(0, 2))
        mol = Chem.MolFromSmiles("CCO")  # Ethanol, 1 HBD
        passes, reasons = filter_obj.check_molecule(mol)
        assert passes

    def test_multiple_criteria(self):
        """Test filter with multiple criteria."""
        filter_obj = PropertyFilter(logp_range=(-1, 5), mw_range=(0, 100), hbd_range=(0, 2))
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        passes, reasons = filter_obj.check_molecule(mol)
        assert passes


class TestFilterGroup:
    """Test filter group."""

    def test_init(self):
        """Test filter group initialization."""
        filters = [PAINSFilter(), ReactiveFilter()]
        group = FilterGroup(filters)
        assert group is not None
        assert len(group.filters) == 2

    def test_single_filter(self):
        """Test filter group with single filter."""
        group = FilterGroup([PAINSFilter()])
        mol = Chem.MolFromSmiles("CCO")  # Ethanol - clean
        passes, details = group.check_molecule(mol)
        assert passes

    def test_multiple_filters_pass(self):
        """Test filter group with multiple filters, molecule passes all."""
        group = FilterGroup([PAINSFilter(), ReactiveFilter()])
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene - clean
        passes, details = group.check_molecule(mol)
        assert passes

    def test_multiple_filters_fail(self):
        """Test filter group with multiple filters, molecule fails one."""
        group = FilterGroup([PAINSFilter(), ReactiveFilter()])
        mol = Chem.MolFromSmiles("CC=O")  # Acetaldehyde - reactive
        passes, details = group.check_molecule(mol)
        assert not passes
        assert "ReactiveFilter" in details

    def test_filter_molecules_list(self):
        """Test filtering a list of molecules."""
        group = FilterGroup([ReactiveFilter()])
        molecules = [
            Chem.MolFromSmiles("CCO"),  # Ethanol - clean
            Chem.MolFromSmiles("CC=O"),  # Acetaldehyde - reactive
            Chem.MolFromSmiles("c1ccccc1"),  # Benzene - clean
        ]
        passed, details = group.filter_molecules(molecules)
        assert len(passed) == 2  # Ethanol and benzene should pass
        assert len(details) == 1  # Acetaldehyde should fail

    def test_combined_filters(self):
        """Test combining structure and property filters."""
        prop_filter = PropertyFilter(mw_range=(0, 100))
        group = FilterGroup([ReactiveFilter(), prop_filter])

        mol = Chem.MolFromSmiles("CCO")  # Ethanol - passes both
        passes, details = group.check_molecule(mol)
        assert passes
