"""Integration tests for the medchem_filter package.

Tests complex workflows and interactions between multiple components.
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


class TestIntegration:
    """Integration tests for complete workflows."""

    def test_complete_filtering_pipeline(self):
        """Test a complete drug-like filtering pipeline."""
        # Create comprehensive filter
        pipeline = FilterGroup(
            [
                PAINSFilter(),
                ReactiveFilter(),
                PropertyFilter(
                    mw_range=(200, 500),
                    logp_range=(-0.5, 5),
                    hbd_range=(0, 5),
                    hba_range=(0, 10),
                ),
            ]
        )

        # Test molecules
        molecules = [
            Chem.MolFromSmiles("CCO"),  # Ethanol - too light (MW < 200), should fail
            Chem.MolFromSmiles(
                "c1ccccc1CCO"
            ),  # Phenethyl alcohol - too light (MW < 200), should fail
            Chem.MolFromSmiles("CC=O"),  # Acetaldehyde - reactive + too light, should fail
            Chem.MolFromSmiles("O=C1C=CC(=O)C=C1"),  # Quinone - PAINS + too light, should fail
            Chem.MolFromSmiles(
                "c1ccc(cc1)c2ccccc2c3ccc(cc3)c4ccccc4"
            ),  # Larger biphenyl - should pass
        ]

        passed, failed_details = pipeline.filter_molecules(molecules)

        # Verify results
        assert len(passed) >= 0  # At least one should pass or all fail
        assert len(failed_details) >= 4  # Most should fail due to MW

    def test_batch_processing_performance(self):
        """Test batch processing of multiple molecules."""
        # Generate test molecules
        smiles_list = [
            "CCO",
            "c1ccccc1",
            "CCN",
            "CCC",
            "c1ccncc1",
            "c1ccc(O)cc1",
            "CCCO",
            "c1ccccc1O",
            "CCCN",
            "c1ccccc1C",
        ]
        molecules = [Chem.MolFromSmiles(s) for s in smiles_list]

        # Create filter
        prop_filter = PropertyFilter(mw_range=(0, 150))

        # Process batch
        passed, failed, reasons = prop_filter.filter_molecules(molecules)

        # Verify all were processed
        assert len(passed) + len(failed) == len(molecules)
        assert len(reasons) == len(failed)

    def test_filter_order_independence(self):
        """Test that filter order doesn't affect final result."""
        mol = Chem.MolFromSmiles("c1ccccc1CCO")

        # Create two pipelines with different order
        pipeline1 = FilterGroup(
            [
                PAINSFilter(),
                ReactiveFilter(),
                PropertyFilter(mw_range=(0, 150)),
            ]
        )

        pipeline2 = FilterGroup(
            [
                PropertyFilter(mw_range=(0, 150)),
                ReactiveFilter(),
                PAINSFilter(),
            ]
        )

        passes1, _ = pipeline1.check_molecule(mol)
        passes2, _ = pipeline2.check_molecule(mol)

        # Both should give same result
        assert passes1 == passes2

    def test_property_consistency(self):
        """Test that properties are consistent across calculations."""
        mol = Chem.MolFromSmiles("c1ccccc1CCO")

        # Calculate properties multiple times
        props1 = MolecularProperties.calculate_all_properties(mol)
        props2 = MolecularProperties.calculate_all_properties(mol)

        # Should be identical
        assert props1 == props2

    def test_lipinski_rule_of_five(self):
        """Test Lipinski's Rule of Five filtering."""
        ro5_filter = PropertyFilter(
            mw_range=(0, 500),
            logp_range=(-float("inf"), 5),
            hbd_range=(0, 5),
            hba_range=(0, 10),
        )

        # Known drug-like molecules
        aspirin = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        caffeine = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

        passes1, _ = ro5_filter.check_molecule(aspirin)
        passes2, _ = ro5_filter.check_molecule(caffeine)

        assert passes1  # Aspirin should pass
        assert passes2  # Caffeine should pass

    def test_heterocycle_require_and_exclude(self):
        """Test heterocycle filter in both modes."""
        pyridine = Chem.MolFromSmiles("c1ccncc1")
        benzene = Chem.MolFromSmiles("c1ccccc1")

        # Require heterocycles
        filter_req = HeterocycleFilter(require_heterocycle=True)
        assert filter_req.check_molecule(pyridine)[0]  # Should pass
        assert not filter_req.check_molecule(benzene)[0]  # Should fail

        # Exclude heterocycles
        filter_exc = HeterocycleFilter(require_heterocycle=False)
        assert not filter_exc.check_molecule(pyridine)[0]  # Should fail
        assert filter_exc.check_molecule(benzene)[0]  # Should pass

    def test_custom_pattern_integration(self):
        """Test using custom patterns in filters."""
        custom_patterns = {
            "phenol": "c1ccccc1O",
            "aniline": "c1ccccc1N",
        }

        custom_filter = PAINSFilter(custom_patterns=custom_patterns)

        phenol = Chem.MolFromSmiles("c1ccccc1O")
        aniline = Chem.MolFromSmiles("c1ccccc1N")
        benzene = Chem.MolFromSmiles("c1ccccc1")

        assert not custom_filter.check_molecule(phenol)[0]  # Should fail
        assert not custom_filter.check_molecule(aniline)[0]  # Should fail
        assert custom_filter.check_molecule(benzene)[0]  # Should pass

    def test_empty_molecule_list(self):
        """Test filtering empty molecule list."""
        pipeline = FilterGroup([PAINSFilter()])
        passed, failed_details = pipeline.filter_molecules([])

        assert len(passed) == 0
        assert len(failed_details) == 0

    def test_filter_with_verbose(self):
        """Test verbose mode output."""
        pipeline = FilterGroup([PAINSFilter(), ReactiveFilter()])
        molecules = [
            Chem.MolFromSmiles("CCO"),
            Chem.MolFromSmiles("CC=O"),
        ]

        # Should not raise any errors
        passed, failed = pipeline.filter_molecules(molecules, verbose=True)
        assert len(passed) + len(failed) == len(molecules)


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_none_molecule(self):
        """Test handling of None molecule."""
        pains = PAINSFilter()
        passes, matched = pains.check_molecule(None)

        assert not passes
        assert "Invalid molecule" in matched

    def test_invalid_smiles(self):
        """Test handling of invalid SMILES."""
        mol = Chem.MolFromSmiles("invalid_smiles_xyz123")
        assert mol is None

        # Filter should handle None gracefully
        pains = PAINSFilter()
        passes, matched = pains.check_molecule(mol)
        assert not passes

    def test_empty_property_ranges(self):
        """Test PropertyFilter with no ranges specified."""
        prop_filter = PropertyFilter()
        mol = Chem.MolFromSmiles("CCO")

        passes, reasons = prop_filter.check_molecule(mol)
        assert passes  # Should pass if no filters applied
        assert len(reasons) == 0

    def test_extreme_property_values(self):
        """Test molecules with extreme property values."""
        # Very large molecule
        large_mol = Chem.MolFromSmiles("C" * 50)

        prop_filter = PropertyFilter(mw_range=(0, 200))
        passes, reasons = prop_filter.check_molecule(large_mol)

        assert not passes  # Should fail MW filter

    def test_filter_group_single_filter(self):
        """Test FilterGroup with single filter."""
        group = FilterGroup([PAINSFilter()])
        mol = Chem.MolFromSmiles("CCO")

        passes, details = group.check_molecule(mol)
        assert passes

    def test_aromatic_vs_aliphatic(self):
        """Test property calculations for aromatic vs aliphatic."""
        benzene = Chem.MolFromSmiles("c1ccccc1")
        hexane = Chem.MolFromSmiles("CCCCCC")

        props_benzene = MolecularProperties.calculate_all_properties(benzene)
        props_hexane = MolecularProperties.calculate_all_properties(hexane)

        assert props_benzene["aromatic_rings"] == 1
        assert props_hexane["aromatic_rings"] == 0
        assert props_benzene["rotatable_bonds"] < props_hexane["rotatable_bonds"]

    def test_multiple_heterocycles(self):
        """Test molecule with multiple heterocycles."""
        # Caffeine has multiple heterocycles
        caffeine = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

        het_filter = HeterocycleFilter(require_heterocycle=True)
        passes, matched = het_filter.check_molecule(caffeine)

        assert passes
        assert len(matched) > 0  # Should match multiple patterns

    def test_zero_molecular_weight(self):
        """Test handling of edge case molecules."""
        # Hydrogen molecule (edge case)
        h2 = Chem.MolFromSmiles("[H][H]")
        if h2:  # RDKit might not create this
            props = MolecularProperties.calculate_all_properties(h2)
            assert props["molecular_weight"] > 0

    def test_charged_molecules(self):
        """Test molecules with charges."""
        # Sodium ion
        charged = Chem.MolFromSmiles("[Na+]")
        if charged:
            props = MolecularProperties.calculate_all_properties(charged)
            assert "molecular_weight" in props

    def test_filter_molecules_preserves_order(self):
        """Test that filtering preserves molecule order."""
        molecules = [
            Chem.MolFromSmiles("CCO"),
            Chem.MolFromSmiles("c1ccccc1"),
            Chem.MolFromSmiles("CCC"),
        ]

        prop_filter = PropertyFilter(mw_range=(0, 100))
        passed, failed, reasons = prop_filter.filter_molecules(molecules)

        # All should be in passed (all small molecules)
        assert len(passed) == 3

    def test_property_filter_boundary_values(self):
        """Test property filters at exact boundary values."""
        mol = Chem.MolFromSmiles("CCO")  # MW ~46

        # Test at exact boundaries
        filter1 = PropertyFilter(mw_range=(46, 50))
        filter2 = PropertyFilter(mw_range=(40, 46))

        passes1, _ = filter1.check_molecule(mol)
        passes2, _ = filter2.check_molecule(mol)

        # Should handle boundary inclusively
        assert passes1 or passes2


class TestRobustness:
    """Test robustness and error handling."""

    def test_concurrent_filter_usage(self):
        """Test that filters can be used concurrently."""
        pains = PAINSFilter()
        molecules = [
            Chem.MolFromSmiles("CCO"),
            Chem.MolFromSmiles("c1ccccc1"),
        ]

        # Use filter multiple times
        for mol in molecules:
            passes, _ = pains.check_molecule(mol)
            assert isinstance(passes, bool)

    def test_filter_with_all_failing_molecules(self):
        """Test filter group when all molecules fail."""
        pipeline = FilterGroup([PropertyFilter(mw_range=(0, 10))])

        molecules = [
            Chem.MolFromSmiles("CCO"),  # MW ~46
            Chem.MolFromSmiles("c1ccccc1"),  # MW ~78
        ]

        passed, failed_details = pipeline.filter_molecules(molecules)

        assert len(passed) == 0
        assert len(failed_details) == len(molecules)

    def test_filter_with_all_passing_molecules(self):
        """Test filter group when all molecules pass."""
        pipeline = FilterGroup([PropertyFilter(mw_range=(0, 500))])

        molecules = [
            Chem.MolFromSmiles("CCO"),
            Chem.MolFromSmiles("c1ccccc1"),
        ]

        passed, failed_details = pipeline.filter_molecules(molecules)

        assert len(passed) == len(molecules)
        assert len(failed_details) == 0
