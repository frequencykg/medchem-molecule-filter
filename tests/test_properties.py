"""Tests for molecular property calculations."""

from rdkit import Chem

from medchem_filter.properties import MolecularProperties


class TestMolecularProperties:
    """Test molecular property calculations."""

    def test_calculate_logp(self):
        """Test logP calculation."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        logp = MolecularProperties.calculate_logp(mol)
        assert isinstance(logp, float)
        assert -1.0 < logp < 1.0  # Ethanol should have low logP

    def test_calculate_tpsa(self):
        """Test TPSA calculation."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        tpsa = MolecularProperties.calculate_tpsa(mol)
        assert isinstance(tpsa, float)
        assert tpsa > 0  # Should have some polar surface area

    def test_calculate_hbd(self):
        """Test HBD calculation."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol - 1 HBD
        hbd = MolecularProperties.calculate_hbd(mol)
        assert hbd == 1

        mol = Chem.MolFromSmiles("O=C(O)C")  # Acetic acid - 1 HBD
        hbd = MolecularProperties.calculate_hbd(mol)
        assert hbd == 1

    def test_calculate_hba(self):
        """Test HBA calculation."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol - 1 HBA
        hba = MolecularProperties.calculate_hba(mol)
        assert hba == 1

        mol = Chem.MolFromSmiles("CCOC")  # Diethyl ether - 1 HBA
        hba = MolecularProperties.calculate_hba(mol)
        assert hba == 1

    def test_calculate_molecular_weight(self):
        """Test molecular weight calculation."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        mw = MolecularProperties.calculate_molecular_weight(mol)
        assert isinstance(mw, float)
        assert 45 < mw < 47  # Ethanol MW is ~46

    def test_calculate_rotatable_bonds(self):
        """Test rotatable bonds calculation."""
        mol = Chem.MolFromSmiles("CCCC")  # Butane - 1 rotatable bond (central C-C)
        rot = MolecularProperties.calculate_rotatable_bonds(mol)
        assert rot == 1

        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene - 0 rotatable bonds
        rot = MolecularProperties.calculate_rotatable_bonds(mol)
        assert rot == 0

    def test_calculate_aromatic_rings(self):
        """Test aromatic rings calculation."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene - 1 aromatic ring
        ar = MolecularProperties.calculate_aromatic_rings(mol)
        assert ar == 1

        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")  # Naphthalene - 2 aromatic rings
        ar = MolecularProperties.calculate_aromatic_rings(mol)
        assert ar == 2

    def test_calculate_all_properties(self):
        """Test calculating all properties at once."""
        mol = Chem.MolFromSmiles("CCO")
        props = MolecularProperties.calculate_all_properties(mol)

        assert "logP" in props
        assert "tpsa" in props
        assert "hbd" in props
        assert "hba" in props
        assert "molecular_weight" in props
        assert "rotatable_bonds" in props
        assert "aromatic_rings" in props

        assert isinstance(props["logP"], float)
        assert isinstance(props["tpsa"], float)
        assert isinstance(props["hbd"], int)
        assert isinstance(props["hba"], int)
        assert isinstance(props["molecular_weight"], float)
        assert isinstance(props["rotatable_bonds"], int)
        assert isinstance(props["aromatic_rings"], int)
