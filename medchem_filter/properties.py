"""
Molecular property calculators for filtering based on calculated properties.
"""

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski


class MolecularProperties:
    """Calculate various molecular properties for filtering."""

    @staticmethod
    def calculate_logp(mol: Chem.Mol) -> float:
        """Calculate the partition coefficient (logP).

        Args:
            mol: RDKit molecule object

        Returns:
            LogP value
        """
        return Crippen.MolLogP(mol)

    @staticmethod
    def calculate_tpsa(mol: Chem.Mol) -> float:
        """Calculate the Topological Polar Surface Area (TPSA).

        Args:
            mol: RDKit molecule object

        Returns:
            TPSA value in Ų
        """
        return Descriptors.TPSA(mol)

    @staticmethod
    def calculate_hbd(mol: Chem.Mol) -> int:
        """Calculate the number of Hydrogen Bond Donors (HBD).

        Args:
            mol: RDKit molecule object

        Returns:
            Number of HBDs
        """
        return Lipinski.NumHDonors(mol)

    @staticmethod
    def calculate_hba(mol: Chem.Mol) -> int:
        """Calculate the number of Hydrogen Bond Acceptors (HBA).

        Args:
            mol: RDKit molecule object

        Returns:
            Number of HBAs
        """
        return Lipinski.NumHAcceptors(mol)

    @staticmethod
    def calculate_molecular_weight(mol: Chem.Mol) -> float:
        """Calculate the molecular weight.

        Args:
            mol: RDKit molecule object

        Returns:
            Molecular weight in g/mol
        """
        return Descriptors.MolWt(mol)

    @staticmethod
    def calculate_rotatable_bonds(mol: Chem.Mol) -> int:
        """Calculate the number of rotatable bonds.

        Args:
            mol: RDKit molecule object

        Returns:
            Number of rotatable bonds
        """
        return Lipinski.NumRotatableBonds(mol)

    @staticmethod
    def calculate_aromatic_rings(mol: Chem.Mol) -> int:
        """Calculate the number of aromatic rings.

        Args:
            mol: RDKit molecule object

        Returns:
            Number of aromatic rings
        """
        return Lipinski.NumAromaticRings(mol)

    @staticmethod
    def calculate_all_properties(mol: Chem.Mol) -> dict:
        """Calculate all available molecular properties.

        Args:
            mol: RDKit molecule object

        Returns:
            Dictionary containing all calculated properties
        """
        return {
            "logP": MolecularProperties.calculate_logp(mol),
            "tpsa": MolecularProperties.calculate_tpsa(mol),
            "hbd": MolecularProperties.calculate_hbd(mol),
            "hba": MolecularProperties.calculate_hba(mol),
            "molecular_weight": MolecularProperties.calculate_molecular_weight(mol),
            "rotatable_bonds": MolecularProperties.calculate_rotatable_bonds(mol),
            "aromatic_rings": MolecularProperties.calculate_aromatic_rings(mol),
        }
