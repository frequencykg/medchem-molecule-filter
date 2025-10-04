"""
medchem_filter: A Python library for filtering small molecules using medicinal chemistry rules.

This library provides easy-to-apply filter groups for filtering lists of small molecules,
including PAINS filters, reactive group filters, heterocycle filters, and property-based
filters (logD, TPSA, HBD/HBA, etc.).
"""

__version__ = "0.1.0"

from .filters import (
    PAINSFilter,
    ReactiveFilter,
    HeterocycleFilter,
    PropertyFilter,
    FilterGroup,
)
from .properties import MolecularProperties

__all__ = [
    "PAINSFilter",
    "ReactiveFilter",
    "HeterocycleFilter",
    "PropertyFilter",
    "FilterGroup",
    "MolecularProperties",
]
