"""
medchem_filter: A Python library for filtering small molecules using medicinal chemistry rules.

This library provides easy-to-apply filter groups for filtering lists of small molecules,
including PAINS filters, reactive group filters, heterocycle filters, and property-based
filters (logD, TPSA, HBD/HBA, etc.).
"""

__version__ = "0.1.0"

from .filters import (
    FilterGroup,
    HeterocycleFilter,
    PAINSFilter,
    PropertyFilter,
    ReactiveFilter,
)
from .properties import MolecularProperties

# Pandas utilities - optional import
try:
    from .pandas_utils import (  # noqa: F401
        add_mol_column,
        apply_filters_to_dataframe,
        calculate_properties_for_dataframe,
        filter_by_properties,
        filter_dataframe,
    )

    _PANDAS_AVAILABLE = True
except ImportError:
    _PANDAS_AVAILABLE = False

__all__ = [
    "PAINSFilter",
    "ReactiveFilter",
    "HeterocycleFilter",
    "PropertyFilter",
    "FilterGroup",
    "MolecularProperties",
]

# Add pandas utilities to __all__ if available
if _PANDAS_AVAILABLE:
    __all__.extend(
        [
            "add_mol_column",
            "filter_dataframe",
            "apply_filters_to_dataframe",
            "calculate_properties_for_dataframe",
            "filter_by_properties",
        ]
    )
