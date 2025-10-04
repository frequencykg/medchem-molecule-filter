medchem-molecule-filter Documentation
======================================

A Python library for filtering small molecules using medicinal chemistry rules and filters.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   api
   filter_rules

Features
--------

* **PAINS Filters**: Identify Pan-Assay Interference Compounds
* **Reactive Group Filters**: Detect chemically reactive functional groups
* **Heterocycle Filters**: Identify heterocyclic scaffolds
* **Property Filters**: Filter by calculated molecular properties
* **Filter Groups**: Combine multiple filters for complex pipelines

Quick Example
-------------

.. code-block:: python

   from rdkit import Chem
   from medchem_filter import PAINSFilter, PropertyFilter, FilterGroup

   # Create a molecule
   mol = Chem.MolFromSmiles("c1ccccc1CCO")

   # Apply filters
   pains_filter = PAINSFilter()
   passes, matched = pains_filter.check_molecule(mol)
   
   # Combine filters
   filter_group = FilterGroup([
       PAINSFilter(),
       PropertyFilter(mw_range=(0, 500))
   ])
   
   passes, details = filter_group.check_molecule(mol)

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
