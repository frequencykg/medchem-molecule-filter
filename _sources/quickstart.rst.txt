Quick Start
===========

Basic Usage
-----------

Filtering with PAINS
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from rdkit import Chem
   from medchem_filter import PAINSFilter

   # Create filter
   pains_filter = PAINSFilter()

   # Check a molecule
   mol = Chem.MolFromSmiles("c1ccccc1")
   passes, matched_patterns = pains_filter.check_molecule(mol)
   
   if passes:
       print("Molecule passed PAINS filter")
   else:
       print(f"Molecule failed: {matched_patterns}")

Property-based Filtering
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from rdkit import Chem
   from medchem_filter import PropertyFilter

   # Create filter with property ranges
   prop_filter = PropertyFilter(
       logp_range=(-1, 5),
       mw_range=(0, 500),
       hbd_range=(0, 5),
       hba_range=(0, 10)
   )

   mol = Chem.MolFromSmiles("CCO")
   passes, reasons = prop_filter.check_molecule(mol)

Combining Multiple Filters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from rdkit import Chem
   from medchem_filter import FilterGroup, PAINSFilter, ReactiveFilter, PropertyFilter

   # Create a comprehensive filter pipeline
   filter_group = FilterGroup([
       PAINSFilter(),
       ReactiveFilter(),
       PropertyFilter(
           mw_range=(200, 500),
           logp_range=(0, 5)
       )
   ])

   # Check single molecule
   mol = Chem.MolFromSmiles("c1ccccc1CCO")
   passes, details = filter_group.check_molecule(mol)
   
   print(f"Passes all filters: {passes}")
   if not passes:
       for filter_name, reasons in details.items():
           print(f"Failed {filter_name}: {reasons}")

Filtering Multiple Molecules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from rdkit import Chem
   from medchem_filter import FilterGroup, PAINSFilter, ReactiveFilter

   # Create molecules
   molecules = [
       Chem.MolFromSmiles("CCO"),           # Ethanol
       Chem.MolFromSmiles("c1ccccc1"),      # Benzene
       Chem.MolFromSmiles("CC=O"),          # Acetaldehyde (reactive)
   ]

   # Apply filters
   filter_group = FilterGroup([PAINSFilter(), ReactiveFilter()])
   passed, failed_details = filter_group.filter_molecules(
       molecules, 
       verbose=True
   )

   print(f"Passed: {len(passed)} molecules")
   print(f"Failed: {len(failed_details)} molecules")

Calculating Properties
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from rdkit import Chem
   from medchem_filter import MolecularProperties

   mol = Chem.MolFromSmiles("c1ccccc1CCO")
   
   # Calculate individual properties
   logp = MolecularProperties.calculate_logp(mol)
   tpsa = MolecularProperties.calculate_tpsa(mol)
   hbd = MolecularProperties.calculate_hbd(mol)
   
   # Calculate all properties at once
   props = MolecularProperties.calculate_all_properties(mol)
   print(props)
   # {
   #     'logP': 1.5,
   #     'tpsa': 20.23,
   #     'hbd': 1,
   #     'hba': 1,
   #     'molecular_weight': 122.16,
   #     'rotatable_bonds': 2,
   #     'aromatic_rings': 1
   # }

Working with Heterocycles
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from rdkit import Chem
   from medchem_filter import HeterocycleFilter

   # Exclude molecules with heterocycles
   filter_exclude = HeterocycleFilter(require_heterocycle=False)
   mol = Chem.MolFromSmiles("c1ccncc1")  # Pyridine
   passes, matched = filter_exclude.check_molecule(mol)
   # passes = False (contains heterocycle)

   # Require molecules with heterocycles
   filter_require = HeterocycleFilter(require_heterocycle=True)
   passes, matched = filter_require.check_molecule(mol)
   # passes = True (contains required heterocycle)

Custom Filter Patterns
~~~~~~~~~~~~~~~~~~~~~~

You can provide custom SMARTS patterns to any filter:

.. code-block:: python

   from medchem_filter import PAINSFilter

   # Custom PAINS patterns
   custom_patterns = {
       "my_pattern_1": "c1ccccc1N=N",
       "my_pattern_2": "C(=O)O[N+](=O)[O-]",
   }

   filter = PAINSFilter(custom_patterns=custom_patterns)
