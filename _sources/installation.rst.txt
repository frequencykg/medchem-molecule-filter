Installation
============

Requirements
------------

* Python >= 3.8
* RDKit >= 2022.9.1

Install from PyPI
-----------------

.. code-block:: bash

   pip install medchem-molecule-filter

Install from source
-------------------

.. code-block:: bash

   git clone https://github.com/frequencykg/medchem-molecule-filter.git
   cd medchem-molecule-filter
   pip install -e .

Development Installation
------------------------

For development, install with additional dependencies:

.. code-block:: bash

   pip install -e ".[dev,docs]"

This will install:

* pytest and pytest-cov for testing
* black and flake8 for code formatting
* mypy for type checking
* sphinx and related packages for documentation
