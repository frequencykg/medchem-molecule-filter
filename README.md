# medchem-molecule-filter

[![Tests](https://github.com/frequencykg/medchem-molecule-filter/workflows/Tests/badge.svg)](https://github.com/frequencykg/medchem-molecule-filter/actions)
[![Documentation](https://github.com/frequencykg/medchem-molecule-filter/workflows/Documentation/badge.svg)](https://frequencykg.github.io/medchem-molecule-filter)
[![PyPI version](https://badge.fury.io/py/medchem-molecule-filter.svg)](https://badge.fury.io/py/medchem-molecule-filter)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python library for filtering small molecules using medicinal chemistry rules and filters.

## Features

- **PAINS Filters**: Identify Pan-Assay Interference Compounds that frequently show up as false positives in high-throughput screens
- **Reactive Group Filters**: Detect chemically reactive functional groups that may cause issues in biological assays
- **Heterocycle Filters**: Identify or require the presence of heterocyclic scaffolds common in medicinal chemistry
- **Property Filters**: Filter by calculated molecular properties:
  - LogP (partition coefficient)
  - TPSA (Topological Polar Surface Area)
  - Molecular Weight
  - HBD/HBA (Hydrogen Bond Donors/Acceptors)
  - Rotatable Bonds
  - Aromatic Rings
- **Filter Groups**: Combine multiple filters for complex filtering pipelines

## Installation

```bash
pip install medchem-molecule-filter
```

### With Pandas Support

For pandas DataFrame integration:
```bash
pip install medchem-molecule-filter[pandas]
```

### From Source

```bash
git clone https://github.com/frequencykg/medchem-molecule-filter.git
cd medchem-molecule-filter
pip install -e ".[dev,docs,pandas]"
```

## Quick Start

```python
from rdkit import Chem
from medchem_filter import PAINSFilter, ReactiveFilter, PropertyFilter, FilterGroup

# Create a molecule
mol = Chem.MolFromSmiles("c1ccccc1CCO")  # Phenethyl alcohol

# Apply individual filters
pains_filter = PAINSFilter()
passes, matched_patterns = pains_filter.check_molecule(mol)
print(f"Passes PAINS filter: {passes}")

# Apply property filters
prop_filter = PropertyFilter(
    logp_range=(-1, 5),
    mw_range=(0, 500),
    hbd_range=(0, 5),
    hba_range=(0, 10)
)
passes, reasons = prop_filter.check_molecule(mol)
print(f"Passes property filter: {passes}")

# Combine multiple filters
filter_group = FilterGroup([
    PAINSFilter(),
    ReactiveFilter(),
    PropertyFilter(mw_range=(0, 500))
])

molecules = [
    Chem.MolFromSmiles("CCO"),
    Chem.MolFromSmiles("c1ccccc1"),
    Chem.MolFromSmiles("CC=O")
]

passed, failed_details = filter_group.filter_molecules(molecules, verbose=True)
print(f"Passed: {len(passed)} molecules")
```

## Pandas Integration

Work with DataFrames containing SMILES strings:

```python
import pandas as pd
from medchem_filter import PAINSFilter, filter_dataframe, calculate_properties_for_dataframe

# Create a DataFrame with SMILES
df = pd.DataFrame({
    "compound_id": ["CMPD001", "CMPD002", "CMPD003"],
    "smiles": ["CCO", "c1ccccc1", "O=C1C=CC(=O)C=C1"]
})

# Filter the DataFrame
filtered_df = filter_dataframe(df, PAINSFilter(), keep_failures=False)

# Calculate properties for all compounds
df_with_props = calculate_properties_for_dataframe(df)
print(df_with_props[["compound_id", "smiles", "molecular_weight", "logP"]])

# Apply multiple filters with detailed results
from medchem_filter import apply_filters_to_dataframe, ReactiveFilter

filters = [PAINSFilter(), ReactiveFilter()]
result = apply_filters_to_dataframe(df, filters, detailed=True)
```

See `examples/pandas_integration.py` for more comprehensive examples.

## Filter Types

### PAINS Filter

Identifies Pan-Assay Interference Compounds using SMARTS patterns for common problematic structures:
- Quinones and derivatives
- Catechols
- Michael acceptors
- Rhodanines
- Phenolic Mannich bases
- And more...

```python
from medchem_filter import PAINSFilter

filter = PAINSFilter()
passes, matched = filter.check_molecule(mol)
```

### Reactive Filter

Detects chemically reactive functional groups:
- Acyl halides
- Aldehydes
- Epoxides
- Isocyanates
- Michael acceptors
- Peroxides
- And more...

```python
from medchem_filter import ReactiveFilter

filter = ReactiveFilter()
passes, matched = filter.check_molecule(mol)
```

### Heterocycle Filter

Identifies heterocyclic scaffolds. Can be used to either exclude or require heterocycles:

```python
from medchem_filter import HeterocycleFilter

# Exclude heterocycles
filter = HeterocycleFilter(require_heterocycle=False)
passes, matched = filter.check_molecule(mol)

# Require heterocycles
filter = HeterocycleFilter(require_heterocycle=True)
passes, matched = filter.check_molecule(mol)
```

### Property Filter

Filter by calculated molecular properties with customizable ranges:

```python
from medchem_filter import PropertyFilter

filter = PropertyFilter(
    logp_range=(-1, 5),          # Lipophilicity
    tpsa_range=(0, 140),         # Topological Polar Surface Area (Ų)
    mw_range=(0, 500),           # Molecular Weight (g/mol)
    hbd_range=(0, 5),            # Hydrogen Bond Donors
    hba_range=(0, 10),           # Hydrogen Bond Acceptors
    rotatable_bonds_max=10,      # Maximum rotatable bonds
    aromatic_rings_range=(1, 3)  # Number of aromatic rings
)
passes, reasons = filter.check_molecule(mol)
```

### Filter Groups

Combine multiple filters for complex pipelines:

```python
from medchem_filter import FilterGroup, PAINSFilter, ReactiveFilter, PropertyFilter

# Create a comprehensive filter
filter_group = FilterGroup([
    PAINSFilter(),
    ReactiveFilter(),
    PropertyFilter(
        mw_range=(200, 500),
        logp_range=(0, 5),
        hbd_range=(0, 5),
        hba_range=(0, 10)
    )
])

# Check single molecule
passes, details = filter_group.check_molecule(mol)

# Filter multiple molecules
passed, failed_details = filter_group.filter_molecules(molecules)
```

## Molecular Properties

Calculate molecular properties directly:

```python
from medchem_filter import MolecularProperties

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
```

## Documentation

Full documentation is available at: https://frequencykg.github.io/medchem-molecule-filter

## Filter Rules Reference

### PAINS Patterns

The PAINS filter includes patterns for:
- Quinones (quinone_A, quinone_B)
- Catechols (catechol, catechol_A)
- Michael acceptors (acyl_hydrazone, chalcone)
- Rhodanines (rhodanine, thiazolidinone, ene_rhodanine)
- Alkylidene barbiturates
- Hydroxyphenyl hydrazones
- Phenolic Mannich bases
- Curcumin-like structures
- Beta-lactams
- Styrene-like structures
- Nitro aromatics
- Azo compounds
- Thioureas
- Isothiocyanates

### Reactive Patterns

The reactive filter includes patterns for:
- Electrophiles (acyl_halide, aldehyde, alkyl_halide, anhydride, epoxide, isocyanate, sulfonyl_halide)
- Michael acceptors (alpha_beta_unsaturated_carbonyl, vinyl_sulfone, vinyl_sulfonamide)
- Peroxides (peroxide, peroxy_acid)
- Diazo compounds (diazo, diazonium)
- Acid halides (sulfonyl_chloride, phosphoryl_halide)
- Metal chelators (phosphonate, hydroxamic_acid)
- Thiols
- Azides
- Hydrazines

### Heterocycle Patterns

The heterocycle filter includes patterns for:
- Five-membered rings: pyrrole, furan, thiophene, imidazole, pyrazole, oxazole, isoxazole, thiazole
- Six-membered rings: pyridine, pyran, pyrimidine, pyrazine, pyridazine, oxazine, triazine
- Fused systems: indole, benzofuran, benzothiophene, benzimidazole, quinoline, isoquinoline, quinazoline, quinoxaline, purine
- Saturated rings: piperidine, piperazine, morpholine, thiomorpholine, pyrrolidine
- Seven-membered rings: azepine, oxepine

## Requirements

- Python >= 3.8
- RDKit >= 2022.9.1

## Development

Install development dependencies:

```bash
pip install -e ".[dev,docs]"
```

### Running Tests

Run all tests:
```bash
pytest tests/ -v
```

Run tests with coverage:
```bash
pytest tests/ --cov=medchem_filter --cov-report=html
```

### Code Quality

Format code with Black:
```bash
black medchem_filter/ tests/ --line-length 100
```

Lint with flake8:
```bash
flake8 medchem_filter/ tests/ --max-line-length 100 --extend-ignore E203,W503
```

### Building Documentation

```bash
cd docs
make html
```

The built documentation will be in `docs/_build/html/`.

## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass (`pytest tests/ -v`)
6. Format code (`black medchem_filter/ tests/`)
7. Commit your changes (`git commit -m 'feat: add amazing feature'`)
8. Push to the branch (`git push origin feature/amazing-feature`)
9. Open a Pull Request

See `.github/copilot-instructions.md` for detailed development guidelines.

## Publishing

### Setup for PyPI

1. Create a PyPI account at https://pypi.org/account/register/
2. Generate an API token in Account Settings
3. Add token to GitHub Secrets as `PYPI_API_TOKEN`
4. Create a release on GitHub to trigger automatic publishing

### Setup for GitHub Pages

1. Go to repository Settings → Pages
2. Select "gh-pages" branch as source
3. Documentation will be automatically built and deployed on push to main

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use this library in your research, please cite:

```
medchem-molecule-filter
https://github.com/frequencykg/medchem-molecule-filter
```
