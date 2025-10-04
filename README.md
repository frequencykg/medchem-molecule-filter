# medchem-molecule-filter

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

Run tests:

```bash
pytest tests/ -v
```

Run tests with coverage:

```bash
pytest tests/ --cov=medchem_filter --cov-report=html
```

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
