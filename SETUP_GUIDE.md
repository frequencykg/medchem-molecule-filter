# medchem-molecule-filter - Setup Guide

This document provides instructions for completing the setup of the medchem-molecule-filter library.

## Current Status ✅

The library is fully implemented with:

- ✅ Complete Python package structure
- ✅ PAINS, Reactive, Heterocycle, and Property filters
- ✅ Comprehensive test suite (31 tests, all passing)
- ✅ Full Sphinx documentation
- ✅ GitHub Actions workflows
- ✅ Example usage scripts
- ✅ README with detailed usage instructions

## Next Steps for Publishing

### 1. Enable GitHub Pages

To publish the documentation:

1. Go to repository Settings → Pages
2. Under "Build and deployment":
   - Source: Deploy from a branch
   - Branch: `gh-pages` (will be created automatically by the workflow)
   - Folder: `/ (root)`
3. Save the settings

The documentation will be automatically built and deployed when you push to the main branch.

### 2. Publish to PyPI

To publish the package to PyPI:

1. Create a PyPI account at https://pypi.org/account/register/
2. Generate an API token:
   - Go to Account Settings → API tokens
   - Create a new token with scope for the project
3. Add the token to GitHub Secrets:
   - Go to repository Settings → Secrets and variables → Actions
   - Add a new secret named `PYPI_API_TOKEN`
   - Paste your PyPI API token as the value
4. Create a release:
   - Go to repository Releases → Create a new release
   - Create a new tag (e.g., `v0.1.0`)
   - Title: `v0.1.0 - Initial Release`
   - Description: Describe the features
   - Publish the release

The GitHub Action will automatically build and publish the package to PyPI.

### 3. Test the Package

Before publishing to PyPI, you can test the package locally:

```bash
# Build the package
python -m pip install build
python -m build

# Install the built package
pip install dist/medchem_molecule_filter-0.1.0-py3-none-any.whl

# Test it
python examples/basic_usage.py
```

### 4. Build Documentation Locally

To build and view the documentation locally:

```bash
# Install documentation dependencies
pip install -e ".[docs]"

# Build documentation
cd docs
make html

# Open docs/_build/html/index.html in a browser
```

## Package Features Summary

### Filters Available

1. **PAINSFilter**: Identifies Pan-Assay Interference Compounds
   - Quinones, catechols, rhodanines, and more
   - 20+ PAINS patterns

2. **ReactiveFilter**: Detects reactive functional groups
   - Electrophiles, Michael acceptors, peroxides
   - 25+ reactive patterns

3. **HeterocycleFilter**: Identifies heterocyclic scaffolds
   - Pyridine, furan, thiophene, indole, and more
   - 30+ heterocycle patterns
   - Can require or exclude heterocycles

4. **PropertyFilter**: Filters by calculated properties
   - LogP, TPSA, molecular weight
   - HBD/HBA counts
   - Rotatable bonds, aromatic rings

5. **FilterGroup**: Combine multiple filters in a pipeline

### Molecular Properties Calculator

Calculate properties for any molecule:
- LogP (partition coefficient)
- TPSA (Topological Polar Surface Area)
- Hydrogen Bond Donors/Acceptors
- Molecular Weight
- Rotatable Bonds
- Aromatic Rings

## Usage Examples

See `examples/basic_usage.py` for comprehensive examples including:
- Basic PAINS filtering
- Property-based filtering
- Combined filter pipelines
- Heterocycle filtering
- Property calculation

## Testing

Run tests with:
```bash
pytest tests/ -v
```

Run tests with coverage:
```bash
pytest tests/ --cov=medchem_filter --cov-report=html
```

## Documentation

Full documentation includes:
- Installation guide
- Quick start tutorial
- Complete API reference
- Filter rules reference with all SMARTS patterns
- Usage examples

## Dependencies

- Python >= 3.8
- RDKit >= 2022.9.1

## License

MIT License - See LICENSE file

## Contributing

See CONTRIBUTING.md for development guidelines.

## Support

For questions or issues:
- Open an issue on GitHub
- Check the documentation at https://frequencykg.github.io/medchem-molecule-filter

---

## Package Files Overview

```
medchem-molecule-filter/
├── medchem_filter/          # Main package
│   ├── __init__.py          # Package initialization
│   ├── filters.py           # Filter implementations
│   ├── properties.py        # Property calculators
│   └── data/                # Filter pattern data
│       ├── pains_patterns.py
│       ├── reactive_patterns.py
│       └── heterocycle_patterns.py
├── tests/                   # Test suite
│   ├── test_filters.py
│   └── test_properties.py
├── docs/                    # Sphinx documentation
│   ├── conf.py
│   ├── index.rst
│   ├── installation.rst
│   ├── quickstart.rst
│   ├── api.rst
│   └── filter_rules.rst
├── examples/                # Usage examples
│   └── basic_usage.py
├── .github/workflows/       # CI/CD workflows
│   ├── tests.yml
│   ├── docs.yml
│   └── publish.yml
├── pyproject.toml           # Package configuration
├── README.md                # Main documentation
├── LICENSE                  # MIT License
├── CONTRIBUTING.md          # Contribution guidelines
└── MANIFEST.in              # Package manifest
```
