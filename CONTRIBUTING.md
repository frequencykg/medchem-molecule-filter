# Contributing to medchem-molecule-filter

Thank you for your interest in contributing to medchem-molecule-filter!

## Development Setup

1. Clone the repository:
```bash
git clone https://github.com/frequencykg/medchem-molecule-filter.git
cd medchem-molecule-filter
```

2. Install development dependencies:
```bash
pip install -e ".[dev,docs]"
```

## Running Tests

Run all tests:
```bash
pytest tests/ -v
```

Run tests with coverage:
```bash
pytest tests/ --cov=medchem_filter --cov-report=html
```

## Code Style

We use Black for code formatting and flake8 for linting:

```bash
# Format code
black medchem_filter/ tests/

# Check linting
flake8 medchem_filter/ tests/
```

## Building Documentation

Build the documentation locally:

```bash
cd docs
make html
```

The built documentation will be in `docs/_build/html/`.

## Adding New Filters

When adding new filter patterns:

1. Add patterns to the appropriate file in `medchem_filter/data/`
2. Update the documentation in `docs/filter_rules.rst`
3. Add tests in `tests/`
4. Update the README with examples

## Submitting Changes

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Submit a pull request

## Questions?

Feel free to open an issue for discussion before making large changes.
