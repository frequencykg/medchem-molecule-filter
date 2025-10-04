# GitHub Copilot Instructions for medchem-molecule-filter

## Project Overview

medchem-molecule-filter is a Python library for filtering small molecules using medicinal chemistry rules. It provides filters for PAINS, reactive groups, heterocycles, and molecular properties.

## Development Guidelines

### Code Quality Standards

**Always ensure:**
- All tests pass before committing (`pytest tests/ -v`)
- Code follows PEP 8 style guidelines
- New functions have docstrings with type hints
- New features include corresponding tests
- Documentation is updated for API changes

### Testing Requirements

**For every code change:**
1. Run the full test suite: `pytest tests/ -v`
2. Check test coverage: `pytest tests/ --cov=medchem_filter --cov-report=term-missing`
3. Aim for >90% code coverage
4. Add tests for:
   - New filter patterns
   - New properties
   - Edge cases (invalid molecules, empty inputs)
   - Error conditions

**Test file locations:**
- `tests/test_filters.py` - Filter functionality tests
- `tests/test_properties.py` - Property calculation tests

### Code Style and Formatting

**Before committing:**
```bash
# Format code with Black
black medchem_filter/ tests/ --line-length 100

# Check linting with flake8
flake8 medchem_filter/ tests/ --max-line-length 100 --extend-ignore E203,W503

# Type checking (optional but recommended)
mypy medchem_filter/ --ignore-missing-imports
```

### Documentation Standards

**Always update documentation when:**
- Adding new filters or filter patterns
- Adding new molecular properties
- Changing public API
- Adding new features

**Documentation locations:**
- `README.md` - Main user documentation with examples
- `docs/` - Sphinx documentation (API reference, tutorials)
- Docstrings - In-code documentation for all public functions/classes

**Docstring format (Google style):**
```python
def example_function(mol: Chem.Mol, threshold: float = 0.5) -> bool:
    """Short one-line description.
    
    Longer description if needed, explaining what the function does,
    when to use it, and any important considerations.
    
    Args:
        mol: RDKit molecule object to process
        threshold: Threshold value for filtering (default: 0.5)
        
    Returns:
        True if molecule passes filter, False otherwise
        
    Raises:
        ValueError: If mol is None or invalid
        
    Example:
        >>> mol = Chem.MolFromSmiles("CCO")
        >>> result = example_function(mol)
        >>> print(result)
        True
    """
```

## Development Workflow

### Setting Up Development Environment

```bash
# Clone repository
git clone https://github.com/frequencykg/medchem-molecule-filter.git
cd medchem-molecule-filter

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode with all dependencies
pip install -e ".[dev,docs]"
```

### Making Changes

1. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**
   - Write code following the style guide
   - Add comprehensive tests
   - Update documentation

3. **Test your changes**
   ```bash
   # Run all tests
   pytest tests/ -v
   
   # Check coverage
   pytest tests/ --cov=medchem_filter --cov-report=html
   
   # Format code
   black medchem_filter/ tests/
   
   # Check linting
   flake8 medchem_filter/ tests/
   ```

4. **Commit and push**
   ```bash
   git add .
   git commit -m "feat: add description of feature"
   git push origin feature/your-feature-name
   ```

### Commit Message Format

Follow conventional commits:
- `feat:` - New feature
- `fix:` - Bug fix
- `docs:` - Documentation changes
- `test:` - Adding or updating tests
- `refactor:` - Code refactoring
- `style:` - Code style changes (formatting)
- `chore:` - Maintenance tasks

Examples:
- `feat: add new PAINS filter pattern for quinones`
- `fix: correct TPSA calculation for aromatic rings`
- `docs: update README with new filter examples`
- `test: add edge cases for PropertyFilter`

## Package Structure

```
medchem-molecule-filter/
├── medchem_filter/          # Main package
│   ├── __init__.py          # Package initialization, exports
│   ├── filters.py           # Filter classes (PAINS, Reactive, etc.)
│   ├── properties.py        # Molecular property calculators
│   └── data/                # Filter pattern data
│       ├── pains_patterns.py
│       ├── reactive_patterns.py
│       └── heterocycle_patterns.py
├── tests/                   # Test suite
│   ├── test_filters.py
│   └── test_properties.py
├── docs/                    # Sphinx documentation
├── examples/                # Example scripts
├── .github/                 # GitHub configurations
│   ├── workflows/           # CI/CD workflows
│   └── copilot-instructions.md  # This file
└── pyproject.toml          # Package configuration
```

## Adding New Features

### Adding a New Filter Pattern

1. **Add pattern to data file**
   ```python
   # In medchem_filter/data/pains_patterns.py
   PAINS_PATTERNS = {
       # ... existing patterns ...
       "new_pattern": "SMARTS_pattern_here",  # Add description
   }
   ```

2. **Add test**
   ```python
   # In tests/test_filters.py
   def test_new_pattern_detection(self):
       """Test detection of new pattern."""
       filter_obj = PAINSFilter()
       mol = Chem.MolFromSmiles("test_smiles")
       passes, matched = filter_obj.check_molecule(mol)
       assert not passes  # Should fail if pattern matches
       assert "new_pattern" in matched
   ```

3. **Update documentation**
   - Add pattern to `docs/filter_rules.rst`
   - Add example to README if significant

### Adding a New Property Calculator

1. **Add method to MolecularProperties class**
   ```python
   # In medchem_filter/properties.py
   @staticmethod
   def calculate_new_property(mol: Chem.Mol) -> float:
       """Calculate new molecular property.
       
       Args:
           mol: RDKit molecule object
           
       Returns:
           Calculated property value
       """
       # Implementation
       return value
   ```

2. **Add to calculate_all_properties**
   ```python
   def calculate_all_properties(mol: Chem.Mol) -> dict:
       return {
           # ... existing properties ...
           'new_property': MolecularProperties.calculate_new_property(mol),
       }
   ```

3. **Add tests**
   ```python
   # In tests/test_properties.py
   def test_calculate_new_property(self):
       """Test new property calculation."""
       mol = Chem.MolFromSmiles("CCO")
       result = MolecularProperties.calculate_new_property(mol)
       assert isinstance(result, float)
       assert result > 0  # Or appropriate validation
   ```

## Common Tasks

### Running Tests Locally
```bash
# All tests with verbose output
pytest tests/ -v

# Specific test file
pytest tests/test_filters.py -v

# Specific test
pytest tests/test_filters.py::TestPAINSFilter::test_quinone_detection -v

# With coverage report
pytest tests/ --cov=medchem_filter --cov-report=html
# Open htmlcov/index.html to view coverage

# Run tests matching a pattern
pytest tests/ -k "property" -v
```

### Building Documentation
```bash
cd docs
make html
# Open docs/_build/html/index.html

# Clean and rebuild
make clean
make html
```

### Building Package
```bash
# Install build tool
pip install build

# Build package
python -m build

# Check distribution
twine check dist/*
```

## Debugging Tips

### Common Issues

1. **Import errors**
   - Ensure package is installed: `pip install -e .`
   - Check PYTHONPATH includes project root

2. **Test failures**
   - Run with `-v` flag for verbose output
   - Use `pytest --pdb` to drop into debugger on failure
   - Check that RDKit is installed correctly

3. **SMARTS pattern errors**
   - Validate SMARTS with RDKit: `Chem.MolFromSmarts("pattern")`
   - Use RDKit documentation for SMARTS syntax

4. **Property calculation errors**
   - Ensure molecule is valid: `mol is not None`
   - Some properties require explicit hydrogens: `Chem.AddHs(mol)`

### Using Copilot Effectively

**Good prompts:**
- "Add a test for filtering molecules with quinone patterns"
- "Create a new filter for metal-containing compounds"
- "Add documentation for the PropertyFilter class"

**Include context:**
- Mention specific files: "In test_filters.py, add a test for..."
- Reference existing patterns: "Similar to PAINSFilter, create..."
- Specify requirements: "The test should check both passing and failing cases"

## CI/CD Pipeline

### GitHub Actions Workflows

**tests.yml** - Runs on every push and PR
- Tests on Python 3.8, 3.9, 3.10, 3.11, 3.12
- Runs linting (flake8)
- Runs tests with coverage
- Uploads coverage to Codecov

**docs.yml** - Runs on push to main
- Builds Sphinx documentation
- Deploys to GitHub Pages

**publish.yml** - Runs on release
- Builds package
- Publishes to PyPI

### Pre-commit Checks (Local)

Before pushing, ensure:
```bash
# Format code
black medchem_filter/ tests/

# Run linting
flake8 medchem_filter/ tests/

# Run all tests
pytest tests/ -v

# Check types (optional)
mypy medchem_filter/ --ignore-missing-imports
```

## Dependencies

### Core Dependencies
- **rdkit** (>=2022.9.1) - Chemistry toolkit for molecular operations

### Development Dependencies
- **pytest** (>=7.0.0) - Testing framework
- **pytest-cov** (>=4.0.0) - Coverage reporting
- **black** (>=22.0.0) - Code formatting
- **flake8** (>=5.0.0) - Linting
- **mypy** (>=0.990) - Type checking

### Documentation Dependencies
- **sphinx** (>=5.0.0) - Documentation generator
- **sphinx-rtd-theme** (>=1.0.0) - ReadTheDocs theme
- **sphinx-autodoc-typehints** (>=1.19.0) - Type hint support

## Best Practices

### Performance
- Compile SMARTS patterns once (done in `__init__`)
- Use batch processing for multiple molecules
- Cache property calculations when possible

### Error Handling
- Always check if molecule is None
- Provide meaningful error messages
- Use appropriate exception types (ValueError, TypeError)

### Testing
- Test both positive and negative cases
- Use realistic test molecules
- Test edge cases (empty inputs, invalid SMILES)
- Mock external dependencies when needed

### Documentation
- Keep README examples simple and working
- Include expected output in examples
- Document all public APIs
- Add comments for complex SMARTS patterns

## Resources

- [RDKit Documentation](https://www.rdkit.org/docs/)
- [SMARTS Tutorial](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
- [Lipinski's Rule of Five](https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five)
- [PAINS Filters](https://en.wikipedia.org/wiki/Pan-assay_interference_compounds)

## Contact

For questions or issues:
- Open an issue on GitHub
- Check existing documentation
- Review test files for examples
