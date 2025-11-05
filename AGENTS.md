# AGENTS.md - Guidelines for AI Coding Agents

## Build/Test Commands
- **Install**: `pip install -e .`
- **Test all**: `pytest`
- **Test single**: `pytest path/to/test_file.py::TestClass::test_method`
- **Lint**: `flake8`
- **Format**: `black .`
- **Type check**: `mypy`

## Code Style Guidelines

### Imports
- Group imports: standard library first, then third-party, then local
- Use absolute imports for local modules
- Example: `from typing import Optional`

### Formatting
- Line length: 100 characters (Black enforced)
- Use Black for consistent formatting

### Types
- Use type hints for all function parameters and return values
- Use `Optional[T]` for nullable types
- Use `Union[T1, T2]` for multiple possible types

### Naming Conventions
- Functions/variables: `snake_case`
- Classes: `PascalCase`
- Constants: `UPPER_CASE`
- Private members: `_prefix`

### Classes
- Use `__slots__` for memory efficiency when appropriate
- Initialize all attributes in `__init__`

### Error Handling
- Use try/except blocks with specific exception types
- Raise descriptive exceptions with context

### Performance
- Use `array.array` for large numeric datasets
- Use `sys.intern()` for repeated strings
- Prefer sets for O(1) membership tests

### Python Version
- Target Python 3.9+ compatibility