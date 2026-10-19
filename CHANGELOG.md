# Change Log

All notable changes to this project will be documented in this file.

## Unreleased

### Added

- #194: Add `sympde.api` as the supported flat public API. It contains the objects used by PSYDAC, IGA-Python, GeLaTo and the SymPDE documentation.
- #194: Add example scripts to the documentation, which are run by the documentation workflow.
- #106: Accept the connectivity of a `Domain` as an iterable of `Interface` objects, in addition to a `Connectivity` object.
- [DEVELOPER] Add a `CHANGELOG.md` file listing changes between different versions.
- [DEVELOPER] Add Python version 3.14 to the testing workflow.
- [DEVELOPER] Add separate documentation dependencies and a strict Sphinx workflow that checks pull requests and publishes the documentation to GitHub Pages.
- [DEVELOPER] Add an `AUTHORS` file listing the project maintainers and contributors.
- [DEVELOPER] #106: Add the public helpers `new_basic()` and `restore_*_argument()` in `sympde.core.basic` for storing non-expression arguments in SymPy objects.

### Changed

- Convert the README to Markdown and expand the installation, testing, and documentation-build instructions.
- #194: Use explicit internal imports and load optional dependencies only when needed.
- #106: Require SymPy 1.12 or newer and support maintained SymPy releases below 1.15.
- #106: Update SymPDE's symbolic objects for modern SymPy's structural arguments, traversal, substitution, and reconstruction requirements.
- #106: Order the parameters of `Domain` like its `args` (`name`, `interiors`, `boundaries`, `mapping`, `logical_domain`, `connectivity`) and allow passing them positionally. `dim` is now keyword-only.
- #106: Allow passing `ornt`, `mapping` and `logical_domain` of an `Interface` positionally, in this order.
- #106: Store only `axis` and `weight` in the `args` of the 1D forms (`Mass`, `Stiffness`, `Advection`, `AdvectionT`, `Bilaplacian`). Their name is now a class attribute.
- [DEVELOPER] Modernize the package metadata and consolidate the pytest and coverage configuration in `pyproject.toml`.
- [DEVELOPER] Modernize the Sphinx configuration, API-documentation generation, and API docstring markup.
- [DEVELOPER] Test the minimum supported SymPy version in a dedicated CI job.

### Fixed

- #140: Fix the rendering of the `Domain.join()` API documentation.
- #166: Restore warning-free documentation builds, including bibliography and mathematical notation support.
- #106: Preserve symbol assumptions and free-symbol information with modern SymPy.
- #106: Keep same-named functions from different function spaces distinct during symbolic simplification.
- #106: Fix form-linearity checks involving integral sums with modern SymPy.
- #106: Raise a `ValueError` when a `Domain` is created without `interiors` or `dim`, instead of failing with a `TypeError` when only a connectivity is given.

### Removed

- #194: Remove the re-exports of public objects from the package and subpackage initializers. Replace imports such as `from sympde.topology import Square` with `from sympde.api import Square`. `from sympde import __version__` keeps working.
- #194: Remove the `sympde.old_sympy_utilities` module. Use `sympy.utilities.iterables.is_sequence` and Python's `metaclass` keyword instead of its `is_sequence` and `with_metaclass`.
- [DEVELOPER] Remove obsolete CI configuration, the standalone pytest configuration, and the legacy test runner script.

### Deprecated

## [0.20.0] - 2026-09-07

### Added

- #185: Store interface orientations when exporting domain connectivity, allowing them to be restored when loading a domain.
- #185: Add consistency checks when constructing mapped domains.

### Changed

- #185: Require an explicit orientation when constructing multipatch interfaces, including through `Interface`, `Domain.join()`, and `Domain.from_file()`. This is a breaking change for callers and domain files that omitted orientation information.
- #185: Clean up expression evaluation and expand its API documentation.
- [DEVELOPER] Run tests in parallel with `pytest-xdist`, while keeping the HDF5 gallery tests serial.

## [0.19.3] - 2026-07-17

### Fixed

- #187: Avoid redundant symbolic expansion in `is_linear_expression()`, fixing severe slowdowns for some expressions.

### Changed

- #187: Refactor and document `is_linear_expression()` and remove its unused `integral` argument.

## [0.19.2] - 2025-04-02

### Changed

- #177: Drop support for Python 3.8 and require Python 3.9 or newer.
- #177: Add testing and installation support for Python 3.13.

## [0.19.1] - 2025-03-05

### Added

- #173: Add `Domain.subdomains` and `Domain.mappings`, which consistently return tuples for both single-patch and multipatch domains.
- #173: Allow `Domain.join()` to be called with a single patch.
- [DEVELOPER] Add terminal-expression tests for exact Navier–Stokes solutions.

### Fixed

- #170: Mark domain coordinates as real SymPy symbols, avoiding code-generation problems for complex expressions.
- #174: Fix README markup so the project description renders correctly on PyPI.

### Changed

- #172: Expand the source-installation and virtual-environment instructions.

## [0.19.0] - 2024-08-14

### Added

- #154: Add support and continuous-integration testing for Python 3.12.
- #158: Accept NumPy scalar values as parameters of analytical mappings.
- #155: Allow `Domain.join()` connectivity entries to reference patch objects directly instead of only their indices.

### Changed

- #155: Document `Domain.join()` and make joined-domain construction less dependent on patch ordering.
- #164: Document how SymPDE branches can be tested against the Psydac test suite.
- #165: Replace the obsolete build badge with the GitHub Actions CI badge.

### Fixed

- #167: Fix README markup rejected by the PyPI uploader.

### Removed

- #165: Remove the obsolete Travis CI configuration.

[0.20.0]: https://github.com/pyccel/sympde/compare/v0.19.3...v0.20.0
[0.19.3]: https://github.com/pyccel/sympde/compare/v0.19.2...v0.19.3
[0.19.2]: https://github.com/pyccel/sympde/compare/v0.19.1...v0.19.2
[0.19.1]: https://github.com/pyccel/sympde/compare/v0.19.0...v0.19.1
[0.19.0]: https://github.com/pyccel/sympde/compare/v0.18.3...v0.19.0
