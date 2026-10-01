# Change Log

All notable changes to this project will be documented in this file.

## Unreleased

### Added

### Changed

### Deprecated

### Deprecated

### Removed

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
