# Change Log

All notable changes to this project will be documented in this file.

## Unreleased

### Added

- #203: Add the logical and physical dimensions `ldim` and `pdim` to the domain classes. `pdim` is the physical dimension of the mapping, or `ldim` without a mapping; for a boundary or an interface, they are the dimensions of its patches.
- #188: Support multipatch interfaces between faces with different normal axes and with any edge or face orientation, several interfaces between the same pair of patches, and periodic and self-interfaces.
- #188: Add `Domain.patches`, `Domain.interface_map`, `Domain.boundary_sides` and `Interface.axis_map` to inspect the topology of multipatch domains.
- #188: Add `Domain.vertices`, `Domain.edges` (in 3D), `Domain.indexed_interfaces`, `Domain.indexed_boundary_sides` and `Domain.patch_index()`, which describe the topology of a multipatch domain with patch indices, for discretization codes. The vertices and edges include unshared ones, record their patch incidences and whether they lie on the exterior boundary, and the edges the direction of each incidence. Their data classes are defined in the new module `sympde.topology.incidence`.
- #188: Add the module `sympde.utilities.plotting`, with `plot_domain()` for plotting 2D and 3D domains with optional topology annotations, and `format_topology()` and `print_topology()` for a text summary of their topology. `plot_domain()` and `print_topology()` are available in `sympde.api`.
- #188: Add the module `sympde.topology.multipatch_gallery`, which builds a collection of multipatch domains (previously in PSYDAC), and the analytical mappings `TransposedPolarMapping`, `TransposedTargetMapping` and `TransposedCzarnyMapping`, whose angle and radius coordinates are swapped. All analytical mappings, including `TorusMapping` and `TwistedTargetMapping`, are now available in `sympde.api`.
- #188: Document the multipatch topology in a new page of the user guide.
- #194: Add `sympde.api` as the supported flat public API. It contains the objects used by PSYDAC, IGA-Python, GeLaTo and the SymPDE documentation.
- #194: Add example scripts to the documentation, which are run by the documentation workflow.
- #194: Draw expression trees in the documentation with Graphviz, generated from the current code at build time.
- #106: Accept the connectivity of a `Domain` as an iterable of `Interface` objects, in addition to a `Connectivity` object.
- [DEVELOPER] Add a `CHANGELOG.md` file listing changes between different versions.
- [DEVELOPER] Add Python version 3.14 to the testing workflow.
- [DEVELOPER] Add separate documentation dependencies and a strict Sphinx workflow that checks pull requests and publishes the documentation to GitHub Pages.
- [DEVELOPER] Add an `AUTHORS` file listing the project maintainers and contributors.
- [DEVELOPER] Add the public helpers `new_basic()` and `restore_*_argument()` in `sympde.core.basic` for storing non-expression arguments in SymPy objects.

### Changed

- #188: Add an optional `patch` argument to `Domain.get_boundary()`, which raises a `ValueError` when several patches of a multipatch domain have a boundary side with the given axis and extremity, instead of returning one of them.
- #203: Replace the `dim` property and the `dim` argument of the domain classes (`Domain`, `InteriorDomain`, `NCube`, `Boundary`, `Interface`, ...) with `ldim`. This is a breaking change: use `domain.ldim` and `Domain(name, ldim=...)`. The `dim` argument of the mappings, a shorthand for `ldim = pdim = dim`, is kept.
- #203: Build the physical coordinates of a mapped domain from its physical dimension, so that a surface in 3D has the coordinates `[x, y, z]`, and map all of them in `LogicalExpr` and `TerminalExpr`.
- #188: Drop support for Python 3.9 and require Python 3.10 or newer.
- #143: Replace the `ornt` argument of `Interface` and `Boundary.join()` with a required `orientation` argument: `None` in 1D, `+1` or `-1` in 2D, and `(flag, sign1, sign2)` in 3D. Use `interface.orientation` instead of `interface.ornt`.
- #188: Accept the connectivity of `Domain.join()` only as `(minus, plus, orientation[, name])` tuples, and name repeated interfaces between the same pair of patches `A|B`, `A|B#2`, and so on. `interfaces=` is accepted as an alias of `connectivity=`. For mapped patches, the names are built from the logical patch names, and a mapped interface and its logical interface have the same name.
- #200: Store `ldim` and `pdim` instead of `dim` in the topology of domain files. This is a breaking change: `Domain.from_file()` rejects files with `dim`.
- #190: Convert the README to Markdown and expand the installation, testing, and documentation-build instructions.
- #194: Use explicit internal imports and load optional dependencies only when needed.
- #194: Remove the leading underscore from `partial_derivatives`, `logical_partial_derivatives` (`sympde.topology.derivatives`) and `split_test_function` (`sympde.expr.evaluation`), and add them to `sympde.api`.
- #194: List only public names in the `__all__` of every module.
- #106: Require SymPy 1.12 or newer and support maintained SymPy releases below 1.15.
- #106: Update SymPDE's symbolic objects for modern SymPy's structural arguments, traversal, substitution, and reconstruction requirements.
- #106: Order the parameters of `Domain` like its `args` (`name`, `interiors`, `boundaries`, `mapping`, `logical_domain`, `connectivity`) and allow passing them positionally. `dim` is now keyword-only.
- #106: Allow passing `mapping` and `logical_domain` of an `Interface` positionally, in this order.
- #106: Store only `axis` and `weight` in the `args` of the 1D forms (`Mass`, `Stiffness`, `Advection`, `AdvectionT`, `Bilaplacian`). Their name is now a class attribute.
- [DEVELOPER] Write the HDF5 files of the unit tests to `tmp_path`, so that the tests can run in parallel without the `xdist_group` marker.
- [DEVELOPER] Modernize the package metadata and consolidate the pytest and coverage configuration in `pyproject.toml`.
- [DEVELOPER] Modernize the Sphinx configuration, API-documentation generation, and API docstring markup.
- [DEVELOPER] Remove the leading underscore from the module-level names used by other SymPDE modules, and rename `_coeffs_registery` to `coefficient_types`.
- [DEVELOPER] Test the minimum supported SymPy version in a dedicated CI job.

### Fixed

- #188: Fix `Domain.get_interface()`, which always failed. It returns the interfaces that join two patches.
- #188: Raise a `ValueError` in `Domain.from_file()` when an interface references an unknown patch, instead of a `KeyError`.
- #188: Fix `split()`, which split a `Line` domain with wrong bounds and an invalid interface.
- #203: Check in the `Domain` constructor that all interiors have the same physical dimension, as `Domain.join()` does.
- #188: Reject interfaces in `Domain.join()` that join a face with itself or use a face that is already part of another interface.
- #188: Reject patches in `Domain.join()` of which only some have a mapping, instead of dropping the mapping of the domain.
- #140: Reject patches with different physical dimensions in `Domain.join()`, instead of creating an inconsistent domain.
- #144: Restore the interface names and exterior boundaries of a domain read with `Domain.from_file()` as they were exported, and read domains whose patches are `InteriorDomain` objects without a `dtype` instead of failing with an `AssertionError`.
- #200: Keep the physical dimension of the mappings in `Domain.export()` and `Domain.from_file()`.
- #166: Restore warning-free documentation builds, including bibliography and mathematical notation support.
- #106: Preserve symbol assumptions and free-symbol information with modern SymPy.
- #106: Keep same-named functions from different function spaces distinct during symbolic simplification.
- #106: Fix form-linearity checks involving integral sums with modern SymPy.
- #106: Raise a `ValueError` when a `Domain` is created without `interiors` or `dim`, instead of failing with a `TypeError` when only a connectivity is given.
- #191: Fix the rendering of the `Domain.join()` API documentation.

### Removed

- #188: Remove the `Domain.connectivity` property, which was the same as the new `Domain.interface_map`. Use `Domain.interface_map`, a read-only mapping from interface names to interfaces. The `connectivity` argument of `Domain` and `Domain.join()` is unchanged.
- #188: Remove `Domain.subdomains`, which was the same as the new `Domain.patches`. Use `Domain.patches`.
- #188: Remove `Domain.set_interfaces()`, the item assignment of `Connectivity` and `Connectivity.patches`. The interfaces of a domain are read-only.
- #188: Remove `Domain.corners`, `Domain.get_shared_corners()` and the classes `CornerBoundary` and `CornerInterface`. Use `Domain.vertices`, which contains all vertices of a domain with their patch indices, and tells whether they lie on the exterior boundary.
- #73: Remove the `Edge` class and the `Interface.axis` property. Use `Boundary`, and `interface.minus.axis` and `interface.plus.axis`.
- #188: Remove `plot_domain()`, `plot_2d()`, `plot_3d()`, `plot_2d_single_patch()` and `plot_3d_single_patch()` from `sympde.utilities.utils`. Use `plot_domain()` from `sympde.api` instead.
- #194: Remove the re-exports of public objects from the package and subpackage initializers. Replace imports such as `from sympde.topology import Square` with `from sympde.api import Square`. `from sympde import __version__` keeps working.
- #194: Remove the `sympde.old_sympy_utilities` module. Use `sympy.utilities.iterables.is_sequence` and Python's `metaclass` keyword instead of its `is_sequence` and `with_metaclass`.
- #194: Remove the outdated Jupyter notebooks and the Binder badge. The example scripts in the documentation replace them.
- #194: Remove the unused `_is_op_test_function` from `sympde.calculus.core`.
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
