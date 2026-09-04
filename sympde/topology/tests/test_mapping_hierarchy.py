# coding: utf-8
"""
Work-package 01: the new mapping base classes (SymbolicMapping, DefinedMapping,
StructuralMapping) introduced alongside the legacy hierarchy.

These checks lock the DefinedMapping interface and the abstractness of
DefinedMapping / StructuralMapping without exercising sympy object construction
(that is covered by later work-packages).
"""
import inspect

import pytest

from sympde.topology import SymbolicMapping, DefinedMapping, StructuralMapping
from sympde.topology import mapping as _mapping_mod
from sympde.topology.mapping import BasicCallableMapping


def test_new_classes_are_exported():
    # the star-import from the package must expose the very objects defined in
    # the module, and the module must advertise them in __all__
    assert SymbolicMapping   is _mapping_mod.SymbolicMapping
    assert DefinedMapping    is _mapping_mod.DefinedMapping
    assert StructuralMapping is _mapping_mod.StructuralMapping
    for name in ('SymbolicMapping', 'DefinedMapping', 'StructuralMapping'):
        assert name in _mapping_mod.__all__


def test_hierarchy_links():
    assert issubclass(DefinedMapping,    SymbolicMapping)
    assert issubclass(StructuralMapping, SymbolicMapping)


def test_defined_mapping_interface_covers_basic_callable_mapping():
    # DefinedMapping is meant to supersede BasicCallableMapping: every method of
    # the old interface must still be required by the new one
    assert set(BasicCallableMapping.__abstractmethods__).issubset(
        DefinedMapping.__abstractmethods__)


def test_symbolic_mapping_is_concrete():
    assert not inspect.isabstract(SymbolicMapping)


@pytest.mark.parametrize('cls', [DefinedMapping, StructuralMapping])
def test_abstract_classes_reject_instantiation(cls):
    assert inspect.isabstract(cls)
    with pytest.raises(TypeError, match='abstract'):
        cls('F')


def test_concrete_subclasses_lose_abstractness_without_instantiation():
    class ConcreteDefined(DefinedMapping):
        def __call__(self, *eta):    return eta
        def jacobian(self, *eta):    return None
        def jacobian_inv(self, *eta): return None
        def metric(self, *eta):      return None
        def metric_det(self, *eta):  return None
        ldim = 2
        pdim = 2

    class ConcreteStructural(StructuralMapping):
        ldim = 2
        pdim = 2

    assert not inspect.isabstract(ConcreteDefined)
    assert not inspect.isabstract(ConcreteStructural)
    assert ConcreteDefined.__abstractmethods__    == frozenset()
    assert ConcreteStructural.__abstractmethods__ == frozenset()
