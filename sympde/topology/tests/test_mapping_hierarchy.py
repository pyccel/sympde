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


# -- work-package 02a: the symbolic `jacobian` property moved to `jacobian_symbol`

def test_jacobian_symbol_holds_the_symbolic_jacobian():
    from sympde.topology import IdentityMapping
    F = IdentityMapping('F', dim=2)
    assert F.jacobian_symbol is F._jacobian
    assert type(F.jacobian_symbol).__name__ == 'JacobianSymbol'


def test_legacy_jacobian_property_is_a_deprecated_alias():
    from sympde.topology import Mapping
    # On a bare (undefined) Mapping, `.jacobian` is still the deprecated
    # symbolic property. On an AnalyticMapping it is shadowed by the numeric
    # point-evaluation method (see test_analytic_mapping_is_directly_point_evaluable).
    F = Mapping('F', dim=2)
    with pytest.warns(DeprecationWarning, match='jacobian_symbol'):
        legacy = F.jacobian
    assert legacy is F.jacobian_symbol


# -- work-package 02b: AnalyticMapping is a concrete, point-evaluable DefinedMapping

def test_defined_mapping_inherits_basic_callable_mapping():
    # D5: the point-eval interface is now stated once, on BasicCallableMapping
    assert issubclass(DefinedMapping, BasicCallableMapping)
    assert DefinedMapping.__abstractmethods__ == BasicCallableMapping.__abstractmethods__


def test_analytic_mapping_hierarchy():
    from sympde.topology import AnalyticMapping, Mapping
    assert issubclass(AnalyticMapping, Mapping)
    assert issubclass(AnalyticMapping, DefinedMapping)
    assert issubclass(AnalyticMapping, BasicCallableMapping)
    assert not inspect.isabstract(AnalyticMapping)


def test_analytical_gallery_reparented_onto_analytic_mapping():
    from sympde.topology import AnalyticMapping, Mapping
    from sympde.topology import (IdentityMapping, AffineMapping, PolarMapping,
                                 TargetMapping, CzarnyMapping, CollelaMapping2D,
                                 TorusMapping)
    for cls in (IdentityMapping, AffineMapping, PolarMapping, TargetMapping,
                CzarnyMapping, CollelaMapping2D, TorusMapping):
        assert issubclass(cls, AnalyticMapping)
        assert issubclass(cls, Mapping)          # isinstance(_, Mapping) still holds


def test_analytic_mapping_is_directly_point_evaluable():
    from sympde.topology import IdentityMapping
    F = IdentityMapping('F', dim=2)
    cm = F.get_callable_mapping()
    # point evaluation goes through the same callable mapping as before
    assert F(0.3, 0.4) == cm(0.3, 0.4)
    J = F.jacobian(0.3, 0.4)
    assert type(J).__name__ != 'JacobianSymbol'      # numeric, not symbolic
    assert list(map(list, J)) == list(map(list, cm.jacobian(0.3, 0.4)))
    assert F.metric_det(0.3, 0.4) == cm.metric_det(0.3, 0.4)


def test_analytic_mapping_still_callable_on_a_domain():
    from sympde.topology import IdentityMapping, Square
    F = IdentityMapping('F', dim=2)
    mapped = F(Square('D'))
    assert type(mapped).__name__ in ('Domain', 'MappedDomain')


def test_analytic_mapping_without_expressions_rejects_point_call():
    from sympde.topology import AnalyticMapping
    F = AnalyticMapping('F', dim=2)          # no _expressions
    with pytest.raises(ValueError, match='analytical expressions'):
        F(0.1, 0.2)


# -- work-package 02b: post-review fixes

def test_analytic_mapping_jacobian_still_shadows_symbolic_alias_but_jacobian_symbol_is_unaffected():
    # code-review fix #1, corrected: AnalyticMapping.jacobian (the numeric
    # method) shadows the deprecated symbolic Mapping.jacobian property with no
    # warning. But investigating the flagged psydac call sites
    # (compute_boundary_jacobian / compute_normal_vector) showed
    # SymbolicExpr() has never supported translating a JacobianSymbol -- even
    # on a bare, pre-refactor Mapping, SymbolicExpr(mapping.jacobian) already
    # raised NotImplementedError. So this is not a regression: those two
    # psydac functions were already non-functional dead code (zero callers in
    # the workspace) before this refactor. What we *can* guarantee here is
    # that jacobian_symbol keeps returning the same JacobianSymbol object
    # regardless of shadowing, so the eventual WP05 rename is behaviour
    # preserving.
    from sympde.topology import PolarMapping
    F = PolarMapping('F', rmin=0, rmax=1, c1=0, c2=0)
    assert type(F.jacobian_symbol).__name__ == 'JacobianSymbol'
    assert type(F.jacobian) is not type(F.jacobian_symbol)   # method vs symbol


def test_interface_mapping_copy_preserves_analytic_mapping_type():
    # code-review fix #2: Mapping.copy() used to hardcode Mapping(...), so
    # InterfaceMapping (which copies both legs) silently downgraded an
    # AnalyticMapping subclass to a plain, non-point-evaluable Mapping.
    import numpy as np
    from sympde.topology import IdentityMapping, InterfaceMapping
    itf = InterfaceMapping(IdentityMapping('F1', dim=2), IdentityMapping('F2', dim=2))
    assert type(itf.minus) is IdentityMapping
    assert type(itf.plus)  is IdentityMapping
    expected = itf.minus.get_callable_mapping().jacobian(0.3, 0.4)
    assert np.array_equal(itf.minus.jacobian(0.3, 0.4), expected)
    assert itf.minus(0.3, 0.4) == (0.3, 0.4)


def test_copy_preserves_user_set_callable_mapping():
    # post-commit finding A: copy() used to write the name-mangled
    # `__callable_map` instead of `_callable_map`, silently dropping a
    # user-supplied callable mapping (set via set_callable_mapping) on copy --
    # observable e.g. through InterfaceMapping, which always copies its legs.
    from sympde.topology import AnalyticMapping
    from sympde.topology.mapping import BasicCallableMapping

    class Custom(BasicCallableMapping):
        def __call__(self, *eta):     return eta
        def jacobian(self, *eta):     return None
        def jacobian_inv(self, *eta): return None
        def metric(self, *eta):       return None
        def metric_det(self, *eta):   return None
        ldim = 2
        pdim = 2

    F  = AnalyticMapping('F', dim=2)   # no _expressions
    cm = Custom()
    F.set_callable_mapping(cm)
    assert F.copy()._callable_map is cm


# -- work-package 03: InverseMapping / InterfaceMapping / MultiPatchMapping
# are now also StructuralMapping subclasses (additive: they keep Mapping too)

def test_structural_subclasses_hierarchy():
    from sympde.topology import (StructuralMapping, InverseMapping,
                                 InterfaceMapping, MultiPatchMapping, Mapping)
    for cls in (InverseMapping, InterfaceMapping, MultiPatchMapping):
        assert issubclass(cls, StructuralMapping)
        assert issubclass(cls, Mapping)              # no regression
        assert not inspect.isabstract(cls)
        # own ldim/pdim must win over StructuralMapping's abstract ones
        assert cls.ldim is not StructuralMapping.ldim
        assert cls.pdim is not StructuralMapping.pdim


def test_interface_mapping_rejects_point_call_but_stays_domain_callable():
    from sympde.topology import IdentityMapping, InterfaceMapping, Square
    itf = InterfaceMapping(IdentityMapping('F1', dim=2), IdentityMapping('F2', dim=2))
    assert itf.ldim == 2
    with pytest.raises(TypeError, match='StructuralMapping'):
        itf(0.3, 0.4)
    mapped = itf(Square('D'))
    assert type(mapped).__name__ in ('Domain', 'MappedDomain')


def test_multi_patch_mapping_rejects_point_call():
    from sympde.topology import MultiPatchMapping, IdentityMapping
    mp = MultiPatchMapping({'p1': IdentityMapping('F1', dim=2),
                            'p2': IdentityMapping('F2', dim=2)})
    assert mp.ldim == 2
    with pytest.raises(TypeError, match='StructuralMapping'):
        mp(0.3, 0.4)
    # NOTE: mp(some_domain) (the domain-call path) is not exercised here: it
    # hits a pre-existing bug independent of WP03 -- MultiPatchMapping.__new__
    # uses Basic.__new__(cls, dic) and never sets self._name, so any domain
    # call that reaches Mapping.name raises AttributeError. Reproduced on the
    # pre-WP03 baseline too (git stash), so it is not a regression from this
    # slice. Flagged for a future fix; out of scope here.


def test_structural_mapping_still_abstract():
    # WP01 invariant must not regress
    from sympde.topology import StructuralMapping
    assert inspect.isabstract(StructuralMapping)
    with pytest.raises(TypeError, match='abstract'):
        StructuralMapping('F')


# -- post-review fixes (found by /code-review on the WP03 diff)

def test_structural_mapping_call_accepts_domain_keyword():
    # StructuralMapping.__call__(self, *args) initially dropped the `domain`
    # keyword that Mapping.__call__(self, domain) supported -- restored via
    # explicit positional/keyword normalization before delegating.
    from sympde.topology import IdentityMapping, InterfaceMapping, Square
    itf = InterfaceMapping(IdentityMapping('F1', dim=2), IdentityMapping('F2', dim=2))
    positional = itf(Square('D'))
    keyword    = itf(domain=Square('D'))
    assert type(positional) is type(keyword)


def test_structural_mapping_call_delegates_via_super():
    # DRY fix: the domain-call branch now delegates to Mapping.__call__ (via
    # super()) instead of duplicating its body.
    from sympde.topology import StructuralMapping, Mapping
    import inspect as _inspect
    src = _inspect.getsource(StructuralMapping.__call__)
    assert 'super().__call__' in src


def test_ldim_pdim_alias_mapping_descriptor():
    # DRY fix: InverseMapping/InterfaceMapping alias Mapping's ldim/pdim
    # property objects rather than re-typing their bodies.
    from sympde.topology import InverseMapping, InterfaceMapping, Mapping
    assert InverseMapping.ldim   is Mapping.ldim
    assert InverseMapping.pdim   is Mapping.pdim
    assert InterfaceMapping.ldim is Mapping.ldim
    assert InterfaceMapping.pdim is Mapping.pdim


def test_interface_mapping_copy_is_broken_pre_existing():
    # post-review finding: InterfaceMapping inherits Mapping.copy(), which
    # calls type(self)(self.name, ldim=..., pdim=..., evaluate=False) -- but
    # InterfaceMapping.__new__(cls, minus, plus) doesn't accept those kwargs.
    # Confirmed pre-existing (reproduces on the pre-WP03 baseline too, and
    # traces back to the WP02b post-review copy() fix that generalized
    # Mapping(...) to type(self)(...)); documenting the current behaviour
    # rather than silently fixing it -- .copy() is never called on an
    # InterfaceMapping instance itself anywhere in sympde or psydac today.
    from sympde.topology import IdentityMapping, InterfaceMapping
    itf = InterfaceMapping(IdentityMapping('F1', dim=2), IdentityMapping('F2', dim=2))
    with pytest.raises(TypeError, match='ldim'):
        itf.copy()


def test_mapped_multipatch_domain_with_interface_is_broken_pre_existing():
    # post-review finding: MappedDomain.__new__ rebuilds each patch Interface
    # as Interface(e.name, mapping(e.minus), mapping(e.plus)) without
    # forwarding e.ornt (or mapping/logical_domain). Since commit 9783335
    # ("Rigorously check arguments of Interface constructor") made
    # Interface.__new__ validate ornt strictly, this now raises instead of
    # silently building an under-specified Interface. Confirmed pre-existing
    # (reproduces identically on the pre-WP03 baseline via git stash) --
    # unrelated to StructuralMapping/Mapping base-list changes, but it lives
    # in mapping.py's MappedDomain.__new__, a function this file's tests
    # exercise elsewhere. Documenting rather than fixing: the fix belongs to
    # whoever owns the Interface/ornt validation, not this work-package.
    from sympde.topology import Domain, Square, IdentityMapping
    patches = [Square('D1'), Square('D2')]
    domain  = Domain.join(patches, [((0, 0, 1), (1, 0, -1), 1)], 'domain')
    F = IdentityMapping('F', dim=2)
    with pytest.raises(AssertionError, match='ornt'):
        F(domain)
