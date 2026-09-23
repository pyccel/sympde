# coding: utf-8
"""
Work-package 01: the new mapping base classes (SymbolicMapping, DefinedMapping,
StructuralMapping) introduced alongside the legacy hierarchy.

These checks lock the DefinedMapping interface and the abstractness of
DefinedMapping / StructuralMapping without exercising sympy object construction
(that is covered by later work-packages).
"""
import inspect
import warnings

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
    # DefinedMapping supersedes BasicCallableMapping: every name of the old
    # interface is still exposed. 06d-4c: ldim / pdim / __call__ are now
    # provided concretely by SymbolicMapping, so only the numeric point-eval
    # methods stay abstract until a concrete subclass implements them.
    for name in BasicCallableMapping.__abstractmethods__:
        assert hasattr(DefinedMapping, name)
    assert {'jacobian', 'jacobian_inv', 'metric', 'metric_det'}.issubset(
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
    # `.jacobian` (the deprecated symbolic property) lives on `Mapping` only --
    # `SymbolicMapping` deliberately does not carry it. Since WP D8,
    # `AnalyticMapping` no longer inherits `Mapping`, so nothing is "shadowed"
    # any more: `.jacobian` is simply confined to `Mapping` and its legacy
    # subclasses, while `AnalyticMapping.jacobian` is a distinct, unrelated
    # numeric method (see test_analytic_mapping_is_directly_point_evaluable).
    with pytest.warns(DeprecationWarning, match='SymbolicMapping'):   # 06d-4c: ctor
        F = Mapping('F', dim=2)
    with pytest.warns(DeprecationWarning, match='jacobian_symbol'):
        legacy = F.jacobian
    assert legacy is F.jacobian_symbol


# -- work-package 02b: AnalyticMapping is a concrete, point-evaluable DefinedMapping

def test_defined_mapping_inherits_basic_callable_mapping():
    # D5: the point-eval interface is now stated once, on BasicCallableMapping.
    # 06d-4c: SymbolicMapping provides ldim / pdim / __call__ concretely, so
    # DefinedMapping only still-abstracts the numeric point-eval methods.
    assert issubclass(DefinedMapping, BasicCallableMapping)
    assert DefinedMapping.__abstractmethods__ == frozenset(
        {'jacobian', 'jacobian_inv', 'metric', 'metric_det'})


def test_analytic_mapping_hierarchy():
    from sympde.topology import AnalyticMapping, Mapping
    # WP D8: AnalyticMapping no longer inherits the deprecated Mapping.
    assert not issubclass(AnalyticMapping, Mapping)
    assert issubclass(AnalyticMapping, DefinedMapping)
    assert issubclass(AnalyticMapping, BasicCallableMapping)
    assert not inspect.isabstract(AnalyticMapping)


def test_analytical_gallery_reparented_onto_analytic_mapping():
    from sympde.topology import AnalyticMapping, Mapping, SymbolicMapping
    from sympde.topology import (IdentityMapping, AffineMapping, PolarMapping,
                                 TargetMapping, CzarnyMapping, CollelaMapping2D,
                                 TorusMapping)
    for cls in (IdentityMapping, AffineMapping, PolarMapping, TargetMapping,
                CzarnyMapping, CollelaMapping2D, TorusMapping):
        assert issubclass(cls, AnalyticMapping)
        # WP D8: no longer under the deprecated Mapping; use SymbolicMapping
        # as the "any symbolic mapping" check (same as 06d-4a for structural).
        assert not issubclass(cls, Mapping)
        assert issubclass(cls, SymbolicMapping)


# -- D8: AnalyticMapping no longer inherits the deprecated Mapping

def test_analytic_mapping_no_longer_inherits_deprecated_mapping():
    from sympde.topology import AnalyticMapping, Mapping, DefinedMapping, SymbolicMapping
    from sympde.topology import PolarMapping
    assert not issubclass(AnalyticMapping, Mapping)
    assert Mapping not in PolarMapping.__mro__
    assert issubclass(AnalyticMapping, DefinedMapping)
    assert issubclass(AnalyticMapping, SymbolicMapping)
    assert not inspect.isabstract(AnalyticMapping)
    # no C3 duplication of SymbolicMapping in the MRO
    assert [c.__name__ for c in AnalyticMapping.__mro__].count('SymbolicMapping') == 1
    assert issubclass(Mapping, SymbolicMapping)


def test_analytic_expression_machinery_survives_the_move():
    from sympy import ImmutableDenseMatrix, eye
    from sympde.topology import IdentityMapping, PolarMapping, TorusSurfaceMapping, AnalyticMapping

    F = IdentityMapping('F', dim=2)
    assert F.jacobian_expr == ImmutableDenseMatrix(eye(2))
    assert F.jacobian_inv_expr == ImmutableDenseMatrix(eye(2))
    assert F.metric_det_expr == 1

    # no numeric constants given: all four stay symbolic Constants
    from sympy import Symbol, cos, sin
    from sympde.core.basic import Constant
    P0 = PolarMapping('P0', dim=2)
    assert set(a.name for a in P0.constants) == {'c1', 'c2', 'rmin', 'rmax'}
    c1, c2, rmin, rmax = (Constant(n) for n in ('c1', 'c2', 'rmin', 'rmax'))
    # the logical coordinates on the mapping are real-valued Symbols
    x1, x2 = (Symbol(n, real=True) for n in ('x1', 'x2'))
    expected = (c1 + (rmin*(1 - x1) + rmax*x1)*cos(x2),
                c2 + (rmin*(1 - x1) + rmax*x1)*sin(x2))
    assert P0.expressions == expected

    # pdim != ldim: the _inv_jac branch is None
    T = TorusSurfaceMapping('T', ldim=2, pdim=3, R0=1., a=0.3)
    assert T.jacobian_expr.shape == (3, 2)
    assert T.jacobian_inv_expr is None
    assert T.metric_expr.shape == (2, 2)

    # early-return branch of the new __new__: no _expressions on the class
    bare = AnalyticMapping('bare', dim=2)
    assert bare.is_analytical is False


def test_legacy_mapping_subclass_still_expands_expressions():
    from sympy import ImmutableDenseMatrix
    from sympde.topology import Mapping, AnalyticMapping

    class M(Mapping):
        _expressions = {'x': '2*x1', 'y': '3*x2'}
        _ldim = 2
        _pdim = 2

    with warnings.catch_warnings():
        warnings.simplefilter('error', DeprecationWarning)
        m = M('m')          # `cls is Mapping` gate: subclass construction is quiet

    assert m.jacobian_expr == ImmutableDenseMatrix([[2, 0], [0, 3]])
    assert m.metric_det_expr == 36
    assert not isinstance(m, AnalyticMapping)


def test_analytic_mapping_construction_is_warning_free():
    from sympde.topology import IdentityMapping, PolarMapping, TorusSurfaceMapping
    with warnings.catch_warnings():
        warnings.simplefilter('error', DeprecationWarning)
        F = IdentityMapping('F', dim=2)
        P = PolarMapping('P', dim=2, c1=0., c2=0., rmin=.3, rmax=1.)
        T = TorusSurfaceMapping('T', ldim=2, pdim=3, R0=1., a=0.3)
        F.copy()
        F.func(*F.args)


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
    # WP13/D2: AnalyticMapping.set_callable_mapping now raises (it is always
    # its own callable), so this is retargeted to a plain SymbolicMapping.
    from sympde.topology import SymbolicMapping
    from sympde.topology.mapping import BasicCallableMapping

    class Custom(BasicCallableMapping):
        def __call__(self, *eta):     return eta
        def jacobian(self, *eta):     return None
        def jacobian_inv(self, *eta): return None
        def metric(self, *eta):       return None
        def metric_det(self, *eta):   return None
        ldim = 2
        pdim = 2

    F  = SymbolicMapping('F', dim=2)
    cm = Custom()
    with pytest.warns(DeprecationWarning):
        F.set_callable_mapping(cm)
    # copy() assigns `_callable_map` directly, not through the deprecated
    # setter, so it must not warn. Sitting outside the `pytest.warns` block
    # above does NOT assert that -- pytest.warns says nothing about warnings
    # raised elsewhere -- so pin it explicitly, as the guards' test does.
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        G = F.copy()
    assert G._callable_map is cm


# -- work-package 03 / 06d-4a: InverseMapping / InterfaceMapping /
# MultiPatchMapping are PURE StructuralMapping subclasses (WP03 added
# StructuralMapping alongside Mapping; 06d-4a severed Mapping).

def test_structural_subclasses_hierarchy():
    from sympde.topology import (StructuralMapping, InverseMapping,
                                 InterfaceMapping, MultiPatchMapping, Mapping,
                                 SymbolicMapping)
    for cls in (InverseMapping, InterfaceMapping, MultiPatchMapping):
        assert issubclass(cls, StructuralMapping)
        assert issubclass(cls, SymbolicMapping)      # still a symbolic mapping
        assert not issubclass(cls, Mapping)          # 06d-4a: Mapping severed
        assert not inspect.isabstract(cls)
        # own concrete ldim/pdim, not StructuralMapping's abstract ones
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
    from sympde.topology import MultiPatchMapping, IdentityMapping, Square
    mp = MultiPatchMapping({'p1': IdentityMapping('F1', dim=2),
                            'p2': IdentityMapping('F2', dim=2)})
    assert mp.ldim == 2
    with pytest.raises(TypeError, match='StructuralMapping'):
        mp(0.3, 0.4)
    # 06d-4a fixed the parked bug: MultiPatchMapping.__new__ now sets _name, so
    # `mp.name`, `mp == mp` and the domain-call path all work.
    assert mp.name == 'F1|F2'
    assert mp == mp
    mapped = mp(Square('D'))
    assert type(mapped).__name__ in ('Domain', 'MappedDomain')


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


def test_structural_mapping_call_returns_mapped_domain():
    # 06d-4a: StructuralMapping.__call__ owns the domain-call body directly
    # (Mapping is no longer in the MRO to super() into). Observable result is
    # unchanged: a single BasicDomain -> a mapped domain.
    from sympde.topology import IdentityMapping, InterfaceMapping, Square
    itf = InterfaceMapping(IdentityMapping('F1', dim=2), IdentityMapping('F2', dim=2))
    assert type(itf(Square('D'))).__name__ in ('Domain', 'MappedDomain')


def test_structural_ldim_pdim_are_concrete():
    # 06d-4a: the `ldim = Mapping.ldim` MRO work-around is gone -- each
    # structural class has its own concrete ldim/pdim returning stored values.
    from sympde.topology import (InverseMapping, InterfaceMapping,
                                 StructuralMapping, IdentityMapping)
    F = IdentityMapping('F', dim=2)
    itf = InterfaceMapping(IdentityMapping('A', dim=2), IdentityMapping('B', dim=2))
    assert InverseMapping(F).ldim == 2 and InverseMapping(F).pdim == 2
    assert itf.ldim == 2 and itf.pdim == 2
    assert InterfaceMapping.ldim is not StructuralMapping.ldim
    assert InverseMapping.pdim   is not StructuralMapping.pdim


def test_interface_mapping_copy_reconstructs():
    # 06d-4a fixed the parked bug: InterfaceMapping used to inherit
    # Mapping.copy() (type(self)(self.name, ldim=...) -- wrong signature).
    # It now has its own copy() that rebuilds from the stored legs.
    from sympde.topology import IdentityMapping, InterfaceMapping
    itf = InterfaceMapping(IdentityMapping('F1', dim=2), IdentityMapping('F2', dim=2))
    c   = itf.copy()
    assert isinstance(c, InterfaceMapping)
    assert c == itf
    assert c.minus.name == itf.minus.name and c.plus.name == itf.plus.name


def test_structural_mappings_are_not_analytic():
    # 06d-4a: the structural object itself carries no analytic `_expressions`.
    # 06d-4c: the callable-mapping API lives on SymbolicMapping now (psydac
    # attaches spline callables to undefined mappings), so a structural mapping
    # inherits get_callable_mapping -- but it has nothing to build / return.
    from sympde.topology import IdentityMapping, InterfaceMapping, MultiPatchMapping
    itf = InterfaceMapping(IdentityMapping('F1', dim=2), IdentityMapping('F2', dim=2))
    mp  = MultiPatchMapping({'p': IdentityMapping('F', dim=2)})
    for m in (itf, mp):
        assert getattr(m, '_expressions', None) is None
        with pytest.raises((ValueError, TypeError, AttributeError)):
            m.get_callable_mapping()


def test_structural_mappings_reject_set_callable_mapping():
    # WP09: a StructuralMapping is symbolic and not point-evaluable (__call__
    # already rejects point calls) -- set_callable_mapping (inherited
    # unguarded from SymbolicMapping otherwise) must not be able to silently
    # make get_callable_mapping() start returning something. Mirrors
    # DiscreteMapping.set_callable_mapping's guard (WP07d).
    from sympde.topology import IdentityMapping, InterfaceMapping, MultiPatchMapping
    itf = InterfaceMapping(IdentityMapping('F1', dim=2), IdentityMapping('F2', dim=2))
    mp  = MultiPatchMapping({'p': IdentityMapping('F', dim=2)})
    for m in (itf, mp):
        with pytest.raises(TypeError):
            m.set_callable_mapping(object())


def test_analytic_mapping_set_callable_mapping_raises():
    # WP13/D2: mirrors DiscreteMapping (WP07d) / StructuralMapping (WP09) --
    # an AnalyticMapping is already its own callable mapping.
    from sympde.topology import IdentityMapping, PolarMapping
    for F in (IdentityMapping('F', dim=2),
              PolarMapping('F', dim=2, rmin=0., rmax=1., c1=0., c2=0.)):
        with pytest.raises(TypeError, match='to_defined_mapping'):
            F.set_callable_mapping(object())


def test_analytic_mapping_get_callable_mapping_is_always_self():
    from sympde.topology import IdentityMapping, PolarMapping
    for F in (IdentityMapping('F', dim=2),
              PolarMapping('F', dim=2, rmin=0., rmax=1., c1=0., c2=0.)):
        assert F.get_callable_mapping() is F


def test_symbolic_mapping_set_callable_mapping_is_deprecated_but_works():
    # D3-b: SymbolicMapping.set_callable_mapping is deprecated (conventions:
    # deprecated symbols keep working until the explicit remove-aliases
    # work-package) -- plain SymbolicMapping.set_callable_mapping/
    # get_callable_mapping stays unguarded, but now warns.
    from sympde.topology import SymbolicMapping
    from sympde.topology.mapping import BasicCallableMapping

    class Custom(BasicCallableMapping):
        def __call__(self, *eta):     return eta
        def jacobian(self, *eta):     return None
        def jacobian_inv(self, *eta): return None
        def metric(self, *eta):       return None
        def metric_det(self, *eta):   return None
        ldim = 2
        pdim = 2

    F  = SymbolicMapping('F', dim=2)
    cm = Custom()
    with pytest.warns(DeprecationWarning):
        F.set_callable_mapping(cm)
    assert F.get_callable_mapping() is cm


def test_set_callable_mapping_warning_names_discrete_mapping():
    # D3-b: the deprecation warning must point callers at the replacement.
    from sympde.topology import SymbolicMapping
    from sympde.topology.mapping import BasicCallableMapping

    class Custom(BasicCallableMapping):
        def __call__(self, *eta):     return eta
        def jacobian(self, *eta):     return None
        def jacobian_inv(self, *eta): return None
        def metric(self, *eta):       return None
        def metric_det(self, *eta):   return None
        ldim = 2
        pdim = 2

    F = SymbolicMapping('F', dim=2)
    with pytest.warns(DeprecationWarning, match='DiscreteMapping'):
        F.set_callable_mapping(Custom())


def test_set_callable_mapping_guards_do_not_warn():
    # D3-b: the three overriding guards (AnalyticMapping, DiscreteMapping,
    # StructuralMapping subclasses) raise TypeError without calling super(),
    # so they must not start emitting the new DeprecationWarning.
    from sympde.topology import (IdentityMapping, InterfaceMapping,
                                 MultiPatchMapping)
    from sympde.topology.mapping import DiscreteMapping, BasicCallableMapping

    class Custom(BasicCallableMapping):
        def __call__(self, *eta):     return eta
        def jacobian(self, *eta):     return None
        def jacobian_inv(self, *eta): return None
        def metric(self, *eta):       return None
        def metric_det(self, *eta):   return None
        ldim = 2
        pdim = 2

    F1  = IdentityMapping('F1', dim=2)
    F2  = IdentityMapping('F2', dim=2)
    G   = DiscreteMapping(Custom(), 'G')
    itf = InterfaceMapping(F1, F2)
    mp  = MultiPatchMapping({'p1': F1, 'p2': F2})

    with warnings.catch_warnings():
        warnings.simplefilter('error')
        for mapping in (F1, G, itf, mp):
            with pytest.raises(TypeError):
                mapping.set_callable_mapping(Custom())


def test_symbolicexpr_lowers_indexed_structural_mapping_to_coordinate():
    # 06d-4a-1: SymbolicExpr.eval's Indexed branch must recognise a severed
    # structural mapping (SymbolicMapping, not Mapping) so `itf[i]` lowers to a
    # coordinate symbol, not the `{base.name}_{i}` fallback (an invalid
    # 'F1|F2_0').
    from sympde.topology import IdentityMapping, InterfaceMapping
    from sympde.topology.mapping import SymbolicExpr
    from sympy import Symbol
    itf = InterfaceMapping(IdentityMapping('F1', dim=2), IdentityMapping('F2', dim=2))
    assert SymbolicExpr(itf[0]) == Symbol('x')
    assert SymbolicExpr(itf[1]) == Symbol('y')


def test_empty_multipatch_mapping_still_constructs():
    # 06d-4a-1: MultiPatchMapping({}) constructed pre-WP06d-4a via Basic.__new__;
    # the _name/_coordinates derivation must not StopIteration on an empty dict.
    from sympde.topology import MultiPatchMapping
    mp = MultiPatchMapping({})
    assert mp.name == ''


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


# -- work-package 06a: Mapping re-parented onto SymbolicMapping

def test_mapping_is_a_symbolic_mapping():
    from sympde.topology import (Mapping, SymbolicMapping, IdentityMapping,
                                 InterfaceMapping)
    assert issubclass(Mapping, SymbolicMapping)
    # every branch of the hierarchy is now isinstance(_, SymbolicMapping):
    with pytest.warns(DeprecationWarning):                                # 06d-4c
        bare = Mapping('F', dim=2)
    assert isinstance(bare, SymbolicMapping)                              # undefined
    assert isinstance(IdentityMapping('G', dim=2), SymbolicMapping)       # analytic
    itf = InterfaceMapping(IdentityMapping('A', dim=2), IdentityMapping('B', dim=2))
    assert isinstance(itf, SymbolicMapping)                               # structural
    # SymbolicMapping appears exactly once in the MRO (no C3 duplication)
    assert [c.__name__ for c in Mapping.__mro__].count('SymbolicMapping') == 1


# -- work-package 06c: AnalyticMapping is its own callable mapping

def test_analytic_mapping_is_its_own_callable_mapping():
    from sympde.topology import IdentityMapping, Square
    F = IdentityMapping('F', dim=2)
    assert F.get_callable_mapping() is F
    # point evaluation (positional) and domain call (positional or `domain=`)
    assert F(0.3, 0.4) == (0.3, 0.4)
    assert type(F(Square('D'))).__name__ in ('Domain', 'MappedDomain')
    assert type(F(domain=Square('D'))).__name__ in ('Domain', 'MappedDomain')


def test_analytic_mapping_with_symbolic_constants_rejects_point_call():
    from sympde.topology import PolarMapping
    P = PolarMapping('P', dim=2)          # rmin/rmax/c1/c2 left symbolic
    with pytest.raises(ValueError, match='symbolic constants'):
        P(0.5, 0.5)


def test_analytic_mapping_rejects_attached_callable():
    # WP13/D2: an AnalyticMapping is already its own callable mapping, so
    # attaching a different one (which used to silently "win" for point
    # evaluation while is_analytical stayed True for assembly) now raises.
    # (Historically -- pre-WP13 -- this test asserted the opposite: that the
    # attached callable "won". Kept as its own regression test, alongside the
    # more general test_analytic_mapping_set_callable_mapping_raises above,
    # because it exercises a real BasicCallableMapping implementation rather
    # than a bare object().)
    from sympde.topology import PolarMapping
    from sympde.topology.mapping import BasicCallableMapping

    class Const(BasicCallableMapping):
        ldim = pdim = 2
        def __call__(self, *eta):      return (7.0, 8.0)
        def jacobian(self, *eta):      return [[1.0, 0.0], [0.0, 1.0]]
        def jacobian_inv(self, *eta):  return [[1.0, 0.0], [0.0, 1.0]]
        def metric(self, *eta):        return [[1.0, 0.0], [0.0, 1.0]]
        def metric_det(self, *eta):    return 1.0

    F = PolarMapping('F', dim=2, rmin=0.0, rmax=1.0, c1=0.0, c2=0.0)
    with pytest.raises(TypeError, match='to_defined_mapping'):
        F.set_callable_mapping(Const())


# -- work-package 06d-2: CallableMapping deleted, BasicMapping folded away

def test_callable_mapping_is_removed():
    import sympde.topology.callable_mapping as cm
    with pytest.raises(AttributeError, match='removed'):
        cm.CallableMapping


def test_bare_mapping_with_expressions_rejects_get_callable_mapping():
    # a Mapping subclass carrying _expressions but not parented under
    # AnalyticMapping is now a mistake -- there is no CallableMapping to build.
    from sympde.topology import Mapping

    class M(Mapping):
        _expressions = {'x': '2*x1', 'y': '3*x2'}
        _ldim = 2
        _pdim = 2

    with pytest.raises(TypeError, match='AnalyticMapping'):
        M('m').get_callable_mapping()


def test_undefined_mapping_still_valueerrors_on_get_callable_mapping():
    # 06d-4c: get_callable_mapping() moved onto SymbolicMapping (an undefined
    # mapping with no _expressions and no attached callable still ValueErrors).
    from sympde.topology import SymbolicMapping
    with pytest.raises(ValueError):
        SymbolicMapping('F', dim=2).get_callable_mapping()


def test_basicmapping_alias_removed():
    # 06d-2 folded BasicMapping into SymbolicMapping (keeping a deprecated
    # alias); 06d-4b drops the alias; 06d-4c relocates SymbolicMapping out of
    # sympde.core.basic into sympde.topology.mapping.
    with pytest.raises(ImportError):
        from sympde.core.basic import BasicMapping  # noqa: F401
    with pytest.raises(ImportError):
        from sympde.core.basic import SymbolicMapping  # noqa: F401

    from sympde.topology import Mapping, SymbolicMapping, IdentityMapping
    assert issubclass(Mapping, SymbolicMapping)
    assert isinstance(IdentityMapping('G', dim=2), SymbolicMapping)
    # SymbolicMapping appears exactly once in the MRO (BasicMapping is gone)
    assert [c.__name__ for c in Mapping.__mro__].count('SymbolicMapping') == 1


# -- work-package 06d-4c: SymbolicMapping is the undefined-mapping constructor,
# `Mapping` is a DeprecationWarning shell over it

def test_symbolicmapping_is_the_undefined_mapping_constructor():
    from sympde.topology import SymbolicMapping, Square
    F = SymbolicMapping('F', dim=2)
    assert F.name == 'F' and F.ldim == 2 and F.pdim == 2
    assert F.is_analytical is False
    assert type(F(Square('D'))).__name__ in ('Domain', 'MappedDomain')


def test_bare_mapping_construction_is_deprecated():
    import warnings
    from sympde.topology import Mapping, IdentityMapping
    with pytest.warns(DeprecationWarning, match='SymbolicMapping'):
        Mapping('M', dim=2)
    with warnings.catch_warnings():          # analytic subclasses do NOT warn
        warnings.simplefilter('error', DeprecationWarning)
        IdentityMapping('G', dim=2)


def test_bare_mapping_rebuild_does_not_warn():
    # /code-review finding 3 (06d-4c-1a): sympy's internal `func(*args)` rebuild
    # (Basic.rebuild, cse, pickling, deepcopy) re-enters Mapping.__new__ with
    # `name` already a Symbol -- it must stay quiet, or a downstream running
    # `filterwarnings=error` breaks on expressions it merely stored / copied.
    import warnings, pickle, copy
    from sympde.topology import Mapping

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        M = Mapping('F_rebuild', dim=2)

    with warnings.catch_warnings():
        warnings.simplefilter('error', DeprecationWarning)
        M.func(*M.args)                       # the exact rebuild path
        copy.deepcopy(M)
        pickle.loads(pickle.dumps(M))


def test_symbolicmapping_honours_injected_jacobian():
    # /code-review finding 2 (06d-4c-1a): the `jacobian=` constructor hook that
    # Mapping.__new__ honoured must survive the move into SymbolicMapping.__new__,
    # and must not leak into the analytic-constants dict of the Mapping shell.
    from sympy import ImmutableDenseMatrix, eye
    from sympde.topology import SymbolicMapping, Mapping

    J = ImmutableDenseMatrix(eye(2))

    F = SymbolicMapping('F_injjac', dim=2, jacobian=J)
    assert F.jacobian_symbol is J

    with pytest.warns(DeprecationWarning):          # bare Mapping still deprecated
        G = Mapping('G_injjac', dim=2, jacobian=J)
    # forwarded through the shell to SymbolicMapping.__new__ (not swallowed by
    # `**kwargs` and replaced by the default JacobianSymbol).
    assert G.jacobian_symbol is J


def test_interface_mapping_from_bare_mapping_legs_is_warning_free_on_copy():
    # /code-review finding 2 (06d-4c-1b): SymbolicMapping.copy() re-invokes the
    # constructor; for a bare Mapping leg that re-tripped the 06d-4c-1a
    # `isinstance(name, str)` warning gate. InterfaceMapping.__new__ copies both
    # legs, so building one from bare-Mapping legs must not warn -- while an
    # explicit Mapping(...) call still does.
    import warnings
    from sympde.topology import Mapping, InterfaceMapping

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        a = Mapping('a_bare', dim=2)
        b = Mapping('b_bare', dim=2)

    with warnings.catch_warnings():
        warnings.simplefilter('error', DeprecationWarning)
        InterfaceMapping(a, b)                  # copies a, b internally -- silent

    with pytest.warns(DeprecationWarning):      # explicit construction still warns
        Mapping('c_bare', dim=2)


# -- WP07: DiscreteMapping -- a DefinedMapping backed by an external callable

class _FakeCallable(BasicCallableMapping):
    """ Minimal BasicCallableMapping for the DiscreteMapping tests. """
    def __init__(self, ldim=2, pdim=2, name=None):
        self._l, self._p, self._n = ldim, pdim, name
    def __call__(self, *e):      return tuple(e)
    def jacobian(self, *e):      return [[1, 0], [0, 1]]
    def jacobian_inv(self, *e):  return [[1, 0], [0, 1]]
    def metric(self, *e):        return [[1, 0], [0, 1]]
    def metric_det(self, *e):    return 1.0
    @property
    def ldim(self): return self._l
    @property
    def pdim(self): return self._p
    @property
    def name(self): return self._n


def test_discrete_mapping_is_a_concrete_defined_mapping():
    from sympde.topology import DiscreteMapping, DefinedMapping

    G = DiscreteMapping(_FakeCallable(), name='D', dim=2)
    assert isinstance(G, DefinedMapping)
    assert isinstance(G, SymbolicMapping)
    assert not inspect.isabstract(DiscreteMapping)
    assert G.name == 'D' and G.ldim == 2 and G.pdim == 2
    assert G.is_analytical is False


def test_discrete_mapping_delegates_point_evaluation():
    from sympde.topology import DiscreteMapping

    f = _FakeCallable()
    G = DiscreteMapping(f, name='D', dim=2)
    assert G.get_callable_mapping() is f
    assert G.jacobian(0.1, 0.2) == [[1, 0], [0, 1]]
    assert G(0.3, 0.4) == (0.3, 0.4)          # point call, delegated


def test_discrete_mapping_is_callable_on_a_domain():
    from sympde.topology import DiscreteMapping, Square

    G  = DiscreteMapping(_FakeCallable(), name='D', dim=2)
    Om = G(Square('S'))
    assert type(Om).__name__ in ('Domain', 'MappedDomain')
    assert Om.mapping is G
    assert Om.logical_domain == Square('S')


def test_discrete_mapping_copy_round_trips():
    from sympde.topology import DiscreteMapping

    f = _FakeCallable()
    G = DiscreteMapping(f, name='D', dim=2)
    c = G.copy()
    assert type(c) is DiscreteMapping
    assert c == G
    assert c.get_callable_mapping() is f


def test_discrete_mapping_name_inference_and_error():
    from sympde.topology import DiscreteMapping

    assert DiscreteMapping(_FakeCallable(name='H')).name == 'H'   # inferred
    with pytest.raises(ValueError):
        DiscreteMapping(_FakeCallable())                          # no name anywhere


def test_discrete_mapping_rejects_non_callable():
    from sympde.topology import DiscreteMapping
    with pytest.raises(TypeError):
        DiscreteMapping(object(), name='D', dim=2)


def test_defined_mapping_factory_sugar():
    from sympde.topology import DefinedMapping, DiscreteMapping

    G = DefinedMapping(_FakeCallable(), name='D', dim=2)
    assert type(G) is DiscreteMapping
    # a plain abstract call still fails
    with pytest.raises(TypeError, match='abstract'):
        DefinedMapping('F')


def test_discrete_mapping_usable_as_interface_leg():
    from sympde.topology import DiscreteMapping, InterfaceMapping

    a = DiscreteMapping(_FakeCallable(), name='a', dim=2)
    b = DiscreteMapping(_FakeCallable(), name='b', dim=2)
    itf = InterfaceMapping(a, b)
    assert itf.is_analytical is False


# -- WP07-1: post-/code-review fixes for DiscreteMapping

def test_discrete_mapping_identity_includes_the_wrapped_callable():
    # F1: two same-named DiscreteMappings over different geometries must be
    # distinct, else MappedDomain's @cacheit conflates them.
    from sympde.topology import DiscreteMapping, Square

    fa, fb = _FakeCallable(), _FakeCallable()
    Ga = DiscreteMapping(fa, name='G', dim=2)
    Gb = DiscreteMapping(fb, name='G', dim=2)
    assert Ga != Gb
    assert hash(Ga) != hash(Gb)
    assert len({Ga, Gb}) == 2

    Da, Db = Ga(Square('S1')), Gb(Square('S1'))
    assert Da is not Db
    assert Da.mapping.get_callable_mapping() is fa
    assert Db.mapping.get_callable_mapping() is fb

    # same callable + name + dims -> still equal
    assert DiscreteMapping(fa, name='G', dim=2) == Ga


def test_discrete_mapping_dims_come_from_the_callable():
    # F4: a dim= / ldim= / pdim= argument may only confirm the callable's dims.
    from sympde.topology import DiscreteMapping

    surf = _FakeCallable(ldim=2, pdim=3)
    G = DiscreteMapping(surf, name='S', ldim=2, pdim=3)
    assert G.ldim == 2 and G.pdim == 3

    with pytest.raises(ValueError):
        DiscreteMapping(surf, name='S2', dim=2)          # 2 != pdim 3
    with pytest.raises(ValueError):
        DiscreteMapping(_FakeCallable(), name='S3', pdim=5)

    assert DiscreteMapping(_FakeCallable(), name='D', dim=2).pdim == 2   # consistent, ok


def test_discrete_mapping_copy_preserves_interface_tags():
    # F2: copy() must carry _is_plus / _is_minus (and coordinates), not just
    # (callable, name, ldim, pdim).
    from sympde.topology import DiscreteMapping, InterfaceMapping

    itf = InterfaceMapping(DiscreteMapping(_FakeCallable(), 'a', dim=2),
                           DiscreteMapping(_FakeCallable(), 'b', dim=2))
    assert itf.minus.is_minus is True
    assert itf.minus.copy().is_minus is True
    assert itf.plus.copy().is_plus is True


def test_discrete_mapping_reconstruction_fails_clearly():
    # F3: sympy's func(*args) / pickle / deepcopy shape must raise an actionable
    # error, not a confusing "got Symbol".
    from sympy import Symbol, Tuple
    from sympde.topology import DiscreteMapping

    with pytest.raises(TypeError, match='reconstructed'):
        DiscreteMapping(Symbol('D'), Tuple(2))


def test_discrete_mapping_get_callable_guards_none():
    # F3: defensive -- a callable-less DiscreteMapping raises, not AttributeError.
    from sympde.topology import DiscreteMapping

    G = DiscreteMapping(_FakeCallable(), name='D', dim=2)
    G._callable_map = None
    with pytest.raises(ValueError):
        G.get_callable_mapping()


def test_discrete_mapping_has_callable_mapping():
    # WP07e: a non-raising predicate for get_callable_mapping()'s guard, so
    # callers don't need to catch ValueError as control flow.
    from sympde.topology import DiscreteMapping

    G = DiscreteMapping(_FakeCallable(), name='D', dim=2)
    assert G.has_callable_mapping() is True

    G._callable_map = None
    assert G.has_callable_mapping() is False


def test_discrete_mapping_set_callable_mapping_raises():
    # WP07d: the wrapped callable is fixed at construction (part of identity,
    # see _hashable_content) -- set_callable_mapping must not be able to swap
    # it, unlike the base SymbolicMapping.set_callable_mapping it would
    # otherwise inherit unguarded.
    from sympde.topology import DiscreteMapping

    G = DiscreteMapping(_FakeCallable(), name='D', dim=2)
    with pytest.raises(TypeError):
        G.set_callable_mapping(_FakeCallable())


# -- WP07-2: second /code-review round for DiscreteMapping

def test_discrete_mapping_jacobian_expr_uses_the_final_identity():
    # G1: __new__ must attach _callable_map before building _jac / _metric, so
    # `G[i]` and the Indexed(G, i) baked into jacobian_expr hash identically.
    from sympy import Symbol
    from sympde.topology import DiscreteMapping

    G = DiscreteMapping(_FakeCallable(), name='D', dim=2)
    assert G[0] in G.jacobian_expr.free_symbols
    assert G[1] in G.metric_expr.free_symbols
    assert G.jacobian_expr.subs(G[0], Symbol('v')) != G.jacobian_expr

    Ga = DiscreteMapping(_FakeCallable(), name='S', dim=2)
    Gb = DiscreteMapping(_FakeCallable(), name='S', dim=2)
    assert all(s.base is Ga for s in Ga.jacobian_expr.free_symbols
               if hasattr(s, 'base'))
    assert all(s.base is Gb for s in Gb.jacobian_expr.free_symbols
               if hasattr(s, 'base'))


def test_discrete_mapping_delegators_guard_none_callable():
    # G2: every delegating point-eval method routes through
    # get_callable_mapping(), so a None _callable_map gives a ValueError, not
    # AttributeError / TypeError.
    from sympde.topology import DiscreteMapping

    G = DiscreteMapping(_FakeCallable(), name='D', dim=2)
    G._callable_map = None
    for call in (lambda: G.jacobian(0., 0.),
                 lambda: G.jacobian_inv(0., 0.),
                 lambda: G.metric(0., 0.),
                 lambda: G.metric_det(0., 0.),
                 lambda: G(0.1, 0.2)):
        with pytest.raises(ValueError):
            call()


def test_discrete_mapping_name_rules():
    # G3: an explicit name wins and must be non-empty; a falsy name is not
    # silently replaced by the callable's.
    from sympde.topology import DiscreteMapping

    with pytest.raises(ValueError):
        DiscreteMapping(_FakeCallable(name='H'), name='')
    with pytest.raises(ValueError):
        DiscreteMapping(_FakeCallable(name=None))
    assert DiscreteMapping(_FakeCallable(name='H')).name == 'H'
    assert DiscreteMapping(_FakeCallable(name='H'), name='K').name == 'K'


def test_discrete_mapping_rejects_unexpected_kwargs():
    from sympde.topology import DiscreteMapping
    with pytest.raises(TypeError):
        DiscreteMapping(_FakeCallable(), name='D', dim=2, bogus=1)


# -- WP07b-1: DiscreteMapping survives sympy's func(*args) rebuild

def test_discrete_mapping_func_is_the_identity_rebuild():
    # sympy walks an expression tree calling `node.func(*node.args)`
    # (Basic.rebuild, cse, and sympy.simplify's replace-reducer). The base
    # `Basic.func` is `type(self)`, so it would call
    # `DiscreteMapping(Symbol(name), Tuple(pdim))` -> TypeError. A DiscreteMapping
    # is a symbolic leaf (the callable is not in `.args`), so `func(*args)` must
    # be the identity -- this is what lets a bilinear form be discretised on a
    # DiscreteMapping-carried domain (Jacobian(M).inv() -> simplify -> func).
    import sympy
    from sympde.topology import DiscreteMapping
    from sympde.topology.mapping import Jacobian

    G = DiscreteMapping(_FakeCallable(), name='D', dim=2)

    assert G.func(*G.args) is G
    assert G.func(sympy.Symbol('x'), (99,), foo='bar') is G
    assert sympy.simplify(G[0]**2 / G[1]) == G[0]**2 / G[1]     # no crash
    Jacobian(G).inv()                                            # the assembly path; no crash
