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
    # `SymbolicMapping` deliberately does not carry it. On an AnalyticMapping it
    # is shadowed by the numeric point-evaluation method
    # (see test_analytic_mapping_is_directly_point_evaluable).
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


def test_analytic_mapping_honours_explicitly_attached_callable():
    # 06c amendment 2: an AnalyticMapping's point-eval methods must agree with
    # get_callable_mapping() -- a callable attached via set_callable_mapping()
    # wins over the mapping's own lambdified expressions.
    from sympde.topology import PolarMapping
    from sympde.topology.mapping import BasicCallableMapping

    class Const(BasicCallableMapping):          # deliberately "wrong" values,
        ldim = pdim = 2                         # so analytic vs attached differ
        def __call__(self, *eta):      return (7.0, 8.0)
        def jacobian(self, *eta):      return [[1.0, 0.0], [0.0, 1.0]]
        def jacobian_inv(self, *eta):  return [[1.0, 0.0], [0.0, 1.0]]
        def metric(self, *eta):        return [[1.0, 0.0], [0.0, 1.0]]
        def metric_det(self, *eta):    return 1.0

    F = PolarMapping('F', dim=2, rmin=0.0, rmax=1.0, c1=0.0, c2=0.0)
    F.set_callable_mapping(Const())
    assert isinstance(F.get_callable_mapping(), Const)
    assert F(0.5, 0.5) == (7.0, 8.0)            # not the analytic value
    assert F.jacobian(0.5, 0.5) == [[1.0, 0.0], [0.0, 1.0]]
    assert F.metric_det(0.5, 0.5) == 1.0


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
