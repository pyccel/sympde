# coding: utf-8

from sympy import symbols
from sympy import Tuple
from sympy import Matrix
from sympy import srepr

from sympde.core import Constant
from sympde.calculus import grad, dot, inner
from sympde.topology import Domain, element_of
from sympde.topology import get_index_derivatives_atom
from sympde.topology import get_max_partial_derivatives
from sympde.topology import ScalarFunctionSpace
from sympde.topology import (dx, dy, dz)
from sympde.topology import SymbolicMapping


def indices_as_str(a):
    a = dict(sorted(a.items()))
    code = ''
    for k,n in list(a.items()):
        code += k*n
    return code



# ...
def test_partial_derivatives_1():
    print('============ test_partial_derivatives_1 ==============')

    # ...
    domain = Domain('Omega', dim=2)
    M      = SymbolicMapping('M', dim=2)

    mapped_domain = M(domain)

    x,y = mapped_domain.coordinates

    V = ScalarFunctionSpace('V', mapped_domain)

    F,u,v,w = [element_of(V, name=i) for i in ['F', 'u', 'v', 'w']]
    uvw = Tuple(u,v,w)

    alpha = Constant('alpha')
    beta = Constant('beta')
    # ...

    assert(dx(x**2) == 2*x)
    assert(dy(x**2) == 0)
    assert(dz(x**2) == 0)

    assert(dx(y**2) == 0)
    assert(dy(y**2) == 2*y)
    assert(dz(y**2) == 0)

    assert(dx(x*F) == F + x*dx(F))
    assert(dx(uvw) == Matrix([[dx(u), dx(v), dx(w)]]))
    assert(dx(uvw) + dy(uvw) == Matrix([[dx(u) + dy(u),
                                         dx(v) + dy(v),
                                         dx(w) + dy(w)]]))

    expected = Matrix([[alpha*dx(u) + beta*dy(u),
                        alpha*dx(v) + beta*dy(v),
                        alpha*dx(w) + beta*dy(w)]])
    assert(alpha * dx(uvw) + beta * dy(uvw) == expected)
    # ...

#    expr = alpha * dx(uvw) + beta * dy(uvw)
#    print(expr)

#    print('> ', srepr(expr))
#    print('')
# ...

# ...
def test_partial_derivatives_2():
    print('============ test_partial_derivatives_2 ==============')

    # ...
    domain = Domain('Omega', dim=2)
    M      = SymbolicMapping('M', dim=2)

    mapped_domain = M(domain)

    V = ScalarFunctionSpace('V', mapped_domain)
    F = element_of(V, name='F')

    alpha = Constant('alpha')
    beta = Constant('beta')
    # ...

    # ...
    expr = alpha * dx(F)

    indices = get_index_derivatives_atom(expr, F)[0]
    assert(indices_as_str(indices) == 'x')
    # ...

    # ...
    expr = dy(dx(F))

    indices = get_index_derivatives_atom(expr, F)[0]
    assert(indices_as_str(indices) == 'xy')
    # ...

    # ...
    expr = alpha * dx(dy(dx(F)))

    indices = get_index_derivatives_atom(expr, F)[0]
    assert(indices_as_str(indices) == 'xxy')
    # ...

    # ...
    expr = alpha * dx(dx(F)) + beta * dy(F) + dx(dy(F))

    indices = get_index_derivatives_atom(expr, F)
    indices = [indices_as_str(i) for i in indices]
    assert(sorted(indices) == ['xx', 'xy', 'y'])
    # ...

    # ...
    expr = alpha * dx(dx(F)) + beta * dy(F) + dx(dy(F))

    d = get_max_partial_derivatives(expr, F)
    assert(indices_as_str(d) == 'xxy')

    d = get_max_partial_derivatives(expr)
    assert(indices_as_str(d) == 'xxy')
    # ...
# ...


# ...
def test_logical_derivative_through_symbolic_mapping_index():
    # /code-review finding 1 (WP06d-4c-1a): `_DifferentialOperator.eval`'s
    # logical chain-rule branch matched `expr.atoms(Mapping)`, which no longer
    # catches a bare `SymbolicMapping` (the constructor WP06d-4c recommends and
    # migrated every test to). The chain-rule term through the mapping component
    # `M[i]` was then silently dropped and `dx1(M[0]**2)` collapsed to 0.
    from sympde.topology import dx1, dx2

    M = SymbolicMapping('M', dim=2)

    assert dx1(M[0]**2) == 2 * M[0] * dx1(M[0])
    assert dx2(M[1]**2) == 2 * M[1] * dx2(M[1])

    N = SymbolicMapping('N', dim=3)
    assert dx1(N[2]**3) == 3 * N[2]**2 * dx1(N[2])
# ...


# ...
def test_logical_derivative_through_interface_mapping_component():
    # /code-review finding 1 (WP06d-4c-1b): the 06d-4c-1a widening to
    # `expr.atoms(SymbolicMapping)` now also selects structural mappings. The
    # indexable ones (InterfaceMapping, InverseMapping) must still be
    # chain-ruled through -- the `_shape` filter added in 4c-1b keeps them.
    from sympde.topology import dx1, InterfaceMapping, IdentityMapping

    itf = InterfaceMapping(IdentityMapping('A', dim=2), IdentityMapping('B', dim=2))
    assert dx1(itf[0]**2) == 2 * itf[0] * dx1(itf[0])


def test_multipatch_mapping_is_excluded_from_the_chain_rule_branch():
    # WP06d-4c-1b: the chain-rule branch indexes the mapping (`M[i]`), so it
    # filters `expr.atoms(SymbolicMapping)` to indexable mappings via
    # `getattr(m, '_shape', None) is not None`. MultiPatchMapping is built via
    # Basic.__new__ and has no `_shape`; without the filter `M[i]` raised
    # AttributeError. (It cannot actually be built into a scalar Expr reaching
    # that branch -- its `.args` is a raw dict, which breaks Expr.is_number
    # first -- so the filter is a defensive guard; lock the discriminator here.)
    from sympde.topology import (MultiPatchMapping, InterfaceMapping,
                                 IdentityMapping)

    mp  = MultiPatchMapping({'p': IdentityMapping('F', dim=2)})
    itf = InterfaceMapping(IdentityMapping('A', dim=2), IdentityMapping('B', dim=2))
    assert getattr(mp,  '_shape', None) is None
    assert getattr(itf, '_shape', None) is not None
# ...


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
