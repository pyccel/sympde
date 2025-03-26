# coding: utf-8

import pytest

from sympde.topology import Domain
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import element_of
from sympde.calculus import grad, dot, inner, cross, curl, div

#from sympde.exterior import d, wedge, ip, delta, jp
from sympde.exterior import ld
from sympde.exterior import DifferentialForm
from sympde.exterior import ExteriorCalculusExpr, augmented_expression
from sympde.calculus.errors import ArgumentTypeError

#==============================================================================

def test_compiler_3d_1():

    domain = Domain('Omega', dim=3)

    H1    = ScalarFunctionSpace('V0', domain, kind='H1')
    Hcurl = VectorFunctionSpace('V1', domain, kind='Hcurl')
    Hdiv  = VectorFunctionSpace('V2', domain, kind='Hdiv')
    L2    = ScalarFunctionSpace('V3', domain, kind='L2')
    V     = VectorFunctionSpace('V', domain)

    X = H1 * Hcurl * Hdiv * L2

    v0, v1, v2, v3 = element_of(X, name='v0, v1, v2, v3')

    beta = element_of(V, 'beta')

#    # ...
#    expr = grad(v0)
#    expected = d(DifferentialForm('v0', index=0, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr) == expected)
#    # ...
#
#    # ...
#    expr = curl(v1)
#    expected = d(DifferentialForm('v1', index=1, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr) == expected)
#    # ...
#
#    # ...
#    expr = div(v2)
#    expected = d(DifferentialForm('v2', index=2, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr) == expected)
#    # ...
#
#    # ...
#    expr = grad(v0)
#    expected = - delta(DifferentialForm('v0', index=3, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr, tests=[v0]) == expected)
#    # ...
#
#    # ...
#    expr = curl(v1)
#    expected = delta(DifferentialForm('v1', index=2, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr, tests=[v1]) == expected)
#    # ...
#
#    # ...
#    expr = div(v2)
#    expected = -delta(DifferentialForm('v2', index=1, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr, tests=[v2]) == expected)
#    # ...
#
#    # ...
#    expr = dot(beta, v1)
#    expected = ip(beta, DifferentialForm('v1', index=1, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr) == expected)
#    # ...
#
#    # ...
#    expr = cross(beta, v2)
#    expected = ip(beta, DifferentialForm('v2', index=2, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr) == expected)
#    # ...
#
#    # ...
#    expr = beta*v3
#    expected = ip(beta, DifferentialForm('v3', index=3, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr) == expected)
#    # ...
#
#    # ...
#    expr = dot(beta, v1)
#    expected = jp(beta, DifferentialForm('v1', index=3, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr, tests=[v1]) == expected)
#    # ...
#
#    # ...
#    expr = cross(beta, v2)
#    expected = -jp(beta, DifferentialForm('v2', index=2, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr, tests=[v2]) == expected)
#    # ...
#
#    # ...
#    expr = beta*v3
#    expected = jp(beta, DifferentialForm('v3', index=1, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr, tests=[v3]) == expected)
#    # ...

    # ..................................................
    #     LIE DERIVATIVES
    # ..................................................
    # ...
    expr = dot(beta, grad(v0))
    expected = ld(beta, DifferentialForm('v0', index=0, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr) == expected)
    # ...

    # ...
    expr = grad(dot(beta,v1)) + cross(curl(v1), beta)
    expected = ld(beta, DifferentialForm('v1', index=1, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr) == expected)
    # ...

    # ...
    expr = curl(cross(v2, beta)) + div(v2)*beta
    expected = ld(beta, DifferentialForm('v2', index=2, dim=domain.dim))

    #assert(ExteriorCalculusExpr(expr) == expected)
    # ...

    # ...
    expr = div(beta*v3)
    expected = ld(beta, DifferentialForm('v3', index=3, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr) == expected)
    # ...

    # ... TODO commutativity of LieDerivative
#    expr = grad(dot(beta,v1)) + cross(curl(v1), beta) + div(beta*v3)
#    expected  = ld(beta, DifferentialForm('v1', index=1, dim=domain.dim))
#    expected += ld(beta, DifferentialForm('v3', index=3, dim=domain.dim))
#
#    assert(ExteriorCalculusExpr(expr) == expected)
    # ...

#    print(expr)
#    print('')
#    print(ExteriorCalculusExpr(expr))
    # ..................................................

#==============================================================================
def test_compiler_3d_2():

    domain = Domain('Omega', dim=3)

    H1    = ScalarFunctionSpace('V0', domain, kind='H1')
    Hcurl = VectorFunctionSpace('V1', domain, kind='Hcurl')
    Hdiv  = VectorFunctionSpace('V2', domain, kind='Hdiv')
    L2    = ScalarFunctionSpace('V3', domain, kind='L2')
    V     = VectorFunctionSpace('V', domain)

    X = H1 * Hcurl * Hdiv * L2

    v = element_of(X, name='v0, v1, v2, v3')
    u = element_of(X, name='u0, u1, u2, u3')

#    # ... Dot operator
#    expr = dot(u1, v1)
#    print(ExteriorCalculusExpr(expr, tests=[v1]))
#
#    expr = dot(u2, v2)
#    print(ExteriorCalculusExpr(expr, tests=[v2]))
#
#    expr = dot(grad(v0), u1)
#    print(ExteriorCalculusExpr(expr, tests=[v0]))
#
#    expr = dot(grad(u0), v1)
#    print(ExteriorCalculusExpr(expr, tests=[v1]))
#
#    expr = dot(curl(u1), v2)
#    print(ExteriorCalculusExpr(expr, tests=[v2]))
#
#    expr = dot(curl(v1), u2)
#    print(ExteriorCalculusExpr(expr, tests=[v1]))
#    # ...

    # ... Mul operator
    expr = u[0] * v[0]
    print(ExteriorCalculusExpr(expr, tests=[v[0]]))

    expr = u[0] * div(v[2])
    print(ExteriorCalculusExpr(expr, tests=[v[2]]))

    expr = v[0] * div(u[2])
    print(ExteriorCalculusExpr(expr, tests=[v[0]]))
    # ...

#    # ... Add operator
#    expr = dot(u1, v1) + dot(u2, v2)
#    print(ExteriorCalculusExpr(expr, tests=[v1,v2]))
#
#    expr = dot(u1, v1) + dot(grad(v0), u1)
#    print(ExteriorCalculusExpr(expr, tests=[v0,v1]))
#
#    expr = dot(u1, v1) + dot(grad(v0), u1) + dot(grad(u0), v1)
#    print(ExteriorCalculusExpr(expr, tests=[v0,v1]))
#    # ...


#==============================================================================
def test_compiler_3d_poisson():

    domain = Domain('Omega', dim=3)

    H1    = ScalarFunctionSpace('V0', domain, kind='H1')
    Hcurl = VectorFunctionSpace('V1', domain, kind='Hcurl')
    Hdiv  = VectorFunctionSpace('V2', domain, kind='Hdiv')
    L2    = ScalarFunctionSpace('V3', domain, kind='L2')
    V     = VectorFunctionSpace('V', domain)

    X = Hdiv * L2

    sigma, u = element_of(X, name='sigma, u')
    tau,   v = element_of(X, name='tau,   v')

    expr = dot(sigma, tau) + div(tau)*u + div(sigma)*v
    print(ExteriorCalculusExpr(expr, tests=[tau,v]))

#==============================================================================
def test_compiler_3d_stokes():

    domain = Domain('Omega', dim=3)

    # ...
    #    by setting the space type, we cannot evaluate grad of Hdiv function, then
    #    ArgumentTypeError will be raised.
    #    In order to avoid this problem, we need first to declare our space as an
    #    undefined type.
    Hdiv  = VectorFunctionSpace('V2', domain, kind='Hdiv')
    L2    = ScalarFunctionSpace('V3', domain, kind='L2')

    X = Hdiv * L2

    u, p = element_of(X, name='u, p')
    v, q = element_of(X, name='v, q')

    with pytest.raises(ArgumentTypeError):
        expr = inner(grad(u), grad(v)) - div(v)*p + q*div(u)
    # ...

    # ...
    Hdiv  = VectorFunctionSpace('V2', domain)
    L2    = ScalarFunctionSpace('V3', domain)

    X = Hdiv * L2

    u, p = element_of(X, name='u, p')
    v, q = element_of(X, name='v, q')

    expr = inner(grad(u), grad(v)) - div(v)*p + q*div(u)
    atoms = {u: DifferentialForm('u', index=2, dim=domain.dim),
             v: DifferentialForm('v', index=2, dim=domain.dim),
             p: DifferentialForm('p', index=3, dim=domain.dim),
             q: DifferentialForm('q', index=3, dim=domain.dim)}
    newexpr = ExteriorCalculusExpr(expr, tests=[v,q], atoms=atoms)
    print('===== BEFORE =====')
    print(newexpr)

    newexpr = augmented_expression(newexpr, tests=[v,q], atoms=atoms, weak=False)
    print('===== AFTER  =====')
    print(newexpr)
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

#test_compiler_3d_1()
#test_compiler_3d_2()
test_compiler_3d_poisson()
#test_compiler_3d_stokes()
