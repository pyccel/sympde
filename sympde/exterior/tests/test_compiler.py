# coding: utf-8

import pytest


from sympy import symbols
from sympy import Tuple
from sympy import Matrix
from sympy import srepr

from sympde import Constant

from sympde.topology import Domain
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import ProductSpace
from sympde.topology import H1Space, HcurlSpace, HdivSpace, L2Space, UndefinedSpace
from sympde.topology import TestFunction
from sympde.topology import Field
from sympde.calculus import grad, dot, inner, cross, rot, curl, div

from sympde.exterior import d, wedge, ip, delta, jp
from sympde.exterior import DifferentialForm
from sympde.exterior import ExteriorCalculusExpr
from sympde.calculus.errors import ArgumentTypeError


#==============================================================================
def test_compiler_3d_1():

    domain = Domain('Omega', dim=3)

    H1    =       FunctionSpace('V0', domain, kind='H1')
    Hcurl = VectorFunctionSpace('V1', domain, kind='Hcurl')
    Hdiv  = VectorFunctionSpace('V2', domain, kind='Hdiv')
    L2    =       FunctionSpace('V3', domain, kind='L2')
    V     = VectorFunctionSpace('V', domain)

    X = H1 * Hcurl * Hdiv * L2

    v0, v1, v2, v3 = TestFunction(X, ['v0', 'v1', 'v2', 'v3'])

    beta = Field(V, 'beta')

    # ...
    expr = grad(v0)
    expected = d(DifferentialForm('v0', index=0, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr) == expected)
    # ...

    # ...
    expr = curl(v1)
    expected = d(DifferentialForm('v1', index=1, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr) == expected)
    # ...

    # ...
    expr = div(v2)
    expected = d(DifferentialForm('v2', index=2, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr) == expected)
    # ...

    # ...
    expr = grad(v0)
    expected = - delta(DifferentialForm('v0', index=3, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr, tests=[v0]) == expected)
    # ...

    # ...
    expr = curl(v1)
    expected = delta(DifferentialForm('v1', index=2, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr, tests=[v1]) == expected)
    # ...

    # ...
    expr = div(v2)
    expected = -delta(DifferentialForm('v2', index=1, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr, tests=[v2]) == expected)
    # ...

    # ...
    expr = dot(beta, v1)
    expected = ip(beta, DifferentialForm('v1', index=1, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr) == expected)
    # ...

    # ...
    expr = cross(beta, v2)
    expected = ip(beta, DifferentialForm('v2', index=2, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr) == expected)
    # ...

    # ...
    expr = beta*v3
    expected = ip(beta, DifferentialForm('v3', index=3, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr) == expected)
    # ...

    # ...
    expr = dot(beta, v1)
    expected = jp(beta, DifferentialForm('v1', index=3, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr, tests=[v1]) == expected)
    # ...

    # ...
    expr = cross(beta, v2)
    expected = -jp(beta, DifferentialForm('v2', index=2, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr, tests=[v2]) == expected)
    # ...

    # ...
    expr = beta*v3
    expected = jp(beta, DifferentialForm('v3', index=1, dim=domain.dim))

    assert(ExteriorCalculusExpr(expr, tests=[v3]) == expected)
    # ...

#==============================================================================
def test_compiler_3d_2():

    domain = Domain('Omega', dim=3)

    H1    =       FunctionSpace('V0', domain, kind='H1')
    Hcurl = VectorFunctionSpace('V1', domain, kind='Hcurl')
    Hdiv  = VectorFunctionSpace('V2', domain, kind='Hdiv')
    L2    =       FunctionSpace('V3', domain, kind='L2')
    V     = VectorFunctionSpace('V', domain)

    X = H1 * Hcurl * Hdiv * L2

    v0, v1, v2, v3 = TestFunction(X, ['v0', 'v1', 'v2', 'v3'])
    u0, u1, u2, u3 = TestFunction(X, ['u0', 'u1', 'u2', 'u3'])

    beta = Field(V, 'beta')

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
    expr = u0*v0
    print(ExteriorCalculusExpr(expr, tests=[v0]))

    expr = u0*div(v2)
    print(ExteriorCalculusExpr(expr, tests=[v2]))

    expr = v0*div(u2)
    print(ExteriorCalculusExpr(expr, tests=[v0]))
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

    H1    =       FunctionSpace('V0', domain, kind='H1')
    Hcurl = VectorFunctionSpace('V1', domain, kind='Hcurl')
    Hdiv  = VectorFunctionSpace('V2', domain, kind='Hdiv')
    L2    =       FunctionSpace('V3', domain, kind='L2')
    V     = VectorFunctionSpace('V', domain)

    X = Hdiv * L2

    sigma, u = TestFunction(X, ['sigma', 'u'])
    tau,   v = TestFunction(X, [  'tau', 'v'])

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
    L2    =       FunctionSpace('V3', domain, kind='L2')

    X = Hdiv * L2

    u,p = TestFunction(X, ['u', 'p'])
    v,q = TestFunction(X, ['v', 'q'])

    with pytest.raises(ArgumentTypeError):
        expr = inner(grad(u), grad(v)) - div(v)*p + q*div(u)
    # ...

    # ...
    Hdiv  = VectorFunctionSpace('V2', domain)
    L2    =       FunctionSpace('V3', domain)

    X = Hdiv * L2

    u,p = TestFunction(X, ['u', 'p'])
    v,q = TestFunction(X, ['v', 'q'])

    expr = inner(grad(u), grad(v)) - div(v)*p + q*div(u)
    atoms = {u: DifferentialForm('u', index=2, dim=domain.dim),
             v: DifferentialForm('v', index=2, dim=domain.dim),
             p: DifferentialForm('p', index=3, dim=domain.dim),
             q: DifferentialForm('q', index=3, dim=domain.dim)}
    print(ExteriorCalculusExpr(expr, tests=[v,q], atoms=atoms))
    # ...

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()

#test_compiler_3d_1()
#test_compiler_3d_2()
#test_compiler_3d_poisson()
#test_compiler_3d_stokes()
