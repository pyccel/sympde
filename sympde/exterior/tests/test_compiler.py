# coding: utf-8


from sympy import symbols
from sympy import Tuple
from sympy import Matrix
from sympy import srepr

from sympde import Constant

from sympde.topology import Domain
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import ProductSpace
from sympde.topology import H1Space, HcurlSpace, HdivSpace, L2Space, UndefinedSpace
from sympde.topology import TestFunction, ScalarTestFunction, VectorTestFunction
from sympde.topology import Field, ScalarField, VectorField
from sympde.calculus import grad, dot, inner, cross, rot, curl, div

from sympde.exterior import d, wedge, ip, delta, jp
from sympde.exterior import DifferentialForm
from sympde.exterior import ExteriorCalculusExpr


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

    beta = Field(V, 'beta')

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()

test_compiler_3d_1()
#test_compiler_3d_2()
