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

from sympde.exterior import d, wedge, ip
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
    print(ExteriorCalculusExpr(expr))
    # ...

    # ...
    expr = curl(v1)
    print(ExteriorCalculusExpr(expr))
    # ...

    # ...
    expr = div(v2)
    print(ExteriorCalculusExpr(expr))
    # ...

    # ... TODO to be improved
    expr = grad(v0)
    print(ExteriorCalculusExpr(expr, index=3))
    # ...

    # ... TODO to be improved
    expr = curl(v1)
    print(ExteriorCalculusExpr(expr, index=2))
    # ...

    # ... TODO to be improved
    expr = div(v2)
    print(ExteriorCalculusExpr(expr, index=1))
    # ...

    # ...
    expr = dot(beta, v1)
    print(ExteriorCalculusExpr(expr))
    # ...

    # ...
    expr = cross(beta, v2)
    print(ExteriorCalculusExpr(expr))
    # ...

    # ...
    expr = beta*v3
    print(ExteriorCalculusExpr(expr))
    # ...

    # ... TODO to be improved
    expr = dot(beta, v1)
    print(ExteriorCalculusExpr(expr, index=2))
    # ...

    # ... TODO to be improved
    expr = cross(beta, v2)
    print(ExteriorCalculusExpr(expr, index=1))
    # ...

    # ... TODO to be improved
    expr = beta*v3
    print(ExteriorCalculusExpr(expr, index=0))
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

#test_compiler_3d_1()
#test_compiler_3d_2()
