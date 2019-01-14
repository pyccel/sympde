# coding: utf-8

import pytest

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import IndexedBase
from sympy import Matrix
from sympy import Function
from sympy import pi, cos, sin
from sympy import srepr
from sympy.physics.quantum import TensorProduct

from sympde.core import Constant
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import laplace, hessian, bracket
from sympde.topology import Domain
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import Field, VectorField


#==============================================================================
def test_calculus_2d():
    domain = Domain('Omega', dim=2)

    V = FunctionSpace('V', domain)
    W = VectorFunctionSpace('W', domain)

    alpha = Constant('alpha')
    beta  = Constant('beta')

    f,g = [Field(i, space=V) for i in ['f','g']]
    F,G = [VectorField(W, i) for i in ['F','G']]

    # ... gradient properties
    assert( grad(f+g) == grad(f) + grad(g) )
    assert( grad(alpha*f) == alpha*grad(f) )
    assert( grad(f*g) == f*grad(g) + g*grad(f) )
    assert( grad(f/g) == -f*grad(g)/g**2 + grad(f)/g )
    # ...

    # ... curl properties
    assert( curl(f+g) == curl(f) + curl(g) )
    assert( curl(alpha*f) == alpha*curl(f) )
    # ...

    # ... laplace properties
    assert( laplace(f+g) == laplace(f) + laplace(g) )
    assert( laplace(alpha*f) == alpha*laplace(f) )
    # ...

    # ... divergence properties
    assert( div(F+G) == div(F) + div(G) )
    assert( div(alpha*F) == alpha*div(F) )
    # ...

    # ... rot properties
    assert( rot(F+G) == rot(F) + rot(G) )
    assert( rot(alpha*F) == alpha*rot(F) )
    # ...


#    expr = grad(F/G)
#    print(expr)

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()

test_calculus_2d()
