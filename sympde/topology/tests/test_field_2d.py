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
from sympde.calculus import grad, dot, inner, cross, rot, curl, div
from sympde.calculus import laplace, hessian
from sympde.topology import Domain
from sympde.topology import dx, dy, dz
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import VectorField


DIM = 2
domain = Domain('Omega', dim=DIM)

# ...
def test_field_2d_1():
    print('============ test_field_2d_1 =============')

    W = VectorFunctionSpace('W', domain)
    x,y = W.coordinates

    F = VectorField(W, 'F')

    assert( dx(F) == Matrix([[dx(F[0]), dx(F[1])]]) )

    # TODO not working yet => check it for VectorTestFunction also
#    print(dx(x*F))

    expr = inner(grad(F), grad(F))
    print(expr)


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
