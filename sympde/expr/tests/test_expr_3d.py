# coding: utf-8

# TODO - split the asserts between algebraic and weak formulations ones
#      - add assert for grad in vector case
# TODO: - __call__ examples are not working anymore

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

from sympde.core import dx, dy, dz
from sympde.core import Constant
from sympde.core import Field
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import laplace
from sympde.core import hessian
from sympde.core import FunctionSpace, VectorFunctionSpace
from sympde.core import ProductSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import BilinearForm, LinearForm, Integral
from sympde.core import atomize
from sympde.core import evaluate
from sympde.core import tensorize
from sympde.core import Mass, Stiffness, Advection, AdvectionT
from sympde.core import Unknown
from sympde.core import FormCall
from sympde.core import Domain, Boundary, NormalVector, TangentVector
from sympde.core import UnionBoundary, ComplementBoundary
from sympde.core import Trace, trace_0, trace_1
from sympde.core import Equation, DirichletBC
from sympde.core import Projection
from sympde.core import Norm

from sympde.core.errors import UnconsistentError
from sympde.core.errors import UnconsistentLhsError
from sympde.core.errors import UnconsistentRhsError
from sympde.core.errors import UnconsistentBCError


DIM = 3
domain = Domain('Omega', dim=DIM)


# ...
def test_tensorize_3d():
    print('============ test_tensorize_3d =============')

    V = FunctionSpace('V', domain)
    U = FunctionSpace('U', domain)
    W1 = VectorFunctionSpace('W1', domain)
    T1 = VectorFunctionSpace('T1', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')
    w1 = VectorTestFunction(W1, name='w1')
    t1 = VectorTestFunction(T1, name='t1')

    x,y,z = domain.coordinates

    alpha = Constant('alpha')

    # ...
    expr = dot(grad(v), grad(u))
    a = BilinearForm((v,u), expr, name='a')
    print(a)
    print(tensorize(a))
    print('')
    # ...

    # ...
    expr = x*dx(v)*dx(u) + y*dy(v)*dy(u)
    a = BilinearForm((v,u), expr, name='a')
    print(a)
    print(tensorize(a))
    print('')
    # ...

    # ...
    expr = sin(x)*dx(v)*dx(u)
    a = BilinearForm((v,u), expr, name='a')
    print(a)
    print(tensorize(a))
    print('')
    # ...

    # ...
    expr = dot(curl(w1), curl(t1)) + div(w1)*div(t1)
    a = BilinearForm((w1, t1), expr, name='a')
    print(a)
    print(tensorize(a))
    print('')
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

# .....................................................
if __name__ == '__main__':

    test_tensorize_3d()

