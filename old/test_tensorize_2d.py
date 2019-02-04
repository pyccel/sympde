# coding: utf-8

# TODO - tensorize is not working due to the last changes on calculus

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
from sympde.calculus import laplace, hessian, bracket
from sympde.topology import (dx, dy, dz)
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import Field, VectorField
from sympde.topology import ProductSpace
from sympde.topology import TestFunction
from sympde.topology import VectorTestFunction
from sympde.topology import Unknown
from sympde.topology import InteriorDomain, Union
from sympde.topology import Boundary, NormalVector, TangentVector
from sympde.topology import Domain
from sympde.topology import Trace, trace_0, trace_1
from sympde.topology import Mapping
from sympde.topology import Square

from sympde.expr import BilinearForm, LinearForm, Integral
from sympde.expr import atomize
from sympde.expr import evaluate
from sympde.expr import tensorize
from sympde.expr import Mass, Stiffness, Advection, AdvectionT
from sympde.expr import Projection
from sympde.expr import Norm
from sympde.expr import FormCall

from sympde.expr.errors import UnconsistentError
from sympde.expr.errors import UnconsistentLhsError
from sympde.expr.errors import UnconsistentRhsError
from sympde.expr.errors import UnconsistentBCError

DIM = 2
VERBOSE = False
VERBOSE = True

##==============================================================================
#def test_tensorize_2d():
#    domain = Domain('Omega', dim=DIM)
#
#    V = FunctionSpace('V', domain)
#    U = FunctionSpace('U', domain)
#    W1 = VectorFunctionSpace('W1', domain)
#    T1 = VectorFunctionSpace('T1', domain)
#
#    v = TestFunction(V, name='v')
#    u = TestFunction(U, name='u')
#    w1 = VectorTestFunction(W1, name='w1')
#    t1 = VectorTestFunction(T1, name='t1')
#
#    x,y = domain.coordinates
#
#    alpha = Constant('alpha')
#
#    # ...
#    expr = dot(grad(v), grad(u))
#    a = BilinearForm((v,u), expr, name='a')
#    print(a)
#    print(tensorize(a))
#    print('')
#    # ...
#
#    # ...
#    expr = x*dx(v)*dx(u) + y*dy(v)*dy(u)
#    a = BilinearForm((v,u), expr, name='a')
#    print(a)
#    print(tensorize(a))
#    print('')
#    # ...
#
#    # ...
#    expr = sin(x)*dx(v)*dx(u)
#    a = BilinearForm((v,u), expr, name='a')
#    print(a)
#    print(tensorize(a))
#    print('')
#    # ...
#
#    # ...
##    expr = rot(w1)*rot(t1) + div(w1)*div(t1)
#    expr = rot(w1)*rot(t1) #+ div(w1)*div(t1)
#    a = BilinearForm((w1, t1), expr, name='a')
#    print(a)
#    print(tensorize(a))
#    print('')
#    # ...
#
##==============================================================================
#def test_tensorize_2d_stokes():
#    domain = Domain('Omega', dim=DIM)
#
#    # ... abstract model
#    V = VectorFunctionSpace('V', domain)
#    W = FunctionSpace('W', domain)
#
#    v = VectorTestFunction(V, name='v')
#    u = VectorTestFunction(V, name='u')
#    p = TestFunction(W, name='p')
#    q = TestFunction(W, name='q')
#
#    a = BilinearForm((v,u), inner(grad(v), grad(u)), name='a')
#    b = BilinearForm((v,p), div(v)*p, name='b')
#    A = BilinearForm(((v,q),(u,p)), a(v,u) - b(v,p) + b(u,q), name='A')
#    # ...
#
#    print(A)
#    print(tensorize(A))
#    print('')


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
