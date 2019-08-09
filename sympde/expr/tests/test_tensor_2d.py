# coding: utf-8

# TODO: - add assert to every test

import pytest

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import IndexedBase
from sympy import Matrix
from sympy import Function
from sympy import pi, cos, sin, exp
from sympy import srepr
from sympy import expand
from sympy.physics.quantum import TensorProduct

from sympde.core import Constant
from sympde.calculus import grad, dot, inner, cross, rot, curl, div
from sympde.calculus import laplace, hessian, bracket, convect
from sympde.topology import (dx, dy, dz)
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import Field, TestFunction
from sympde.topology import ProductSpace
from sympde.topology import InteriorDomain, Union
from sympde.topology import Boundary, NormalVector, TangentVector
from sympde.topology import Domain
from sympde.topology import Trace, trace_0, trace_1
from sympde.topology import Mapping
from sympde.topology import Square
from sympde.topology import ElementDomain
from sympde.topology import Area

from sympde.expr.expr import LinearExpr, BilinearExpr
from sympde.expr.expr import LinearForm, BilinearForm
from sympde.expr.expr import DomainIntegral, BoundaryIntegral
from sympde.expr.expr import Functional, Norm
from sympde.expr.expr import linearize
from sympde.expr.evaluation import TensorExpr


#==============================================================================
def test_tensorize_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = ScalarFunctionSpace('V', domain)

    u,v = [TestFunction(V, name=i) for i in ['u', 'v']]

    # ...
#    a = BilinearForm((u,v), u*v)
    a = BilinearForm((u,v), mu*u*v + dot(grad(u),grad(v)))
#    a = BilinearForm((u,v), dot(grad(u),grad(v)))
#    a = BilinearForm((u,v), dx(u)*v)
#    a = BilinearForm((u,v), laplace(u)*laplace(v))

    expr = TensorExpr(a)
    print(expr)
    # ...

#==============================================================================
def test_tensorize_2d_2():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,v = [TestFunction(V, name=i) for i in ['u', 'v']]

    # ...
#    a = BilinearForm((u,v), dot(u,v))
    a = BilinearForm((u,v), curl(u)*curl(v) + div(u)*div(v))

    expr = TensorExpr(a)
    print(expr)
    # ...

#==============================================================================
def test_tensorize_2d_1_mapping():

    DIM = 2

    M = Mapping('Map', DIM)

    domain = Domain('Omega', dim=DIM)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = ScalarFunctionSpace('V', domain)

    u,v = [TestFunction(V, name=i) for i in ['u', 'v']]

    # ...
#    a = BilinearForm((u,v), u*v)
#    a = BilinearForm((u,v), mu*u*v + dot(grad(u),grad(v)))
    a = BilinearForm((u,v), dot(grad(u),grad(v)))
#    a = BilinearForm((u,v), dx(u)*v)
#    a = BilinearForm((u,v), laplace(u)*laplace(v))

    expr = TensorExpr(a, mapping=M)
    print(expr)
    # ...

#==============================================================================
def test_tensorize_2d_2_mapping():

    DIM = 2
    M = Mapping('M', DIM)
    domain = Domain('Omega', dim=DIM)

    V = VectorFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    c = Constant('c')

    a = BilinearForm((u,v), c * div(v) * div(u) + curl(v) * curl(u))
    expr = TensorExpr(a, mapping=M)
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

#test_tensorize_2d_1()
#test_tensorize_2d_2()
#test_tensorize_2d_1_mapping()
#test_tensorize_2d_2_mapping()
