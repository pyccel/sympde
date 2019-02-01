# coding: utf-8

import pytest

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import IndexedBase
from sympy import Matrix
from sympy import Function
from sympy import pi, cos, sin, exp
from sympy import srepr
from sympy.physics.quantum import TensorProduct

from sympde.core import Constant
from sympde.calculus import grad, dot, inner, cross, rot, curl, div
from sympde.calculus import laplace, hessian, bracket, convect
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
from sympde.topology import ElementDomain
from sympde.topology import Area

from sympde.expr import LinearExpr, BilinearExpr
from sympde.expr.expr import LinearForm, BilinearForm
from sympde.expr.expr import DomainIntegral, BoundaryIntegral
from sympde.expr.expr import Call
from sympde.expr.expr import TerminalExpr
from sympde.expr import atomize
from sympde.expr import evaluate
from sympde.expr import Mass, Stiffness, Advection, AdvectionT
from sympde.expr import Projection
from sympde.expr import Norm
from sympde.expr import FormCall
from sympde.expr import is_linear_form, is_bilinear_form
from sympde.expr import linearize

from sympde.expr.errors import UnconsistentError
from sympde.expr.errors import UnconsistentLinearExpressionError
from sympde.expr.errors import UnconsistentLhsError
from sympde.expr.errors import UnconsistentRhsError
from sympde.expr.errors import UnconsistentBCError

DIM = 2
VERBOSE = False
VERBOSE = True

#==============================================================================
def test_linear_expr_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    l = LinearExpr(v, x*y*v)
    print(l)
    print(l.expr)
    print(l(v1))
    print(l(v1+v2))
    print('')
    # ...

    # ...
    l = LinearExpr((v1,v2), x*v1 + y*v2)
    print(l)
    print(l.expr)
    print(l(u1, u2))
    print(l(u1+v1, u2+v2))
    print('')
    # ...

#==============================================================================
def test_bilinear_expr_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    a = BilinearExpr((u,v), u*v)
    print(a)
    print(a.expr)
    print(a(u1,v1))
    print(a(u1+u2,v1+v2))
    print('')
    # ...

    # ...
    a1 = BilinearExpr((u,v), u*v)
    a2 = BilinearExpr((u,v), dot(grad(u),grad(v)))
    print(a1(u1,v1) + a2(u2,v2))
    print('')
    # ...

#==============================================================================
def test_integral_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    a = BilinearExpr((u,v), u*v)
    print(a)
    a = DomainIntegral(a)
    print(a)
    print('')
    # ...

    # ...
    l = LinearExpr(v, x*y*v)
    print(l)
    l = DomainIntegral(l)
    print(l)
    # ...

#==============================================================================
def test_bilinear_form_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    a = BilinearForm((u,v), u*v)
    print(a)
    print(a.domain)
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), u*v + dot(grad(u), grad(v)))
    print(a)
    print(a.domain)
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), v*trace_1(grad(u), B1))
    print(a)
    print(a.domain)
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), u*v + v*trace_1(grad(u), B1))
    print(a)
    print(a.domain)
    print('')
    # ...

#==============================================================================
def test_linear_form_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    l = LinearForm(v, x*y*v)
    print(l)
    print(l.domain)
    print('')
    # ...

    # ...
    g = Tuple(x**2, y**2)
    l = LinearForm(v, v*trace_1(g, B1))
    print(l)
    print(l.domain)
    print('')
    # ...

    # ...
    g = Tuple(x**2, y**2)
    l = LinearForm(v, v*trace_1(g, B1) + x*y*v)
    print(l)
    print('')
    # ...


#==============================================================================
def test_call_linear_expr_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    l = LinearExpr(v, x*y*v)
    print(Call(l, v1))
    print(Call(l, v1, evaluate=True))
    print(l(v1))
    print(l(v1, evaluate=True))
    print('')
    # ...

#==============================================================================
def test_call_bilinear_expr_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    a = BilinearExpr((u,v), u*v)
    print(Call(a, (u1,v1)))
    print(Call(a, (u1,v1), evaluate=True))
    print(a(u1,v1))
    print(a(u1,v1, evaluate=True))
    print('')
    # ...

#==============================================================================
def test_call_linear_form_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    l = LinearForm(v, x*y*v)
    print(Call(l, v1))
    print(Call(l, v1, evaluate=True))
    print(l(v1))
    print(l(v1, evaluate=True))
    print('')
    # ...

#==============================================================================
def test_call_bilinear_form_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    a = BilinearForm((u,v), u*v)
    print(Call(a, (u1,v1)))
    print(Call(a, (u1,v1), evaluate=True))
    print(a(u1,v1))
    print(a(u1,v1, evaluate=True))
    print('')
    # ...

#==============================================================================
def test_terminal_expr_bilinear_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    a = BilinearForm((u,v), u*v)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), dot(grad(u),grad(v)))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), u*v + dot(grad(u),grad(v)))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), u*v + dot(grad(u),grad(v)) + v*trace_1(grad(u), B1) )
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm(((u1,u2),(v1,v2)), u1*v1 + u2*v2)
    print(TerminalExpr(a))
    print('')
    # ...

#==============================================================================
def test_terminal_expr_linear_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    l = LinearForm(v, x*y*v)
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    l = LinearForm(v, x*y*v + v)
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    g = Tuple(x**2, y**2)
    l = LinearForm(v, v*trace_1(g, B1))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    g = Tuple(x**2, y**2)
    l = LinearForm(v, v*trace_1(g, B1) + x*y*v)
    print(TerminalExpr(l))
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

#test_linear_expr_2d_1()
#print('================')
#test_bilinear_expr_2d_1()
#print('================')
#test_integral_2d_1()
#print('================')
#test_bilinear_form_2d_1()
#print('================')
#test_linear_form_2d_1()
#print('================')
#test_call_linear_expr_2d_1()
#print('================')
#test_call_bilinear_expr_2d_1()
#print('================')
#test_call_linear_form_2d_1()
#print('================')
#test_call_bilinear_form_2d_1()
#print('================')
test_terminal_expr_linear_2d_1()
print('================')
test_terminal_expr_bilinear_2d_1()
