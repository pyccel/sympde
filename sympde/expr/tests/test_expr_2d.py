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

from sympde.expr.expr import LinearExpr, BilinearExpr
from sympde.expr.expr import LinearForm, BilinearForm
from sympde.expr.expr import DomainIntegral, BoundaryIntegral
from sympde.expr.expr import Functional, Norm
from sympde.expr.expr import linearize
from sympde.expr.evaluation import TerminalExpr

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
    # TODO
#    print(l(v1+v2))
    print('')
    # ...

    # ...
    l = LinearExpr((v1,v2), x*v1 + y*v2)
    print(l)
    print(l.expr)
    print(l(u1, u2))
    # TODO
#    print(l(u1+v1, u2+v2))
    print('')
    # ...

#==============================================================================
def test_linear_expr_2d_2():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    g = Tuple(x,y)
    l = LinearExpr(v, dot(g, v))
    print(l)
    print(l.expr)
    print(l(v1))
    # TODO
#    print(l(v1+v2))
    print('')
    # ...

    # ...
    g1 = Tuple(x,0)
    g2 = Tuple(0,y)
    l = LinearExpr((v1,v2), dot(g1, v1) + dot(g2, v2))
    print(l)
    print(l.expr)
    print(l(u1, u2))
    # TODO
#    print(l(u1+v1, u2+v2))
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
    # TODO
#    print(a(u1+u2,v1+v2))
    print('')
    # ...

    # ...
    a1 = BilinearExpr((u,v), u*v)
    a2 = BilinearExpr((u,v), dot(grad(u),grad(v)))
    print(a1(u1,v1) + a2(u2,v2))
    print('')
    # ...

#==============================================================================
def test_bilinear_expr_2d_2():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    a = BilinearExpr((u,v), dot(u,v))
    print(a)
    print(a.expr)
    print(a(u1,v1))
    # TODO
#    print(a(u1+u2,v1+v2))
    print('')
    # ...

    # ...
    a1 = BilinearExpr((u,v), dot(u,v))
    a2 = BilinearExpr((u,v), inner(grad(u),grad(v)))
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
def test_linear_form_2d_2():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    g = Tuple(x,y)
    l = LinearForm(v, dot(g, v))
    print(l)
    print(l.domain)
    print('')
    # ...

    # TODO
#    # ...
#    g = Tuple(x**2, y**2)
#    l = LinearForm(v, v*trace_1(g, B1))
#    print(l)
#    print(l.domain)
#    print('')
#    # ...

    # TODO
#    # ...
#    g = Tuple(x**2, y**2)
#    l = LinearForm(v, v*trace_1(g, B1) + x*y*v)
#    print(l)
#    print('')
#    # ...

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
def test_bilinear_form_2d_2():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    a = BilinearForm((u,v), dot(u,v))
    print(a)
    print(a.domain)
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), dot(u,v) + inner(grad(u), grad(v)))
    print(a)
    print(a.domain)
    print('')
    # ...

    # TODO
#    # ...
#    a = BilinearForm((u,v), v*trace_1(grad(u), B1))
#    print(a)
#    print(a.domain)
#    print('')
#    # ...
#
#    # ...
#    a = BilinearForm((u,v), u*v + v*trace_1(grad(u), B1))
#    print(a)
#    print(a.domain)
#    print('')
#    # ...


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
    print(l)
    print(l(v1))
    print('')
    # ...

#==============================================================================
def test_call_linear_expr_2d_2():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    g = Tuple(x,y)
    l = LinearExpr(v, dot(g, v))
    print(l(v1))
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
    print(a(u1,v1))
    print('')
    # ...

#==============================================================================
def test_call_bilinear_expr_2d_2():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    a = BilinearExpr((u,v), dot(u,v))
    print(a(u1,v1))
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
    print(l(v1))
    print('')
    # ...

    # ...
    l1 = LinearForm(v1, x*y*v1)
    l = LinearForm(v, l1(v))
    print(l)
    print('')
    # ...

    # ...
    g = Tuple(x,y)
    l1 = LinearForm(v1, x*y*v1)
    l2 = LinearForm(v2, dot(grad(v2), g))

    l = LinearForm(v, l1(v) + l2(v))
    print(l)
    print('')
    # ...

#==============================================================================
def test_call_linear_form_2d_2():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    g = Tuple(x,y)
    l = LinearForm(v, dot(g, v))
    print(l(v1))
    print('')
    # ...

    # ...
    g = Tuple(x,y)
    l1 = LinearForm(v1, dot(g, v1))
    l = LinearForm(v, l1(v))
    print(l)
    print('')
    # ...

    # ...
    g1 = Tuple(x,0)
    g2 = Tuple(0,y)
    l1 = LinearForm(v1, dot(v1, g1))
    l2 = LinearForm(v2, dot(v2, g2))

    l = LinearForm(v, l1(v) + l2(v))
    print(l)
    print('')
    # ...

#==============================================================================
def test_call_linear_form_2d_3():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    V = FunctionSpace('V', domain)

    pn = Field('pn', V)
    wn = Field('wn', V)

    dp    = TestFunction(V, name='dp')
    dw    = TestFunction(V, name='dw')
    tau   = TestFunction(V, name='tau')
    sigma = TestFunction(V, name='sigma')

    Re    = Constant('Re', real=True)
    dt    = Constant('dt', real=True)
    alpha = Constant('alpha', real=True)

    l1 = LinearForm(tau, bracket(pn, wn)*tau - 1./Re * dot(grad(tau), grad(wn)))

    l = LinearForm((tau, sigma), dt*l1(tau))
    print(l)

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
    print(a(u1,v1))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), u1*v1)
    a = BilinearForm((u,v), a1(u,v))
    print(a)
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), u1*v1)
    a2 = BilinearForm((u2,v2), dot(grad(u2), grad(v2)))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))
    print(a)
    print('')
    # ...

#==============================================================================
def test_call_bilinear_form_2d_2():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    a = BilinearForm((u,v), dot(u,v))
    print(a(u1,v1))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), dot(u1,v1))
    a = BilinearForm((u,v), a1(u,v))
    print(a)
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), dot(u1,v1))
    a2 = BilinearForm((u2,v2), inner(grad(u2), grad(v2)))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))
    print(a)
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

    # ...
    l1 = LinearForm(v1, x*y*v1)
    l = LinearForm(v, l1(v))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    g = Tuple(x,y)
    l1 = LinearForm(v1, x*y*v1)
    l2 = LinearForm(v2, dot(grad(v2), g))

    l = LinearForm(v, l1(v) + l2(v))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    l1 = LinearForm(v1, x*y*v1)
    l2 = LinearForm(v1, v1)
    l = LinearForm(v, l1(v) + kappa*l2(v))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    g = Tuple(x**2, y**2)
    l1 = LinearForm(v1, x*y*v1)
    l2 = LinearForm(v1, v1)
    l3 = LinearForm(v, v*trace_1(g, B1))
    l = LinearForm(v, l1(v) + kappa*l2(v) + mu*l3(v))
    print(TerminalExpr(l))
    print('')
    # ...

#==============================================================================
def test_terminal_expr_linear_2d_2():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    g = Tuple(x,y)
    l = LinearForm(v, dot(g, v))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    g = Tuple(x,y)
    l = LinearForm(v, dot(g, v) + div(v))
    print(TerminalExpr(l))
    print('')
    # ...

    # TODO
#    # ...
#    g = Tuple(x**2, y**2)
#    l = LinearForm(v, v*trace_1(g, B1))
#    print(TerminalExpr(l))
#    print('')
#    # ...
#
#    # ...
#    g = Tuple(x**2, y**2)
#    l = LinearForm(v, v*trace_1(g, B1) + x*y*v)
#    print(TerminalExpr(l))
#    print('')
#    # ...
#
#    # ...
#    l1 = LinearForm(v1, x*y*v1)
#    l = LinearForm(v, l1(v))
#    print(TerminalExpr(l))
#    print('')
#    # ...
#
#    # ...
#    g = Tuple(x,y)
#    l1 = LinearForm(v1, x*y*v1)
#    l2 = LinearForm(v2, dot(grad(v2), g))
#
#    l = LinearForm(v, l1(v) + l2(v))
#    print(TerminalExpr(l))
#    print('')
#    # ...
#
#    # ...
#    l1 = LinearForm(v1, x*y*v1)
#    l2 = LinearForm(v1, v1)
#    l = LinearForm(v, l1(v) + kappa*l2(v))
#    print(TerminalExpr(l))
#    print('')
#    # ...
#
#    # ...
#    g = Tuple(x**2, y**2)
#    l1 = LinearForm(v1, x*y*v1)
#    l2 = LinearForm(v1, v1)
#    l3 = LinearForm(v, v*trace_1(g, B1))
#    l = LinearForm(v, l1(v) + kappa*l2(v) + mu*l3(v))
#    print(TerminalExpr(l))
#    print('')
#    # ...

#==============================================================================
def test_terminal_expr_linear_2d_3():

    domain = Square()

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    l = LinearForm(v, trace_1(grad(v), domain.boundary))
    print(TerminalExpr(l))
    print('')
    # ...

#==============================================================================
def test_terminal_expr_linear_2d_4():

    D1 = InteriorDomain('D1', dim=2)
    D2 = InteriorDomain('D2', dim=2)
    domain = Union(D1, D2)

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

#==============================================================================
def test_terminal_expr_linear_2d_5(boundary=['Gamma_1', 'Gamma_3']):

    # ... abstract model
    domain = Square()

    V = FunctionSpace('V', domain)

    B_neumann = [domain.get_boundary(i) for i in boundary]
    if len(B_neumann) == 1:
        B_neumann = B_neumann[0]

    else:
        B_neumann = Union(*B_neumann)

    x,y = domain.coordinates

    F = Field('F', V)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    expr = dot(grad(v), grad(u))
    a = BilinearForm((v,u), expr)

    solution = cos(0.5*pi*x)*cos(0.5*pi*y)
    f        = (1./2.)*pi**2*solution

    expr = f*v
    l0 = LinearForm(v, expr)

    expr = v*trace_1(grad(solution), B_neumann)
    l_B_neumann = LinearForm(v, expr)

    expr = l0(v) + l_B_neumann(v)
    l = LinearForm(v, expr)

    print(TerminalExpr(l))
    print('')
    # ...

#==============================================================================
def test_terminal_expr_bilinear_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    eps   = Constant('eps', real=True)

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

    # ...
    a1 = BilinearForm((u1,v1), u1*v1)
    a = BilinearForm((u,v), a1(u,v))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), u1*v1)
    a2 = BilinearForm((u2,v2), dot(grad(u2), grad(v2)))
    a = BilinearForm((u,v), a1(u,v) + a2(u,v))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), u1*v1)
    a2 = BilinearForm((u2,v2), dot(grad(u2), grad(v2)))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), u1*v1)
    a2 = BilinearForm((u2,v2), dot(grad(u2), grad(v2)))
    a3 = BilinearForm((u,v), v*trace_1(grad(u), B1))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v) + mu*a3(u,v))
    print(TerminalExpr(a))
    print('')
    # ...

    # ... Poisson with Nitsch method
    a0 = BilinearForm((u,v), dot(grad(u),grad(v)))
    a_B1 = BilinearForm((u,v), - kappa * u*trace_1(grad(v), B1)
                               - v*trace_1(grad(u), B1)
                               + trace_0(u, B1) * trace_0(v, B1) / eps)
    a = BilinearForm((u,v), a0(u,v) + a_B1(u,v))
    print(TerminalExpr(a))
    print('')
    # ...

#==============================================================================
def test_terminal_expr_bilinear_2d_2():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    a = BilinearForm((u,v), dot(u,v))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), inner(grad(u),grad(v)))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), dot(u,v) + inner(grad(u),grad(v)))
    print(TerminalExpr(a))
    print('')
    # ...

    # TODO
#    # ...
#    a = BilinearForm((u,v), u*v + dot(grad(u),grad(v)) + v*trace_1(grad(u), B1) )
#    print(TerminalExpr(a))
#    print('')
#    # ...
#
#    # ...
#    a = BilinearForm(((u1,u2),(v1,v2)), u1*v1 + u2*v2)
#    print(TerminalExpr(a))
#    print('')
#    # ...
#
#    # ...
#    a1 = BilinearForm((u1,v1), u1*v1)
#    a = BilinearForm((u,v), a1(u,v))
#    print(TerminalExpr(a))
#    print('')
#    # ...
#
#    # ...
#    a1 = BilinearForm((u1,v1), u1*v1)
#    a2 = BilinearForm((u2,v2), dot(grad(u2), grad(v2)))
#    a = BilinearForm((u,v), a1(u,v) + a2(u,v))
#    print(TerminalExpr(a))
#    print('')
#    # ...
#
#    # ...
#    a1 = BilinearForm((u1,v1), u1*v1)
#    a2 = BilinearForm((u2,v2), dot(grad(u2), grad(v2)))
#    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))
#    print(TerminalExpr(a))
#    print('')
#    # ...
#
#    # ...
#    a1 = BilinearForm((u1,v1), u1*v1)
#    a2 = BilinearForm((u2,v2), dot(grad(u2), grad(v2)))
#    a3 = BilinearForm((u,v), v*trace_1(grad(u), B1))
#    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v) + mu*a3(u,v))
#    print(TerminalExpr(a))
#    print('')
#    # ...

#==============================================================================
# stokes
def test_terminal_expr_bilinear_2d_3():

    domain = Domain('Omega', dim=2)

    x,y = domain.coordinates

    V = VectorFunctionSpace('V', domain)
    W = FunctionSpace('W', domain)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')
    p = TestFunction(W, name='p')
    q = TestFunction(W, name='q')

    a = BilinearForm((u,v), inner(grad(v), grad(u)))
    b = BilinearForm((v,p), div(v)*p)
    A = BilinearForm(((u,p),(v,q)), a(v,u) - b(v,p) + b(u,q))

    print(TerminalExpr(A))
    print('')


#==============================================================================
def test_linearize_expr_2d_1():
    domain = Domain('Omega', dim=DIM)
    x,y = domain.coordinates

    V1 = FunctionSpace('V1', domain)
    W1 = VectorFunctionSpace('W1', domain)

    v1 = TestFunction(V1, name='v1')
    w1 = VectorTestFunction(W1, name='w1')

    alpha = Constant('alpha')

    F = Field('F', space=V1)
    G = VectorField(W1, 'G')

    # ...
    l = LinearExpr(v1, F**2*v1)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearExpr(v1, dot(grad(F), grad(F))*v1)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearExpr(v1, exp(-F)*v1)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearExpr(v1, cos(F)*v1)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearExpr(v1, cos(F**2)*v1)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearExpr(v1, F**2*dot(grad(F), grad(v1)))
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearExpr(w1, dot(rot(G), grad(G))*w1)
    a = linearize(l, G, trials='u1')
    print(a)
    # ...

#==============================================================================
def test_linearize_expr_2d_2():
    domain = Domain('Omega', dim=DIM)
    x,y = domain.coordinates

    V1 = FunctionSpace('V1', domain)

    v1 = TestFunction(V1, name='v1')

    alpha = Constant('alpha')

    F = Field('F', space=V1)
    G = Field('G', space=V1)

    # ...
    l1 = LinearExpr(v1, F**2*v1)
    l  = LinearExpr(v1, l1(v1))

    a = linearize(l, F, trials='u1')
    print(a)

    expected = linearize(l1, F, trials='u1')
    assert( linearize(l, F, trials='u1') == expected )
    # ...

#==============================================================================
def test_linearize_form_2d_1():
    domain = Domain('Omega', dim=DIM)
    x,y = domain.coordinates

    V1 = FunctionSpace('V1', domain)
    W1 = VectorFunctionSpace('W1', domain)

    v1 = TestFunction(V1, name='v1')
    w1 = VectorTestFunction(W1, name='w1')

    alpha = Constant('alpha')

    F = Field('F', space=V1)
    G = VectorField(W1, 'G')

    # ...
    l = LinearForm(v1, F**2*v1)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, dot(grad(F), grad(F))*v1)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, exp(-F)*v1)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, cos(F)*v1)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, cos(F**2)*v1)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, F**2*dot(grad(F), grad(v1)))
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(w1, dot(rot(G), grad(G))*w1)
    a = linearize(l, G, trials='u1')
    print(a)
    # ...

#==============================================================================
def test_linearize_form_2d_2():
    domain = Domain('Omega', dim=DIM)
    x,y = domain.coordinates

    V1 = FunctionSpace('V1', domain)

    v1 = TestFunction(V1, name='v1')

    alpha = Constant('alpha')

    F = Field('F', space=V1)
    G = Field('G', space=V1)

    # ...
    l1 = LinearForm(v1, F**2*v1)
    l  = LinearForm(v1, l1(v1))

    a = linearize(l, F, trials='u1')
    print(a)

    expected = linearize(l1, F, trials='u1')
    assert( linearize(l, F, trials='u1') == expected )
    # ...

#==============================================================================
def test_linearize_form_2d_3():
    """steady Euler equation."""
    domain = Domain('Omega', dim=DIM)
    x,y = domain.coordinates

    U = VectorFunctionSpace('U', domain)
    W =       FunctionSpace('W', domain)

    v   = VectorTestFunction(U, name='v')
    phi =       TestFunction(W, name='phi')
    q   =       TestFunction(W, name='q')

    U_0   = VectorField(U, name='U_0')
    Rho_0 =       Field('Rho_0', W)
    P_0   =       Field('P_0', W)

    # ...
    expr = div(Rho_0*U_0) * phi
    l1 = LinearForm(phi, expr)

    expr = Rho_0*dot(convect(U_0, grad(U_0)), v) + dot(grad(P_0), v)
    l2 = LinearForm(v, expr)

    expr = dot(U_0, grad(P_0)) * q + P_0 * div(U_0) * q
    l3 = LinearForm(q, expr)
    # ...

    a1 = linearize(l1, [Rho_0, U_0], trials=['d_rho', 'd_u'])
    print(a1)
    print('')

    a2 = linearize(l2, [Rho_0, U_0, P_0], trials=['d_rho', 'd_u', 'd_p'])
    print(a2)
    print('')

    a3 = linearize(l3, [P_0, U_0], trials=['d_p', 'd_u'])
    print(a3)
    print('')

    l = LinearForm((phi, v, q), l1(phi) + l2(v) + l3(q))
    a = linearize(l, [Rho_0, U_0, P_0], trials=['d_rho', 'd_u', 'd_p'])
    print(a)

#==============================================================================
def test_area_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    mu    = Constant('mu'   , is_real=True)

    e = ElementDomain(domain)
    area = Area(e)

    V = FunctionSpace('V', domain)

    u,v = [TestFunction(V, name=i) for i in ['u', 'v']]

    # ...
    a = BilinearForm((v,u), area * u * v)
    print(TerminalExpr(a))
    # ...

#==============================================================================
def test_stabilization_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    b1 = 1.
    b2 = 0.
    b = Tuple(b1, b2)

    # right hand side
    f = x*y

    e = ElementDomain()
    area = Area(e)

    V = FunctionSpace('V', domain)

    u,v = [TestFunction(V, name=i) for i in ['u', 'v']]

    # ...
    expr = kappa * dot(grad(u), grad(v)) + dot(b, grad(u)) * v
    a = BilinearForm((v,u), expr)
    # ...

    # ...
    expr = f * v
    l = LinearForm(v, expr)
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * dot(b, grad(v))
    s1 = BilinearForm((v,u), expr)

    expr = - f * dot(b, grad(v))
    l1 = LinearForm(v, expr)
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * ( dot(b, grad(v)) - kappa * laplace(v))
    s2 = BilinearForm((v,u), expr)

    expr = - f * ( dot(b, grad(v)) - kappa * laplace(v))
    l2 = LinearForm(v, expr)
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * ( dot(b, grad(v)) + kappa * laplace(v))
    s3 = BilinearForm((v,u), expr)

    expr = - f * ( dot(b, grad(v)) + kappa * laplace(v))
    l3 = LinearForm(v, expr)
    # ...

    # ...
    expr = a(v,u) + mu*area*s1(v,u)
    a1 = BilinearForm((v,u), expr)
    # ...

    # ...
    expr = a(v,u) + mu*area*s2(v,u)
    a2 = BilinearForm((v,u), expr)
    # ...

    # ...
    expr = a(v,u) + mu*area*s3(v,u)
    a3 = BilinearForm((v,u), expr)
    # ...

    print(a1)
    print(TerminalExpr(a1))
    print('')

    print(a2)
    print(TerminalExpr(a2))
    print('')

    print(a3)
    print(TerminalExpr(a3))
    print('')

#==============================================================================
def test_user_function_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    # right hand side
    f = Function('f')

    V = FunctionSpace('V', domain)

    u,v = [TestFunction(V, name=i) for i in ['u', 'v']]

    # ...
    expr = dot(grad(u), grad(v)) + f(x,y) * u * v
    a = BilinearForm((v,u), expr)

    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    expr = f(x,y) * v
    l = LinearForm(v, expr)

    print(l)
    print(TerminalExpr(l))
    print('')
    # ...

#==============================================================================
def test_functional_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = FunctionSpace('V', domain)
    F = Field('F', space=V)

    # ...
    expr = x*y
    a = Functional(expr, domain)

    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    expr = F - cos(2*pi*x)*cos(3*pi*y)
    expr = dot(grad(expr), grad(expr))
    a = Functional(expr, domain)

    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

#==============================================================================
def test_norm_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    V = FunctionSpace('V', domain)
    F = Field('F', space=V)

    # ...
    expr = x*y
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u))
    print('> h1 norm = ', TerminalExpr(h1_norm_u))
    print('')
    # ...

    # ...
    expr = sin(pi*x)*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u))
    print('> h1 norm = ', TerminalExpr(h1_norm_u))
    print('')
    # ...

    # ...
    expr = F-x*y
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u))
    print('> h1 norm = ', TerminalExpr(h1_norm_u))
    print('')
    # ...

    # ...
    expr = F-sin(pi*x)*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u))
    print('> h1 norm = ', TerminalExpr(h1_norm_u))
    print('')
    # ...

    # ...
    expr = F-sin(0.5*pi*(1.-x))*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u))
    print('> h1 norm = ', TerminalExpr(h1_norm_u))
    print('')
    # ...

    # ...
    expr = F-cos(0.5*pi*x)*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u))
    print('> h1 norm = ', TerminalExpr(h1_norm_u))
    print('')
    # ...

#==============================================================================
def test_norm_2d_2():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    V = VectorFunctionSpace('V', domain)
    F = VectorField(V, 'F')

    # ...
    f = Tuple(sin(pi*x)*sin(pi*y), sin(pi*x)*sin(pi*y))
    expr = Matrix([F[0]-f[0], F[1]-f[1]])
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u))
    print('> h1 norm = ', TerminalExpr(h1_norm_u))
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
