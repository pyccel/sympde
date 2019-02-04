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
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()

#test_linear_expr_2d_1()
#test_linear_expr_2d_2()
#print('================')
#test_bilinear_expr_2d_1()
#test_bilinear_expr_2d_2()
#print('================')
#test_integral_2d_1()
#print('================')
#test_linear_form_2d_1()
#test_linear_form_2d_2()
#print('================')
#test_bilinear_form_2d_1()
#test_bilinear_form_2d_2()
#print('================')
#test_call_linear_expr_2d_1()
#test_call_linear_expr_2d_2()
#print('================')
#test_call_bilinear_expr_2d_1()
#test_call_bilinear_expr_2d_2()
#print('================')
#test_call_linear_form_2d_1()
#test_call_linear_form_2d_2()
#test_call_linear_form_2d_3()
#print('================')
#test_call_bilinear_form_2d_1()
#test_call_bilinear_form_2d_2()
#print('================')
#test_terminal_expr_linear_2d_1()
#test_terminal_expr_linear_2d_2()
#print('================')
test_terminal_expr_bilinear_2d_1()
#test_terminal_expr_bilinear_2d_2()
#test_terminal_expr_bilinear_2d_3()
