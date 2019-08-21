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
from sympy.physics.quantum import TensorProduct

from sympde.core import Constant
from sympde.calculus import grad, dot, inner, cross, rot, curl, div
from sympde.calculus import laplace, hessian, bracket, convect
from sympde.calculus import jump, avg, Dn, minus, plus
from sympde.topology import (dx, dy, dz)
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import ProductSpace
from sympde.topology import element_of, elements_of
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
from sympde.expr.expr import integral, is_linear_form, is_bilinear_form
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

    V = ScalarFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

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

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

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

    V = ScalarFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

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

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

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
def test_linear_form_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    nn    = NormalVector('nn')

    V = ScalarFunctionSpace('V', domain)
    W = VectorFunctionSpace('W', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]
    # ...
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    l = LinearForm(v, int_0(x*y*v))

    assert(l.domain == domain)
    assert(l(v1) == int_0(x*y*v1))
    # ...

    # ...
    g  = Tuple(x**2, y**2)
    l  = LinearForm(v, int_1(v*dot(g, nn)))
    print(l)

    assert(l.domain == B1)
    assert(l(v1) == int_1(v1*trace_1(g, B1)))
    # ...

    # ...
    g = Tuple(x**2, y**2)
    l = LinearForm(v, int_1(v*dot(g, nn)) + int_0(x*y*v))

    assert(len(l.domain.args) == 2)
    for i in l.domain.args:
        assert(i in [domain, B1])

    assert(l(v1) == int_1(v1*trace_1(g, B1)) + int_0(x*y*v1))
    # ...

    # ...
    l1 = LinearForm(v1, int_0(x*y*v1))
    l = LinearForm(v, l1(v))

    assert(l.domain == domain)
    assert(l(u1) == int_0(x*y*u1))
    # ...

    # ...
    g = Tuple(x,y)
    l1 = LinearForm(u1, int_0(x*y*u1))
    l2 = LinearForm(u2, int_0(dot(grad(u2), g)))

    l = LinearForm(v, l1(v) + l2(v))

    assert(l.domain == domain)
    assert(l(v1) == int_0(x*y*v1) + int_0(dot(grad(v1), g)))
    # ...

    # ...
    pn, wn = [element_of(V, name=i) for i in ['pn', 'wn']]

    tau   = element_of(V, name='tau')
    sigma = element_of(V, name='sigma')

    Re    = Constant('Re', real=True)
    dt    = Constant('dt', real=True)
    alpha = Constant('alpha', real=True)

    l1 = LinearForm(tau, int_0(bracket(pn, wn)*tau - 1./Re * dot(grad(tau), grad(wn))))



    l = LinearForm((tau, sigma), dt*l1(tau))

    assert(l.domain == domain)
    assert(l(u1,u2) == int_0(-1.0*dt*dot(grad(u1), grad(wn))/Re) +
                       int_0(dt*u1*bracket(pn, wn)))
    # ...

#==============================================================================
def test_linear_form_2d_2():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    g = Tuple(x,y)
    l = LinearForm(v, int_0(dot(g, v)))

    assert(l.domain == domain)
    assert(l(v1) == int_0(dot(g, v1)))
    # ...

    # ...
    g = Tuple(x,y)
    l1 = LinearForm(v1, int_0(dot(g, v1)))
    l = LinearForm(v, l1(v))

    assert(l.domain == domain)
    assert(l(u1) == int_0(dot(g, u1)))
    # ...

    # ...
    g1 = Tuple(x,0)
    g2 = Tuple(0,y)
    l1 = LinearForm(v1, int_0(dot(v1, g1)))
    l2 = LinearForm(v2, int_0(dot(v2, g2)))

    l = LinearForm(v, l1(v) + l2(v))

    assert(l.domain == domain)
    assert(l(u) == int_0(dot(u, g1)) + int_0(dot(u, g2)))
    # ...

#==============================================================================
def test_bilinear_form_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    nn    = NormalVector('nn')

    V = ScalarFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    a = BilinearForm((u,v), int_0(u*v))

    assert(a.domain == domain)
    assert(a(u1,v1) == int_0(u1*v1))
    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v + dot(grad(u), grad(v))))

    assert(a.domain == domain)
    assert(a(u1,v1) == int_0(u1*v1) + int_0(dot(grad(u1), grad(v1))))
    # ...

    # ...
    a = BilinearForm((u,v), int_1(v*dot(grad(u), nn)))

    assert(a.domain == B1)
    assert(a(u1,v1) == int_1(v1*trace_1(grad(u1), B1)))
    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v) + int_1(v*dot(grad(u), nn)))

    # TODO a.domain are not ordered
    assert(len(a.domain.args) == 2)
    for i in a.domain.args:
        assert(i in [domain, B1])

    assert(a(u1,v1) == int_0(u1*v1) + int_1(v1*trace_1(grad(u1), B1)))
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a = BilinearForm((u,v), a1(u,v))

    assert(a.domain == domain)
    assert(a(u2,v2) == int_0(u2*v2))
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))

    assert(a.domain == domain)
    assert(a(u,v) == int_0(u*v) + int_0(kappa*dot(grad(u), grad(v))))
    # ...

#==============================================================================
def test_bilinear_form_2d_2():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    a = BilinearForm((u,v), int_0(dot(u,v)))

    assert(a.domain == domain)
    assert(a(u1,v1) == int_0(dot(u1,v1)))
    # ...

    # ...
    a = BilinearForm((u,v), int_0(dot(u,v) + inner(grad(u), grad(v))))

    assert(a.domain == domain)
    assert(a(u1,v1) == int_0(dot(u1,v1)) + int_0(inner(grad(u1), grad(v1))))
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(dot(u1,v1)))
    a = BilinearForm((u,v), a1(u,v))

    assert(a.domain == domain)
    assert(a(u2,v2) == int_0(dot(u2,v2)))
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(dot(u1,v1)))
    a2 = BilinearForm((u2,v2), int_0(inner(grad(u2), grad(v2))))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))

    assert(a.domain == domain)
    assert(a(u,v) == int_0(dot(u,v)) + int_0(kappa*inner(grad(u), grad(v))))
    # ...

#==============================================================================
def test_terminal_expr_linear_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    nn    = NormalVector('nn')

    V = ScalarFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    l = LinearForm(v, int_0(x*y*v))

    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    l = LinearForm(v, int_0(x*y*v + v))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    g = Tuple(x**2, y**2)
    l = LinearForm(v, int_1(v*dot(g, nn)))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    g = Tuple(x**2, y**2)
    l = LinearForm(v, int_1(v*dot(g, nn)) + int_0(x*y*v))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    l1 = LinearForm(v1, int_0(x*y*v1))
    l = LinearForm(v, l1(v))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    g = Tuple(x,y)
    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v2, int_0(dot(grad(v2), g)))

    l = LinearForm(v, l1(v) + l2(v))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v1, int_0(v1))
    l = LinearForm(v, l1(v) + kappa*l2(v))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    g = Tuple(x**2, y**2)
    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v1, int_0(v1))
    l3 = LinearForm(v, int_1(v*dot(g, nn)))
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

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    g = Tuple(x,y)
    l = LinearForm(v, int_0(dot(g, v)))
    print(TerminalExpr(l))
    print('')
    # ...

    # ...
    g = Tuple(x,y)
    l = LinearForm(v, int_0(dot(g, v) + div(v)))
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
    B      = domain.boundary

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    nn    = NormalVector('nn')

    V = ScalarFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B, expr)

    l = LinearForm(v, int_1(dot(grad(v), nn)))
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

    V = ScalarFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    int_0 = lambda expr: integral(domain , expr)

    l = LinearForm(v, int_0(x*y*v))
    print(TerminalExpr(l))
    print('')
    # ...
#==============================================================================
def test_terminal_expr_linear_2d_5(boundary=[r'\Gamma_1', r'\Gamma_3']):

    # ... abstract model
    domain = Square()

    V = ScalarFunctionSpace('V', domain)

    B_neumann = [domain.get_boundary(i) for i in boundary]
    if len(B_neumann) == 1:
        B_neumann = B_neumann[0]

    else:
        B_neumann = Union(*B_neumann)

    x,y = domain.coordinates
    nn    = NormalVector('nn')

    F = element_of(V, name='F')

    v = element_of(V, name='v')
    u = element_of(V, name='u')

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B_neumann , expr)

    expr = dot(grad(v), grad(u))
    a = BilinearForm((v,u), int_0(expr))

    solution = cos(0.5*pi*x)*cos(0.5*pi*y)
    f        = (1./2.)*pi**2*solution

    expr = f*v
    l0 = LinearForm(v, int_0(expr))

    expr = v*dot(grad(solution), nn)
    l_B_neumann = LinearForm(v, int_1(expr))

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
    nn    = NormalVector('nn')

    V = ScalarFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    a = BilinearForm((u,v), int_0(u*v))
    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), int_0(dot(grad(u),grad(v))))
    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v + dot(grad(u),grad(v))))
    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v + dot(grad(u),grad(v))) + int_1(v*dot(grad(u), nn)) )
    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm(((u1,u2),(v1,v2)), int_0(u1*v1 + u2*v2))
    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a = BilinearForm((u,v), a1(u,v))
    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a = BilinearForm((u,v), a1(u,v) + a2(u,v))
    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))
    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a3 = BilinearForm((u,v), int_1(v*dot(grad(u), nn)))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v) + mu*a3(u,v))
    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ... Poisson with Nitsch method
    a0 = BilinearForm((u,v), int_0(dot(grad(u),grad(v))))
    a_B1 = BilinearForm((u,v), int_1(- kappa * u*dot(grad(v), nn)
                               - v*dot(grad(u), nn)
                               + u*v / eps))
    a = BilinearForm((u,v), a0(u,v) + a_B1(u,v))
    print(a)
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
    nn    = NormalVector('nn')

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    a = BilinearForm((u,v), int_0(dot(u,v)))
    print(TerminalExpr(a))
    print('')

    # ...
    a = BilinearForm((u,v), int_0(inner(grad(u),grad(v))))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), int_0(dot(u,v) + inner(grad(u),grad(v))))
    print(TerminalExpr(a))
    print('')
    # ...

def test_terminal_expr_bilinear_2d_3():

    domain = Square()

    V = ScalarFunctionSpace('V', domain)

    B = domain.boundary

    v = element_of(V, name='v')
    u = element_of(V, name='u')

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    nn    = NormalVector('nn')

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B, expr)

    # nitsche
    a0  = BilinearForm((u,v), int_0(dot(grad(v),grad(u))))

    a_B = BilinearForm((u,v), int_1(-u*dot(grad(v), nn) \
                              -v*dot(grad(u), nn) \
                              +kappa*u*v))

    a = BilinearForm((u,v), a0(u,v) + a_B(u,v))

    print(TerminalExpr(a))
    print('')

    a = BilinearForm((u,v), int_0(u*v + dot(grad(u),grad(v))) + int_1(v*dot(grad(u), nn)) )
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u,v), int_0(u*v))
    a = BilinearForm((u,v), a1(u,v))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u,v), int_0(u*v))
    a2 = BilinearForm((u,v), int_0(dot(grad(u), grad(v))))
    a = BilinearForm((u,v), a1(u,v) + a2(u,v))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u,v), int_0(u*v))
    a2 = BilinearForm((u,v), int_0(dot(grad(u), grad(v))))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u,v), int_0(u*v))
    a2 = BilinearForm((u,v), int_0(dot(grad(u), grad(v))))
    a3 = BilinearForm((u,v), int_1(v*dot(grad(u), nn)))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v) + mu*a3(u,v))
    print(TerminalExpr(a))
    print('')
#    # ...

#==============================================================================


def test_terminal_expr_bilinear_2d_4():

    domain = Domain('Omega', dim=2)

    x,y = domain.coordinates

    V = VectorFunctionSpace('V', domain)
    W = ScalarFunctionSpace('W', domain)

    v = element_of(V, name='v')
    u = element_of(V, name='u')
    p = element_of(W, name='p')
    q = element_of(W, name='q')

    int_0 = lambda expr: integral(domain , expr)

    # stokes
    a = BilinearForm((u,v), int_0(inner(grad(v), grad(u))))
    b = BilinearForm((v,p), int_0(div(v)*p))
    a = BilinearForm(((u,p),(v,q)), a(v,u) - b(v,p) + b(u,q))

    print(TerminalExpr(a))
    print('')

#==============================================================================
def test_linearize_expr_2d_1():
    domain = Domain('Omega', dim=DIM)
    x,y = domain.coordinates

    V1 = ScalarFunctionSpace('V1', domain)
    W1 = VectorFunctionSpace('W1', domain)

    v1 = element_of(V1, name='v1')
    w1 = element_of(W1, name='w1')

    alpha = Constant('alpha')

    F = element_of(V1, name='F')
    G = element_of(W1, 'G')


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

    V1 = ScalarFunctionSpace('V1', domain)

    v1 = element_of(V1, name='v1')

    alpha = Constant('alpha')

    F = element_of(V1, name='F')
    G = element_of(V1, name='G')

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

    V1 = ScalarFunctionSpace('V1', domain)
    W1 = VectorFunctionSpace('W1', domain)

    v1 = element_of(V1, name='v1')
    w1 = element_of(W1, name='w1')

    alpha = Constant('alpha')

    F = element_of(V1, name='F')
    G = element_of(W1, 'G')

    int_0 = lambda expr: integral(domain , expr)

    # ...
    l = LinearForm(v1, int_0(F**2*v1))
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, int_0(dot(grad(F), grad(F))*v1))
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, int_0(exp(-F)*v1))
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, int_0(cos(F)*v1))
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, int_0(cos(F**2)*v1))
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, int_0(F**2*dot(grad(F), grad(v1))))
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(w1, int_0(dot(rot(G), grad(G))*w1))
    a = linearize(l, G, trials='u1')
    print(a)
    # ...

#==============================================================================
def test_linearize_form_2d_2():
    domain = Domain('Omega', dim=DIM)
    x,y = domain.coordinates

    V1 = ScalarFunctionSpace('V1', domain)

    v1 = element_of(V1, name='v1')

    alpha = Constant('alpha')

    F = element_of(V1, name='F')
    G = element_of(V1, name='G')

    int_0 = lambda expr: integral(domain , expr)

    # ...
    l1 = LinearForm(v1, int_0(F**2*v1))
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
    W = ScalarFunctionSpace('W', domain)

    v   = element_of(U, name='v')
    phi =       element_of(W, name='phi')
    q   =       element_of(W, name='q')

    U_0   = element_of(U, name='U_0')
    Rho_0 = element_of(W, name='Rho_0')
    P_0   = element_of(W, name='P_0')

    int_0 = lambda expr: integral(domain , expr)

    # ...
    expr = div(Rho_0*U_0) * phi
    l1 = LinearForm(phi, int_0(expr))

    expr = Rho_0*dot(convect(U_0, grad(U_0)), v) + dot(grad(P_0), v)
    l2 = LinearForm(v, int_0(expr))

    expr = dot(U_0, grad(P_0)) * q + P_0 * div(U_0) * q
    l3 = LinearForm(q, int_0(expr))
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
def test_linearize_form_2d_4():
    domain = Domain('Omega', dim=DIM)
    Gamma_N = Boundary(r'\Gamma_N', domain)

    x,y = domain.coordinates

    V = ScalarFunctionSpace('V', domain)

    v  = element_of(V, name='v')
    du = element_of(V, name='du')

    u = element_of(V, name='u')

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(Gamma_N, expr)

    g = Tuple(cos(pi*x)*sin(pi*y),
              sin(pi*x)*cos(pi*y))

    expr = dot(grad(v), grad(u)) - 4.*exp(-u)*v # + v*trace_1(g, Gamma_N)

    l = LinearForm(v, int_0(expr) )

    # linearising l around u, using du
    a = linearize(l, u, trials=du)
    # ...
    print(a)
    # ...

#==============================================================================
def test_area_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    mu    = Constant('mu'   , is_real=True)

    e = ElementDomain(domain)
    area = Area(e)

    V = ScalarFunctionSpace('V', domain)

    u,v = [element_of(V, name=i) for i in ['u', 'v']]

    int_0 = lambda expr: integral(domain , expr)

    # ...
    a = BilinearForm((v,u), int_0(area * u * v))
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

    V = ScalarFunctionSpace('V', domain)

    u,v = [element_of(V, name=i) for i in ['u', 'v']]

    int_0 = lambda expr: integral(domain , expr)

    # ...
    expr = kappa * dot(grad(u), grad(v)) + dot(b, grad(u)) * v
    a = BilinearForm((v,u), int_0(expr))
    # ...

    # ...
    expr = f * v
    l = LinearForm(v, int_0(expr))
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * dot(b, grad(v))
    s1 = BilinearForm((v,u), int_0(expr))

    expr = - f * dot(b, grad(v))
    l1 = LinearForm(v, int_0(expr))
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * ( dot(b, grad(v)) - kappa * laplace(v))
    s2 = BilinearForm((v,u), int_0(expr))

    expr = - f * ( dot(b, grad(v)) - kappa * laplace(v))
    l2 = LinearForm(v, int_0(expr))
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * ( dot(b, grad(v)) + kappa * laplace(v))
    s3 = BilinearForm((v,u), int_0(expr))

    expr = - f * ( dot(b, grad(v)) + kappa * laplace(v))
    l3 = LinearForm(v, int_0(expr))
    # ...

    # ...
    expr = a(v,u) + mu*area*s1(v,u)
    a1 = BilinearForm((v,u), int_0(expr))
    # ...

    # ...
    expr = a(v,u) + mu*area*s2(v,u)
    a2 = BilinearForm((v,u), int_0(expr))
    # ...

    # ...
    expr = a(v,u) + mu*area*s3(v,u)
    a3 = BilinearForm((v,u), int_0(expr))
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

    V = ScalarFunctionSpace('V', domain)

    u,v = [element_of(V, name=i) for i in ['u', 'v']]

    int_0 = lambda expr: integral(domain , expr)

    # ...
    expr = dot(grad(u), grad(v)) + f(x,y) * u * v
    a = BilinearForm((v,u), int_0(expr))

    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    expr = f(x,y) * v
    l = LinearForm(v, int_0(expr))

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

    V = ScalarFunctionSpace('V', domain)
    F = element_of(V, name='F')

    int_0 = lambda expr: integral(domain , expr)

    # ...
    expr = x*y
    a = Functional(int_0(expr), domain)

    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

    # ...
    expr = F - cos(2*pi*x)*cos(3*pi*y)
    expr = dot(grad(expr), grad(expr))
    a = Functional(int_0(expr), domain)

    print(a)
    print(TerminalExpr(a))
    print('')
    # ...

#==============================================================================
def test_norm_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    V = ScalarFunctionSpace('V', domain)
    F = element_of(V, name='F')

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
    F = element_of(V, 'F')

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
# this test checks some properties of bilinear forms
def test_bilinear_form_2d_3():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = ScalarFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    # ...
    a = BilinearForm((u,v), int_0(u*v))
    assert(a.is_symmetric)
    # ...

    # ...
    a = BilinearForm((u,v), int_0(dot(grad(u), grad(v))))
    assert(a.is_symmetric)
    # ...

#==============================================================================
# this test checks some properties of bilinear forms
def test_bilinear_form_2d_4():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    V = VectorFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)
    # ...
    a = BilinearForm((u,v), int_0(dot(u,v)))
    assert(a.is_symmetric)
    # ...

    # ...
    a = BilinearForm((u,v), int_0(inner(grad(u), grad(v))))
    assert(a.is_symmetric)
    # ...

def test_linearity_linear_form_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    nn    = NormalVector('nn')

    V = ScalarFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    l = LinearForm(v, int_0(x*y*v))
    assert is_linear_form(l)

    # ...
    l = LinearForm(v, int_0(x*y*v + v))
    assert is_linear_form(l)

    # ...

    # ...
    g = Tuple(x**2, y**2)
    l = LinearForm(v, int_1(v*dot(g, nn)))
    assert is_linear_form(l)

    # ...

    # ...
    g = Tuple(x**2, y**2)
    l = LinearForm(v, int_1(v*dot(g, nn)) + int_0(x*y*v))
    assert is_linear_form(l)

    # ...

    # ...
    l1 = LinearForm(v1, int_0(x*y*v1))
    l = LinearForm(v, l1(v))
    assert is_linear_form(l)


    # ...
    g = Tuple(x,y)
    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v2, int_0(dot(grad(v2), g)))

    l = LinearForm(v, l1(v) + l2(v))
    assert is_linear_form(l)

    # ...

    # ...
    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v1, int_0(v1))
    l = LinearForm(v, l1(v) + kappa*l2(v))
    assert is_linear_form(l)

    # ...

    # ...
    g = Tuple(x**2, y**2)
    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v1, int_0(v1))
    l3 = LinearForm(v, int_1(v*dot(g, nn)))
    l = LinearForm(v, l1(v) + kappa*l2(v) + mu*l3(v))
    assert is_linear_form(l)

#==============================================================================
def test_linearity_bilinear_form_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    eps   = Constant('eps', real=True)
    nn    = NormalVector('nn')

    V = ScalarFunctionSpace('V', domain)

    u,u1,u2 = [element_of(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [element_of(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    a = BilinearForm((u,v), int_0(u*v))
    assert is_bilinear_form(a)


    # ...
    a = BilinearForm((u,v), int_0(dot(grad(u),grad(v))))
    assert is_bilinear_form(a)

    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v + dot(grad(u),grad(v))))
    assert is_bilinear_form(a)

    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v + dot(grad(u),grad(v))) + int_1(v*dot(grad(u), nn)) )
    assert is_bilinear_form(a)


    # ...
    a = BilinearForm(((u1,u2),(v1,v2)), int_0(u1*v1 + u2*v2))
    assert is_bilinear_form(a)

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a = BilinearForm((u,v), a1(u,v))
    assert is_bilinear_form(a)


    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a = BilinearForm((u,v), a1(u,v) + a2(u,v))
    assert is_bilinear_form(a)


    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))
    assert is_bilinear_form(a)

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a3 = BilinearForm((u,v), int_1(v*dot(grad(u), nn)))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v) + mu*a3(u,v))
    assert is_bilinear_form(a)


    # ... Poisson with Nitsch method
    a0 = BilinearForm((u,v), int_0(dot(grad(u),grad(v))))
    a_B1 = BilinearForm((u,v), int_1(- kappa * u*dot(grad(v), nn)
                               - v*dot(grad(u), nn)
                               + u*v / eps))
    a = BilinearForm((u,v), a0(u,v) + a_B1(u,v))
    assert is_bilinear_form(a)

    # ...

#==============================================================================
def test_interface_2d_1():

    # ...
    def two_patches():

        from sympde.topology import InteriorDomain
        from sympde.topology import Connectivity

        A = Square('A')
        B = Square('B')

        A = A.interior
        B = B.interior

        connectivity = Connectivity()

        bnd_A_1 = Boundary(r'\Gamma_1', A, axis=0, ext=-1)
        bnd_A_2 = Boundary(r'\Gamma_2', A, axis=0, ext=1)
        bnd_A_3 = Boundary(r'\Gamma_3', A, axis=1, ext=-1)
        bnd_A_4 = Boundary(r'\Gamma_4', A, axis=1, ext=1)

        bnd_B_1 = Boundary(r'\Gamma_1', B, axis=0, ext=-1)
        bnd_B_2 = Boundary(r'\Gamma_2', B, axis=0, ext=1)
        bnd_B_3 = Boundary(r'\Gamma_3', B, axis=1, ext=-1)
        bnd_B_4 = Boundary(r'\Gamma_4', B, axis=1, ext=1)

        connectivity['I'] = (bnd_A_2, bnd_B_1)

        Omega = Domain('Omega',
                       interiors=[A, B],
                       boundaries=[bnd_A_1, bnd_A_2, bnd_A_3, bnd_A_4, bnd_B_1, bnd_B_2, bnd_B_3, bnd_B_4],
                       connectivity=connectivity)

        return Omega
    # ...

    # create a domain with an interface
    domain = two_patches()
    interfaces = domain.interfaces
#    interfaces = Union(*interfaces) # TODO move to domain.connectivity.interfaces

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u, v')

    print(integral(interfaces, u*v))

    expr  = integral(domain, dot(grad(v),grad(u)))
    expr += integral(interfaces, - avg(Dn(u)) * jump(v)
                                 + avg(Dn(v)) * jump(u))
    a  = BilinearForm((u,v), expr)
    print(a)

#==============================================================================
def test_interface_integral_1():

    # ...
    A = Square('A')
    B = Square('B')

    domain = A.join(B, name = 'domain',
                bnd_minus = A.get_boundary(axis=0, ext=1),
                bnd_plus  = B.get_boundary(axis=0, ext=-1))
    # ...

    x,y = domain.coordinates

    V = ScalarFunctionSpace('V', domain, kind=None)
    assert(V.is_broken)

    u, v = elements_of(V, names='u, v')

    # ...
    I = domain.interfaces
    # ...

#    expr = minus(Dn(u))
#    print(expr)
#    import sys; sys.exit(0)

    # ... bilinear forms
#    a = BilinearForm((u,v), integral(domain, u*v))
#    a = BilinearForm((u,v), integral(domain, dot(grad(u),grad(v))))
#    a = BilinearForm((u,v), integral(I, jump(u) * jump(v)))
#    a = BilinearForm((u,v), integral(I, jump(Dn(u)) * jump(v)))

#    a = BilinearForm((u,v), integral(domain, dot(grad(u),grad(v)))
#                          + integral(I,      jump(u) * jump(v)))

    # Nitsch
    kappa = Constant('kappa')
    expr_I = ( - jump(u) * jump(Dn(v))
               + kappa * jump(u) * jump(v)
               + plus(Dn(u)) * minus(v)
               + minus(Dn(u)) * plus(v) )
    a = BilinearForm((u,v), integral(domain, dot(grad(u),grad(v)))
                          + integral(I,      expr_I))

#    # TODO BUG
#    bnd_A = A.get_boundary(axis=0, ext=1)
#
#    a = BilinearForm((u,v), integral(domain, dot(grad(u),grad(v)))
#                          + integral(I,      jump(u) * jump(v))
#                          + integral(bnd_A,      dx(u)*v))

    expr = TerminalExpr(a)
    print(expr)
    # ...

#    # ... linear forms
#    b = LinearForm(v, integral(domain, sin(x+y)*v)
#                    + integral(I, cos(x+y) * jump(v)))
#
#    expr = TerminalExpr(b)
#    print(expr)
#    # ...

#==============================================================================
def test_interface_integral_2():

    # ...
    A = Square('A')
    B = Square('B')

    domain = A.join(B, name = 'domain',
                bnd_minus = A.get_boundary(axis=0, ext=1),
                bnd_plus  = B.get_boundary(axis=0, ext=-1))
    # ...

    x,y = domain.coordinates

    V = ScalarFunctionSpace('V', domain, kind=None)
    assert(V.is_broken)

    u, u1, u2, u3 = elements_of(V, names='u, u1, u2, u3')
    v, v1, v2, v3 = elements_of(V, names='v, v1, v2, v3')

    # ...
    I = domain.interfaces

    a = BilinearForm((u,v), integral(domain, dot(grad(u),grad(v))))
    b = BilinearForm((u,v), integral(I, jump(u) * jump(v)))

    A = BilinearForm(((u1,u2),(v1,v2)), a(u1,v1) + a(u2,v2) + b(u1,v1) + b(u2,v2) + b(u1, v2) )
    B = BilinearForm(((u1,u2,u3),(v1,v2,v3)), a(u1,v1) + a(u2,v2) + a(u3,v3) + b(u1,v1) + b(u2,v2) + b(u1, v2) )

    print(TerminalExpr(A))
    print(TerminalExpr(B))
    # ...

    # ... linear forms
    b = LinearForm(v, integral(I, jump(v)))

    b = LinearForm((v1,v2), b(v1) + b(v2) )
    expr = TerminalExpr(b)
    print(expr)
    # ...

#==============================================================================
def test_interface_integral_3():

    # ...
    A = Square('A')
    B = Square('B')
    C = Square('C')

    AB = A.join(B, name = 'AB',
               bnd_minus = A.get_boundary(axis=0, ext=1),
               bnd_plus  = B.get_boundary(axis=0, ext=-1))

    domain = AB.join(C, name = 'domain',
               bnd_minus = B.get_boundary(axis=0, ext=1),
               bnd_plus  = C.get_boundary(axis=0, ext=-1))
    # ...

    x,y = domain.coordinates

    V = ScalarFunctionSpace('V', domain, kind=None)
    assert(V.is_broken)

    u, v = elements_of(V, names='u, v')

    # ...
    I = domain.interfaces
#    print(I)
#    print(integral(I, jump(u) * jump(v)))

#    a = BilinearForm((u,v), integral(domain, u*v))
#    a = BilinearForm((u,v), integral(domain, dot(grad(u),grad(v))))
#    a = BilinearForm((u,v), integral(I, jump(u) * jump(v)))

    a = BilinearForm((u,v), integral(domain, dot(grad(u),grad(v)))
                          + integral(I,      jump(u) * jump(v)))

    expr = TerminalExpr(a)
    print(expr)
    # ...

    # ... linear forms
    b = LinearForm(v, integral(domain, sin(x+y)*v)
                    + integral(I, cos(x+y) * jump(v)))

    expr = TerminalExpr(b)
    print(expr)
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
