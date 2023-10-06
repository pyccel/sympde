# coding: utf-8

# TODO: - add assert to every test

import pytest

from sympy.core.containers import Tuple
from sympy import Function
from sympy import pi, cos, sin, exp
from sympy import ImmutableDenseMatrix as Matrix

from sympde.core     import Constant
from sympde.calculus import grad, dot, inner, rot, div
from sympde.calculus import laplace, bracket, convect
from sympde.calculus import jump, avg, Dn, minus, plus

from sympde.topology import dx1, dx2, dx3
from sympde.topology import dx, dy, dz
from sympde.topology import Mapping
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import element_of, elements_of
from sympde.topology import InteriorDomain, Union
from sympde.topology import Boundary, NormalVector
from sympde.topology import Domain
#from sympde.topology import trace_1  # TODO [YG, 27.01.2021]: fix trace
from sympde.topology import Square
from sympde.topology import ElementDomain
from sympde.topology import Area

from sympde.expr.expr import LinearExpr
from sympde.expr.expr import LinearForm, BilinearForm
from sympde.expr.expr import integral
from sympde.expr.expr import Functional, Norm
from sympde.expr.expr import linearize
from sympde.expr.evaluation import TerminalExpr

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

    assert(l.domain == domain.interior)
    assert(l(v1) == int_0(x*y*v1))
    # ...

    # ...
    g  = Tuple(x**2, y**2)
    l  = LinearForm(v, int_1(v*dot(g, nn)))
    print(l)

    assert(l.domain == B1)
#    assert(l(v1) == int_1(v1*trace_1(g, B1))) # TODO [YG, 27.01.2021]: fix trace
    assert(l(v1) == int_1(v1*dot(g, nn)))
    # ...

    # ...
    g = Tuple(x**2, y**2)
    l = LinearForm(v, int_1(v*dot(g, nn)) + int_0(x*y*v))

    assert(len(l.domain.args) == 2)
    for i in l.domain.args:
        assert(i in [domain.interior, B1])

#    assert(l(v1) == int_1(v1*trace_1(g, B1)) + int_0(x*y*v1)) # TODO [YG, 27.01.2021]: fix trace
    assert(l(v1) == int_1(v1*dot(g, nn)) + int_0(x*y*v1))
    # ...

    # ...
    l1 = LinearForm(v1, int_0(x*y*v1))
    l = LinearForm(v, l1(v))

    assert(l.domain == domain.interior)
    assert(l(u1) == int_0(x*y*u1))
    # ...

    # ...
    g = Tuple(x,y)
    l1 = LinearForm(u1, int_0(x*y*u1))
    l2 = LinearForm(u2, int_0(dot(grad(u2), g)))

    l = LinearForm(v, l1(v) + l2(v))

    assert(l.domain == domain.interior)
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

    assert(l.domain == domain.interior)
    assert(l(u1,u2).expand() == int_0(-1.0*dt*dot(grad(u1), grad(wn))/Re) +
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

    g = Matrix((x,y))
    l = LinearForm(v, int_0(dot(g, v)))

    assert(l.domain == domain.interior)
    assert(l(v1) == int_0(dot(g, v1)))
    # ...

    # ...
    g = Matrix((x,y))
    l1 = LinearForm(v1, int_0(dot(g, v1)))
    l = LinearForm(v, l1(v))

    assert(l.domain == domain.interior)
    assert(l(u1) == int_0(dot(g, u1)))
    # ...

    # ...
    g1 = Matrix((x,0))
    g2 = Matrix((0,y))
    l1 = LinearForm(v1, int_0(dot(v1, g1)))
    l2 = LinearForm(v2, int_0(dot(v2, g2)))

    l = LinearForm(v, l1(v) + l2(v))

    assert(l.domain == domain.interior)
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

    assert(a.domain == domain.interior)
    assert(a(u1,v1) == int_0(u1*v1))
    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v + dot(grad(u), grad(v))))

    assert(a.domain == domain.interior)
    assert(a(u1,v1) == int_0(u1*v1) + int_0(dot(grad(u1), grad(v1))))
    # ...

    # ...
    a = BilinearForm((u,v), int_1(v*dot(grad(u), nn)))

    assert(a.domain == B1)
#    assert(a(u1,v1) == int_1(v1*trace_1(grad(u1), B1))) # TODO [YG, 27.01.2021]: fix trace
    assert a(u1,v1) == int_1(v1*dot(grad(u1), nn))
    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v) + int_1(v*dot(grad(u), nn)))

    # TODO a.domain are not ordered
    assert(len(a.domain.args) == 2)
    for i in a.domain.args:
        assert(i in [domain.interior, B1])

#    assert(a(u1,v1) == int_0(u1*v1) + int_1(v1*trace_1(grad(u1), B1))) # TODO [YG, 27.01.2021]: fix trace
    assert a(u1,v1) == int_0(u1*v1) + int_1(v1*dot(grad(u1), nn))
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a = BilinearForm((u,v), a1(u,v))

    assert(a.domain == domain.interior)
    assert(a(u2,v2) == int_0(u2*v2))
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))

    assert(a.domain == domain.interior)
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

    assert(a.domain == domain.interior)
    assert(a(u1,v1) == int_0(dot(u1,v1)))
    # ...

    # ...
    a = BilinearForm((u,v), int_0(dot(u,v) + inner(grad(u), grad(v))))

    assert(a.domain == domain.interior)
    assert(a(u1,v1) == int_0(dot(u1,v1)) + int_0(inner(grad(u1), grad(v1))))
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(dot(u1,v1)))
    a = BilinearForm((u,v), a1(u,v))

    assert(a.domain == domain.interior)
    assert(a(u2,v2) == int_0(dot(u2,v2)))
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(dot(u1,v1)))
    a2 = BilinearForm((u2,v2), int_0(inner(grad(u2), grad(v2))))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))

    assert(a.domain == domain.interior)
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

    print(TerminalExpr(l, domain))
    print('')
    # ...

    # ...
    l = LinearForm(v, int_0(x*y*v + v))
    print(TerminalExpr(l, domain))
    print('')
    # ...

    # ...
    g = Matrix((x**2, y**2))
    l = LinearForm(v, int_1(v*dot(g, nn)))
    print(TerminalExpr(l, domain))
    print('')
    # ...

    # ...
    g = Matrix((x**2, y**2))
    l = LinearForm(v, int_1(v*dot(g, nn)) + int_0(x*y*v))
    print(TerminalExpr(l, domain))
    print('')
    # ...

    # ...
    l1 = LinearForm(v1, int_0(x*y*v1))
    l = LinearForm(v, l1(v))
    print(TerminalExpr(l, domain))
    print('')
    # ...

    # ...
    g = Matrix((x,y))
    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v2, int_0(dot(grad(v2), g)))

    l = LinearForm(v, l1(v) + l2(v))
    print(TerminalExpr(l, domain))
    print('')
    # ...

    # ...
    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v1, int_0(v1))
    l = LinearForm(v, l1(v) + kappa*l2(v))
    print(TerminalExpr(l, domain))
    print('')
    # ...

    # ...
    g = Matrix((x**2, y**2))
    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v1, int_0(v1))
    l3 = LinearForm(v, int_1(v*dot(g, nn)))
    l = LinearForm(v, l1(v) + kappa*l2(v) + mu*l3(v))
    print(TerminalExpr(l, domain))
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

    g = Matrix((x,y))
    l = LinearForm(v, int_0(dot(g, v)))
    print(TerminalExpr(l, domain))
    print('')
    # ...

    # ...
    g = Matrix((x,y))
    l = LinearForm(v, int_0(dot(g, v) + div(v)))
    print(TerminalExpr(l, domain))
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
    print(TerminalExpr(l, domain))
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
    print(TerminalExpr(l, domain))
    print('')
    # ...
#==============================================================================
def test_terminal_expr_linear_2d_5():

    # ... abstract model
    domain = Square()

    V      = ScalarFunctionSpace('V', domain)

    B_neumann = [domain.get_boundary(axis=0, ext=-1), 
                domain.get_boundary(axis=1, ext=-1)]

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

    print(TerminalExpr(l, domain))
    print('')
    # ...

#==============================================================================
def test_terminal_expr_bilinear_2d_1():

    domain = Domain('Omega', dim=2)
    B1     = Boundary(r'\Gamma_1', domain)

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
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), int_0(dot(grad(u),grad(v))))
    print(a)
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v + dot(grad(u),grad(v))))
    print(a)
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v + dot(grad(u),grad(v))) + int_1(v*dot(grad(u), nn)) )
    print(a)
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a = BilinearForm(((u1,u2),(v1,v2)), int_0(u1*v1 + u2*v2))
    print(a)
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a = BilinearForm((u,v), a1(u,v))
    print(a)
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a = BilinearForm((u,v), a1(u,v) + a2(u,v))
    print(a)
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))
    print(a)
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a3 = BilinearForm((u,v), int_1(v*dot(grad(u), nn)))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v) + mu*a3(u,v))
    print(a)
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ... Poisson with Nitsch method
    a0 = BilinearForm((u,v), int_0(dot(grad(u),grad(v))))
    a_B1 = BilinearForm((u,v), int_1(- kappa * u*dot(grad(v), nn)
                               - v*dot(grad(u), nn)
                               + u*v / eps))
    a = BilinearForm((u,v), a0(u,v) + a_B1(u,v))
    print(a)
    print(TerminalExpr(a, domain))
    print('')
    # ...

#==============================================================================
def test_terminal_expr_bilinear_2d_2():

    domain = Domain('Omega', dim=2)
    B1     = Boundary(r'\Gamma_1', domain)

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
    print(TerminalExpr(a, domain))
    print('')

    # ...
    a = BilinearForm((u,v), int_0(inner(grad(u),grad(v))))
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), int_0(dot(u,v) + inner(grad(u),grad(v))))
    print(TerminalExpr(a, domain))
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

    print(TerminalExpr(a, domain))
    print('')

    a = BilinearForm((u,v), int_0(u*v + dot(grad(u),grad(v))) + int_1(v*dot(grad(u), nn)) )
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a = BilinearForm((u,v), int_0(u*v))
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u,v), int_0(u*v))
    a = BilinearForm((u,v), a1(u,v))
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u,v), int_0(u*v))
    a2 = BilinearForm((u,v), int_0(dot(grad(u), grad(v))))
    a = BilinearForm((u,v), a1(u,v) + a2(u,v))
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u,v), int_0(u*v))
    a2 = BilinearForm((u,v), int_0(dot(grad(u), grad(v))))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    a1 = BilinearForm((u,v), int_0(u*v))
    a2 = BilinearForm((u,v), int_0(dot(grad(u), grad(v))))
    a3 = BilinearForm((u,v), int_1(v*dot(grad(u), nn)))
    a = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v) + mu*a3(u,v))
    print(TerminalExpr(a, domain))
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

    print(TerminalExpr(a, domain))
    print('')

#==============================================================================
def test_terminal_expr_bilinear_3d_1():

    domain = Domain('Omega', dim=3)
    M      = Mapping('M', dim=3)

    mapped_domain = M(domain)

    V  = ScalarFunctionSpace('V', domain)
    VM = ScalarFunctionSpace('VM', mapped_domain)

    u,v   = elements_of(V, names='u,v')
    um,vm = elements_of(VM, names='u,v')

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(mapped_domain , expr)

    J   = M.det_jacobian
    det = dx1(M[0])*dx2(M[1])*dx3(M[2]) - dx1(M[0])*dx2(M[2])*dx3(M[1]) - dx1(M[1])*dx2(M[0])*dx3(M[2])\
        + dx1(M[1])*dx2(M[2])*dx3(M[0]) + dx1(M[2])*dx2(M[0])*dx3(M[1]) - dx1(M[2])*dx2(M[1])*dx3(M[0])


    a1 = BilinearForm((u,v), int_0(dot(grad(u),grad(v))))
    a2 = BilinearForm((um,vm), int_1(dot(grad(um),grad(vm))))
    a3 = BilinearForm((u,v), int_0(J*dot(grad(u),grad(v))))

    e1 = TerminalExpr(a1, domain)
    e2 = TerminalExpr(a2, domain)
    e3 = TerminalExpr(a3, domain)

    assert e1[0].expr          == dx1(u)*dx1(v) + dx2(u)*dx2(v) + dx3(u)*dx3(v)
    assert e2[0].expr          == dx(um)*dx(vm) + dy(um)*dy(vm) + dz(um)*dz(vm)
    assert e3[0].expr.factor() == (dx1(u)*dx1(v) + dx2(u)*dx2(v) + dx3(u)*dx3(v))*det

#==============================================================================
@pytest.mark.skip(reason="New linearize() function does not accept 'LinearExpr' objects")
def test_linearize_expr_2d_1():
    domain = Domain('Omega', dim=2)
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
@pytest.mark.skip(reason="New linearize() function does not accept 'LinearExpr' objects")
def test_linearize_expr_2d_2():
    domain = Domain('Omega', dim=2)
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
    domain = Domain('Omega', dim=2)

    V = ScalarFunctionSpace('V', domain)
    W = VectorFunctionSpace('W', domain)

    v, F, u = elements_of(V, names='v, F, u')
    w, G, m = elements_of(W, names='w, G, m')

    int_0 = lambda expr: integral(domain , expr)

    # ...
    l  = LinearForm(v, int_0(F**2*v))
    a  = linearize(l, F, trials=u)
    assert a(u, v) == int_0(2 * F * u * v)
    # ...

    # ...
    l = LinearForm(v, int_0(dot(grad(F), grad(F))*v))
    a = linearize(l, F, trials=u)
    assert a(u, v) == int_0(2 * dot(grad(F), grad(u)) * v)
    # ...

    # ...
    l = LinearForm(v, int_0(exp(-F)*v))
    a = linearize(l, F, trials=u)
    assert a(u, v) == int_0(-exp(-F) * u * v)
    # ...

    # ...
    l = LinearForm(v, int_0(cos(F)*v))
    a = linearize(l, F, trials=u)
    assert a(u, v) == int_0(-sin(F) * u * v)
    # ...

    # ...
    l = LinearForm(v, int_0(cos(F**2)*v))
    a = linearize(l, F, trials=u)
    assert a(u, v) == int_0(-2 * F *sin(F**2) * u * v)
    # ...

    # ...
    l = LinearForm(v, int_0(F**2*dot(grad(F), grad(v))))
    a = linearize(l, F, trials=u)
    assert a(u, v) == int_0(2 * F * u * dot(grad(F), grad(v)) + F**2 * dot(grad(u), grad(v)))
    # ...

    # ...
    l = LinearForm(w, int_0(dot(rot(G), grad(G))*w))
    a = linearize(l, G, trials=m)
    assert a(m, w) == int_0((dot(rot(m), grad(G)) + dot(rot(G), grad(m))) * w)
    # ...

#==============================================================================
def test_linearize_form_2d_2():
    domain = Domain('Omega', dim=2)

    V = ScalarFunctionSpace('V', domain)

    v, F, u = elements_of(V, names='v, F, u')

    int_0 = lambda expr: integral(domain , expr)

    # ...
    l1 = LinearForm(v, int_0(F**2*v))
    l  = LinearForm(v, l1(v))

    a = linearize(l, F, trials=u)

    expected = linearize(l1, F, trials=u)
    assert a == expected
    # ...

#==============================================================================
def test_linearize_form_2d_3():
    """steady Euler equation."""
    domain = Domain('Omega', dim=2)

    U = VectorFunctionSpace('U', domain)
    W = ScalarFunctionSpace('W', domain)

    # Test functions
    v   = element_of(U, name='v')
    phi = element_of(W, name='phi')
    q   = element_of(W, name='q')

    # Steady-state fields
    U_0   = element_of(U, name='U_0')
    Rho_0 = element_of(W, name='Rho_0')
    P_0   = element_of(W, name='P_0')

    # Trial functions (displacements from steady-state)
    d_u   = element_of(U, name='d_u')
    d_rho = element_of(W, name='d_rho')
    d_p   = element_of(W, name='d_p')

    # Shortcut
    int_0 = lambda expr: integral(domain , expr)

    # The Euler equations are a system of three non-linear equations; for each of
    # them we create a linear form in the test functions (phi, v, q) respectively.
    e1 = div(Rho_0 * U_0)
    l1 = LinearForm(phi, int_0(e1 * phi))

    e2 = Rho_0 * convect(U_0, U_0) + grad(P_0)
    l2 = LinearForm(v, int_0(dot(e2, v)))

    e3 = div(P_0 * U_0)
    l3 = LinearForm(q, int_0(e3 * q))
    # ...

    # Linearize l1, l2 and l3 separately
    a1 = linearize(l1, fields=[Rho_0, U_0     ], trials=[d_rho, d_u     ])
    a2 = linearize(l2, fields=[Rho_0, U_0, P_0], trials=[d_rho, d_u, d_p])
    a3 = linearize(l3, fields=[       U_0, P_0], trials=[       d_u, d_p])

    # Check individual bilinear forms
    d_e1 = div(U_0 * d_rho + Rho_0 * d_u)
    d_e2 = d_rho * convect(U_0, U_0) + \
           Rho_0 * convect(d_u, U_0) + \
           Rho_0 * convect(U_0, d_u) + grad(d_p)
    d_e3 = div(d_p * U_0 + P_0 * d_u)

    assert a1([d_rho, d_u     ], phi) == int_0(d_e1 * phi)
    assert a2([d_rho, d_u, d_p], v  ) == int_0(dot(d_e2, v))
    assert a3([       d_u, d_p], q  ) == int_0(d_e3 * q)

    # Linearize linear form of system: l = l1 + l2 + l3
    l = LinearForm((phi, v, q), l1(phi) + l2(v) + l3(q))
    a = linearize(l, fields=[Rho_0, U_0, P_0], trials=[d_rho, d_u, d_p])

    # Check composite linear form
    assert a([d_rho, d_u, d_p], [phi, v, q]) == \
            int_0(d_e1 * phi + dot(d_e2, v) + d_e3 * q)

#==============================================================================
def test_linearize_form_2d_4():
    domain  = Domain('Omega', dim=2)
    Gamma_N = Boundary(r'\Gamma_N', domain)

    x,y = domain.coordinates

    V = ScalarFunctionSpace('V', domain)

    v  = element_of(V, name='v')
    u  = element_of(V, name='u')
    du = element_of(V, name='du')

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(Gamma_N, expr)

#    g = Matrix((cos(pi*x)*sin(pi*y),
#              sin(pi*x)*cos(pi*y)))

    expr = dot(grad(v), grad(u)) - 4.*exp(-u)*v # + v*trace_1(g, Gamma_N)

    l = LinearForm(v, int_0(expr) )

    # linearising l around u, using du
    a = linearize(l, u, trials=du)

    assert a(du, v) == int_0(dot(grad(v), grad(du)) + 4.*exp(-u) * du * v)

#==============================================================================
def test_area_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    mu    = Constant('mu'   , is_real=True)

    e = ElementDomain()
    area = Area(e)

    V = ScalarFunctionSpace('V', domain)

    u,v = [element_of(V, name=i) for i in ['u', 'v']]

    int_0 = lambda expr: integral(domain , expr)

    # ...
    a = BilinearForm((v,u), int_0(area * u * v))
    print(TerminalExpr(a, domain))
    # ...

#==============================================================================
def test_stabilization_2d_1():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    b1 = 1.
    b2 = 0.
    b = Matrix((b1, b2))

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
    print(TerminalExpr(a1, domain))
    print('')

    print(a2)
    print(TerminalExpr(a2, domain))
    print('')

    print(a3)
    print(TerminalExpr(a3, domain))
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
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    expr = f(x,y) * v
    l = LinearForm(v, int_0(expr))

    print(l)
    print(TerminalExpr(l, domain))
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
    print(TerminalExpr(a, domain))
    print('')
    # ...

    # ...
    expr = F - cos(2*pi*x)*cos(3*pi*y)
    expr = dot(grad(expr), grad(expr))
    a = Functional(int_0(expr), domain)

    print(a)
    print(TerminalExpr(a, domain))
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

    print('> l2 norm = ', TerminalExpr(l2_norm_u, domain))
    print('> h1 norm = ', TerminalExpr(h1_norm_u, domain))
    print('')
    # ...

    # ...
    expr = sin(pi*x)*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u, domain))
    print('> h1 norm = ', TerminalExpr(h1_norm_u, domain))
    print('')
    # ...

    # ...
    expr = F-x*y
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u, domain))
    print('> h1 norm = ', TerminalExpr(h1_norm_u, domain))
    print('')
    # ...

    # ...
    expr = F-sin(pi*x)*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u, domain))
    print('> h1 norm = ', TerminalExpr(h1_norm_u, domain))
    print('')
    # ...

    # ...
    expr = F-sin(0.5*pi*(1.-x))*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u, domain))
    print('> h1 norm = ', TerminalExpr(h1_norm_u, domain))
    print('')
    # ...

    # ...
    expr = F-cos(0.5*pi*x)*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u, domain))
    print('> h1 norm = ', TerminalExpr(h1_norm_u, domain))
    print('')
    # ...

#==============================================================================
def test_norm_2d_2():

    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    V = VectorFunctionSpace('V', domain)
    F = element_of(V, 'F')

    # ...
    f = Matrix((sin(pi*x)*sin(pi*y), sin(pi*x)*sin(pi*y)))
    expr = Matrix([F[0]-f[0], F[1]-f[1]])
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', TerminalExpr(l2_norm_u, domain))
    print('> h1 norm = ', TerminalExpr(h1_norm_u, domain))
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

#==============================================================================
def test_linearity_linear_form_2d_1():

    from sympde.expr.errors import UnconsistentLinearExpressionError

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    nn    = NormalVector('nn')

    V = ScalarFunctionSpace('V', domain)
    W = VectorFunctionSpace('W', domain)

    v, v1, v2 = elements_of(V, names='v, v1, v2')
    w = element_of(W, name='w')

    # ...
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    # The following integral expressions are linear, hence it must be possible
    # to create LinearForm objects from them

    _ = LinearForm(v, int_0(x * y * v))
    _ = LinearForm(v, int_0(x * y * v + v))

    g = Matrix((x**2, y**2))
    _ = LinearForm(v, int_0(v * dot(g, nn)))

    g = Matrix((x**2, y**2))
    _ = LinearForm(v, int_1(v*dot(g, nn)) + int_0(x*y*v))

    l1 = LinearForm(v1, int_0(x * y * v1))
    _  = LinearForm(v, l1(v))

    g  = Matrix((x,y))
    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v2, int_0(dot(grad(v2), g)))
    _  = LinearForm(v, l1(v) + l2(v))

    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v1, int_0(v1))
    _  = LinearForm(v, l1(v) + kappa*l2(v))

    g  = Matrix((x**2, y**2))
    l1 = LinearForm(v1, int_0(x*y*v1))
    l2 = LinearForm(v1, int_0(v1))
    l3 = LinearForm(v, int_1(v*dot(g, nn)))
    _  = LinearForm(v, l1(v) + kappa*l2(v) + mu*l3(v))

    l1 = LinearForm(w, int_0(5 * w * x))
    l2 = LinearForm(w, int_1(w * y))
    _  = LinearForm(w, l1(w) + kappa*l2(w))

    l0 = LinearForm(v, int_0(x * y * v))
    l1 = LinearForm(w, int_0(w[0] * y))
    l2 = LinearForm(w, int_1(w[1] * x))
    l3 = LinearForm(w, kappa * l1(w) + mu * l2(w))
    _  = LinearForm((v, w), l0(v) + l3(w))

    # The following integral expressions are not linear, hence LinearForm must
    # raise an exception

    with pytest.raises(UnconsistentLinearExpressionError):
        _ = LinearForm(v, int_0(x * y * v + 1))

    with pytest.raises(UnconsistentLinearExpressionError):
        _ = LinearForm(v, int_0(x * v**2))

    with pytest.raises(UnconsistentLinearExpressionError):
        _ = LinearForm(v, int_0(x * y * v) + int_1(v * exp(v)))

    with pytest.raises(UnconsistentLinearExpressionError):
        _ = LinearForm(w, int_0(w * w))

    with pytest.raises(UnconsistentLinearExpressionError):
        _ = LinearForm(w, int_0(w[0] * w[1]))

    with pytest.raises(UnconsistentLinearExpressionError):
        _ = LinearForm((v, w), int_0(x * w[0]) + int_1(v * w[1]))

#==============================================================================
def test_linearity_bilinear_form_2d_1():

    from sympde.expr.errors import UnconsistentLinearExpressionError

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    eps   = Constant('eps', real=True)
    nn    = NormalVector('nn')

    V = ScalarFunctionSpace('V', domain)

    u, u1, u2 = elements_of(V, names='u, u1, u2')
    v, v1, v2 = elements_of(V, names='v, v1, v2')

    # ...

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    # The following integral expressions are bilinear, hence it must be possible
    # to create BilinearForm objects from them

    _ = BilinearForm((u,v), int_0(u*v))
    _ = BilinearForm((u,v), int_0(dot(grad(u),grad(v))))
    _ = BilinearForm((u,v), int_0(u*v + dot(grad(u),grad(v))))
    _ = BilinearForm((u,v), int_0(u*v + dot(grad(u),grad(v))) + int_1(v*dot(grad(u), nn)) )
    _ = BilinearForm(((u1,u2),(v1,v2)), int_0(u1*v1 + u2*v2))

    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    _  = BilinearForm((u,v), a1(u,v))

    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    _  = BilinearForm((u,v), a1(u,v) + a2(u,v))


    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    _  = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v))

    a1 = BilinearForm((u1,v1), int_0(u1*v1))
    a2 = BilinearForm((u2,v2), int_0(dot(grad(u2), grad(v2))))
    a3 = BilinearForm((u,v), int_1(v*dot(grad(u), nn)))
    _  = BilinearForm((u,v), a1(u,v) + kappa*a2(u,v) + mu*a3(u,v))

    # ... Poisson with Nitsch method
    a0  = BilinearForm((u,v), int_0(dot(grad(u),grad(v))))
    a_B1 = BilinearForm((u,v), int_1(- kappa * u*dot(grad(v), nn)
                               - v*dot(grad(u), nn)
                               + u*v / eps))
    _ = BilinearForm((u,v), a0(u,v) + a_B1(u,v))
    # ...

    # The following integral expressions are not bilinear, hence BilinearForm must
    # raise an exception

    with pytest.raises(UnconsistentLinearExpressionError):
        _ = BilinearForm((u, v), int_0(x * y * dot(grad(u), grad(v)) + 1))

    with pytest.raises(UnconsistentLinearExpressionError):
        _ = BilinearForm((u, v), int_0(x * dot(grad(u), grad(v**2))))

    with pytest.raises(UnconsistentLinearExpressionError):
        _ = BilinearForm((u, v), int_0(u * v) + int_1(v * exp(u)))

#==============================================================================
def test_interface_2d_1():

    # ...
    def two_patches():

        from sympde.topology import InteriorDomain
        from sympde.topology import Connectivity, Interface

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

        connectivity['I'] = Interface('I', bnd_A_2, bnd_B_1)

        Omega = Domain('Omega',
                       interiors=[A, B],
                       boundaries=[bnd_A_1, bnd_A_2, bnd_A_3, bnd_A_4, bnd_B_1, bnd_B_2, bnd_B_3, bnd_B_4],
                       connectivity=connectivity)

        return Omega
    # ...

    # create a domain with an interface
    domain = two_patches()
    interfaces = domain.interfaces

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

    domains = [A, B]
    connectivity = [((0, 0, 1),(1, 0, -1))]
    domain = Domain.join(domains, connectivity, 'domain')
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

    expr = TerminalExpr(a, domain)
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


    domains = [A, B]
    connectivity = [((0, 0, 1),(1, 0, -1))]
    domain = Domain.join(domains, connectivity, 'domain')
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

    print(TerminalExpr(A, domain))
    print(TerminalExpr(B, domain))
    # ...

    # ... linear forms
    b = LinearForm(v, integral(I, jump(v)))

    b = LinearForm((v1,v2), b(v1) + b(v2) )
    expr = TerminalExpr(b, domain)
    print(expr)
    # ...

#==============================================================================
def test_interface_integral_3():

    # ...
    A = Square('A')
    B = Square('B')
    C = Square('C')


    domains = [A, B, C]
    connectivity = [((0, 0, 1),(1, 0, -1)),((1,0,1),(2,0,-1))]
    domain = Domain.join(domains, connectivity, 'domain')

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

    expr = TerminalExpr(a, domain)
    print(expr)
    # ...

    # ... linear forms
    b = LinearForm(v, integral(domain, sin(x+y)*v)
                    + integral(I, cos(x+y) * jump(v)))

    expr = TerminalExpr(b, domain)
    print(expr)
    # ...

#==============================================================================
def test_interface_integral_4():

    # ...
    A = Square('A')
    B = Square('B')

    domains = [A, B]
    connectivity = [((0, 0, 1),(1, 0, -1))]
    domain = Domain.join(domains, connectivity, 'AB')

    x,y = domain.coordinates
    assert all([x.name == 'x1', y.name == 'x2'])

    V1 = ScalarFunctionSpace('V1', domain, kind=None)
    V2 = VectorFunctionSpace('V2', domain, kind=None)
    assert(V1.is_broken)

    u1, v1 = elements_of(V1, names='u1, v1')
    u2, v2 = elements_of(V2, names='u2, v2')

    # ...
    I = domain.interfaces

    a1 = BilinearForm((u1,v1), integral(I, plus(u1)*plus(v1)))

    expr = TerminalExpr(a1, domain)

    assert len(expr) == 1
    assert expr[0].target == B.get_boundary(axis=0, ext=-1)
    assert expr[0].expr   == u1*v1

    a2 = BilinearForm((u1,v1), integral(I, minus(u1)*minus(v1)))

    expr = TerminalExpr(a2, domain)

    assert len(expr) == 1
    assert expr[0].target == A.get_boundary(axis=0, ext=1)
    assert expr[0].expr   == u1*v1


    a1 = BilinearForm((u2,v2), integral(I, dot(plus(u2),plus(v2))))

    expr = TerminalExpr(a1, domain)

    assert len(expr) == 1
    assert expr[0].target == B.get_boundary(axis=0, ext=-1)
    assert expr[0].expr[0,0]   == u2[0]*v2[0]
    assert expr[0].expr[1,1]   == u2[1]*v2[1]

    a2 = BilinearForm((u2,v2), integral(I, dot(minus(u2),minus(v2))))

    expr = TerminalExpr(a2, domain)

    assert len(expr) == 1
    assert expr[0].target == A.get_boundary(axis=0, ext=1)
    assert expr[0].expr[0,0]   == u2[0]*v2[0]
    assert expr[0].expr[1,1]   == u2[1]*v2[1]

    # ...
#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
