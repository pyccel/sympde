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
from sympy import pi, cos, sin, exp
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
def test_boundary_2d_1():
    domain = Domain('Omega', dim=DIM)

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
    W1 = VectorFunctionSpace('W1', domain)
    W2 = VectorFunctionSpace('W2', domain)
    T1 = VectorFunctionSpace('T1', domain)
    T2 = VectorFunctionSpace('T2', domain)

    v1 = TestFunction(V1, name='v1')
    v2 = TestFunction(V2, name='v2')
    u1 = TestFunction(U1, name='u1')
    u2 = TestFunction(U2, name='u2')
    w1 = VectorTestFunction(W1, name='w1')
    w2 = VectorTestFunction(W2, name='w2')
    t1 = VectorTestFunction(T1, name='t1')
    t2 = VectorTestFunction(T2, name='t2')

    x,y = V1.coordinates

    alpha = Constant('alpha')

    B1 = Boundary(r'\Gamma_1', domain)
    B2 = Boundary(r'\Gamma_2', domain)
    B3 = Boundary(r'\Gamma_3', domain)

    # ...
    with pytest.raises(UnconsistentError):
        expr = dot(grad(v1), grad(u1)) + v1*trace_0(u1, B1)
        a = BilinearForm((v1,u1), expr, name='a')
    # ...

    # ...
    with pytest.raises(UnconsistentError):
        expr = v1*trace_0(u1, B3) + v1*trace_1(grad(u1), B3) + u1*trace_0(v1, B2)
        a1 = BilinearForm((v1, u1), expr, name='a1')
    # ...

    # ...
    expr = dot(grad(v1), grad(u1))
    a_0 = BilinearForm((v1,u1), expr, name='a_0')

    expr = v1*trace_0(u1, B1)
    a_bnd = BilinearForm((v1, u1), expr, name='a_bnd')

    expr = a_0(v1,u1) + a_bnd(v1,u1)
    a = BilinearForm((v1,u1), expr, name='a')
    print(a)
    print(evaluate(a, verbose=True))
    print('')
#    import sys; sys.exit(0)
    # ...


    # ...
    expr = v1*trace_0(u1, B1) + v1*trace_1(grad(u1), B1)
    a1 = BilinearForm((v1, u1), expr, name='a1')

    expr = u1*trace_1(grad(v1), B2)
    a2 = BilinearForm((v1, u1), expr, name='a2')

    expr = a1(v2, u2) + a2(v2, u2)
    # as expected, we can define the form call, but we cannot create a Bilinear
    # form out of it.
    # TODO add assert on exception type
#    a = BilinearForm((v2, u2), expr, name='a')

    print(expr)
    print('')
    # ...

    # ...
    expr = v1*trace_0(u1, B1)
    a0 = BilinearForm((v1, u1), expr, name='a0')

    expr = v1*trace_1(grad(u1), B1)
    a1 = BilinearForm((v1, u1), expr, name='a1')

    expr = v1*trace_1(grad(u1), B2)
    a2 = BilinearForm((v1, u1), expr, name='a2')

    expr = dot(grad(u1), grad(v1))
    a3 = BilinearForm((v1, u1), expr, name='a3')

    expr = u1*v1
    a4 = BilinearForm((v1, u1), expr, name='a4')

    # TODO Mul not treated yet
#    expr = a0(v2, u2) + a1(v2, u2) + alpha * a2(v2, u2) + a3(v2, u2) + alpha*a4(v2, u2)
#    a = BilinearForm((v2, u2), expr, name='a')
##    print(expr)
#    print(evaluate(expr, verbose=True))
#    print('')
#
#    print(evaluate(a, verbose=True))
#    print('')


#    expr = a(v2, u2) + a1(v2, u2)
#    b = BilinearForm((v2, u2), expr, name='b')
#    print(b)
#    print(evaluate(b, verbose=True))
    # ...

    # ...
    g = Tuple(x**2, y**2)
    expr = v1*trace_1(g, B1)
    l1 = LinearForm(v1, expr, name='l1')
    print(l1)
#    print(atomize(l1))
#    print(evaluate(l1))
    print('')
    # ...

#==============================================================================
def test_boundary_2d_2():
    Omega_1 = InteriorDomain('Omega_1', dim=2)

    B1 = Boundary('B1', Omega_1)
    B2 = Boundary('B2', Omega_1)
    B3 = Boundary('B3', Omega_1)

    domain = Domain('Omega', interiors=[Omega_1],
                    boundaries=[B1, B2, B3])

    V = FunctionSpace('V', domain)
    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    x,y = V.coordinates

    alpha = Constant('alpha')

    # ...
    print('==== l0 ====')
    l0 = LinearForm(v, x*y*v, name='l0')

    print(evaluate(l0, verbose=VERBOSE))
    print('')
    # ...

    # ...
    print('==== l1 ====')
    g = Tuple(x**2, y**2)
    l1 = LinearForm(v, v*trace_1(g, domain.boundary))

    print(evaluate(l1, verbose=VERBOSE))
    print('')
    # ...

    # ...
    print('==== l2 ====')
    B_neumann = Union(B1, B2)
    g = Tuple(x**2, y**2)
    l2 = LinearForm(v, v*trace_1(g, B_neumann), name='l2')

    print(evaluate(l2, verbose=VERBOSE))
    print('')
    # ...

    # ...
    print('==== l3 ====')
    l3 = LinearForm(v, l2(v))

    assert(l3(v).__str__ == l2(v).__str__)

    print(evaluate(l3, verbose=VERBOSE))
    print('')
    # ...

    # ...
    print('==== l4 ====')
    l4 = LinearForm(v, l0(v) + l2(v))

    print(evaluate(l4, verbose=VERBOSE))
    print('')
    # ...

#    # ...
#    print('==== a1 ====')
#    a1 = BilinearForm((v, u), v*trace_0(u, domain.boundary))
#
#    print(evaluate(a1, verbose=VERBOSE))
#    print('')
#    # ...


#==============================================================================
def test_calls_2d():
    domain = Domain('Omega', dim=DIM)

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
    W1 = VectorFunctionSpace('W1', domain)
    W2 = VectorFunctionSpace('W2', domain)
    T1 = VectorFunctionSpace('T1', domain)
    T2 = VectorFunctionSpace('T2', domain)

    v1 = TestFunction(V1, name='v1')
    v2 = TestFunction(V2, name='v2')
    u1 = TestFunction(U1, name='u1')
    u2 = TestFunction(U2, name='u2')
    w1 = VectorTestFunction(W1, name='w1')
    w2 = VectorTestFunction(W2, name='w2')
    t1 = VectorTestFunction(T1, name='t1')
    t2 = VectorTestFunction(T2, name='t2')

    V = ProductSpace(V1, V2)
    U = ProductSpace(U1, U2)

    x,y = V1.coordinates

    alpha = Constant('alpha')

    F = Field('F', space=V1)

    # ...
    a1 = BilinearForm((v1, u1), u1*v1, name='a1')
    print(a1)
    print(atomize(a1))
    print(evaluate(a1))
    print('')

    expr = a1(v2, u2)
    a = BilinearForm((v2, u2), expr, name='a')
    print(a)
    print(atomize(a))
    print(evaluate(a))
    print('')
    # ...

    # ...
    a = BilinearForm((v1, u1), dot(grad(v1), grad(u1)), name='a')

    print(a)
    print(atomize(a))
    print(evaluate(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((v1, u1), dot(grad(v1), grad(u1)), name='a1')

    expr = a1(v2, u2)
    a = BilinearForm((v2, u2), expr, name='a')
    print(a)
    print(atomize(a))
    print(evaluate(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((v1, u1), u1*v1, name='a1')
    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1), name='a2')

    expr = a1(v2, u2) + a2(v2, u2)
    a = BilinearForm((v2, u2), expr, name='a')
    print(a)
    print(atomize(a))
    print(evaluate(a))
    print('')
    # ...

    # ...
    a1 = BilinearForm((v1, u1), u1*v1, name='a1')
    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1), name='a2')

    expr =  a1(v1, u2)
    print('> before = ', expr)
    expr = expr.subs(u2, u1)
    print('> after  = ', expr)
    print('')

    expr =  a1(v1, u2) + a1(v2, u2)
    print('> before = ', expr)
    expr = expr.subs(u2, u1)
    print('> after  = ', expr)
    print('')
    # ...

    # ...
    a1 = BilinearForm((v1, u1), u1*v1, name='a1')
    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1), name='a2')

    expr =  a1(v1, u2) + a2(v2, u1)
    a = BilinearForm(((v1,v2),(u1,u2)), expr, name='a')
    print(a)
    print(atomize(a))
    print(evaluate(a))
    print('')
    # ...

    # ...
    a = BilinearForm((w1, t1), rot(w1)*rot(t1) + div(w1)*div(t1), name='a')
    print(a)
    print(atomize(a))
    print(evaluate(a))
    # ...

    # ...
    a1 = BilinearForm((v1, u1), u1*v1, name='a1')
    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1), name='a2')
    a3 = BilinearForm((w1, t1), rot(w1)*rot(t1) + div(w1)*div(t1), name='a3')
    a4 = BilinearForm((w1, u1), div(w1)*u1, name='a4')

    expr = a3(w2,t2) + a2(v2,u2) + a4(w2,u2)
    a = BilinearForm(((w2,v2),(t2,u2)), expr, name='a')
    print(a)
    print(atomize(a))
    print(evaluate(a))
    # ...

    # ...
    a1 = BilinearForm((v1, u1), laplace(u1)*laplace(v1), name='a1')
    print(a1)
    print(atomize(a1))
    print(evaluate(a1))
    print('')
    # ...

    # ...
    a1 = BilinearForm((v1, u1), inner(hessian(u1),hessian(v1)), name='a1')
    print('================================')
    print(a1)
    print(atomize(a1))
    print(evaluate(a1))
    print('')
    # ...

    # ...
    l1 = LinearForm(v1, x*y*v1, name='11')

    expr = l1(v2)
    l = LinearForm(v2, expr, name='1')
    print(l)
    print(atomize(l))
    print(evaluate(l))
    # ...

    # ...
    l1 = LinearForm(v1, x*y*v1, name='l1')
    l2 = LinearForm(v2, cos(x+y)*v2, name='l2')

    expr = l1(u1) + l2(u2)
    l = LinearForm((u1,u2), expr, name='1')
    print(l)
    print(atomize(l))
    print(evaluate(l))
    # ...

    # ...
    l1 = LinearForm(v1, x*y*v1, name='l1')
    l2 = LinearForm(v2, cos(x+y)*v2, name='l2')

    expr = l1(u1) + alpha * l2(u2)
    l = LinearForm((u1,u2), expr, name='1')
    print(l)
    print(atomize(l))
    print(evaluate(l))
    # ...

    # ...
    l1 = LinearForm(v1, x*y*v1, name='l1')
    l2 = LinearForm(w1, div(w1), name='l2')

    expr = l1(v2) + l2(w2)
    l = LinearForm((v2,w2), expr, name='1')
    print(l)
    print(atomize(l))
    print(evaluate(l))
    # ...

    # ...
    I1 = Integral(x*y, domain, name='I1')

    print(I1)
    print(atomize(I1))
    print(evaluate(I1))
    # ...

    # ...
    expr = F - cos(2*pi*x)*cos(3*pi*y)
    expr = dot(grad(expr), grad(expr))
    I2 = Integral(expr, domain, name='I2')

    print(I2)
    print(atomize(I2))
    print(evaluate(I2))
    # ...

    # ...
    expr = F - cos(2*pi*x)*cos(3*pi*y)
    expr = dot(grad(expr), grad(expr))
    I2 = Integral(expr, domain, name='I2')

    print(I2)
    print(atomize(I2))
    print(evaluate(I2))
    # ...

    # ... stokes
    V = VectorFunctionSpace('V', domain)
    W = FunctionSpace('W', domain)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')
    p = TestFunction(W, name='p')
    q = TestFunction(W, name='q')

    a = BilinearForm((v,u), inner(grad(v), grad(u)), name='a')
    b = BilinearForm((v,p), div(v)*p, name='b')
    A = BilinearForm(((v,q),(u,p)), a(v,u) - b(v,p) + b(u,q), name='A')

    print(A)
    print(atomize(A))
    print(evaluate(A))
    # ...

#==============================================================================
def test_projection_2d():
    domain = Domain('Omega', dim=DIM)

    V = FunctionSpace('V', domain)
    x,y = domain.coordinates

    alpha = Constant('alpha')

    u = Projection(x**2+alpha*y, V, name='u')


#==============================================================================
def test_norm_2d():
    domain = Domain('Omega', dim=DIM)

    x,y = domain.coordinates
    V = FunctionSpace('V', domain)
    F = Field('F', V)

    # ...
    expr = x*y
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', evaluate(l2_norm_u))
    print('> h1 norm = ', evaluate(h1_norm_u))
    print('')
    # ...

    # ...
    expr = sin(pi*x)*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', evaluate(l2_norm_u))
    print('> h1 norm = ', evaluate(h1_norm_u))
    print('')
    # ...

    # ...
    expr = F-x*y
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', evaluate(l2_norm_u))
    print('> h1 norm = ', evaluate(h1_norm_u))
    print('')
    # ...

    # ...
    expr = F-sin(pi*x)*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', evaluate(l2_norm_u))
    print('> h1 norm = ', evaluate(h1_norm_u))
    print('')
    # ...

    # ...
    expr = F-sin(0.5*pi*(1.-x))*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', evaluate(l2_norm_u))
    print('> h1 norm = ', evaluate(h1_norm_u))
    print('')
    # ...

    # ...
    expr = F-cos(0.5*pi*x)*sin(pi*y)
    l2_norm_u = Norm(expr, domain, kind='l2')
    h1_norm_u = Norm(expr, domain, kind='h1')

    print('> l2 norm = ', evaluate(l2_norm_u))
    print('> h1 norm = ', evaluate(h1_norm_u))
    print('')
    # ...


#==============================================================================
def test_vector_2d_1():
    domain = Domain('Omega', dim=DIM)

    W1 = VectorFunctionSpace('W1', domain)
    T1 = VectorFunctionSpace('T1', domain)

    w1 = VectorTestFunction(W1, name='w1')
    t1 = VectorTestFunction(T1, name='t1')

    x,y = W1.coordinates

    F = VectorField(W1, 'F')

#    # ...
#    l1 = LinearForm(w1, dot(w1, F), name='l1')
#    print(l1)
#    print(atomize(l1))
#    print(evaluate(l1))
#    print('')
#    # ...
#
#    # ...
#    l2 = LinearForm(w1, rot(w1)*rot(F) + div(w1)*div(F), name='l2')
#    print(l2)
#    print(atomize(l2))
#    print(evaluate(l2))
#    print('')
#    # ...

    # ...
    f = Tuple(sin(pi*x)*sin(pi*y), sin(pi*x)*sin(pi*y))
    error = Matrix([F[0]-f[0], F[1]-f[1]])
    l2_norm = Norm(error, domain, kind='l2')
    print(l2_norm)
    print(atomize(l2_norm))
    print(evaluate(l2_norm))
    print('')
    # ...

    # ...
    f = Tuple(sin(pi*x)*sin(pi*y), sin(pi*x)*sin(pi*y))
    error = Matrix([F[0]-f[0], F[1]-f[1]])
    h1_norm = Norm(error, domain, kind='h1')
    print(h1_norm)
    print(atomize(h1_norm))
    print(evaluate(h1_norm))
    print('')
    # ...

#==============================================================================
def test_expr_mapping_2d():

    F = Mapping('F', DIM)
    patch = Domain('Omega', dim=DIM)
    domain = F(patch)

    V = FunctionSpace('V', domain)
    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    x,y = V.coordinates

    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    assert(a.mapping is F)

    l = LinearForm(v, x*y*v)
    assert(l.mapping is F)

#==============================================================================
def test_system_2d():

    domain = Square()

    V = FunctionSpace('V', domain)
    x,y = V.coordinates

    u,v,p,q = [TestFunction(V, name=i) for i in ['u','v','p','q']]

    a1,a2,b1,b2 = [Constant(i, real=True) for i in ['a1','a2','b1','b2']]

    # ...
    a = BilinearForm((v,u), dot(grad(u), grad(v)))
    m = BilinearForm((v,u), u*v)

    expr = a(p,u) + a1*m(p,u) + b1*m(p,v)  + a(q,v) + a2*m(q,u) +  b2*m(q,v)
    b = BilinearForm(((p,q), (u,v)), expr)

    print(evaluate(b, verbose=True))
    # ...

    # ...
    f1 = x*y
    f2 = x+y
    l1 = LinearForm(p, f1*p)
    l2 = LinearForm(q, f2*q)

    expr = l1(p) + l2(q)
    l = LinearForm((p,q), expr)

    print(evaluate(l, verbose=True))
    # ...

#==============================================================================
def test_curldiv_2d():

    domain = Square()

    W1 = VectorFunctionSpace('W1', domain)
    T1 = VectorFunctionSpace('T1', domain)

    w1 = VectorTestFunction(W1, name='w1')
    t1 = VectorTestFunction(T1, name='t1')

    mu = Constant('mu')

    # ...
    a = BilinearForm((w1, t1), rot(w1)*rot(t1) + mu*div(w1)*div(t1), name='a')
    print(a)
    print(atomize(a))
    print(evaluate(a))
    # ...


#==============================================================================
def test_calls_2d_2():

    domain = Square()

    V = FunctionSpace('V', domain)
    x,y = V.coordinates

    u,v = [TestFunction(V, name=i) for i in ['u','v']]
    Un = Field('Un', V)

    # ...
    a = BilinearForm((v,u), dot(grad(u), grad(v)))

    expr = a(v, Un)
    print(evaluate(expr, verbose=True))
    # ...

    # ...
    l = LinearForm(v, a(v, Un))

    print(evaluate(l, verbose=True))
    # ...

#==============================================================================
def test_calls_2d_3():

    domain = Square()

    V = FunctionSpace('V', domain)

    x,y = domain.coordinates

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

    # ...
    l = LinearForm((tau, sigma), dt*l1(tau))

    print(evaluate(l, verbose=True))
    # ...

#==============================================================================
def test_evaluation_2d_1():
    domain = Domain('Omega', dim=2)
    B_neumann = Boundary(r'\Gamma_1', domain)

    V = FunctionSpace('V', domain)
    W = VectorFunctionSpace('W', domain)

    p,q = [TestFunction(V, name=i) for i in ['p', 'q']]
    u,v = [VectorTestFunction(W, name=i) for i in ['u', 'v']]

    alpha = Constant('alpha')

    x,y = V.coordinates
    F = Field('F', space=V)

    a1 = BilinearForm((p, q), dot(grad(p), grad(q)))
    m  = BilinearForm((p, q), p*q)
    a2 = BilinearForm((p, q), a1(p,q) + alpha*m(p,q))
    a3 = BilinearForm((u, v), rot(u)*rot(v) + alpha*div(u)*div(v))

    a11 = BilinearForm((v,u), inner(grad(v), grad(u)))
    a12 = BilinearForm((v,p), div(v)*p)
    a4  = BilinearForm(((v,q),(u,p)), a11(v,u) - a12(v,p) + a12(u,q))

    l0 = LinearForm(p, F*p)
    l_neu = LinearForm(p, p*trace_1(grad(F), B_neumann))
    l = LinearForm(p, l0(p) + l_neu(p))

    # ...
    print(a1)
    print(evaluate(a1))
    print('')
    # ...

    # ...
    print(a2)
    print(evaluate(a2))
    print('')
    # ...

    # ...
    print(a3)
    print(evaluate(a3))
    print('')
    # ...

    # ...
    print(a4)
    print(evaluate(a4))
    print('')
    # ...

    # ...
    print(l)
    print(evaluate(l))
    print('')
    # ...

#==============================================================================
def test_evaluation_2d_2():
    domain = Square()
    x,y = domain.coordinates

    f0 = Tuple(2*pi**2*sin(pi*x)*sin(pi*y),
              2*pi**2*sin(pi*x)*sin(pi*y))

    f1 = cos(pi*x)*cos(pi*y)

    W = VectorFunctionSpace('W', domain)
    V = FunctionSpace('V', domain)
    X = ProductSpace(W, V)

    # TODO improve: naming are not given the same way
    F = VectorField(W, name='F')
    G = Field('G', V)

    u,v = [VectorTestFunction(W, name=i) for i in ['u', 'v']]
    p,q = [      TestFunction(V, name=i) for i in ['p', 'q']]

    a0 = BilinearForm((v,u), inner(grad(v), grad(u)))
    a1 = BilinearForm((q,p), p*q)
    a  = BilinearForm(((v,q),(u,p)), a0(v,u) + a1(q,p))

    l0 = LinearForm(v, dot(f0, v))
    l1 = LinearForm(q, f1*q)
    l  = LinearForm((v,q), l0(v) + l1(q))

    # ...
    print(a)
    print(evaluate(a))
    print('')
    # ...

    # ...
    print(l)
    print(evaluate(l))
    print('')
    # ...

#==============================================================================
def test_linearity_2d_1():
    domain = Square()
    x,y = domain.coordinates

    alpha = Constant('alpha')

    f1 = x*y
    f2 = x+y
    f = Tuple(f1, f2)

    V = FunctionSpace('V', domain)

    # TODO improve: naming are not given the same way
    G = Field('G', V)

    p,q = [TestFunction(V, name=i) for i in ['p', 'q']]

    #####################################
    # linear expressions
    #####################################
    # ...
    expr = p
    assert(is_linear_form(expr, p))
    # ...

    # ...
    expr = alpha*p
    assert(is_linear_form(expr, p))
    # ...

    # ...
    expr = dx(p)
    assert(is_linear_form(expr, p))
    # ...

    # ...
    expr = dot(grad(p), f)
    assert(is_linear_form(expr, p))
    # ...

    # ...
    expr = laplace(p)
    assert(is_linear_form(expr, p))
    # ...

    # ...
    expr = alpha*p + dot(grad(p), f) + dx(p) + laplace(p)
    assert(is_linear_form(expr, p))
    # ...
    #####################################

    #####################################
    # nonlinear expressions
    #####################################
    # ...
    with pytest.raises(UnconsistentLinearExpressionError):
        expr = p**2
        is_linear_form(expr, p)
    # ...

    # ...
    with pytest.raises(UnconsistentLinearExpressionError):
        expr = dot(grad(p), grad(p))
        is_linear_form(expr, p)
    # ...
    #####################################

#    expr = dot(grad(p), grad(p))
#    print(is_linear_form(expr, p))


#==============================================================================
def test_bilinearity_2d_1():
    domain = Square()
    x,y = domain.coordinates

    alpha = Constant('alpha')
    beta = Constant('beta')

    f1 = x*y
    f2 = x+y
    f = Tuple(f1, f2)

    V = FunctionSpace('V', domain)

    # TODO improve: naming are not given the same way
    G = Field('G', V)

    p,q = [TestFunction(V, name=i) for i in ['p', 'q']]

    #####################################
    # linear expressions
    #####################################
    # ...
    expr = p*q
    assert(is_bilinear_form(expr, (p,q)))
    # ...

    # ...
    expr = dot(grad(p), grad(q))
    assert(is_bilinear_form(expr, (p,q)))
    # ...

    # ...
    expr = alpha*dot(grad(p), grad(q)) + beta*p*q + laplace(p)*laplace(q)
    assert(is_bilinear_form(expr, (p,q)))
    # ...
    #####################################

    #####################################
    # nonlinear expressions
    #####################################
    # ...
    with pytest.raises(UnconsistentLinearExpressionError):
        expr = alpha*dot(grad(p**2), grad(q)) + beta*p*q
        is_bilinear_form(expr, (p,q))
    # ...
    #####################################

#    expr = p*q
#    print(is_bilinear_form(expr, (p,q)))

#==============================================================================
def test_linearity_2d_2():
    domain = Domain('Omega', dim=DIM)

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
    W1 = VectorFunctionSpace('W1', domain)
    W2 = VectorFunctionSpace('W2', domain)
    T1 = VectorFunctionSpace('T1', domain)
    T2 = VectorFunctionSpace('T2', domain)

    v1 = TestFunction(V1, name='v1')
    v2 = TestFunction(V2, name='v2')
    u1 = TestFunction(U1, name='u1')
    u2 = TestFunction(U2, name='u2')
    w1 = VectorTestFunction(W1, name='w1')
    w2 = VectorTestFunction(W2, name='w2')
    t1 = VectorTestFunction(T1, name='t1')
    t2 = VectorTestFunction(T2, name='t2')

    V = ProductSpace(V1, V2)
    U = ProductSpace(U1, U2)

    x,y = V1.coordinates

    alpha = Constant('alpha')

    F = Field('F', space=V1)

    # ...
    l1 = LinearForm(v1, x*y*v1, check=True)

    l = LinearForm(v2, l1(v2), check=True)
    # ...

    # ...
    l1 = LinearForm(v1, x*y*v1, check=True)
    l2 = LinearForm(v2, cos(x+y)*v2, check=True)

    l = LinearForm((u1,u2), l1(u1) + l2(u2), check=True)
    # ...

    # ...
    l1 = LinearForm(v1, x*y*v1, check=True)
    l2 = LinearForm(v2, cos(x+y)*v2, check=True)

    l = LinearForm((u1,u2), l1(u1) + alpha * l2(u2), check=True)
    # ...

    # ...
    l1 = LinearForm(v1, x*y*v1, check=True)
    l2 = LinearForm(w1, div(w1), check=True)

    l = LinearForm((v2,w2), l1(v2) + l2(w2), check=True)
    # ...

    ################################
    #    non bilinear forms
    ################################
    # ...
    with pytest.raises(UnconsistentLinearExpressionError):
        l = LinearForm(v1, x*y*v1**2, check=True)
    # ...

    # ...
    with pytest.raises(UnconsistentLinearExpressionError):
        l = LinearForm(v1, x*y, check=True)
    # ...
    ################################

#==============================================================================
def test_bilinearity_2d_2():
    domain = Domain('Omega', dim=DIM)

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
    W1 = VectorFunctionSpace('W1', domain)
    W2 = VectorFunctionSpace('W2', domain)
    T1 = VectorFunctionSpace('T1', domain)
    T2 = VectorFunctionSpace('T2', domain)

    v1 = TestFunction(V1, name='v1')
    v2 = TestFunction(V2, name='v2')
    u1 = TestFunction(U1, name='u1')
    u2 = TestFunction(U2, name='u2')
    w1 = VectorTestFunction(W1, name='w1')
    w2 = VectorTestFunction(W2, name='w2')
    t1 = VectorTestFunction(T1, name='t1')
    t2 = VectorTestFunction(T2, name='t2')

    V = ProductSpace(V1, V2)
    U = ProductSpace(U1, U2)

    x,y = V1.coordinates

    alpha = Constant('alpha')

    F = Field('F', space=V1)

    # ...
    a1 = BilinearForm((v1, u1), u1*v1, check=True)
    a  = BilinearForm((v2, u2), a1(v2, u2), check=True)
    # ...

    # ...
    a  = BilinearForm((v1, u1), dot(grad(v1), grad(u1)), check=True)
    # ...

    # ...
    a1 = BilinearForm((v1, u1), dot(grad(v1), grad(u1)), check=True)
    a  = BilinearForm((v2, u2), a1(v2, u2), check=True)
    # ...

    # ...
    a1 = BilinearForm((v1, u1), u1*v1)
    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1), check=True)
    a = BilinearForm((v2, u2), a1(v2, u2) + a2(v2, u2), check=True)
    # ...

    # ...
    a1 = BilinearForm((v1, u1), u1*v1, check=True)
    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1), check=True)
    # ...

    # ...
    a1 = BilinearForm((v1, u1), u1*v1, check=True)
    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1), check=True)
    a  = BilinearForm(((v1,v2),(u1,u2)), a1(v1, u2) + a2(v2, u1), check=True)
    # ...

    # ...
    a = BilinearForm((w1, t1), rot(w1)*rot(t1) + div(w1)*div(t1), check=True)
    # ...

    # ...
    a1 = BilinearForm((v1, u1), u1*v1, check=True)
    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1), check=True)
    a3 = BilinearForm((w1, t1), rot(w1)*rot(t1) + div(w1)*div(t1), check=True)
    a4 = BilinearForm((w1, u1), div(w1)*u1, check=True)

    a = BilinearForm(((w2,v2),(t2,u2)), a3(w2,t2) + a2(v2,u2) + a4(w2,u2), check=True)
    # ...

    # ...
    a1 = BilinearForm((v1, u1), laplace(u1)*laplace(v1), check=True)
    # ...

    # ...
    a1 = BilinearForm((v1, u1), inner(hessian(u1),hessian(v1)), check=True)
    # ...

    # ... stokes
    V = VectorFunctionSpace('V', domain)
    W = FunctionSpace('W', domain)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')
    p = TestFunction(W, name='p')
    q = TestFunction(W, name='q')

    a = BilinearForm((v,u), inner(grad(v), grad(u)), check=True)
    b = BilinearForm((v,p), div(v)*p, check=True)
    A = BilinearForm(((v,q),(u,p)), a(v,u) - b(v,p) + b(u,q), check=True)
    # ...

    ################################
    #    non bilinear forms
    ################################
    # ...
    with pytest.raises(UnconsistentLinearExpressionError):
        a  = BilinearForm((v1, u1), dot(grad(v1), grad(u1)) + v1, check=True)
    # ...

    # ...
    with pytest.raises(UnconsistentLinearExpressionError):
        a  = BilinearForm((v1, u1), v1**2*u1, check=True)
    # ...

    # ...
    with pytest.raises(UnconsistentLinearExpressionError):
        a  = BilinearForm((v1, u1), dot(grad(v1), grad(v1)), check=True)
    # ...
    ################################

#==============================================================================
#def test_nonlinear_2d_1():
#
#    domain = Square()
#
#    V = FunctionSpace('V', domain)
#    x,y = V.coordinates
#
#    u,v = [TestFunction(V, name=i) for i in ['u','v']]
#    Un = Field('Un', V)
#
#    from sympy import diff as sympy_diff
#
#    def diff(L, u):
#        ls = []
#        expr = evaluate(L)
#        for i in expr:
#            e = sympy_diff(i.expr, u)
#            ls.append(e)
#
#        if len(ls) == 1:
#            return ls[0]
#
#        else:
#            return ls
#
#    # ...
#    expr = Un**2 * dot(grad(Un), grad(v))
#    l = LinearForm(v, expr)
#
#    a = diff(l, Un)
#    print(a)
#    # ...

#==============================================================================
def test_linearize_2d_1():
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
    l = LinearForm(v1, F**2*v1, check=True)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, dot(grad(F), grad(F))*v1, check=True)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, exp(-F)*v1, check=True)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, cos(F)*v1, check=True)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, cos(F**2)*v1, check=True)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(v1, F**2*dot(grad(F), grad(v1)), check=True)
    a = linearize(l, F, trials='u1')
    print(a)
    # ...

    # ...
    l = LinearForm(w1, dot(rot(G), grad(G))*w1, check=True)
    a = linearize(l, G, trials='u1')
    print(a)
    # ...

#==============================================================================
def test_linearize_2d_2():
    domain = Domain('Omega', dim=DIM)
    x,y = domain.coordinates

    V1 = FunctionSpace('V1', domain)

    v1 = TestFunction(V1, name='v1')

    alpha = Constant('alpha')

    F = Field('F', space=V1)
    G = Field('G', space=V1)

    # ...
    l1 = LinearForm(v1, F**2*v1, check=True)
    l = LinearForm(v1, l1(v1))

    a = linearize(l, F, trials='u1')
    print(a)

    expected = linearize(l1, F, trials='u1')
    assert( linearize(l, F, trials='u1') == expected )
    # ...

#==============================================================================
def test_linearize_2d_3():
    """steady Euler equation."""
    domain = Domain('Omega', dim=DIM)
    x,y = domain.coordinates

    U = VectorFunctionSpace('U', domain)
    W =       FunctionSpace('W', domain)

#    u   = VectorTestFunction(U, name='u')
#    rho =       TestFunction(W, name='rho')
#    p   =       TestFunction(W, name='p')

    v   = VectorTestFunction(U, name='v')
    phi =       TestFunction(W, name='phi')
    q   =       TestFunction(W, name='q')

    U_0   = VectorField(U, name='U_0')
    Rho_0 =       Field('Rho_0', W)
    P_0   =       Field('P_0', W)

    # ...
    expr = div(Rho_0*U_0) * phi
#    expr = Rho_0*dot(U_0, grad(phi))
    l1 = LinearForm(phi, expr, check=True)

    # TODO
#    expr = Rho_0*convect(U_0, grad(U_0))
#    l2 = LinearForm(phi, expr, check=True)

    expr = dot(U_0, grad(P_0)) * q + P_0 * div(U_0) * q
    l3 = LinearForm(q, expr, check=True)
    # ...

    a1 = linearize(l1, [Rho_0, U_0], trials=['d_rho', 'd_u'])
    print(a1)

    a3 = linearize(l3, [P_0, U_0], trials=['d_p', 'd_u'])
    print(a3)

#    l = LinearForm((phi, q), l1(phi) + l3(q))
#    a = linearize(l, [Rho_0, U_0, P_0], trials=['d_rho', 'd_u', 'd_p'])
#    print(a)

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()

#test_linearize_2d_1()
#test_linearize_2d_2()
#test_linearize_2d_3()

#test_linearity_2d_1()
#test_linearity_2d_2()
#test_bilinearity_2d_1()
#test_bilinearity_2d_2()

#test_boundary_2d_1()
#test_boundary_2d_2()
#test_projection_2d()
#test_norm_2d()
#test_vector_2d_1()
#test_expr_mapping_2d()
#test_system_2d()
#test_curldiv_2d()
#test_calls_2d()
#test_calls_2d_2()
#test_calls_2d_3()
#test_evaluation_2d_1()
#test_evaluation_2d_2()
