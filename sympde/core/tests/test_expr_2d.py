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


DIM = 2
domain = Domain('Omega', dim=DIM)

# ...
def test_boundary_2d():
    print('============ test_boundary_2d =============')

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

    assert(isinstance(B1 + B2, UnionBoundary))
    assert(isinstance(B1 + B2 + B3, UnionBoundary))
    assert(isinstance(-(B1 + B2), ComplementBoundary))

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

    expr = a0(v2, u2) + a1(v2, u2) + alpha * a2(v2, u2) + a3(v2, u2) + alpha*a4(v2, u2)
    a = BilinearForm((v2, u2), expr, name='a')
#    print(expr)
    print(evaluate(expr, verbose=True))
    print('')

    print(evaluate(a, verbose=True))
    print('')


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

# ...
def test_calls_2d():
    print('============ test_calls_2d =============')

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
# ...

# ...
def test_equation_2d():
    print('============ test_equation_2d =============')

    V = FunctionSpace('V', domain)
    U = FunctionSpace('U', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    x,y = domain.coordinates

    alpha = Constant('alpha')

    B1 = Boundary(r'\Gamma_1', domain)
    B2 = Boundary(r'\Gamma_2', domain)
    B3 = Boundary(r'\Gamma_3', domain)

    # ... bilinear/linear forms
    expr = dot(grad(v), grad(u))
    a1 = BilinearForm((v,u), expr, name='a1')

    expr = v*u
    a2 = BilinearForm((v,u), expr, name='a2')

    expr = v*trace_0(u, B1)
    a_B1 = BilinearForm((v, u), expr, name='a_B1')

    expr = x*y*v
    l1 = LinearForm(v, expr, name='l1')

    expr = cos(x+y)*v
    l2 = LinearForm(v, expr, name='l2')

    expr = x*y*trace_0(v, B2)
    l_B2 = LinearForm(v, expr, name='l_B2')
    # ...

    # ...
    with pytest.raises(UnconsistentLhsError):
        equation = Equation(a1, l1(v))

    with pytest.raises(UnconsistentLhsError):
        equation = Equation(l1(v), l1(v))

    with pytest.raises(UnconsistentLhsError):
        equation = Equation(a1(v,u) + alpha*a2(v,u), l1(v))
    # ...

    # ...
    with pytest.raises(UnconsistentRhsError):
        equation = Equation(a1(v,u), l1)

    with pytest.raises(UnconsistentRhsError):
        equation = Equation(a1(v,u), a1(v,u))

    with pytest.raises(UnconsistentRhsError):
        equation = Equation(a1(v,u), l1(v) + l2(v))
    # ...

    # ...
    equation = Equation(a1(v,u), l1(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + alpha*a2(v,u))
    equation = Equation(a(v,u), l1(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    equation = Equation(a(v,u), l1(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    l = LinearForm(v, l1(v) + l2(v))
    equation = Equation(a(v,u), l(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    l = LinearForm(v, l1(v) + alpha*l_B2(v))
    equation = Equation(a(v,u), l(v))
    # ...

    # ... using bc
    equation = Equation(a1(v,u), l1(v), bc=DirichletBC(B1))
    # ...

    # ... using bc
    equation = Equation(a1(v,u), l1(v), bc=DirichletBC(ComplementBoundary(B1)))
    # ...

    # ...
    with pytest.raises(UnconsistentBCError):
        a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
        equation = Equation(a(v,u), l1(v), bc=DirichletBC(B1))
    # ...

    # ...
    with pytest.raises(UnconsistentBCError):
        l = LinearForm(v, l1(v) + alpha*l_B2(v))
        equation = Equation(a1(v,u), l(v), bc=DirichletBC(B2))
    # ...


def test_projection_2d():
    print('============ test_projection_2d =============')

    V = FunctionSpace('V', domain)
    x,y = domain.coordinates

    alpha = Constant('alpha')

    u = Projection(x**2+alpha*y, V, name='u')


def test_norm_2d():
    print('============ test_norm_2d =============')

    x,y = domain.coordinates

    expr = x*y
    l2_norm_u = Norm(expr, domain, kind='l2', name='u')
    h1_norm_u = Norm(expr, domain, kind='h1', name='u')

    print('> l2 norm = ', evaluate(l2_norm_u))
    print('> l2 norm = ', evaluate(h1_norm_u))

# ...
def test_tensorize_2d():
    print('============ test_tensorize_2d =============')

    V = FunctionSpace('V', domain)
    U = FunctionSpace('U', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    x,y = domain.coordinates

    alpha = Constant('alpha')

    # ...
    expr = dot(grad(v), grad(u))
    a = BilinearForm((v,u), expr, name='a')
    print(a)
    print(tensorize(a))
    print('')
#    import sys; sys.exit(0)
    # ...

# .....................................................
if __name__ == '__main__':

#    test_calls_2d()
#    test_boundary_2d()
#    test_equation_2d()
#    test_projection_2d()
#    test_norm_2d()
    test_tensorize_2d()

