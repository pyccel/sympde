# coding: utf-8

# TODO - split the asserts between algebraic and weak formulations ones
#      - add assert for grad in vector case
# TODO: - __call__ examples are not working anymore

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
from sympde.core import FunctionSpace
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
from sympde.core import Trace, trace_0, trace_1

DIM = 2
domain = Domain('Omega', dim=DIM)

# ...
def test_atomize_2d_1():
    print('============ test_atomize_2d_1 =============')

    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    w = TestFunction(V, name='w')
    c = Constant('c')
    F = Field('F', space=V)

    # ...
    assert(atomize(grad(v)) == Tuple(dx(v),
                                      dy(v)))
    assert(atomize(grad(c*v)) == Tuple(c*dx(v),
                                        c*dy(v)))
    assert(atomize(grad(F*v)) == Tuple(F*dx(v) + v*dx(F),
                                        F*dy(v) + v*dy(F)))

    assert(atomize(dot(grad(v), grad(w))) == dx(v)*dx(w) + dy(v)*dy(w))
    # ...

#    expr = grad(F*v)
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
# ...

# ...
def test_evaluate_2d_1():
    print('============ test_evaluate_2d_1 =============')

    V = FunctionSpace('V', domain)
    U = FunctionSpace('U', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    x,y = V.coordinates

    c = Constant('c')
    F = Field('F', space=V)
    f1 = Function('f1')
    f2 = Function('f2')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    bx, by = symbols('bx by')
    b = Tuple(bx, by)

    f = Tuple(f1(x,y), f2(x,y))

    a00 = Constant('a00')
    a10 = Constant('a10')
    a01 = Constant('a01')
    a11 = Constant('a11')
    A = Matrix([[a00, a01], [a10, a11]])

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)) + c*v*u)
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + c*Ni*Nj)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)) + F*v*u)
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + F*Ni*Nj)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(F*v), grad(u)))
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == F*Ni_x*Nj_x + F*Ni_y*Nj_y + Ni*Nj_x*dx(F) + Ni*Nj_y*dy(F))
    # ...

# ...

# ...
def test_atomize_2d_2():
    print('============ test_atomize_2d_2 =============')

    V = FunctionSpace('V', domain, is_block=True, shape=2)

    v = VectorTestFunction(V, name='v')

    assert(atomize(rot(v)) == -dx(v[1]) + dy(v[0]))
    assert(atomize(div(v)) == dx(v[0]) + dy(v[1]))

#    expr = div(v)
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
# ...

# ...
def test_boundary_2d_3():
    print('============ test_boundary_2d_3 =============')

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
    W1 = FunctionSpace('W1', domain, is_block=True, shape=2)
    W2 = FunctionSpace('W2', domain, is_block=True, shape=2)
    T1 = FunctionSpace('T1', domain, is_block=True, shape=2)
    T2 = FunctionSpace('T2', domain, is_block=True, shape=2)

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
    # these two expressions are returning, as expected, an error.
    # TODO add assert on exception type
#    expr = v1*trace_0(u1, B3) + v1*trace_1(grad(u1), B3) + u1*trace_0(v1, B2)
#    expr = dot(grad(u1), grad(v1)) + v1*trace_1(grad(u1), B3)

    expr = v1*trace_0(u1, B1) + v1*trace_1(grad(u1), B1)

    a1 = BilinearForm((v1, u1), expr, name='a1')
    print(a1)
    print(evaluate(a1))
    print('')
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
def test_calls_2d_3():
    print('============ test_calls_2d_3 =============')

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
    W1 = FunctionSpace('W1', domain, is_block=True, shape=2)
    W2 = FunctionSpace('W2', domain, is_block=True, shape=2)
    T1 = FunctionSpace('T1', domain, is_block=True, shape=2)
    T2 = FunctionSpace('T2', domain, is_block=True, shape=2)

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
    I1 = Integral(x*y, domain, coordinates=[x,y], name='I1')

    print(I1)
    print(atomize(I1))
    print(evaluate(I1))
    # ...

    # ...
    expr = F - cos(2*pi*x)*cos(3*pi*y)
    expr = dot(grad(expr), grad(expr))
    I2 = Integral(expr, domain, coordinates=[x,y], name='I2')

    print(I2)
    print(atomize(I2))
    print(evaluate(I2))
    # ...

    # ...
    expr = F - cos(2*pi*x)*cos(3*pi*y)
    expr = dot(grad(expr), grad(expr))
    I2 = Integral(expr, domain, coordinates=[x,y], name='I2')

    print(I2)
    print(atomize(I2))
    print(evaluate(I2))
    # ...

    # ... stokes
    V = FunctionSpace('V', domain, is_block=True, shape=DIM)
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
def test_evaluate_2d_3():
    print('============ test_evaluate_2d_3 =============')

    V1 = FunctionSpace('V1', domain)
    U1 = FunctionSpace('U1', domain)
    V2 = FunctionSpace('V2', domain)
    U2 = FunctionSpace('U2', domain)

    v1 = TestFunction(V1, name='v1')
    u1 = TestFunction(U1, name='u1')
    v2 = TestFunction(V2, name='v2')
    u2 = TestFunction(U2, name='u2')

    V = ProductSpace(V1, V2)
    U = ProductSpace(U1, U2)

    v1v2 = VectorTestFunction(V, name='v1v2')
    u1u2 = VectorTestFunction(U, name='u1u2')

    c = Constant('c')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    basis = {v1v2: 'Ni', u1u2: 'Nj'}

    # ...
    expr = v1*u1 + dx(v2)*dx(u2)
    a = BilinearForm(((v1, v2), (u1, u2)), expr)

    expected = Matrix([[Ni*Nj, 0], [0, Ni_x*Nj_x]])
    assert(evaluate(a, basis=basis) == expected)
    # ...

    # ...
    expr = v1*u1 + dy(v2)*u1 + v1*dx(u2) + dx(v2)*dx(u2)
    a = BilinearForm(((v1, v2), (u1, u2)), expr)

    expected = Matrix([[Ni*Nj, Ni_y*Nj], [Ni*Nj_x, Ni_x*Nj_x]])
    assert(evaluate(a, basis=basis) == expected)
    # ...

#    expr = v1*u1 + dy(v2)*u1 + v1*dx(u2) + dx(v2)*dx(u2)
#    expr = BilinearForm(((v1, v2), (u1, u2)), expr)
#    print('> input         >>> {0}'.format(expr))
#    print('> normal form   >>> {0}'.format(evaluate(expr, basis=basis)))
# ...

# ...
#def test_bilinear_form_2d_10():
#    print('============ test_bilinear_form_2d_10 =============')
#
#    U = FunctionSpace('U', domain)
#    V = FunctionSpace('V', domain)
#
#    u = TestFunction(U, name='u')
#    v = TestFunction(V, name='v')
#
#    u1 = TestFunction(U, name='u1')
#    v1 = TestFunction(V, name='v1')
#
#    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
#    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')
#
#    c1 = Symbol('c1')
#    c2 = Symbol('c2')
#
#    a = BilinearForm((v,u), inner(grad(u), grad(v)))
#    b = BilinearForm((v,u), u*v)
#    adv = BilinearForm((v,u), dx(u)*v)
#
#    # ...
#    expected = Ni*Nj + Ni_x*Nj_x + Ni_y*Nj_y
#    assert(evaluate(a + b, basis={v: 'Nj', u: 'Ni'}) == expected)
#    # ...
#
#    # ...
#    expected = 2*Ni_x*Nj_x + 2*Ni_y*Nj_y
#    assert(evaluate(2 * a, basis={v: 'Nj', u: 'Ni'}) == expected)
#    # ...
#
#    # ...
#    expected = c1*(Ni_x*Nj_x + Ni_y*Nj_y)
#    assert(evaluate(c1*a, basis={v: 'Nj', u: 'Ni'}) == expected)
#    # ...
#
#    # ...
#    expected = Ni*Nj*c2 + c1*(Ni_x*Nj_x + Ni_y*Nj_y)
#    assert(evaluate(c1*a + c2*b, basis={v: 'Nj', u: 'Ni'}) == expected)
#    # ...
#
#    # ...
#    expected = c1*(Ni_x*Nj_x + Ni_y*Nj_y) + c2*(Ni*Nj + Ni_x*Nj)
#    assert(evaluate(c1*a  + c2*(b + adv), basis={v: 'Nj', u: 'Ni'}) == expected)
#    # ...
#
##    expr = c1*a  + c2*(b + adv)
##    print('> input      >>> {0}'.format(expr))
##    print('> evaluated  >>> {0}'.format(evaluate(expr, basis={v: 'Nj', u: 'Ni'}) ))
##    print('')
# ...

# ...
def test_linear_form_2d_10():
    print('============ test_linear_form_2d_10 =============')

    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')

    x,y = V.coordinates
    f = Function('f')
    g = Function('g')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    bx, by = symbols('bx by')
    b = Tuple(bx, by)
    fg = Tuple(f(x,y), g(x,y))

    a = LinearForm(v, cos(2*pi*x)*cos(4*pi*y)*v)

    # ...
    expected = cos(2*pi*x)*cos(4*pi*y)*Ni
    assert(evaluate(LinearForm(v, cos(2*pi*x)*cos(4*pi*y)*v),
                    basis={v: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x,y)*Ni
    assert(evaluate(LinearForm(v, f(x,y)*v),
                    basis={v: 'Ni'}) == expected)
    # ...

    # ...
    expected = bx*Ni_x + by*Ni_y + f(x,y)*Ni
    assert(evaluate(LinearForm(v, dot(b, grad(v)) + f(x,y)*v),
                    basis={v: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x,y)*Ni_x + g(x,y)*Ni_y
    assert(evaluate(LinearForm(v, dot(fg, grad(v))),
                    basis={v: 'Ni'}) == expected)
    # ...

#    expr =
#    print('> input      >>> {0}'.format(expr))
#    print('> evaluated  >>> {0}'.format(evaluate(expr, basis={v: 'Ni'}) ))
#    print('')
# ...

# ...
def test_function_form_2d_10():
    print('============ test_function_form_2d_10 =============')

    V = FunctionSpace('V', domain)

    F = Field('F', space=V)

    x,y = V.coordinates

    f = Function('f')
    g = Function('g')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    bx, by = symbols('bx by')
    b = Tuple(bx, by)
    fg = Tuple(f(x,y), g(x,y))

    # ...
    expected = 4*pi**2*sin(2*pi*x)**2*cos(3*pi*y)**2 + 9*pi**2*sin(3*pi*y)**2*cos(2*pi*x)**2
    e = cos(2*pi*x)*cos(3*pi*y)
    assert(evaluate(Integral(dot(grad(e), grad(e)), coordinates=[x,y])) == expected)
    # ...

    # ...
    expected = F - cos(2*pi*x)*cos(3*pi*y)
    assert(evaluate(Integral(F-cos(2*pi*x)*cos(3*pi*y))) == expected)
    # ...

    # ...
    expected = (F - cos(2*pi*x)*cos(3*pi*y))**2
    assert(evaluate(Integral((F - cos(2*pi*x)*cos(3*pi*y))**2)) == expected)
    # ...

    # ...
    expected = (dx(F) + 2*pi*sin(2*pi*x)*cos(3*pi*y))**2 + (dy(F) + 3*pi*sin(3*pi*y)*cos(2*pi*x))**2
    e = F -cos(2*pi*x)*cos(3*pi*y)
    assert(evaluate(Integral(dot(grad(e), grad(e)))) == expected)
    # ...

#    e = F -cos(2*pi*x)*cos(3*pi*y)
#    expr = Integral(dot(grad(e), grad(e)))
#    print('> input      >>> {0}'.format(expr))
#    print('> evaluated  >>> {0}'.format(evaluate(expr) ))
#    print('')
# ...

# ...
def test_tensorize_2d_1():
    print('============ test_tensorize_2d_1 =============')

    V = FunctionSpace('V', domain)
    V_0 = FunctionSpace('V_0', domain, coordinates=['x'])
    V_1 = FunctionSpace('V_1', domain, coordinates=['y'])

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    v0 = TestFunction(V_0, name='v0')
    u0 = TestFunction(V_0, name='u0')

    v1 = TestFunction(V_1, name='v1')
    u1 = TestFunction(V_1, name='u1')

    c = Constant('c')

    bx = Constant('bx')
    by = Constant('by')
    b = Tuple(bx, by)

    # ...
    expected = TensorProduct(Mass(v1, u1), Mass(v0, u0))
    assert(tensorize(BilinearForm((v,u), u*v)) == expected)
    # ...

    # ...
    expected = TensorProduct(Mass(v1, u1), Stiffness(v0, u0))
    assert(tensorize(BilinearForm((v,u), dx(u)*dx(v))) == expected)
    # ...

    # ...
    expected = TensorProduct(Advection(v1, u1), Mass(v0, u0))
    assert(tensorize(BilinearForm((v,u), dy(u) * v)) == expected)
    # ...

    # ...
    expected =  TensorProduct(Mass(v1,u1), Advection(v0,u0))
    assert(tensorize(BilinearForm((v,u), dx(u) * v)) == expected)
    # ...

    # ...
    expected = TensorProduct(Mass(v1,u1), Stiffness(v0,u0)) + TensorProduct(Stiffness(v1,u1), Mass(v0,u0))
    assert(tensorize(BilinearForm((v,u), dot(grad(v), grad(u)))) == expected)
    # ...

    # ...
    expected = (TensorProduct(Advection(v1,u1), Mass(v0,u0)) +
               TensorProduct(Mass(v1,u1), Advection(v0,u0)) +
               TensorProduct(Mass(v1,u1), Stiffness(v0,u0)) +
               TensorProduct(Stiffness(v1,u1), Mass(v0,u0)))
    assert(tensorize(BilinearForm((v,u), dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v)) == expected)
    # ...

    # ...
    expected = (TensorProduct(bx, Mass(v1,u1), AdvectionT(v0,u0)) +
               TensorProduct(by, AdvectionT(v1,u1), Mass(v0,u0)))
    assert(tensorize(BilinearForm((v,u), dot(b, grad(v)) * u)) == expected)
    # ...

    # ...
    expected = (TensorProduct(bx**2,Mass(v1,u1),Stiffness(v0,u0)) +
                TensorProduct(bx,by,Advection(v1,u1),AdvectionT(v0,u0)) +
                TensorProduct(bx,by,AdvectionT(v1,u1),Advection(v0,u0)) +
                TensorProduct(by**2,Stiffness(v1,u1),Mass(v0,u0)))
    assert(tensorize(BilinearForm((v,u), dot(b, grad(v)) * dot(b, grad(u)))) == expected)
    # ...

#    expr = dot(b, grad(v)) * u
#    expr = BilinearForm((v,u), expr)
#
#    print('> input         >>> {0}'.format(expr))
#    print('> tensorized    >>> {0}'.format(tensorize(expr)))
# ...


# ...
def test_tensorize_2d_2():
    print('============ test_tensorize_2d_2 =============')

    V = FunctionSpace('V', domain, is_block=True, shape=2)
#    V_0 = FunctionSpace('V_0', domain, coordinates=['x'])
#    V_1 = FunctionSpace('V_1', domain, coordinates=['y'])

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')

#    v0 = TestFunction(V_0, name='v0')
#    u0 = TestFunction(V_0, name='u0')
#
#    v1 = TestFunction(V_1, name='v1')
#    u1 = TestFunction(V_1, name='u1')

    c = Constant('c')

    bx = Constant('bx')
    by = Constant('by')
    b = Tuple(bx, by)

#    # ...
#    expected = Mass(v1, u1)*Mass(v0, u0)
#    assert(tensorize(BilinearForm((v,u), div(v) * div(u))) == expected)
#    # ...

#    expr = div(v) * div(u) + rot(v) * rot(u)
#    expr = BilinearForm((v,u), expr)
#
#    print('> input         >>> {0}'.format(expr))
#    print('> tensorized    >>> {0}'.format(tensorize(expr)))
# ...

# ...
def test_unknown_2d_1():
    print('============ test_unknown_2d_1 =============')

    domain = Domain('Omega', dim=DIM)

    v = Unknown('v', domain)
    c = Constant('c')

    # ...
    assert(atomize(grad(v)) == Tuple(dx(v),
                                      dy(v)))
    assert(atomize(grad(c*v)) == Tuple(c*dx(v),
                                        c*dy(v)))
    # ...
# ...

# .....................................................
if __name__ == '__main__':
#    test_atomize_2d_1()
#    test_evaluate_2d_1()
#
#    test_atomize_2d_2()
#
##    test_bilinear_form_2d_10() # TODO not working, since args are the same
#    test_linear_form_2d_10()
#    test_function_form_2d_10()
#
#    test_evaluate_2d_3()
#
#    test_tensorize_2d_1()
#
#    test_tensorize_2d_2()
#
#    test_unknown_2d_1()



    test_calls_2d_3()
    test_boundary_2d_3()
