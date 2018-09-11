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

from sympde.core import dx, dy, dz
from sympde.core import Constant
from sympde.core import Field
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import FunctionSpace
from sympde.core import ProductSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import BilinearForm, LinearForm, Integral
from sympde.core import atomize, normalize, matricize
from sympde.core import evaluate
from sympde.core import tensorize
from sympde.core import inv_normalize
from sympde.core import Mass, Stiffness, Advection, AdvectionT
from sympde.core import Unknown
from sympde.core import FormCall
from sympy.physics.quantum import TensorProduct

# ...
def test_atomize_2d_1():
    print('============ test_atomize_2d_1 =============')

    V = FunctionSpace('V', ldim=2)

    v = TestFunction(V, name='v')
    w = TestFunction(V, name='w')
    c = Constant('c')
    F = Field('F', space=V)

    # ... expressions that can be normalized (valid for a weak formulation)
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
def test_normalize_2d_1():
    print('============ test_normalize_2d_1 =============')

    V = FunctionSpace('V', ldim=2)
    U = FunctionSpace('U', ldim=2)

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
    assert(normalize(grad(v), basis={v: 'Ni'}) == Tuple(Ni_x, Ni_y))
    assert(normalize(grad(c*v), basis={v: 'Ni'}) == Tuple(c*Ni_x, c*Ni_y))
    assert(normalize(dot(b, grad(v)), basis={v: 'Ni'}) == Ni_x*bx + Ni_y*by)
    assert(normalize(dot(b, grad(v)) + c*v, basis={v: 'Ni'}) == Ni_x*bx + Ni_y*by + c*Ni)
    assert(normalize(dot(f, grad(v)), basis={v: 'Ni'}) == Ni_x*f1(x,y) + Ni_y*f2(x,y))
    assert(normalize(dot(Tuple(2, 3), grad(v)), basis={v: 'Ni'}) == 2*Ni_x + 3*Ni_y)
    assert(normalize(grad(F*v), basis={v: 'Ni'}) == Tuple(F*Ni_x + Ni*dx(F),
                                                          F*Ni_y + Ni*dy(F)))
    # TODO debug
#    assert(normalize(A*grad(v), basis={v: 'Ni'}) == 2*Ni_x + 3*Ni_y)

    assert(normalize(dot(grad(v), grad(u)), basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y)
    assert(normalize(dot(grad(v), grad(u)) + c*v*u, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + c*Ni*Nj)
    assert(normalize(dot(grad(F*v), grad(u)), basis={v: 'Ni', u: 'Nj'}) == Nj_x*(F*Ni_x + Ni*dx(F)) + Nj_y*(F*Ni_y + Ni*dy(F)))
    # ...

#    expr = dot(A, grad(v))
#    expr = div(dot(A, grad(v)))
#    print('> input         >>> {0}'.format(expr))

#    print('> normal form   >>> {0}'.format(normalize(expr, basis={v: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={v: 'Ni', u: 'Nj'})))
# ...


# ...
def test_evaluate_2d_1():
    print('============ test_evaluate_2d_1 =============')

    V = FunctionSpace('V', ldim=2)
    U = FunctionSpace('U', ldim=2)

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

    V = FunctionSpace('V', ldim=2, is_block=True, shape=2)

    v = VectorTestFunction(V, name='v')

    assert(atomize(rot(v)) == -dx(v[1]) + dy(v[0]))
    assert(atomize(div(v)) == dx(v[0]) + dy(v[1]))

#    expr = div(v)
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
# ...

# ...
def test_normalize_2d_2():
    print('============ test_normalize_2d_2 =============')

    V = FunctionSpace('V', ldim=2, is_block=True, shape=2)
    U = FunctionSpace('U', ldim=2, is_block=True, shape=2)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(U, name='u')

    Ni = IndexedBase('Ni', shape=2)
    Ni_x = IndexedBase('Ni_x', shape=2)
    Ni_y = IndexedBase('Ni_y', shape=2)

    Nj = IndexedBase('Nj', shape=2)
    Nj_x = IndexedBase('Nj_x', shape=2)
    Nj_y = IndexedBase('Nj_y', shape=2)

    assert(normalize(v[0], basis={v: 'Ni'}) == Ni[0])
    assert(normalize(dx(v[0]), basis={v: 'Ni'}) == Ni_x[0])
    assert(normalize(div(v), basis={v: 'Ni'}) == Ni_x[0] + Ni_y[1])
    assert(normalize(rot(v), basis={v: 'Ni'}) == -Ni_x[1] + Ni_y[0])

    expected = Tuple(Matrix([[dx(v[0]), dx(v[1])]]), Matrix([[dy(v[0]), dy(v[1])]]))
    assert(normalize(grad(v), basis={v: 'Ni'}) == expected)

    assert(normalize(v[0]*u[0], basis={v: 'Ni', u: 'Nj'}) == Ni[0]*Nj[0])
    assert(normalize(v[1]*dx(u[0]), basis={v: 'Ni', u: 'Nj'}) == Ni[1]*Nj_x[0])
    assert(normalize(dy(v[0])*u[1], basis={v: 'Ni', u: 'Nj'}) == Ni_y[0]*Nj[1])
    assert(normalize(dx(v[1])*dy(u[1]), basis={v: 'Ni', u: 'Nj'}) == Ni_x[1]*Nj_y[1])

    expected = (Ni_x[0] + Ni_y[1]) * (Nj_x[0] + Nj_y[1])
    assert(normalize(div(v) * div(u), basis={v: 'Ni', u: 'Nj'}) == expected)

    expected = (-Ni_x[1] + Ni_y[0]) * (-Nj_x[1] + Nj_y[0])
    assert(normalize(rot(v) * rot(u), basis={v: 'Ni', u: 'Nj'}) == expected)

    expected = Ni_x[0]*Nj_x[0] + Ni_x[1]*Nj_x[1] + Ni_y[0]*Nj_y[0] + Ni_y[1]*Nj_y[1]
    assert(normalize(inner(grad(v), grad(u)), basis={v: 'Ni', u: 'Nj'}) == expected)

#    expr = inner(grad(v), grad(u))
#    print('> input         >>> {0}'.format(expr))
#
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={v: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={v: 'Ni', u: 'Nj'})))
# ...

# ...
def test_matricize_2d_2():
    print('============ test_matricize_2d_2 =============')

    V = FunctionSpace('V', ldim=2, is_block=True, shape=2)
    U = FunctionSpace('U', ldim=2, is_block=True, shape=2)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(U, name='u')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    # ...
    expr = v[0]*u[0]
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni*Nj, 0], [0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = v[1]*dx(u[0])
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[0, Ni*Nj_x], [0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dy(v[0])*u[1]
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[0, 0], [Ni_y*Nj, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dx(v[1])*dy(u[1])
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[0, 0], [0, Ni_x*Nj_y]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = div(v) * div(u)
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x, Ni_y*Nj_x], [Ni_x*Nj_y, Ni_y*Nj_y]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = rot(v) * rot(u)
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni_y*Nj_y, -Ni_x*Nj_y], [-Ni_y*Nj_x, Ni_x*Nj_x]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = div(v) * div(u) + rot(v) * rot(u)
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x + Ni_y*Nj_y, -Ni_x*Nj_y + Ni_y*Nj_x],
                       [Ni_x*Nj_y - Ni_y*Nj_x, Ni_x*Nj_x + Ni_y*Nj_y]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = c1 * div(v) * div(u) + rot(v) * rot(u)
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x*c1 + Ni_y*Nj_y, -Ni_x*Nj_y + Ni_y*Nj_x*c1],
                       [Ni_x*Nj_y*c1 - Ni_y*Nj_x, Ni_x*Nj_x + Ni_y*Nj_y*c1]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = c1 * div(v) * div(u) + c2 * rot(v) * rot(u)
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x*c1 + Ni_y*Nj_y*c2, -Ni_x*Nj_y*c2 + Ni_y*Nj_x*c1],
                       [Ni_x*Nj_y*c1 - Ni_y*Nj_x*c2, Ni_x*Nj_x*c2 + Ni_y*Nj_y*c1]])
    assert(matricize(expr) == expected)
    # ...

#    expr = c1 * div(v) * div(u) + rot(v) * rot(u)
#    print('> input         >>> {0}'.format(expr))
#    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
#    print('> matricize     >>> {0}'.format(matricize(expr)))
# ...

# ...
def test_matricize_2d_3():
    print('============ test_matricize_2d_3 =============')

    V = FunctionSpace('V', ldim=2, is_block=True, shape=2)
    U = FunctionSpace('U', ldim=2)

    v = VectorTestFunction(V, name='v')
    u = TestFunction(U, name='u')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    # ...
    expr = v[0]*u
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni*Nj], [0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dot(v, grad(u))
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni*Nj_x], [Ni*Nj_y]])
    assert(matricize(expr) == expected)
    # ...

#    expr = v[0]*u
#    print('> input         >>> {0}'.format(expr))
#    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
#    print('> matricize     >>> {0}'.format(matricize(expr)))

# ...
def test_calls_2d_3():
    print('============ test_calls_2d_3 =============')

    V1 = FunctionSpace('V1', ldim=2)
    V2 = FunctionSpace('V2', ldim=2)
    U1 = FunctionSpace('U1', ldim=2)
    U2 = FunctionSpace('U2', ldim=2)
    W1 = FunctionSpace('W1', ldim=2, is_block=True, shape=2)
    W2 = FunctionSpace('W2', ldim=2, is_block=True, shape=2)
    T1 = FunctionSpace('T1', ldim=2, is_block=True, shape=2)
    T2 = FunctionSpace('T2', ldim=2, is_block=True, shape=2)

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

    v1v2 = VectorTestFunction(V, name='v1v2')
    u1u2 = VectorTestFunction(U, name='u1u2')

    # ...
    a1 = BilinearForm((v1, u1), u1*v1, name='a1')

    expr = a1(v2, u2)
    a = BilinearForm((v2, u2), expr, name='a')
    print(a)
    print(atomize(a))
#    print(evaluate(a))
    print('')
    # ...

#    # ...
#    a1 = BilinearForm((v1, u1), dot(grad(v1), grad(u1)), name='a1')
#
#    expr = a1(v2, u2)
#    a = BilinearForm((v2, u2), expr, name='a')
#    print(a)
#    print(atomize(a))
##    print(evaluate(a))
#    print('')
#    # ...
#
#    # ...
#    a1 = BilinearForm((v1, u1), u1*v1, name='a1')
#    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1), name='a2')
#
#    expr = a1(v2, u2) + a2(v2, u2)
#    a = BilinearForm((v2, u2), expr, name='a')
#    print(a)
#    print(atomize(a))
##    print(evaluate(a))
#    print('')
#    # ...
#
#    # ...
#    a1 = BilinearForm((v1, u1), u1*v1, name='a1')
#    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1), name='a2')
#
#    expr =  a1(v1, u2)
#    print('> before = ', expr)
#    expr = expr.subs(u2, u1)
#    print('> after  = ', expr)
#    print('')
#
#    expr =  a1(v1, u2) + a1(v2, u2)
#    print('> before = ', expr)
#    expr = expr.subs(u2, u1)
#    print('> after  = ', expr)
#    print('')
#    # ...
#
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

#    # ...
#    a1 = BilinearForm((v1, u1), u1*v1, name='a1')
#    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1), name='a2')
#    a3 = BilinearForm((w1, t1), rot(w1)*rot(t1) + div(w1)*div(t1), name='a3')
#    a4 = BilinearForm((w1, u1), div(w1)*u1, name='a4')
#
#    expr = a3(w2,t2) + a2(v2,u2) + a4(w2,u2)
#    a = BilinearForm(((w2,v2),(t2,u2)), expr, name='a')
#    print(a)
#    print(atomize(a))
#    print(evaluate(a3))
##    print(evaluate(a))
#    # ...

    import sys; sys.exit(0)

    # ...
    l1 = LinearForm(v1, x*y*v1, name='11')

    expr = l1(v2)
    # ...

    # ...
    l1 = LinearForm(v1, x*y*v1, name='l1')
    l2 = LinearForm(v2, cos(x+y)*v2, name='l2')

    expr = l1(u1) + l2(u2)
    # ...
# ...

# ...
def test_evaluate_2d_3():
    print('============ test_evaluate_2d_3 =============')

    V1 = FunctionSpace('V1', ldim=2)
    U1 = FunctionSpace('U1', ldim=2)
    V2 = FunctionSpace('V2', ldim=2)
    U2 = FunctionSpace('U2', ldim=2)

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
def test_inv_normalize_2d_1():
    print('============ test_inv_normalize_2d_1 =============')

    U = FunctionSpace('U', ldim=2)
    V = FunctionSpace('V', ldim=2)

    u = TestFunction(U, name='u')
    v = TestFunction(V, name='v')

    Ni, Ni_x, Ni_y, Ni_xx, Ni_xy, Ni_yy = symbols('Ni Ni_x Ni_y Ni_xx Ni_xy Ni_yy')
    Nj, Nj_x, Nj_y, Nj_xx, Nj_xy, Nj_yy = symbols('Nj Nj_x Nj_y Nj_xx Nj_xy Nj_yy')

#    expr = Nj*Ni_x + Nj_x*Ni + Nj*Ni
#    expr = Nj*Ni_xx + Nj_x*Ni + Nj*Ni
#
#    expr = inv_normalize(expr, {Ni: v, Nj: u})
#    print(expr)


# ...
#def test_bilinear_form_2d_10():
#    print('============ test_bilinear_form_2d_10 =============')
#
#    U = FunctionSpace('U', ldim=2)
#    V = FunctionSpace('V', ldim=2)
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

    V = FunctionSpace('V', ldim=2)

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

    V = FunctionSpace('V', ldim=2)

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

    V = FunctionSpace('V', ldim=2)
    V_0 = FunctionSpace('V_0', ldim=1, coordinates=['x'])
    V_1 = FunctionSpace('V_1', ldim=1, coordinates=['y'])

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

    V = FunctionSpace('V', ldim=2, is_block=True, shape=2)
#    V_0 = FunctionSpace('V_0', ldim=1, coordinates=['x'])
#    V_1 = FunctionSpace('V_1', ldim=1, coordinates=['y'])

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

    v = Unknown('v', ldim=2)
    c = Constant('c')

    # ... expressions that can be normalized (valid for a weak formulation)
    assert(atomize(grad(v)) == Tuple(dx(v),
                                      dy(v)))
    assert(atomize(grad(c*v)) == Tuple(c*dx(v),
                                        c*dy(v)))
    # ...
# ...

# .....................................................
if __name__ == '__main__':
#    test_atomize_2d_1()
#    test_normalize_2d_1()
#    test_evaluate_2d_1()
#
#    test_atomize_2d_2()
#    test_normalize_2d_2()
#    test_matricize_2d_2()
#
##    test_bilinear_form_2d_10() # TODO not working, since args are the same
#    test_linear_form_2d_10()
#    test_function_form_2d_10()
#
#    test_matricize_2d_3()
#    test_evaluate_2d_3()
#
#    test_tensorize_2d_1()
#
#    test_inv_normalize_2d_1()
#    test_tensorize_2d_2()
#
    test_calls_2d_3()
#
#    test_unknown_2d_1()
