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
from sympde.core import Mass, Stiffness, Advection, AdvectionT
from sympde.core import Unknown
from sympy.physics.quantum import TensorProduct


# ...
def test_atomize_3d_1():
    print('============ test_atomize_3d_1 =============')

    V = FunctionSpace('V', ldim=3)

    v = TestFunction(V, name='v')
    w = TestFunction(V, name='w')
    c = Constant('c')
    F = Field('F', space=V)

    # ... expressions that can be normalized (valid for a weak formulation)
    assert(atomize(grad(v)) == Tuple(dx(v),
                                      dy(v),
                                      dz(v)))
    assert(atomize(grad(c*v)) == Tuple(c*dx(v),
                                        c*dy(v),
                                        c*dz(v)))
    assert(atomize(grad(F*v)) == Tuple(F*dx(v) + v*dx(F),
                                        F*dy(v) + v*dy(F),
                                        F*dz(v) + v*dz(F)))

    assert(atomize(dot(grad(v), grad(w))) == dx(v)*dx(w) + dy(v)*dy(w) + dz(v)*dz(w))
    # ...

#    expr = grad(F*v)
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
# ...

# ...
def test_normalize_3d_1():
    print('============ test_normalize_3d_1 =============')

    V = FunctionSpace('V', ldim=3)
    U = FunctionSpace('U', ldim=3)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    x,y,z = V.coordinates

    c = Constant('c')
    F = Field('F', space=V)
    f1 = Function('f1')
    f2 = Function('f2')
    f3 = Function('f3')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    bx, by, bz = symbols('bx by bz')
    b = Tuple(bx, by, bz)

    f = Tuple(f1(x,y,z), f2(x,y,z), f3(x,y,z))

    a00 = Constant('a00')
    a10 = Constant('a10')
    a20 = Constant('a20')
    a01 = Constant('a01')
    a11 = Constant('a11')
    a21 = Constant('a21')
    a02 = Constant('a02')
    a12 = Constant('a12')
    a22 = Constant('a22')
    A = Matrix([[a00, a01, a02], [a10, a11, a12], [a20, a21, a22]])

    # ...
    assert(normalize(grad(v), basis={v: 'Ni'}) == Tuple(Ni_x, Ni_y, Ni_z))
    assert(normalize(grad(c*v), basis={v: 'Ni'}) == Tuple(c*Ni_x, c*Ni_y, c*Ni_z))
    assert(normalize(dot(b, grad(v)), basis={v: 'Ni'}) == Ni_x*bx + Ni_y*by + Ni_z*bz)
    assert(normalize(dot(b, grad(v)) + c*v, basis={v: 'Ni'}) == Ni_x*bx + Ni_y*by + Ni_z*bz + c*Ni)
    assert(normalize(dot(f, grad(v)), basis={v: 'Ni'}) == Ni_x*f1(x,y,z) + Ni_y*f2(x,y,z) + Ni_z*f3(x,y,z))
    assert(normalize(dot(Tuple(2, 3, 4), grad(v)), basis={v: 'Ni'}) == 2*Ni_x + 3*Ni_y + 4*Ni_z)
    assert(normalize(grad(F*v), basis={v: 'Ni'}) == Tuple(F*Ni_x + Ni*dx(F),
                                                          F*Ni_y + Ni*dy(F),
                                                          F*Ni_z + Ni*dz(F)))

    assert(normalize(dot(grad(v), grad(u)), basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z)
    assert(normalize(dot(grad(v), grad(u)) + c*v*u, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z + c*Ni*Nj)
    assert(normalize(dot(grad(F*v), grad(u)), basis={v: 'Ni', u: 'Nj'}) == Nj_x*(F*Ni_x + Ni*dx(F)) + Nj_y*(F*Ni_y + Ni*dy(F)) + Nj_z*(F*Ni_z + Ni*dz(F)))
    # ...

#    expr = dot(A, grad(v))
#    expr = div(dot(A, grad(v)))
#    print('> input         >>> {0}'.format(expr))

#    print('> normal form   >>> {0}'.format(normalize(expr, basis={v: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={v: 'Ni', u: 'Nj'})))
# ...

# ...
def test_evaluate_3d_1():
    print('============ test_evaluate_3d_1 =============')

    V = FunctionSpace('V', ldim=3)
    U = FunctionSpace('U', ldim=3)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    x,y,z = V.coordinates

    c = Constant('c')
    F = Field('F', space=V)
    f1 = Function('f1')
    f2 = Function('f2')
    f3 = Function('f3')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    bx, by, bz = symbols('bx by bz')
    b = Tuple(bx, by, bz)

    f = Tuple(f1(x,y,z), f2(x,y,z), f3(x,y,z))

    a00 = Constant('a00')
    a10 = Constant('a10')
    a20 = Constant('a20')
    a01 = Constant('a01')
    a11 = Constant('a11')
    a21 = Constant('a21')
    a02 = Constant('a02')
    a12 = Constant('a12')
    a22 = Constant('a22')
    A = Matrix([[a00, a01, a02], [a10, a11, a12], [a20, a21, a22]])

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)) + c*v*u)
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z + c*Ni*Nj)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)) + F*v*u)
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z + F*Ni*Nj)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(F*v), grad(u)))
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == F*Ni_x*Nj_x + F*Ni_y*Nj_y + F*Ni_z*Nj_z + Ni*Nj_x*dx(F) + Ni*Nj_y*dy(F) + Ni*Nj_z*dz(F))
    # ...

# ...

# ...
def test_atomize_3d_2():
    print('============ test_atomize_3d_2 =============')

    V = FunctionSpace('V', ldim=3, is_vector=True, shape=3)

    v = VectorTestFunction(V, name='v')

    assert(atomize(curl(v)) == Tuple( dy(v[2]) - dz(v[1]),
                                      -dx(v[2]) + dz(v[0]),
                                       dx(v[1]) - dy(v[0])))
    assert(atomize(div(v)) == dx(v[0]) + dy(v[1]) + dz(v[2]))

#    expr = curl(v)
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
# ...

# ...
def test_normalize_3d_2():
    print('============ test_normalize_3d_2 =============')

    V = FunctionSpace('V', ldim=3, is_block=True, shape=3)
    U = FunctionSpace('U', ldim=3, is_block=True, shape=3)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(U, name='u')

    Ni = IndexedBase('Ni', shape=3)
    Ni_x = IndexedBase('Ni_x', shape=3)
    Ni_y = IndexedBase('Ni_y', shape=3)
    Ni_z = IndexedBase('Ni_z', shape=3)

    Nj = IndexedBase('Nj', shape=3)
    Nj_x = IndexedBase('Nj_x', shape=3)
    Nj_y = IndexedBase('Nj_y', shape=3)
    Nj_z = IndexedBase('Nj_z', shape=3)

    assert(normalize(v[0], basis={v: 'Ni'}) == Ni[0])
    assert(normalize(dx(v[0]), basis={v: 'Ni'}) == Ni_x[0])
    assert(normalize(div(v), basis={v: 'Ni'}) == Ni_x[0] + Ni_y[1] + Ni_z[2])
    assert(normalize(curl(v), basis={v: 'Ni'}) == Tuple( Ni_y[2] - Ni_z[1],
                                                        -Ni_x[2] + Ni_z[0],
                                                         Ni_x[1] - Ni_y[0]))

    expected = Tuple(Matrix([[dx(v[0]), dx(v[1]), dx(v[2])]]),
                     Matrix([[dy(v[0]), dy(v[1]), dy(v[2])]]),
                     Matrix([[dz(v[0]), dz(v[1]), dz(v[2])]]))
    assert(normalize(grad(v), basis={v: 'Ni'}) == expected)

    assert(normalize(v[0]*u[0], basis={v: 'Ni', u: 'Nj'}) == Ni[0]*Nj[0])
    assert(normalize(v[1]*dx(u[0]), basis={v: 'Ni', u: 'Nj'}) == Ni[1]*Nj_x[0])
    assert(normalize(dy(v[0])*u[1], basis={v: 'Ni', u: 'Nj'}) == Ni_y[0]*Nj[1])
    assert(normalize(dx(v[1])*dy(u[1]), basis={v: 'Ni', u: 'Nj'}) == Ni_x[1]*Nj_y[1])
    assert(normalize(dx(v[1])*dz(u[2]), basis={v: 'Ni', u: 'Nj'}) == Ni_x[1]*Nj_z[2])

    expected = (Ni_x[0] + Ni_y[1] + Ni_z[2])*(Nj_x[0] + Nj_y[1] + Nj_z[2])
    assert(normalize(div(v) * div(u), basis={v: 'Ni', u: 'Nj'}) == expected)

    expected = ((Ni_x[1] - Ni_y[0])*(Nj_x[1] - Nj_y[0])
             + (-Ni_x[2] + Ni_z[0])*(-Nj_x[2] + Nj_z[0])
             + (Ni_y[2] - Ni_z[1])*(Nj_y[2] - Nj_z[1]))
    assert(normalize(dot(curl(v), curl(u)), basis={v: 'Ni', u: 'Nj'}) == expected)

    expected = (Ni_x[0]*Nj_x[0] + Ni_x[1]*Nj_x[1] + Ni_x[2]*Nj_x[2] +
                Ni_y[0]*Nj_y[0] + Ni_y[1]*Nj_y[1] + Ni_y[2]*Nj_y[2] +
                Ni_z[0]*Nj_z[0] + Ni_z[1]*Nj_z[1] + Ni_z[2]*Nj_z[2])
    assert(normalize(inner(grad(v), grad(u)), basis={v: 'Ni', u: 'Nj'}) == expected)

#    expr = inner(grad(v), grad(u))
#    print('> input         >>> {0}'.format(expr))
#
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={v: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={v: 'Ni', u: 'Nj'})))
# ...

# ...
def test_matricize_3d_2():
    print('============ test_matricize_3d_2 =============')

    V = FunctionSpace('V', ldim=3, is_block=True, shape=3)
    U = FunctionSpace('U', ldim=3, is_block=True, shape=3)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(U, name='u')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    # ...
    expr = v[0]*u[0]
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni*Nj, 0, 0], [0, 0, 0], [0, 0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = v[1]*dx(u[0])
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[0, Ni*Nj_x, 0], [0, 0, 0], [0, 0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dy(v[0])*u[1]
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[0, 0, 0], [Ni_y*Nj, 0, 0], [0, 0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dx(v[1])*dy(u[1])
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[0, 0, 0], [0, Ni_x*Nj_y, 0], [0, 0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = div(v) * div(u)
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x, Ni_y*Nj_x, Ni_z*Nj_x],
                       [Ni_x*Nj_y, Ni_y*Nj_y, Ni_z*Nj_y],
                       [Ni_x*Nj_z, Ni_y*Nj_z, Ni_z*Nj_z]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dot(curl(v), curl(u))
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni_y*Nj_y + Ni_z*Nj_z, -Ni_x*Nj_y, -Ni_x*Nj_z],
                       [-Ni_y*Nj_x, Ni_x*Nj_x + Ni_z*Nj_z, -Ni_y*Nj_z],
                       [-Ni_z*Nj_x, -Ni_z*Nj_y, Ni_x*Nj_x + Ni_y*Nj_y]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = div(v) * div(u) + dot(curl(v), curl(u))
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z, -Ni_x*Nj_y + Ni_y*Nj_x, -Ni_x*Nj_z + Ni_z*Nj_x],
                       [Ni_x*Nj_y - Ni_y*Nj_x, Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z, -Ni_y*Nj_z + Ni_z*Nj_y],
                       [Ni_x*Nj_z - Ni_z*Nj_x, Ni_y*Nj_z - Ni_z*Nj_y, Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = c1 * div(v) * div(u) + dot(curl(v), curl(u))
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected =  Matrix([[Ni_x*Nj_x*c1 + Ni_y*Nj_y + Ni_z*Nj_z,
                         -Ni_x*Nj_y + Ni_y*Nj_x*c1,
                         -Ni_x*Nj_z + Ni_z*Nj_x*c1],
                        [Ni_x*Nj_y*c1 - Ni_y*Nj_x,
                         Ni_x*Nj_x + Ni_y*Nj_y*c1 + Ni_z*Nj_z,
                         -Ni_y*Nj_z + Ni_z*Nj_y*c1],
                        [Ni_x*Nj_z*c1 - Ni_z*Nj_x,
                         Ni_y*Nj_z*c1 - Ni_z*Nj_y,
                         Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z*c1]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = c1 * div(v) * div(u) + c2 * dot(curl(v), curl(u))
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x*c1 + Ni_y*Nj_y*c2 + Ni_z*Nj_z*c2,
                        -Ni_x*Nj_y*c2 + Ni_y*Nj_x*c1,
                        -Ni_x*Nj_z*c2 + Ni_z*Nj_x*c1],
                       [Ni_x*Nj_y*c1 - Ni_y*Nj_x*c2,
                        Ni_x*Nj_x*c2 + Ni_y*Nj_y*c1 + Ni_z*Nj_z*c2,
                        -Ni_y*Nj_z*c2 + Ni_z*Nj_y*c1],
                       [Ni_x*Nj_z*c1 - Ni_z*Nj_x*c2,
                        Ni_y*Nj_z*c1 - Ni_z*Nj_y*c2,
                        Ni_x*Nj_x*c2 + Ni_y*Nj_y*c2 + Ni_z*Nj_z*c1]])
    assert(matricize(expr) == expected)
    # ...

#    expr = c1 * div(v) * div(u) + c2 * dot(curl(v), curl(u))
#    print('> input         >>> {0}'.format(expr))
#    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
#    print('> matricize     >>> {0}'.format(matricize(expr)))
# ...

# ...
#def test_bilinear_form_3d_10():
#    print('============ test_bilinear_form_3d_10 =============')
#
#    U = FunctionSpace('U', ldim=3)
#    V = FunctionSpace('V', ldim=3)
#
#    u = TestFunction(U, name='u')
#    v = TestFunction(V, name='v')
#
#    u1 = TestFunction(U, name='u1')
#    v1 = TestFunction(V, name='v1')
#
#    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
#    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')
#
#    c1 = Symbol('c1')
#    c2 = Symbol('c2')
#
#    a = BilinearForm((v,u), inner(grad(u), grad(v)))
#    b = BilinearForm((v,u), u*v)
#    adv = BilinearForm((v,u), dx(u)*v)
#
#    # ...
#    expected = Ni*Nj + Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z
#    assert(evaluate(a + b, basis={v: 'Nj', u: 'Ni'}) == expected)
#    # ...
#
#    # ...
#    expected = Ni*Nj + Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z
#    assert(evaluate(a + b, basis={v: 'Nj', u: 'Ni'}) == expected)
#    # ...
#
#    # ...
#    expected = 2*Ni_x*Nj_x + 2*Ni_y*Nj_y + 2*Ni_z*Nj_z
#    assert(evaluate(2 * a, basis={v: 'Nj', u: 'Ni'}) == expected)
#    # ...
#
#    # ...
#    expected = c1*(Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z)
#    assert(evaluate(c1*a, basis={v: 'Nj', u: 'Ni'}) == expected)
#    # ...
#
#    # ...
#    expected = Ni*Nj*c2 + c1*(Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z)
#    assert(evaluate(c1*a + c2*b, basis={v: 'Nj', u: 'Ni'}) == expected)
#    # ...
#
#    # ...
#    expected = c1*(Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z) + c2*(Ni*Nj + Ni_x*Nj)
#    assert(evaluate(c1*a  + c2*(b + adv), basis={v: 'Nj', u: 'Ni'}) == expected)
#    # ...
#
#    # ...
#    assert(evaluate(a(u1, v1), basis={v: 'Nj', u: 'Ni'}) == evaluate(a(v1, u1), basis={v: 'Nj', u: 'Ni'}))
#    # ...
#
##    # ... TODO debug
##    expected = Ni_x*Nj
##    assert(evaluate(adv(v1, u1), basis={v: 'Nj', u: 'Ni'}) == expected)
##
##    expected = Nj_x*Ni
##    assert(evaluate(adv(u1, v1), basis={v: 'Nj', u: 'Ni'}) == expected)
##    # ...
#
##    expr = c1*a  + c2*(b + adv)
##    print('> input      >>> {0}'.format(expr))
##    print('> evaluated  >>> {0}'.format(evaluate(expr, basis={v: 'Nj', u: 'Ni'}) ))
##    print('')
# ...

# ...
def test_linear_form_3d_10():
    print('============ test_linear_form_3d_10 =============')

    V = FunctionSpace('V', ldim=3)

    v = TestFunction(V, name='v')
    x,y,z = V.coordinates

    f = Function('f')
    g = Function('g')
    r = Function('r')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    bx, by, bz = symbols('bx by bz')
    b = Tuple(bx, by, bz)
    fgr = Tuple(f(x,y,z), g(x,y,z), r(x,y,z))

    a = LinearForm(v, cos(2*pi*x)*cos(4*pi*y)*cos(5*pi*z)*v)

    # ...
    expected = cos(2*pi*x)*cos(4*pi*y)*cos(5*pi*z)*Ni
    assert(evaluate(LinearForm(v, cos(2*pi*x)*cos(4*pi*y)*cos(5*pi*z)*v),
                    basis={v: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x,y,z)*Ni
    assert(evaluate(LinearForm(v, f(x,y,z)*v),
                    basis={v: 'Ni'}) == expected)
    # ...

    # ...
    expected = bx*Ni_x + by*Ni_y + bz*Ni_z + f(x,y,z)*Ni
    assert(evaluate(LinearForm(v, dot(b, grad(v)) + f(x,y,z)*v),
                    basis={v: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x,y,z)*Ni_x + g(x,y,z)*Ni_y + r(x,y,z)*Ni_z
    assert(evaluate(LinearForm(v, dot(fgr, grad(v))),
                    basis={v: 'Ni'}) == expected)
    # ...

#    expr = c1*a  + c2*(b + adv)
#    print('> input      >>> {0}'.format(expr))
#    print('> evaluated  >>> {0}'.format(evaluate(expr, basis={v: 'Ni'}) ))
#    print('')
# ...

# ...
def test_function_form_3d_10():
    print('============ test_function_form_3d_10 =============')

    V = FunctionSpace('V', ldim=3)

    F = Field('F', space=V)

    x,y,z = V.coordinates

    f = Function('f')
    g = Function('g')
    r = Function('r')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    bx, by, bz = symbols('bx by bz')
    b = Tuple(bx, by, bz)
    fgr = Tuple(f(x,y,z), g(x,y,z), r(x,y,z))

    # ...
    expected = cos(2*pi*x)*cos(3*pi*y)*cos(5*pi*z)
    assert(evaluate(Integral(cos(2*pi*x)*cos(3*pi*y)*cos(5*pi*z), coordinates=[x,y,z])) == expected)
    # ...

    # ...
    expected = x**2 + y**2 + 1
    e = x*y + z
    assert(evaluate(Integral(dot(grad(e), grad(e)), coordinates=[x,y,z])) == expected)
    # ...

    # ...
    expected = F - cos(2*pi*x)*cos(3*pi*y)*cos(5*pi*z)
    assert(evaluate(Integral(F-cos(2*pi*x)*cos(3*pi*y)*cos(5*pi*z))) == expected)
    # ...

    # ...
    expected = (F - x*y - z)**2
    assert(evaluate(Integral((F - x*y - z)**2)) == expected)
    # ...

    # ...
    expected = dx(F)**2 + dy(F)**2 + dz(F)**2
    assert(evaluate(Integral(dot(grad(F), grad(F)))) == expected)
    # ...

    # ...
    expected = (-x + dy(F))**2 + (-y + dx(F))**2 + (dz(F) - 1)**2
    e = F - (x*y + z)
    assert(evaluate(Integral(dot(grad(e), grad(e)), coordinates=[x,y,z])) == expected)
    # ...

    # ... TODO debug. => infinite recursion!!! why?
    #          must be a problem with Mul treatements
#    e = cos(2*pi*x)*cos(3*pi*y)*cos(5*pi*z)
#    e = cos(2*pi*x)*cos(3*pi*y)*z
#    e = x*y*z
    # ...

#    e = F - (x*y + z)
#    expr = Integral(dot(grad(e), grad(e)), coordinates=[x,y,z])
#    print('> input      >>> {0}'.format(expr))
#    print('> evaluated  >>> {0}'.format(evaluate(expr) ))
#    print('')
# ...

# ...
def test_matricize_3d_3():
    print('============ test_matricize_3d_3 =============')

    V = FunctionSpace('V', ldim=3, is_block=True, shape=3)
    U = FunctionSpace('U', ldim=3)

    v = VectorTestFunction(V, name='v')
    u = TestFunction(U, name='u')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    # ...
    expr = v[0]*u
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni*Nj], [0], [0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dot(v, grad(u))
    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
    expected = Matrix([[Ni*Nj_x], [Ni*Nj_y], [Ni*Nj_z]])
    assert(matricize(expr) == expected)
    # ...

#    expr = v[0]*u
#    print('> input         >>> {0}'.format(expr))
#    expr = normalize(expr, basis={v: 'Ni', u: 'Nj'})
#    print('> matricize     >>> {0}'.format(matricize(expr)))

# ...
def test_calls_3d_3():
    print('============ test_calls_3d_3 =============')

    V1 = FunctionSpace('V1', ldim=3)
    V2 = FunctionSpace('V2', ldim=3)
    U1 = FunctionSpace('U1', ldim=3)
    U2 = FunctionSpace('U2', ldim=3)
    W1 = FunctionSpace('W1', ldim=3, is_block=True, shape=3)
    W2 = FunctionSpace('W2', ldim=3, is_block=True, shape=3)
    T1 = FunctionSpace('T1', ldim=3, is_block=True, shape=3)
    T2 = FunctionSpace('T2', ldim=3, is_block=True, shape=3)

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

    x,y,z = V1.coordinates

    v1v2 = VectorTestFunction(V, name='v1v2')
    u1u2 = VectorTestFunction(U, name='u1u2')

    # ...
    a1 = BilinearForm((v1, u1), u1*v1)

    assert(u2*v2 == a1(v2, u2))
    # ...

    # ...
    a1 = BilinearForm((v1, u1), u1*v1)
    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1))

    assert(u2*v2 + dx(u2)*dx(v2) == a1(v2, u2) + a2(v2, u2))
    # ...

    # ...
    a1 = BilinearForm((v1, u1), u1*v1)
    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1))

    expr =  dx(u1u2[0])*dx(v1v2[1]) + u1u2[1]*v1v2[0]
    expected = BilinearForm((v1v2, u1u2), expr)
    result = BilinearForm(((v1,v2), (u1,u2)), a1(v1, u2) + a2(v2, u1))
    # we can only compare the strings, since the spaces are not the same
    assert(str(result) == str(expected))
    # ...

    # ...
    a1 = BilinearForm((v1, u1), u1*v1)
    a2 = BilinearForm((v1, u1), dx(u1)*dx(v1))
    a3 = BilinearForm((w1, t1), dot(curl(w1), curl(t1)) + div(w1)*div(t1))
    a4 = BilinearForm((w1, u1), div(w1)*u1)

    expr = a3(w2,t2) + a2(v2,u2) + a4(w2,u2)
    b = BilinearForm(((w2,v2), (t2,u2)), expr)
#    print(expr)
#    print(b)
    # ...

    # ...
    l1 = LinearForm(v1, x*y*z*v1)
#    print(l1(v2))
    # ...
# ...

# ...
def test_evaluate_3d_3():
    print('============ test_evaluate_3d_3 =============')

    V1 = FunctionSpace('V1', ldim=3)
    U1 = FunctionSpace('U1', ldim=3)
    V2 = FunctionSpace('V2', ldim=3)
    U2 = FunctionSpace('U2', ldim=3)

    v1 = TestFunction(V1, name='v1')
    u1 = TestFunction(U1, name='u1')
    v2 = TestFunction(V2, name='v2')
    u2 = TestFunction(U2, name='u2')

    V = ProductSpace(V1, V2)
    U = ProductSpace(U1, U2)

    v1v2 = VectorTestFunction(V, name='v1v2')
    u1u2 = VectorTestFunction(U, name='u1u2')

    c = Constant('c')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    basis = {v1v2: 'Ni', u1u2: 'Nj'}

    # ...
    expr = v1*u1 + dz(v2)*dz(u2)
    a = BilinearForm(((v1, v2), (u1, u2)), expr)

    expected = Matrix([[Ni*Nj, 0], [0, Ni_z*Nj_z]])
    assert(evaluate(a, basis=basis) == expected)
    # ...

    # ...
    expr = v1*u1 + dx(v2)*dx(u1) + dy(v1)*dy(u2) + dz(v2)*dz(u2)
    a = BilinearForm(((v1, v2), (u1, u2)), expr)

    expected = Matrix([[Ni*Nj, Ni_x*Nj_x], [Ni_y*Nj_y, Ni_z*Nj_z]])
    assert(evaluate(a, basis=basis) == expected)
    # ...

#    expr = v1*u1 + dx(v2)*dx(u1) + dy(v1)*dy(u2) + dz(v2)*dz(u2)
#    expr = BilinearForm(((v1, v2), (u1, u2)), expr)
#    print('> input         >>> {0}'.format(expr))
#    print('> normal form   >>> {0}'.format(evaluate(expr, basis=basis)))
# ...

# ...
def test_tensorize_3d_1():
    print('============ test_tensorize_3d_1 =============')

    V = FunctionSpace('V', ldim=3)
    V_0 = FunctionSpace('V_0', ldim=1, coordinates=['x'])
    V_1 = FunctionSpace('V_1', ldim=1, coordinates=['y'])
    V_2 = FunctionSpace('V_2', ldim=1, coordinates=['z'])

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    v0 = TestFunction(V_0, name='v0')
    u0 = TestFunction(V_0, name='u0')

    v1 = TestFunction(V_1, name='v1')
    u1 = TestFunction(V_1, name='u1')

    v2 = TestFunction(V_2, name='v2')
    u2 = TestFunction(V_2, name='u2')

    c = Constant('c')

    bx = Constant('bx')
    by = Constant('by')
    bz = Constant('bz')
    b = Tuple(bx, by, bz)

    # ...
    expected =  TensorProduct(Mass(v2,u2),Mass(v1,u1),Mass(v0,u0))
    assert(tensorize(BilinearForm((v,u), u*v)) == expected)
    # ...

    # ...
    expected =  TensorProduct(Mass(v2,u2),Mass(v1,u1),Stiffness(v0,u0))
    assert(tensorize(BilinearForm((v,u), dx(u)*dx(v))) == expected)
    # ...

    # ...
    expected = TensorProduct(Mass(v2,u2),Advection(v1,u1),Mass(v0,u0))
    assert(tensorize(BilinearForm((v,u), dy(u) * v)) == expected)
    # ...

    # ...
    expected =  TensorProduct(Mass(v2,u2),Mass(v1,u1),Advection(v0,u0))
    assert(tensorize(BilinearForm((v,u), dx(u) * v)) == expected)
    # ...

    # ...
    expected = (TensorProduct(Mass(v2,u2),Mass(v1,u1),Stiffness(v0,u0)) +
                TensorProduct(Mass(v2,u2),Stiffness(v1,u1),Mass(v0,u0)) +
                TensorProduct(Stiffness(v2,u2),Mass(v1,u1),Mass(v0,u0)))
    assert(tensorize(BilinearForm((v,u), dot(grad(v), grad(u)))) == expected)
    # ...

    # ...
    expected = (TensorProduct(Mass(v2,u2),Advection(v1,u1),Mass(v0,u0)) +
                TensorProduct(Mass(v2,u2),Mass(v1,u1),Advection(v0,u0)) +
                TensorProduct(Mass(v2,u2),Mass(v1,u1),Stiffness(v0,u0)) +
                TensorProduct(Mass(v2,u2),Stiffness(v1,u1),Mass(v0,u0)) +
                TensorProduct(Stiffness(v2,u2),Mass(v1,u1),Mass(v0,u0)))
    assert(tensorize(BilinearForm((v,u), dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v)) == expected)
    # ...

    # ...
    expected = (TensorProduct(bx,Mass(v2,u2),Mass(v1,u1),AdvectionT(v0,u0)) +
                TensorProduct(by,Mass(v2,u2),AdvectionT(v1,u1),Mass(v0,u0)) +
                TensorProduct(bz,AdvectionT(v2,u2),Mass(v1,u1),Mass(v0,u0)))

    assert(tensorize(BilinearForm((v,u), dot(b, grad(v)) * u)) == expected)
    # ...

    # ...
    expected = (TensorProduct(bx**2,Mass(v2,u2),Mass(v1,u1),Stiffness(v0,u0)) +
                TensorProduct(bx,by,Mass(v2,u2),Advection(v1,u1),AdvectionT(v0,u0)) +
                TensorProduct(bx,by,Mass(v2,u2),AdvectionT(v1,u1),Advection(v0,u0)) +
                TensorProduct(bx,bz,Advection(v2,u2),Mass(v1,u1),AdvectionT(v0,u0)) +
                TensorProduct(bx,bz,AdvectionT(v2,u2),Mass(v1,u1),Advection(v0,u0)) +
                TensorProduct(by**2,Mass(v2,u2),Stiffness(v1,u1),Mass(v0,u0)) +
                TensorProduct(by,bz,Advection(v2,u2),AdvectionT(v1,u1),Mass(v0,u0)) +
                TensorProduct(by,bz,AdvectionT(v2,u2),Advection(v1,u1),Mass(v0,u0)) +
                TensorProduct(bz**2,Stiffness(v2,u2),Mass(v1,u1),Mass(v0,u0)))

    assert(tensorize(BilinearForm((v,u), dot(b, grad(v)) * dot(b, grad(u)))) == expected)
    # ...

#    expr = dot(b, grad(v)) * dot(b, grad(u))
#    expr = BilinearForm((v,u), expr)
#
#    print('> input         >>> {0}'.format(expr))
#    print('> tensorized    >>> {0}'.format(tensorize(expr)))
# ...

# ...
def test_unknown_3d_1():
    print('============ test_unknown_3d_1 =============')

    v = Unknown('v', ldim=3)
    c = Constant('c')

    # ... expressions that can be normalized (valid for a weak formulation)
    assert(atomize(grad(v)) == Tuple(dx(v),
                                      dy(v),
                                      dz(v)))
    assert(atomize(grad(c*v)) == Tuple(c*dx(v),
                                        c*dy(v),
                                        c*dz(v)))
    # ...
# ...

# .....................................................
if __name__ == '__main__':
    test_atomize_3d_1()
    test_normalize_3d_1()
    test_evaluate_3d_1()

    test_atomize_3d_2()
    test_normalize_3d_2()
    test_matricize_3d_2()

#    test_bilinear_form_3d_10() # TODO not working, since args are the same
    test_linear_form_3d_10()
    test_function_form_3d_10()

    test_matricize_3d_3()
    test_evaluate_3d_3()

    test_tensorize_3d_1()
    test_calls_3d_3()

    test_unknown_3d_1()
