# coding: utf-8

# TODO split the asserts between algebraic and weak formulations ones
# TODO: - __call__ examples are not working anymore

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Function
from sympy import Matrix
from sympy import pi, cos, sin
from sympy import srepr

from sympde.core import dx, dy, dz
from sympde.core import Constant
from sympde.core import Field
from sympde.core import grad, dot, cross, rot, curl, div
from sympde.core import FunctionSpace
from sympde.core import ProductSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import BilinearForm, LinearForm, FunctionForm
from sympde.core import atomize, normalize
from sympde.core import evaluate
from sympde.core import tensorize
from sympde.core import Mass, Stiffness, Advection, AdvectionT
from sympde.core import Unknown
from sympy.physics.quantum import TensorProduct


# ...
def test_atomize_1d_1():
    print('============ test_atomize_1d_1 =============')

    V = FunctionSpace('V', ldim=1)

    v = TestFunction(V, name='v')
    w = TestFunction(V, name='w')
    c = Constant('c')
    F = Field('F', space=V)
    x = Symbol('x')
    f = Function('f')

    # ...
    assert(atomize(grad(v)) == dx(v))
    assert(atomize(grad(c*v)) == c*dx(v))
    assert(atomize(grad(F*v)) == F*dx(v) + v*dx(F))
    assert(atomize(f(x)*grad(v)) == dx(v)*f(x))

    assert(atomize(dot(grad(v), grad(w))) == dx(v)*dx(w))
    # ...

    # ...
    assert(atomize(grad(v*w)) == w*dx(v) + v*dx(w))
    assert(atomize(div(grad(v*w))) == 2*dx(v)*dx(w) + dx(dx(v))*w + dx(dx(w))*v)
    # ...

#    expr = div(grad(v*w))
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
#    print(expr.is_commutative)
# ...

# ...
def test_normalize_1d_1():
    print('============ test_normalize_1d_1 =============')

    V = FunctionSpace('V', ldim=1)
    U = FunctionSpace('U', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')
    c = Constant('c')
    F = Field('F', space=V)

    Ni, Ni_x = symbols('Ni Ni_x')
    Nj, Nj_x = symbols('Nj Nj_x')

    # ...
    assert(normalize(grad(v), basis={v: 'Ni'}) == Ni_x)
    assert(normalize(grad(c*v), basis={v: 'Ni'}) == c*Ni_x)
    assert(normalize(grad(v) + c*v, basis={v: 'Ni'}) == Ni_x + c*Ni)
    assert(normalize(grad(F*v), basis={v: 'Ni'}) == F*Ni_x + Ni*dx(F))

    assert(normalize(dot(grad(v), grad(u)), basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x)
    assert(normalize(dot(grad(v), grad(u)) + c*v*u, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + c*Ni*Nj)
    assert(normalize(dot(grad(F*v), grad(u)), basis={v: 'Ni', u: 'Nj'}) == Nj_x*(F*Ni_x + Ni*dx(F)))
    # ...

#    expr = dot(grad(F*v), grad(u))
#    print('> input         >>> {0}'.format(expr))

#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni', U: 'Nj'})))
# ...

# ...
def test_evaluate_1d_1():
    print('============ test_evaluate_1d_1 =============')

    V = FunctionSpace('V', ldim=1)
    U = FunctionSpace('U', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')
    c = Constant('c')
    F = Field('F', space=V)

    Ni, Ni_x = symbols('Ni Ni_x')
    Nj, Nj_x = symbols('Nj Nj_x')

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)) + c*v*u)
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + c*Ni*Nj)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)) + F*v*u)
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == Ni_x*Nj_x + F*Ni*Nj)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(F*v), grad(u)))
    assert(evaluate(a, basis={v: 'Ni', u: 'Nj'}) == F*Ni_x*Nj_x + Ni*Nj_x*dx(F))
    # ...
# ...

# ...
def test_calls_1d_3():
    print('============ test_calls_1d_3 =============')

    V1 = FunctionSpace('V1', ldim=1)
    V2 = FunctionSpace('V2', ldim=1)
    U1 = FunctionSpace('U1', ldim=1)
    U2 = FunctionSpace('U2', ldim=1)

    v1 = TestFunction(V1, name='v1')
    v2 = TestFunction(V2, name='v2')
    u1 = TestFunction(U1, name='u1')
    u2 = TestFunction(U2, name='u2')

    V = ProductSpace(V1, V2)
    U = ProductSpace(U1, U2)

    x = V1.coordinates

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
    l1 = LinearForm(v1, x*(1.-x)*v1)
#    print(l1(v2))
    # ...
# ...

# ...
def test_evaluate_1d_3():
    print('============ test_evaluate_1d_3 =============')

    V1 = FunctionSpace('V1', ldim=1)
    U1 = FunctionSpace('U1', ldim=1)
    V2 = FunctionSpace('V2', ldim=1)
    U2 = FunctionSpace('U2', ldim=1)

    v1 = TestFunction(V1, name='v1')
    u1 = TestFunction(U1, name='u1')
    v2 = TestFunction(V2, name='v2')
    u2 = TestFunction(U2, name='u2')

    V = ProductSpace(V1, V2)
    U = ProductSpace(U1, U2)

    v1v2 = VectorTestFunction(V, name='v1v2')
    u1u2 = VectorTestFunction(U, name='u1u2')

    c = Constant('c')

    Ni, Ni_x = symbols('Ni Ni_x')
    Nj, Nj_x = symbols('Nj Nj_x')

    basis = {v1v2: 'Ni', u1u2: 'Nj'}

    # ...
    expr = v1*u1 + dx(v2)*dx(u2)
    a = BilinearForm(((v1, v2), (u1, u2)), expr)

    expected = Matrix([[Ni*Nj, 0], [0, Ni_x*Nj_x]])
    assert(evaluate(a, basis=basis) == expected)
    # ...

    # ...
    expr = v1*u1 + dx(v2)*dx(u2) + dx(v2)*u1 + v1*dx(u2)
    a = BilinearForm(((v1, v2), (u1, u2)), expr)

    expected = Matrix([[Ni*Nj, Ni_x*Nj], [Ni*Nj_x, Ni_x*Nj_x]])
    assert(evaluate(a, basis=basis) == expected)
    # ...

#    expr = v1*u1 + dx(v2)*dx(u2) + dx(v2)*u1 + v1*dx(u2)
#    expr = BilinearForm(((v1, v2), (u1, u2)), expr)
#    print('> input         >>> {0}'.format(expr))
#    print('> normal form   >>> {0}'.format(evaluate(expr, basis=basis)))

# ...

# ...
def test_bilinear_form_1d_10():
    print('============ test_bilinear_form_1d_10 =============')

    U = FunctionSpace('U', ldim=1)
    V = FunctionSpace('V', ldim=1)

    u = TestFunction(U, name='u')
    v = TestFunction(V, name='v')

    u1 = TestFunction(U, name='u1')
    v1 = TestFunction(V, name='v1')

    Ni, Ni_x = symbols('Ni Ni_x')
    Nj, Nj_x = symbols('Nj Nj_x')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    a = BilinearForm((v,u), dot(grad(u), grad(v)))
    b = BilinearForm((v,u), u*v)
    adv = BilinearForm((v,u), dx(u)*v)

    # ...
    expected = Ni*Nj + Ni_x*Nj_x
    assert(evaluate(a + b, basis={v: 'Nj', u: 'Ni'}) == expected)
    # ...

    # ...
    expected = 2*Ni_x*Nj_x
    assert(evaluate(2 * a, basis={v: 'Nj', u: 'Ni'}) == expected)
    # ...

    # ...
    expected = c1*Ni_x*Nj_x
    assert(evaluate(c1*a, basis={v: 'Nj', u: 'Ni'}) == expected)
    # ...

    # ...
    expected = c2*Ni*Nj + c1*Ni_x*Nj_x
    assert(evaluate(c1*a + c2*b, basis={v: 'Nj', u: 'Ni'}) == expected)
    # ...

    # ...
    expected = Ni_x*Nj_x*c1 + c2*(Ni*Nj + Ni_x*Nj)
    assert(evaluate(c1*a  + c2*(b + adv), basis={v: 'Nj', u: 'Ni'}) == expected)
    # ...

#    expr = adv(v1, u1)
#    print('> input      >>> {0}'.format(expr))
#    print('> evaluated  >>> {0}'.format(evaluate(expr, basis={V: 'Nj', U: 'Ni'}) ))
#    print('')
# ...

# ...
def test_linear_form_1d_10():
    print('============ test_linear_form_1d_10 =============')

    V = FunctionSpace('V', ldim=1)

    v = TestFunction(V, name='v')

    x = V.coordinates
    f = Function('f')

    Ni, Ni_x, Ni_xx = symbols('Ni Ni_x Ni_xx')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    # ...
    expected = cos(2*pi*x)*Ni
    assert(evaluate(LinearForm(v, cos(2*pi*x)*v),
                    basis={v: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x)*Ni
    assert(evaluate(LinearForm(v, f(x)*v),
                    basis={v: 'Ni'}) == expected)
    # ...

    # ...
    expected = cos(2*pi*x)*Ni_x
    assert(evaluate(LinearForm(v, cos(2*pi*x)*dx(v)),
                    basis={v: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x)*Ni_xx
    assert(evaluate(LinearForm(v, f(x)*dx(dx(v))),
                    basis={v: 'Ni'}) == expected)
    # ...

#    expr = LinearForm(v, cos(2*pi*x)*dx(v))
#    print('> input      >>> {0}'.format(expr))
#    print('> evaluated  >>> {0}'.format(evaluate(expr, basis={V: 'Ni'}) ))
#    print('')
# ...

# ...
def test_function_form_1d_10():
    print('============ test_function_form_1d_10 =============')

    V = FunctionSpace('V', ldim=1)

    F = Field('F', space=V)

    x = V.coordinates
    f = Function('f')

    Ni, Ni_x, Ni_xx = symbols('Ni Ni_x Ni_xx')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    # ...
    expected = -2*pi*sin(2*pi*x)
    assert(evaluate(FunctionForm(grad(cos(2*pi*x)), coordinates=[x])) == expected)
    # ...

    # ...
    expected = F-cos(2*pi*x)
    assert(evaluate(FunctionForm(F-cos(2*pi*x))) == expected)
    # ...

    # ...
    expected = (F-cos(2*pi*x))**2
    assert(evaluate(FunctionForm((F-cos(2*pi*x))**2)) == expected)
    # ...

    # ...
    expected = dx(F) + 2*pi*sin(2*pi*x)
    assert(evaluate(FunctionForm(grad(F-cos(2*pi*x)))) == expected)
    # ...

    # ...
    expected = (dx(F) + 2*pi*sin(2*pi*x))**2
    assert(evaluate(FunctionForm((grad(F-cos(2*pi*x)))**2)) == expected)
    # ...

#    expr = FunctionForm()
#    print('> input      >>> {0}'.format(expr))
#    print('> evaluated  >>> {0}'.format(evaluate(expr) ))
#    print('')
# ...

# ...
def test_tensorize_1d_1():
    print('============ test_tensorize_1d_1 =============')

    V = FunctionSpace('V', ldim=1)
    V_0 = FunctionSpace('V_0', ldim=1, coordinates=['x'])

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    v0 = TestFunction(V_0, name='v0')
    u0 = TestFunction(V_0, name='u0')

    c = Constant('c')

    # ...
    expected = Mass(v0, u0)
    assert(tensorize(BilinearForm((v,u), u*v)) == expected)
    # ...

    # ...
    expected = Stiffness(v0, u0)
    assert(tensorize(BilinearForm((v,u), dx(u)*dx(v))) == expected)
    # ...

    # ...
    expected = Advection(v0,u0)
    assert(tensorize(BilinearForm((v,u), dx(u) * v)) == expected)
    # ...

    # ...
    expected = Advection(v0,u0) + AdvectionT(v0,u0) + TensorProduct(c ,Stiffness(v0,u0))
    assert(tensorize(BilinearForm((v,u), dx(v) * u + v * dx(u) + c * dx(v)*dx(u))) == expected)
    # ...

#    expr = dx(v) * u + v * dx(u) + c * dx(v)*dx(u)
#    expr = BilinearForm((v,u), expr)
#
#    print('> input         >>> {0}'.format(expr))
#    print('> tensorized    >>> {0}'.format(tensorize(expr)))
# ...

# ...
def test_unknown_1d_1():
    print('============ test_unknown_1d_1 =============')

    v = Unknown('v', ldim=1)
    c = Constant('c')

    # ...
    assert(atomize(grad(v)) == dx(v))
    assert(atomize(grad(c*v)) == c*dx(v))
    # ...
# ...

# .....................................................
if __name__ == '__main__':
    test_atomize_1d_1()
    test_normalize_1d_1()
    test_evaluate_1d_1()

    test_evaluate_1d_3()

    test_bilinear_form_1d_10()
    test_linear_form_1d_10()
    test_function_form_1d_10()

    test_tensorize_1d_1()
    test_calls_1d_3()

    test_unknown_1d_1()