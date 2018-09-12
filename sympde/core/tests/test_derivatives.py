# coding: utf-8

# TODO: - add more tests and asserts

from sympy import symbols
from sympy import Tuple
from sympy import Matrix
from sympy import srepr

from sympde.core import (dx, dy, dz)
from sympde.core import grad, dot, inner
from sympde.core import Field
from sympde.core import get_index_derivatives_atom
from sympde.core import FunctionSpace
from sympde.core import Domain


# ...
def test_partial_derivatives_1():
    print('============ test_partial_derivatives_1 ==============')

    x, y, z = symbols('x y z')
    alpha, beta = symbols('alpha beta')

    domain = Domain('Omega', dim=2)

    V = FunctionSpace('V', domain)
    F = Field('F', space=V)

    u = Field('u', space=V)
    v = Field('v', space=V)
    w = Field('w', space=V)
    uvw = Tuple(u,v,w)

    # ...
    assert(dx(x**2) == 2*x)
    assert(dy(x**2) == 0)
    assert(dz(x**2) == 0)

    assert(dx(x*F) == F + x*dx(F))
    assert(dx(uvw) == Matrix([[dx(u), dx(v), dx(w)]]))
    assert(dx(uvw) + dy(uvw) == Matrix([[dx(u) + dy(u),
                                         dx(v) + dy(v),
                                         dx(w) + dy(w)]]))

    expected = Matrix([[alpha*dx(u) + beta*dy(u),
                        alpha*dx(v) + beta*dy(v),
                        alpha*dx(w) + beta*dy(w)]])
    assert(alpha * dx(uvw) + beta * dy(uvw) == expected)
    # ...

#    expr = alpha * dx(uvw) + beta * dy(uvw)
#    print(expr)

#    print('> ', srepr(expr))
#    print('')
# ...

# ...
def test_partial_derivatives_2():
    print('============ test_partial_derivatives_2 ==============')

    x, y, z = symbols('x y z')
    alpha, beta = symbols('alpha beta')

    domain = Domain('Omega', dim=2)

    V = FunctionSpace('V', domain)
    F = Field('F', space=V)

    expr = alpha * dx(F) + beta * dy(F) + dx(dy(F))
#    expr = alpha * dx(F)
    print('> expr = ', expr)
    indices = get_index_derivatives_atom(expr, F)
    print('> indices = ', indices)

#    print('> ', srepr(expr))
#    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_partial_derivatives_1()
    test_partial_derivatives_2()
