# coding: utf-8

from collections import OrderedDict

from sympy import symbols
from sympy import Tuple
from sympy import Matrix
from sympy import srepr

from sympde.core import (dx, dy, dz)
from sympde.core import grad, dot, inner
from sympde.core import Field, Constant
from sympde.core import get_index_derivatives_atom
from sympde.core import partial_derivative_as_str
from sympde.core import get_max_partial_derivatives
from sympde.core import FunctionSpace
from sympde.core import Domain


def indices_as_str(a):
    code = ''
    for k,n in list(a.items()):
        code += k*n
    return code



# ...
def test_partial_derivatives_1():
    print('============ test_partial_derivatives_1 ==============')

    # ...
    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    V = FunctionSpace('V', domain)

    F = Field('F', space=V)
    u = Field('u', space=V)
    v = Field('v', space=V)
    w = Field('w', space=V)
    uvw = Tuple(u,v,w)

    alpha = Constant('alpha')
    beta = Constant('beta')
    # ...

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

    # ...
    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    V = FunctionSpace('V', domain)
    F = Field('F', space=V)

    alpha = Constant('alpha')
    beta = Constant('beta')
    # ...

    # ...
    expr = alpha * dx(F)

    indices = get_index_derivatives_atom(expr, F)[0]
    assert(indices_as_str(indices) == 'x')
    # ...

    # ...
    expr = dy(dx(F))

    indices = get_index_derivatives_atom(expr, F)[0]
    assert(indices_as_str(indices) == 'xy')
    # ...

    # ...
    expr = alpha * dx(dy(dx(F)))

    indices = get_index_derivatives_atom(expr, F)[0]
    assert(indices_as_str(indices) == 'xxy')
    # ...

    # ...
    expr = alpha * dx(dx(F)) + beta * dy(F) + dx(dy(F))

    indices = get_index_derivatives_atom(expr, F)
    indices = [indices_as_str(i) for i in indices]
    assert(sorted(indices) == ['xx', 'xy', 'y'])
    # ...

    # ...
    expr = alpha * dx(dx(F)) + beta * dy(F) + dx(dy(F))

    d = get_max_partial_derivatives(expr, F)
    assert(indices_as_str(d) == 'xxy')

    d = get_max_partial_derivatives(expr)
    assert(indices_as_str(d) == 'xxy')
    # ...

#    print('> ', srepr(expr))
#    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_partial_derivatives_1()
    test_partial_derivatives_2()
