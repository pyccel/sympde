# coding: utf-8

import numpy as np

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda

from sympde.core import (dx, dy, dz)
from sympde.core import LinearOperator
from sympde.core import Field
from sympde.core import grad, dot, inner
from sympde.core import get_index_derivatives


# ...
def test_0():
    print('============ test_0 ==============')

    x,y, a = symbols('x y a')

    # ...
    expr = x+y
    print('> expr := {0}'.format(expr))

    expr = LinearOperator(expr)
    print('> gelatized := {0}'.format(expr))
    print('')
    # ...

    # ...
    expr = 2*x+y
    print('> expr := {0}'.format(expr))

    expr = LinearOperator(expr)
    print('> gelatized := {0}'.format(expr))
    print('')
    # ...

    # ...
    expr = a*x+y
    print('> expr := {0}'.format(expr))

    expr = LinearOperator(expr)
    print('> gelatized := {0}'.format(expr))
    # ...

    # ...
    expr = 2*a*x+y
    print('> expr := {0}'.format(expr))

    expr = LinearOperator(expr)
    print('> gelatized := {0}'.format(expr))
    # ...
# ...

# ...
def test_1():
    print('============ test_1 ==============')

    u, v, a = symbols('u v a')

    # ...
    expr = u+v
    print('> expr := {0}'.format(expr))

    expr = dx(expr)
    print('> gelatized := {0}'.format(expr))
    print('')
    # ...

    # ...
    expr = 2*u*v
    print('> expr := {0}'.format(expr))

    expr = dx(expr)
    print('> gelatized := {0}'.format(expr))
    print('')
    # ...

    # ... dx should not operate on u^2,
    #     since we consider only linearized weak formulations
    expr = u*u
    print('> expr := {0}'.format(expr))

    expr = dx(expr)
    print('> gelatized := {0}'.format(expr))
    print('')
    # ...
# ...

# ...
def test_2():
    print('============ test_2 ==============')

    u, v = symbols('u v')
    F = Field('F')

    # ...
    expr = F*v*u
    print('> expr := {0}'.format(expr))

    expr = dx(expr)
    print('> gelatized := {0}'.format(expr))
    print('')
    # ...
# ...

# ...
def test_3():
    print('============ test_3 ==============')

    u = symbols('u')

    # ...
    expr = dx(u)
    d = get_index_derivatives(expr)
    assert(d['x'] == 1)
    assert(d['y'] == 0)
    assert(d['z'] == 0)
    # ...

    # ...
    expr = dx(dy(u))
    d = get_index_derivatives(expr)
    assert(d['x'] == 1)
    assert(d['y'] == 1)
    assert(d['z'] == 0)
    # ...

    # ...
    expr = dx(dy(dx(u)))
    d = get_index_derivatives(expr)
    assert(d['x'] == 2)
    assert(d['y'] == 1)
    assert(d['z'] == 0)
    # ...

    print('')
# ...

# ...
def test_poisson():
    print('============ test_poisson ==============')

    u, v = symbols('u v')

    # ...
    expr = inner(grad(v), grad(u))
    print('> expr := {0}'.format(expr))
    # ...
# ...

# .....................................................
if __name__ == '__main__':
    test_0()
    test_1()
    test_2()
    test_3()
    test_poisson()
