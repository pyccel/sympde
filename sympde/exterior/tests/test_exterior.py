# coding: utf-8


from sympy import symbols
from sympy import Tuple
from sympy import Matrix
from sympy import srepr

from sympde import Constant

from sympde.exterior import d, wedge, ip, jp, delta
from sympde.exterior import DifferentialForm
from sympde.exterior import PullBack
from sympde.exterior import infere_index


#==============================================================================
def test_feec_1():

    x, y, z = symbols('x y z')
    a = Constant('a')

    # ...
    u_0 = DifferentialForm('u_0', index=0)
    v_0 = DifferentialForm('v_0', index=0)

    u_1 = DifferentialForm('u_1', index=1)
    v_1 = DifferentialForm('v_1', index=1)

    u_2 = DifferentialForm('u_2', index=2)
    v_2 = DifferentialForm('v_2', index=2)

    u_3 = DifferentialForm('u_3', index=3)
    v_3 = DifferentialForm('v_3', index=3)

#    # TODO not working yet correctly
#    L = PullBack(name='L')

    # ... exterior derivative
    assert(d(d(u_0)) == 0)
    assert(d(u_0+v_0) == d(u_0) + d(v_0))
    assert(d(2*u_0) == 2*d(u_0))
    assert(d(a*u_0+v_0) == a*d(u_0) + d(v_0))
    # ...

    # ... exterior product
    print('> ', wedge(u_0, u_1))
    # ...

#    # ... interior product
#    # TODO must add the notion of Vector
##    print('> ', ip(u_0, u_1))
#    # ...
#
#    # ... pullback
#    expr = wedge(u_0, u_1)
#    print('> ', PullBack(expr))
##    print('> ', L(expr))
#    # ...
#
#    # ... type/index inference
#    print('> ', infere_index(u_0))
#    print('> ', infere_index(d(u_1)))
#    print('> ', infere_index(wedge(u_2, u_2)))
#
#    expr = d(u_1) + u_0
#    print('> ', infere_index(expr))
#
#    # this one will raise an error
#    expr = d(u_1) + u_1
#    print('> ', infere_index(expr))
#    # ...
#
#
##    expr = alpha * dx(uvw) + beta * dy(uvw)
##    print(expr)
#
##    print('> ', srepr(expr))
##    print('')

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
