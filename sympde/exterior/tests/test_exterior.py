# coding: utf-8


from sympy import symbols
from sympy import Tuple
from sympy import Matrix
from sympy import srepr
from sympy import Symbol

from sympde import Constant

from sympde.exterior import d, wedge, ip, jp, delta, hodge
from sympde.exterior import DifferentialForm



#==============================================================================
def test_exterior_1():

    x, y, z = symbols('x y z')
    a = Constant('a')
    n = Symbol('n')

    # ...
    u_0 = DifferentialForm('u_0', index=0, dim=n)
    v_0 = DifferentialForm('v_0', index=0, dim=n)

    u_1 = DifferentialForm('u_1', index=1, dim=n)
    v_1 = DifferentialForm('v_1', index=1, dim=n)

    u_2 = DifferentialForm('u_2', index=2, dim=n)
    v_2 = DifferentialForm('v_2', index=2, dim=n)

    u_3 = DifferentialForm('u_3', index=3, dim=n)
    v_3 = DifferentialForm('v_3', index=3, dim=n)

    u_n = DifferentialForm('u_n', index=n, dim=n)
    v_n = DifferentialForm('v_n', index=n, dim=n)
    # ...

    # ... exterior derivative
    assert(d(d(u_0)) == 0)
    assert(d(u_0+v_0) == d(u_0) + d(v_0))
    assert(d(2*u_0) == 2*d(u_0))
    assert(d(a*u_0+v_0) == a*d(u_0) + d(v_0))
    assert(d(u_n) == 0)
    assert(not(d(u_0) == 0))
    # ...

    # ...
    assert(delta(u_0) == 0)
    assert(delta(u_1+v_1) == delta(u_1) + delta(v_1))
    assert(delta(2*u_1) == 2*delta(u_1))
    assert(delta(a*u_1+v_1) == a*delta(u_1) + delta(v_1))
    assert(not(delta(u_n) == 0))
    # ...

    # ... exterior product
    print('> ', wedge(u_0, u_1))
    # ...

    # ... hodge operator
    print('> ', hodge(u_0))
    # ...

    print(hodge(hodge(u_0)))
    print(hodge(hodge(u_1)))
    print(hodge(hodge(u_2)))
    print(hodge(hodge(u_3)))


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()

#test_exterior_1()
