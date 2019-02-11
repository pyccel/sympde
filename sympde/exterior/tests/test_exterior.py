# coding: utf-8


from sympy import symbols
from sympy import Tuple
from sympy import Matrix
from sympy import srepr

from sympde import Constant

from sympde.exterior import d, wedge, ip, jp, delta, hodge
from sympde.exterior import DifferentialForm
from sympde.exterior import PullBack
from sympde.exterior import infere_type
from sympde.exterior import ZeroFormType, OneFormType, TwoFormType, ThreeFormType
from sympde.exterior import FourFormType, FiveFormType, SixFormType



#==============================================================================
def test_exterior_1():

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

    # ... exterior derivative
    assert(d(d(u_0)) == 0)
    assert(d(u_0+v_0) == d(u_0) + d(v_0))
    assert(d(2*u_0) == 2*d(u_0))
    assert(d(a*u_0+v_0) == a*d(u_0) + d(v_0))
    # ...

    # ... exterior product
    print('> ', wedge(u_0, u_1))
    # ...

    # ... hodge operator
    print('> ', hodge(u_0))
    # ...


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()

#test_exterior_1()
