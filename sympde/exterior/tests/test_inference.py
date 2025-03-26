# coding: utf-8


from sympy import symbols
from sympy import Tuple
from sympy import Matrix
from sympy import srepr
from sympy import Symbol

from sympde import Constant

from sympde.exterior import d, wedge, ip, jp, delta, hodge
from sympde.exterior import DifferentialForm
from sympde.exterior import PullBack
from sympde.exterior import infere_type
from sympde.exterior import ZeroFormType, OneFormType, TwoFormType, ThreeFormType
from sympde.exterior import FourFormType, FiveFormType, SixFormType


#==============================================================================
def test_type_inference_1():

    x, y, z = symbols('x y z')
    a = Constant('a')
#    n = Symbol('n')
    n = 3

    # ...
    u_0 = DifferentialForm('u_0', index=0, dim=n)
    v_0 = DifferentialForm('v_0', index=0, dim=n)

    u_1 = DifferentialForm('u_1', index=1, dim=n)
    v_1 = DifferentialForm('v_1', index=1, dim=n)

    u_2 = DifferentialForm('u_2', index=2, dim=n)
    v_2 = DifferentialForm('v_2', index=2, dim=n)

    u_3 = DifferentialForm('u_3', index=3, dim=n)
    v_3 = DifferentialForm('v_3', index=3, dim=n)

    # ...
    expr = u_0
    assert(isinstance(infere_type(expr), ZeroFormType))
    # ...

    # ...
    expr = d(u_1)
    assert(isinstance(infere_type(expr), TwoFormType))
    # ...

    # ...
    expr = delta(u_1)
    assert(isinstance(infere_type(expr), ZeroFormType))
    # ...

    # ...
    expr = hodge(u_1)
    assert(isinstance(infere_type(expr), TwoFormType))
    # ...

    # ...
    expr = wedge(u_2, u_2)
    assert(isinstance(infere_type(expr), FourFormType))
    # ...

#    print(infere_type(expr))

#    expr = d(u_1) + u_0
#    print('> ', infere_type(expr))
#
#    # this one will raise an error
#    expr = d(u_1) + u_1
#    print('> ', infere_type(expr))
    # ...

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()

#test_type_inference_1()
