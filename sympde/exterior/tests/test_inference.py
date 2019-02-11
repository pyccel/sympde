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
def test_type_inference_1():

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
    expr = wedge(u_2, u_2)
    assert(isinstance(infere_type(expr), FourFormType))
    # ...

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
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()

test_type_inference_1()
