#==============================================================================
def test_constant_1():
    from sympde.core import Constant
    from sympde.core import ScalarConstant


    i = Constant('i', value=5)        # integer constant with value 5
    r = Constant('r', dtype=float)    # undefined real constant
    c = Constant('c', dtype=complex)  # undefined complex constant

#    w = Constant('w', shape=2, dtype=float, value=[3, 7])  # 2-vector of floats with value [3.0, 7.0]
#    n = Constant('n', shape=3, dtype=int)                  # undefined 3-vector of integers
#
#    M = Constant('M', shape=(2, 2), value=[[1, 2], [3, 4]]) # 2x2 matrix of integers with value [[1, 2], [3, 4]]
#    A = Constant('A', shape=(3, 3), dtype=complex)          # undefined 3x3 matrix of complex numbers

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()

test_constant_1()
