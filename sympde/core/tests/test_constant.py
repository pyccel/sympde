from sympde.core import constant
from sympde.core import ScalarConstant
from sympde.core import VectorConstant
from sympde.core import MatrixConstant

from sympy import Tuple

#==============================================================================
def test_constant_1():

    # ...
    i = constant('i', value=5)        # integer constant with value 5
    f = constant('f', dtype=float, value=5)        # float constant with value 5
    r = constant('r', dtype=float)    # undefined real constant
    c = constant('c', dtype=complex)  # undefined complex constant

    assert isinstance(f, ScalarConstant)
    assert i.name == 'i'
    assert i.value == 5
    assert     i.is_integer
    assert not i.is_real
    assert not i.is_complex

    assert isinstance(f, ScalarConstant)
    assert f.name == 'f'
    assert f.value == 5.0
    assert not f.is_integer
    assert     f.is_real
    assert not f.is_complex

    assert isinstance(r, ScalarConstant)
    assert r.name == 'r'
    assert not r.is_integer
    assert     r.is_real
    assert not r.is_complex

    assert isinstance(c, ScalarConstant)
    assert c.name == 'c'
    assert not c.is_integer
    assert not c.is_real
    assert     c.is_complex
    # ...

    # ...
    w = constant('w', shape=2, dtype=float, value=[3, 7])  # 2-vector of floats with value [3.0, 7.0]
    n = constant('n', shape=3, dtype=int)                  # undefined 3-vector of integers

    assert isinstance(w, VectorConstant)
    assert w.name == 'w'
    assert w.value == Tuple(3.0, 7.0)
    assert not w.is_integer
    assert     w.is_real
    assert not w.is_complex

    assert isinstance(n, VectorConstant)
    assert n.name == 'n'
    assert     n.is_integer
    assert not n.is_real
    assert not n.is_complex
    # ...

    # ...
    M = constant('M', shape=(2, 2), value=[[1, 2], [3, 4]]) # 2x2 matrix of integers with value [[1, 2], [3, 4]]
    A = constant('A', shape=(3, 3), dtype=complex)          # undefined 3x3 matrix of complex numbers
    # ...

    assert isinstance(M, MatrixConstant)
    assert M.name == 'M'
    assert M.value == Tuple(Tuple(1, 2), Tuple(3, 4))
    assert     M.is_integer
    assert not M.is_real
    assert not M.is_complex

    assert isinstance(A, MatrixConstant)
    assert A.name == 'A'
    assert not A.is_integer
    assert not A.is_real
    assert     A.is_complex

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================
def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
