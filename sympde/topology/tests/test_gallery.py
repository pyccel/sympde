# coding: utf-8

from sympy.abc import x,y,z
from sympy import Tuple

from sympde.topology import Interval, ProductDomain
from sympde.topology import Domain, Line, Square, Cube

#==============================================================================
def test_interval_1():
    Ix = Interval('Ix', coordinate=x)
    Iy = Interval('Iy', coordinate=y)
    Iz = Interval('Iz', coordinate=z)

    for I,i in zip([Ix, Iy, Iz], [x, y, z]):
        assert(I.coordinates == i)

    D_xy = ProductDomain(Ix, Iy)
    assert(D_xy.dim == 2)

#==============================================================================
def test_gallery():
    # ...
    line = Line()
    assert(line.dim == 1)
    assert(len(line.boundary) == 2)
    # ...

    # ...
    cube = Cube()
    assert(cube.dim == 3)
    assert(len(cube.boundary) == 6)
    # ...

#==============================================================================
def test_square():
    square = Square('square')

    assert(square.dim == 2)
    assert(len(square.boundary) == 4)

    square.export('square.h5')

    # read it again
    D = Domain.from_file('square.h5')

    # TODO BUG not working after naming boundaries \Gamma_1 etc
#    assert( D == square )



#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
