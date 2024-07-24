# coding: utf-8

from sympy.abc import x,y,z
from sympy import Tuple
from sympy import symbols

x1, x2, x3 = symbols('x1, x2, x3')

from sympde.topology import Interval, ProductDomain, InteriorDomain, Domain
from sympde.topology import Line, Square, Cube, NCubeInterior

#==============================================================================
def test_interval():
    Ix = Interval('Ix', coordinate=x)
    Iy = Interval('Iy', coordinate=y)
    Iz = Interval('Iz', coordinate=z)

    for I,i in zip([Ix, Iy, Iz], [x, y, z]):
        assert(I.coordinates == i)

    D_xy = ProductDomain(Ix, Iy)
    assert(D_xy.dim == 2)

#==============================================================================
def test_unit_line():

    # Create 1D domain (Line) from interval [0, 1]
    domain = Line('line')

    assert isinstance(domain, Line)

    # BasicDomain's attributes
    assert domain.dim  == 1
    assert domain.name == 'line'
    assert domain.coordinates == x1

    # Domain's attributes
    assert isinstance(domain.interior, NCubeInterior)
    assert len(domain.boundary) == 2
    assert domain.dtype == {'type': 'Line',
                            'parameters': {'bounds': [0, 1]}}

    # NCube's attributes
    assert domain.min_coords == (0,)
    assert domain.max_coords == (1,)

    # Line's attributes
    assert domain.bounds == (0, 1)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain

#==============================================================================
def test_generic_line():

    # Create 1D domain (Line) from interval [-3, 4]
    domain = Line('line', bounds=(-3, 4))

    assert isinstance(domain, Line)

    # BasicDomain's attributes
    assert domain.dim  == 1
    assert domain.name == 'line'
    assert domain.coordinates == x1

    # Domain's attributes
    assert isinstance(domain.interior, NCubeInterior)
    assert len(domain.boundary) == 2
    assert domain.dtype == {'type': 'Line',
                            'parameters': {'bounds': [-3, 4]}}

    # NCube's attributes
    assert domain.min_coords == (-3,)
    assert domain.max_coords == ( 4,)

    # Line's attributes
    assert domain.bounds == (-3, 4)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain

#==============================================================================
def test_unit_square():

    # Create 2D square domain [0, 1]^2
    domain = Square('square')

    assert isinstance(domain, Square)

    # BasicDomain's attributes
    assert domain.dim  == 2
    assert domain.name == 'square'
    assert domain.coordinates == (x1, x2)

    # Domain's attributes
    assert isinstance(domain.interior, InteriorDomain)

    assert len(domain.boundary) == 4
    assert domain.dtype == {'type': 'Square',
                            'parameters': {'bounds1': [0, 1],
                                           'bounds2': [0, 1]}}

    # NCube's attributes
    assert domain.min_coords == (0, 0)
    assert domain.max_coords == (1, 1)

    # Square's attributes
    assert domain.bounds1 == (0, 1)
    assert domain.bounds2 == (0, 1)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain

#==============================================================================
def test_rectangle():

    # Create 2D rectangular domain [1, 5] X [3, 7]
    domain = Square('rectangle', bounds1=(1, 5), bounds2=(3, 7))

    assert isinstance(domain, Square)

    # BasicDomain's attributes
    assert domain.dim  == 2
    assert domain.name == 'rectangle'
    assert domain.coordinates == (x1, x2)

    # Domain's attributes
    assert isinstance(domain.interior, InteriorDomain)
    assert len(domain.boundary) == 4
    assert domain.dtype == {'type': 'Square',
                            'parameters': {'bounds1': [1, 5],
                                           'bounds2': [3, 7]}}

    # NCube's attributes
    assert domain.min_coords == (1, 3)
    assert domain.max_coords == (5, 7)

    # Square's attributes
    assert domain.bounds1 == (1, 5)
    assert domain.bounds2 == (3, 7)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain

#==============================================================================
def test_unit_cube():

    # Create 3D cube domain [0, 1]^3
    domain = Cube('cube')

    assert isinstance(domain, Cube)

    # Check object attributes
    assert domain.dim  == 3
    assert domain.name == 'cube'
    assert domain.coordinates == (x1, x2, x3)

    # Domain's attributes
    assert isinstance(domain.interior, InteriorDomain)
    assert len(domain.boundary) == 6
    assert domain.dtype == {'type': 'Cube',
                            'parameters': {'bounds1': [0, 1],
                                           'bounds2': [0, 1],
                                           'bounds3': [0, 1]}}

    # NCube's attributes
    assert domain.min_coords == (0, 0, 0)
    assert domain.max_coords == (1, 1, 1)

    # Cube's attributes
    assert domain.bounds1 == (0, 1)
    assert domain.bounds2 == (0, 1)
    assert domain.bounds3 == (0, 1)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain

#==============================================================================
def test_orthogonal_hexahedron():

    # Create 3D orthogonal hexahedron [-1, 1] X [0, 10] X [0, 2]
    domain = Cube('hexahedron', bounds1=(-1, 1), bounds2=(0, 10), bounds3=(0, 2))

    assert isinstance(domain, Cube)

    # Check object attributes
    assert domain.dim  == 3
    assert domain.name == 'hexahedron'
    assert domain.coordinates == (x1, x2, x3)

    # Domain's attributes
    assert isinstance(domain.interior, InteriorDomain)
    assert len(domain.boundary) == 6
    assert domain.dtype == {'type': 'Cube',
                            'parameters': {'bounds1': [-1, 1],
                                           'bounds2': [0, 10],
                                           'bounds3': [0,  2]}}

    # NCube's attributes
    assert domain.min_coords == (-1, 0, 0)
    assert domain.max_coords == (1, 10, 2)

    # Cube's attributes
    assert domain.bounds1 == (-1, 1)
    assert domain.bounds2 == (0, 10)
    assert domain.bounds3 == (0,  2)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

    # Remove output file generated by tests
    import os
    fname = 'domain.h5'
    if os.path.exists(fname):
        os.remove(fname)

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
