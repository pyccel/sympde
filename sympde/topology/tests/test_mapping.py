# coding: utf-8

import numpy as np

from sympy.core.containers import Tuple
from sympy import Matrix
from sympy.tensor import IndexedBase
from sympy import symbols, simplify

from sympde.topology import Mapping, MappedDomain, AffineMapping
from sympde.topology import dx, dy, dz
from sympde.topology import dx1, dx2, dx3
from sympde.topology import Domain

from sympde.topology.mapping import Jacobian, Covariant, Contravariant
# ...
def test_mapping_1d():
    print('============ test_mapping_1d ==============')

    dim = 1

    F = Mapping('F', dim=dim)

    assert(F.name == 'F')

    # ...
    expected = Matrix([[dx1(F[0])]])
    assert(Jacobian(F) == expected)
    # ...

    # ...
    expected = dx1(F[0])
    assert(Jacobian(F).det() == expected)
    # ...
# ...

# ...
def test_mapping_2d():
    print('============ test_mapping_2d ==============')

    dim = 2

    F = Mapping('F', dim=dim)

    a,b = symbols('a b')
    ab = Tuple(a, b)

    assert(F.name == 'F')

    # ...
    expected = Matrix([[dx1(F[0]), dx2(F[0])],
                       [dx1(F[1]), dx2(F[1])]])
    assert(Jacobian(F) == expected)
    # ...

    # ...
    expected = dx1(F[0])*dx2(F[1]) - dx1(F[1])*dx2(F[0])
    assert(Jacobian(F).det() == expected)
    # ...

    # ...
    expected = Tuple(a*dx2(F[1])/(dx1(F[0])*dx2(F[1]) - dx1(F[1])*dx2(F[0]))
                     - b*dx1(F[1])/(dx1(F[0])*dx2(F[1]) - dx1(F[1])*dx2(F[0])),
                     - a*dx2(F[0])/(dx1(F[0])*dx2(F[1]) - dx1(F[1])*dx2(F[0]))
                     + b*dx1(F[0])/(dx1(F[0])*dx2(F[1]) - dx1(F[1])*dx2(F[0])))
    assert(Covariant(F, ab) == expected)
    # ...

    # ...
    expected = Tuple(a*dx1(F[0])/(dx1(F[0])*dx2(F[1]) - dx1(F[1])*dx2(F[0]))
                     + b*dx2(F[0])/(dx1(F[0])*dx2(F[1]) - dx1(F[1])*dx2(F[0])),
                     a*dx1(F[1])/(dx1(F[0])*dx2(F[1]) - dx1(F[1])*dx2(F[0]))
                     + b*dx2(F[1])/(dx1(F[0])*dx2(F[1]) - dx1(F[1])*dx2(F[0])))
    assert(simplify(Contravariant(F, ab)) == simplify(expected))
    # ...

# ...
def test_mapping_3d():
    print('============ test_mapping_3d ==============')

    dim = 3

    F = Mapping('F', dim=dim)

    a,b,c = symbols('a b c')
    abc = Tuple(a, b, c)

    assert(F.name == 'F')

    # ...
    expected = Matrix([[dx1(F[0]), dx2(F[0]), dx3(F[0])],
                       [dx1(F[1]), dx2(F[1]), dx3(F[1])],
                       [dx1(F[2]), dx2(F[2]), dx3(F[2])]])
    assert(Jacobian(F) == expected)
    # ...

    # ...
    expected = (dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) -
                dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) +
                dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0]))
    assert(Jacobian(F).det() == expected)
    # ...

    # ...
    expected = Tuple (a*(dx2(F[1])*dx3(F[2]) - dx2(F[2])*dx3(F[1]))/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + b*(-dx1(F[1])*dx3(F[2]) + dx1(F[2])*dx3(F[1]))/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + c*(dx1(F[1])*dx2(F[2]) - dx1(F[2])*dx2(F[1]))/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])), a*(-dx2(F[0])*dx3(F[2]) + dx2(F[2])*dx3(F[0]))/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + b*(dx1(F[0])*dx3(F[2]) - dx1(F[2])*dx3(F[0]))/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + c*(-dx1(F[0])*dx2(F[2]) + dx1(F[2])*dx2(F[0]))/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])), a*(dx2(F[0])*dx3(F[1]) - dx2(F[1])*dx3(F[0]))/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + b*(-dx1(F[0])*dx3(F[1]) + dx1(F[1])*dx3(F[0]))/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + c*(dx1(F[0])*dx2(F[1]) - dx1(F[1])*dx2(F[0]))/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])))
    cov      = Covariant(F, abc)
    cov      = Matrix(cov)
    expected = Matrix(expected)
    diff     = cov-expected
    diff.simplify()

    assert(diff.dot(diff).is_zero)
    # ...

    # ...
    expected = Tuple (a*dx1(F[0])/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + b*dx2(F[0])/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + c*dx3(F[0])/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])), a*dx1(F[1])/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + b*dx2(F[1])/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + c*dx3(F[1])/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])), a*dx1(F[2])/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + b*dx2(F[2])/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])) + c*dx3(F[2])/(dx1(F[0])*dx2(F[1])*dx3(F[2]) - dx1(F[0])*dx2(F[2])*dx3(F[1]) - dx1(F[1])*dx2(F[0])*dx3(F[2]) + dx1(F[1])*dx2(F[2])*dx3(F[0]) + dx1(F[2])*dx2(F[0])*dx3(F[1]) - dx1(F[2])*dx2(F[1])*dx3(F[0])))
    cov = Contravariant(F, abc)
    assert(simplify(cov) == simplify(expected))
    # ...

# ...

# ...
def test_mapping_2d_2():
    print('============ test_mapping_2d_2 ==============')

    dim   = 2
    F      = Mapping('F', dim=dim)
    domain = Domain('Omega', dim=dim)
    D      = F(domain)

def test_AffineMapping():
    print('============ test_AffineMapping ==============')

    dim   = 2
    alpha = np.pi/2
    c1    = 0.2
    c2    = 1.5

    F =  AffineMapping(
        name='F', dim=2, c1=c1, c2=c2,
        a11=np.cos(alpha), a12=-np.sin(alpha),
        a21=np.sin(alpha), a22=np.cos(alpha),
    )

    domain = Domain('Omega', dim=dim)
    D      = F(domain)

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()

