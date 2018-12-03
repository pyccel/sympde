# coding: utf-8

from sympy.core.containers import Tuple
from sympy import Matrix
from sympy.tensor import IndexedBase
from sympy import symbols, simplify

from sympde.core import Mapping
from sympde.core import Jacobian, DetJacobian, Covariant, Contravariant
from sympde.core import dx, dy, dz
from sympde.core import print_expression
from sympde.core import Domain

# ...
def test_mapping_1d():
    print('============ test_mapping_1d ==============')

    rdim = 1

    domain = Domain('Omega', dim=rdim)
    F = Mapping('F', rdim, domain)

    assert(F.name == 'F')

    # ...
    assert(print_expression(F) == 'F')
    assert(print_expression(F[0]) == 'Fx')
    assert(print_expression(dx(F[0])) == 'Fx_x1')
    # ...

    # ...
    expected = Matrix([[dx(F[0])]])
    assert(F.jacobian == expected)
    # ...

    # ...
    expected = dx(F[0])
    assert(F.det_jacobian == expected)
    # ...
# ...

# ...
def test_mapping_2d():
    print('============ test_mapping_2d ==============')

    rdim = 2

    domain = Domain('Omega', dim=rdim)
    F = Mapping('F', rdim, domain)

    a,b = symbols('a b')
    ab = Tuple(a, b)

    assert(F.name == 'F')

    # ...
    expected = Matrix([[dx(F[0]), dy(F[0])],
                       [dx(F[1]), dy(F[1])]])
    assert(F.jacobian == expected)
    # ...

    # ...
    expected = dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])
    assert(F.det_jacobian == expected)
    # ...

    # ...
    expected = Tuple(a*dy(F[1])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))
                     - b*dx(F[1])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])),
                     - a*dy(F[0])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))
                     + b*dx(F[0])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])))
    assert(Covariant(F, ab) == expected)
    # ...

    # ...
    expected = Tuple(a*dx(F[0])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))
                     + b*dy(F[0])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])),
                     a*dx(F[1])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))
                     + b*dy(F[1])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])))
    assert(simplify(Contravariant(F, ab)) == simplify(expected))
    # ...

# ...
def test_mapping_3d():
    print('============ test_mapping_3d ==============')

    rdim = 3

    domain = Domain('Omega', dim=rdim)
    F = Mapping('F', rdim, domain)

    a,b,c = symbols('a b c')
    abc = Tuple(a, b, c)

    assert(F.name == 'F')

    # ...
    expected = Matrix([[dx(F[0]), dy(F[0]), dz(F[0])],
                       [dx(F[1]), dy(F[1]), dz(F[1])],
                       [dx(F[2]), dy(F[2]), dz(F[2])]])
    assert(F.jacobian == expected)
    # ...

    # ...
    expected = (dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) -
                dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) +
                dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0]))
    assert(F.det_jacobian == expected)
    # ...

    # ...
    expected = Tuple (a*(dy(F[1])*dz(F[2]) - dy(F[2])*dz(F[1]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + b*(-dx(F[1])*dz(F[2]) + dx(F[2])*dz(F[1]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + c*(dx(F[1])*dy(F[2]) - dx(F[2])*dy(F[1]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])), a*(-dy(F[0])*dz(F[2]) + dy(F[2])*dz(F[0]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + b*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + c*(-dx(F[0])*dy(F[2]) + dx(F[2])*dy(F[0]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])), a*(dy(F[0])*dz(F[1]) - dy(F[1])*dz(F[0]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + b*(-dx(F[0])*dz(F[1]) + dx(F[1])*dz(F[0]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + c*(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])))
    cov = Covariant(F, abc)
    assert(simplify(cov) == simplify(expected))
    # ...

    # ...
    expected = Tuple (a*dx(F[0])/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + b*dy(F[0])/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + c*dz(F[0])/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])), a*dx(F[1])/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + b*dy(F[1])/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + c*dz(F[1])/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])), a*dx(F[2])/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + b*dy(F[2])/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])) + c*dz(F[2])/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])))
    cov = Contravariant(F, abc)
    assert(simplify(cov) == simplify(expected))
    # ...

# ...

# .....................................................
if __name__ == '__main__':

    test_mapping_1d()
    test_mapping_2d()
    test_mapping_3d()
