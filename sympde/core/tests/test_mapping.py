# coding: utf-8

from sympy.core.containers import Tuple
from sympy import Matrix
from sympy.tensor import IndexedBase
from sympy import symbols, simplify

from sympde.core import Mapping
from sympde.core import Jacobian, DetJacobian, Covariant, Contravariant
from sympde.core import dx, dy, dz
from sympde.core import print_expression

# ...
def test_1d():
    print('============ test_1d ==============')

    rdim = 1

    F = Mapping('F', rdim=rdim)

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
def test_2d():
    print('============ test_2d ==============')

    rdim = 2

    F = Mapping('F', rdim=rdim)

    a,b = symbols('a b')
    ab = Tuple(a, b)

    assert(F.name == 'F')

    # ...
    expected = Matrix([[dx(F[0]), dx(F[1])], [dy(F[0]), dy(F[1])]])
    assert(F.jacobian == expected)
    # ...

    # ...
    expected = dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])
    assert(F.det_jacobian == expected)
    # ...

    # ...
    expected = Tuple(a*dy(F[1])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])) -
                     b*dy(F[0])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])),
                     -a*dx(F[1])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])) +
                     b*dx(F[0])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])))
    assert(Covariant(F, ab) == expected)
    # ...

    # ...
    expected = Tuple((a*dx(F[0]) + b*dx(F[1]))/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])),
                     (a*dy(F[0]) + b*dy(F[1]))/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])))
    assert(simplify(Contravariant(F, ab)) == simplify(expected))
    # ...

# ...
def test_3d():
    print('============ test_3d ==============')

    rdim = 3

    F = Mapping('F', rdim=rdim)

    a,b,c = symbols('a b c')
    abc = Tuple(a, b, c)

    assert(F.name == 'F')

    # ...
    expected = Matrix([[dx(F[0]), dx(F[1]), dx(F[2])],
                       [dy(F[0]), dy(F[1]), dy(F[2])],
                       [dz(F[0]), dz(F[1]), dz(F[2])]])
    assert(F.jacobian == expected)
    # ...

    # ...
    expected = (dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) -
                dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) +
                dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0]))
    assert(F.det_jacobian == expected)
    # ...

    # ...
    expected = Tuple((a*dy(F[1])*dz(F[2]) - a*dy(F[2])*dz(F[1]) - b*dy(F[0])*dz(F[2]) + b*dy(F[2])*dz(F[0]) + c*dy(F[0])*dz(F[1]) - c*dy(F[1])*dz(F[0]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])), (-a*dx(F[1])*dz(F[2]) + a*dx(F[2])*dz(F[1]) + b*dx(F[0])*dz(F[2]) - b*dx(F[2])*dz(F[0]) - c*dx(F[0])*dz(F[1]) + c*dx(F[1])*dz(F[0]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])), (a*dx(F[1])*dy(F[2]) - a*dx(F[2])*dy(F[1]) - b*dx(F[0])*dy(F[2]) + b*dx(F[2])*dy(F[0]) + c*dx(F[0])*dy(F[1]) - c*dx(F[1])*dy(F[0]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])))
    assert(simplify(Covariant(F, abc)) == simplify(expected))
    # ...

    # ...
    expected = Tuple((a*dx(F[0]) + b*dx(F[1]) + c*dx(F[2]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])), (a*dy(F[0]) + b*dy(F[1]) + c*dy(F[2]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])), (a*dz(F[0]) + b*dz(F[1]) + c*dz(F[2]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])))
    assert(simplify(Contravariant(F, abc)) == simplify(expected))
    # ...

# ...

# .....................................................
if __name__ == '__main__':

    test_1d()
    test_2d()
    test_3d()
