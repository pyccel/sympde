# coding: utf-8

from sympy.core.containers import Tuple
from sympy import Matrix
from sympy.tensor import IndexedBase
from sympy import symbols

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
    assert(Jacobian(F) == expected)
    # ...

    # ...
    expected = dx(F[0])
    assert(DetJacobian(F) == expected)
    # ...

#    print(DetJacobian(F))
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
    expected = Tuple(Matrix([[dx(F[0]), dx(F[1])]]),
                     Matrix([[dy(F[0]), dy(F[1])]]))
    assert(Jacobian(F) == expected)
    # ...

    # ...
    expected = dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])
    assert(DetJacobian(F) == expected)
    # ...

    # ...
    expected = Tuple(a*dy(F[1])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])) -
                     b*dy(F[0])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])),
                     -a*dx(F[1])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])) +
                     b*dx(F[0])/(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0])))
    assert(Covariant(F, ab) == expected)
    # ...

    # ...
    expected = Tuple((a*dx(F[0]) + b*dx(F[1]))/(dx(F[0])*dy(F[1]) -
                                                dx(F[1])*dy(F[0])), (a*dy(F[0])
                                                                     +
                                                                     b*dy(F[1]))/(dx(F[0])*dy(F[1])
                                                                                  -
                                                                                  dx(F[1])*dy(F[0])))
    assert(Contravariant(F, ab) == expected)
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
    expected = Tuple(Matrix([[dx(F[0]), dx(F[1]), dx(F[2])]]),
                     Matrix([[dy(F[0]), dy(F[1]), dy(F[2])]]),
                     Matrix([[dz(F[0]), dz(F[1]), dz(F[2])]]))
    assert(Jacobian(F) == expected)
    # ...

    # ...
    expected = (dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) -
                dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) +
                dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0]))
    assert(DetJacobian(F) == expected)
    # ...

#    # ... TODO why does the assert fail?
#    expected = Tuple(a*(((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0])))*dx(F[0])*dy(F[1]) - ((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*dx(F[2]) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*dx(F[1]))*(-(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*dz(F[0]) + (dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0]))*dy(F[0])))/(((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0])))*(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*dx(F[0])) + b*(-((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0])))*dy(F[0]) - (-(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*dz(F[0]) + (dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0]))*dy(F[0]))*(dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0])))/(((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0])))*(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))) + c*(-(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*dz(F[0]) + (dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0]))*dy(F[0]))/((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0]))), a*(-((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0])))*dx(F[0])*dx(F[1]) + ((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*dx(F[2]) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*dx(F[1]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0]))*dx(F[0]))/(((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0])))*(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*dx(F[0])) + b*(((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0])))*dx(F[0]) + (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0]))*dx(F[0]))/(((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0])))*(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))) - c*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0]))*dx(F[0])/((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0]))), -a*((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*dx(F[2]) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*dx(F[1]))/((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0]))) - b*(dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*dx(F[0])/((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0]))) + c*(dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*dx(F[0])/((dx(F[0])*dy(F[1]) - dx(F[1])*dy(F[0]))*(dx(F[0])*dz(F[2]) - dx(F[2])*dz(F[0])) - (dx(F[0])*dy(F[2]) - dx(F[2])*dy(F[0]))*(dx(F[0])*dz(F[1]) - dx(F[1])*dz(F[0]))))
#    assert(Covariant(F, abc) == expected)
#    # ...

    # ...
    expected = Tuple((a*dx(F[0]) + b*dx(F[1]) + c*dx(F[2]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])), (a*dy(F[0]) + b*dy(F[1]) + c*dy(F[2]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])), (a*dz(F[0]) + b*dz(F[1]) + c*dz(F[2]))/(dx(F[0])*dy(F[1])*dz(F[2]) - dx(F[0])*dy(F[2])*dz(F[1]) - dx(F[1])*dy(F[0])*dz(F[2]) + dx(F[1])*dy(F[2])*dz(F[0]) + dx(F[2])*dy(F[0])*dz(F[1]) - dx(F[2])*dy(F[1])*dz(F[0])))
    assert(Contravariant(F, abc) == expected)
    # ...

#    print(Covariant(F, abc))
#    print(Contravariant(F, abc))

# ...

# .....................................................
if __name__ == '__main__':

    test_1d()
    test_2d()
    test_3d()
