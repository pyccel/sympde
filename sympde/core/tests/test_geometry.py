# coding: utf-8

from sympy.tensor import Indexed

from sympde.core import Line, Square, Cube
from sympde.core import Domain, Boundary, NormalVector, TangentVector

# ...
def test_geometry_1d():
    print('============ test_geometry_1d ==============')

    I = Line(0, 1)
    print(I)
# ...

# ...
def test_geometry_2d():
    print('============ test_geometry_2d ==============')

    I = Square([0, 0], [1, 1])
    print(I)

    D = Domain('Omega', dim=2)
    B1 = Boundary('\Gamma_1', D)
    B2 = Boundary('\Gamma_2', D)

    nn = NormalVector('n')
    tt = TangentVector('t')

    assert(isinstance(nn[0], Indexed))

# ...
def test_geometry_3d():
    print('============ test_geometry_3d ==============')

    I = Cube([0, 0, 0], [1, 1, 1])
    print(I)
# ...

# .....................................................
if __name__ == '__main__':

    test_geometry_1d()
    test_geometry_2d()
    test_geometry_3d()
