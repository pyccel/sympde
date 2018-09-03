# coding: utf-8

from sympde.core import Line, Square, Cube

# ...
def test_1d():
    print('============ test_1d ==============')

    I = Line(0, 1)
    print(I)
# ...

# ...
def test_2d():
    print('============ test_2d ==============')

    I = Square([0, 0], [1, 1])
    print(I)

# ...
def test_3d():
    print('============ test_3d ==============')

    I = Cube([0, 0, 0], [1, 1, 1])
    print(I)
# ...

# .....................................................
if __name__ == '__main__':

    test_1d()
    test_2d()
    test_3d()
