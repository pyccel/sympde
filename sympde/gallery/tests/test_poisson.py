# coding: utf-8

from sympde.gallery import Poisson

# ...
def test_poisson_1d():
    print('============ test_poisson_1d ==============')

    model = Poisson(dim=1)

    assert(str(model.space) == 'V')
# ...

# ...
def test_poisson_2d():
    print('============ test_poisson_2d ==============')

    model = Poisson(dim=2)

    assert(str(model.space) == 'V')
# ...

# ...
def test_poisson_3d():
    print('============ test_poisson_3d ==============')

    model = Poisson(dim=3)

    assert(str(model.space) == 'V')
# ...


# .....................................................
if __name__ == '__main__':

    test_poisson_1d()
    test_poisson_2d()
    test_poisson_3d()
