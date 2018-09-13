# coding: utf-8

from sympde.core import Domain
from sympde.gallery import Poisson

# ...
def test_poisson_1d():
    print('============ test_poisson_1d ==============')

    domain = Domain(r'\Omega', dim=1)
    model = Poisson(domain)

    assert(str(model.space) == 'V')

    model.preview(outputTexFile='test_poisson_1d.tex')
# ...

# ...
def test_poisson_2d():
    print('============ test_poisson_2d ==============')

    domain = Domain(r'\Omega', dim=2)
    model = Poisson(domain)

    assert(str(model.space) == 'V')

    model.preview(outputTexFile='test_poisson_2d.tex')
# ...

# ...
def test_poisson_3d():
    print('============ test_poisson_3d ==============')

    domain = Domain(r'\Omega', dim=3)
    model = Poisson(domain)

    assert(str(model.space) == 'V')

    model.preview(outputTexFile='test_poisson_3d.tex')
# ...


# .....................................................
if __name__ == '__main__':

    test_poisson_1d()
    test_poisson_2d()
    test_poisson_3d()
