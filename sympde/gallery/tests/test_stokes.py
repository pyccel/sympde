# coding: utf-8

from sympde.core import Domain
from sympde.gallery import Stokes

# ...
def test_stokes_2d():
    print('============ test_stokes_2d ==============')

    domain = Domain('Omega', dim=2)
    model = Stokes(domain=domain)

    model.preview(outputTexFile='test_stokes_2d.tex')
# ...

# ...
def test_stokes_3d():
    print('============ test_stokes_3d ==============')

    domain = Domain('Omega', dim=3)
    model = Stokes(domain=domain)

    model.preview(outputTexFile='test_stokes_3d.tex')
# ...


# .....................................................
if __name__ == '__main__':

    test_stokes_2d()
#    test_stokes_3d()
