# coding: utf-8

# TODO not working yet

from sympde.core import Domain
from sympde.gallery import Wave

# ...
def test_wave_1d():
    print('============ test_wave_1d ==============')

    domain = Domain(r'\Omega', dim=1)
    model = Wave(domain=domain)

    assert(str(model.space) == 'V')

    model.preview(outputTexFile='test_wave_1d.tex')
# ...


# .....................................................
if __name__ == '__main__':

    test_wave_1d()
