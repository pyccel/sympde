# coding: utf-8

from sympde.gallery import Stokes

# ...
def test_stokes_2d():
    print('============ test_stokes_2d ==============')

    model = Stokes(dim=2)

    assert(str(model.space) == 'V')

    model.preview(outputTexFile='test_stokes_2d.tex')
# ...

# ...
def test_stokes_3d():
    print('============ test_stokes_3d ==============')

    model = Stokes(dim=3)

    assert(str(model.space) == 'V')

    model.preview(outputTexFile='test_stokes_3d.tex')
# ...


# .....................................................
if __name__ == '__main__':

    test_stokes_2d()
    test_stokes_3d()
