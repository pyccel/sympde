# coding: utf-8

from sympy import sin, cos, pi

from sympde.core import dx, dy, dz
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import FunctionSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import BilinearForm, LinearForm, Integral
from sympde.core import Field
from sympde.printing.latex import latex


def test_latex_1d():
    V = FunctionSpace('V', ldim=1)

    x = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    F = Field('F', space=V)

    assert(latex(grad(v)) == r'\nabla{v}')
    assert(latex(dot(grad(v), grad(u))) == r'\nabla{v} \cdot \nabla{u}')

    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    print(latex(a))
#    assert(latex(a) == r'\int_{0}^{1} \nabla{v} \cdot \nabla{u} dx')

    b = LinearForm(v, sin(pi*x)*v)
    print(latex(b))
#    assert(latex(b) == r'\int_{0}^{1} v \sin{\left (\pi x \right )} dx')

    f = Integral(dx(F)-x)
    print(latex(f))
#    assert(latex(f) == r'\int_{0}^{1} - x + \partial_{x}F dx')

def test_latex_2d_1():
    V = FunctionSpace('V', ldim=2)

    x,y = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    F = Field('F', space=V)

    assert(latex(grad(v)) == r'\nabla{v}')
    assert(latex(dot(grad(v), grad(u))) == r'\nabla{v} \cdot \nabla{u}')

    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    print(latex(a))
#    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1} \nabla{v} \cdot \nabla{u} dxdy')

    b = LinearForm(v, sin(pi*x)*cos(pi*y)*v)
    print(latex(b))
#    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1} v \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} dxdy')

    f = Integral(dx(F)-dy(F)-x*y)
    print(latex(f))
#    assert(latex(f) == r'\int_{0}^{1}\int_{0}^{1} - x y + \partial_{x}F - \partial_{y}F dxdy')

def test_latex_2d_2():
    V = FunctionSpace('V', ldim=2, is_block=True, shape=2)

    x,y = V.coordinates

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')

    F = Field('F', space=V)

    assert(latex(v) == r'\mathbf{v}')
    assert(latex(inner(grad(v), grad(u))) == r'\nabla{\mathbf{v}} : \nabla{\mathbf{u}}')

    a = BilinearForm((v,u), inner(grad(v), grad(u)))
    print(latex(a))
#    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1} \nabla{\mathbf{v}} : \nabla{\mathbf{u}} dxdy')

    b = LinearForm(v, sin(pi*x)*cos(pi*y)*div(v))
    print(latex(b))
#    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1} \nabla \cdot \mathbf{v} \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} dxdy')


def test_latex_3d_1():
    V = FunctionSpace('V', ldim=3)

    x,y,z = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    F = Field('F', space=V)

    assert(latex(grad(v)) == r'\nabla{v}')
    assert(latex(dot(grad(v), grad(u))) == r'\nabla{v} \cdot \nabla{u}')

    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    print(latex(a))
#    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} \nabla{v} \cdot \nabla{u} dxdydz')

    b = LinearForm(v, sin(pi*x)*cos(pi*y)*cos(2*pi*z)*v)
    print(latex(b))
#    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} v \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} \cos{\left (2 \pi z \right )} dxdydz')

    f = Integral(dx(F)-dy(F)+dz(F)-x*y*z)
    print(latex(f))
#    assert(latex(f) == r'\int_{0}^{1}\int_{0}^{1} - x y z + \partial_{x}F - \partial_{y}F + \partial_{z}F dxdy')

def test_latex_3d_2():
    V = FunctionSpace('V', ldim=3, is_block=True, shape=3)

    x,y,z = V.coordinates

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')

    F = Field('F', space=V)

    assert(latex(v) == r'\mathbf{v}')
    assert(latex(inner(grad(v), grad(u))) == r'\nabla{\mathbf{v}} : \nabla{\mathbf{u}}')

    a = BilinearForm((v,u), inner(grad(v), grad(u)))
    print(latex(a))
#    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} \nabla{\mathbf{v}} : \nabla{\mathbf{u}} dxdydz')

    b = LinearForm(v, sin(pi*x)*cos(pi*y)*div(v))
    print(latex(b))
#    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} \nabla \cdot \mathbf{v} \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} dxdydz')


def test_latex_model_2d_1():
#    from sympde.gallery import Poisson
#    model = Poisson(dim=2)
#    model.preview(outputTexFile='poisson_2d.tex')

    from sympde.gallery import Stokes
    model = Stokes(dim=2)
    model.preview(outputTexFile='stokes_2d.tex')

####################
if __name__ == '__main__':
#    test_latex_1d()
#    test_latex_2d_1()
#    test_latex_3d_1()
#
#    test_latex_2d_2()
#    test_latex_3d_2()

    test_latex_model_2d_1()

