# coding: utf-8

from sympy import sin, cos, pi

from sympde.core import dx, dy, dz
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import H1Space
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import BilinearForm, LinearForm, FunctionForm
from sympde.core import Field
from sympde.printing.latex import latex


def test_1d():
    V = H1Space('V', ldim=1)

    x = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    F = Field('F', space=V)

    assert(latex(grad(v)) == r'\nabla{v}')
    assert(latex(dot(grad(v), grad(u))) == r'\nabla{v} \cdot \nabla{u}')

    a = BilinearForm((v,u), dot(grad(v), grad(u)))
#    print(latex(a))
    assert(latex(a) == r'\int_{0}^{1} \nabla{v} \cdot \nabla{u} dx')

    b = LinearForm(v, sin(pi*x)*v)
#    print(latex(b))
    assert(latex(b) == r'\int_{0}^{1} v \sin{\left (\pi x \right )} dx')

    f = FunctionForm(dx(F)-x)
#    print(latex(f))
    assert(latex(f) == r'\int_{0}^{1} - x + \partial_{x}F dx')

def test_2d_1():
    V = H1Space('V', ldim=2)

    x,y = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    F = Field('F', space=V)

    assert(latex(grad(v)) == r'\nabla{v}')
    assert(latex(dot(grad(v), grad(u))) == r'\nabla{v} \cdot \nabla{u}')

    a = BilinearForm((v,u), dot(grad(v), grad(u)))
#    print(latex(a))
    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1} \nabla{v} \cdot \nabla{u} dxdy')

    b = LinearForm(v, sin(pi*x)*cos(pi*y)*v)
#    print(latex(b))
    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1} v \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} dxdy')

    f = FunctionForm(dx(F)-dy(F)-x*y)
#    print(latex(f))
    assert(latex(f) == r'\int_{0}^{1}\int_{0}^{1} - x y + \partial_{x}F - \partial_{y}F dxdy')

def test_2d_2():
    V = H1Space('V', ldim=2, is_block=True, shape=2)

    x,y = V.coordinates

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')

    F = Field('F', space=V)

    assert(latex(v) == r'\mathbf{v}')
    assert(latex(inner(grad(v), grad(u))) == r'\nabla{\mathbf{v}} : \nabla{\mathbf{u}}')

    a = BilinearForm((v,u), inner(grad(v), grad(u)))
#    print(latex(a))
    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1} \nabla{\mathbf{v}} : \nabla{\mathbf{u}} dxdy')

    b = LinearForm(v, sin(pi*x)*cos(pi*y)*div(v))
#    print(latex(b))
    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1} \nabla \cdot \mathbf{v} \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} dxdy')


def test_3d_1():
    V = H1Space('V', ldim=3)

    x,y,z = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    F = Field('F', space=V)

    assert(latex(grad(v)) == r'\nabla{v}')
    assert(latex(dot(grad(v), grad(u))) == r'\nabla{v} \cdot \nabla{u}')

    a = BilinearForm((v,u), dot(grad(v), grad(u)))
#    print(latex(a))
    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} \nabla{v} \cdot \nabla{u} dxdydz')

    b = LinearForm(v, sin(pi*x)*cos(pi*y)*cos(2*pi*z)*v)
#    print(latex(b))
    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} v \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} \cos{\left (2 \pi z \right )} dxdydz')

    f = FunctionForm(dx(F)-dy(F)+dz(F)-x*y*z)
#    print(latex(f))
    assert(latex(f) == r'\int_{0}^{1}\int_{0}^{1} - x y z + \partial_{x}F - \partial_{y}F + \partial_{z}F dxdy')

def test_3d_2():
    V = H1Space('V', ldim=3, is_block=True, shape=3)

    x,y,z = V.coordinates

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')

    F = Field('F', space=V)

    assert(latex(v) == r'\mathbf{v}')
    assert(latex(inner(grad(v), grad(u))) == r'\nabla{\mathbf{v}} : \nabla{\mathbf{u}}')

    a = BilinearForm((v,u), inner(grad(v), grad(u)))
#    print(latex(a))
    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} \nabla{\mathbf{v}} : \nabla{\mathbf{u}} dxdydz')

    b = LinearForm(v, sin(pi*x)*cos(pi*y)*div(v))
#    print(latex(b))
    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} \nabla \cdot \mathbf{v} \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} dxdydz')

####################
if __name__ == '__main__':
    test_1d()
    test_2d_1()
    test_3d_1()

    test_2d_2()
    test_3d_2()
