# coding: utf-8

from sympy import sin, cos, pi

from sympde.core import dx, dy, dz
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import FunctionSpace, VectorFunctionSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import BilinearForm, LinearForm, Integral
from sympde.core import Field
from sympde.core import Domain, Boundary, NormalVector, TangentVector
from sympde.core import Trace, trace_0, trace_1
from sympde.core import tensorize
from sympde.core import VectorField
from sympde.printing.latex import latex


def test_latex_1d():

    DIM = 1
    domain = Domain('Omega', dim=DIM)

    V = FunctionSpace('V', domain)

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

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    V = FunctionSpace('V', domain)

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

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    V = FunctionSpace('V', domain, is_block=True, shape=2)

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

    DIM = 3
    domain = Domain('Omega', dim=DIM)

    V = FunctionSpace('V', domain)

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

    DIM = 3
    domain = Domain('Omega', dim=DIM)

    V = FunctionSpace('V', domain, is_block=True, shape=3)

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
    domain = Domain(r'\Omega', dim=2)
    model = Stokes(domain=domain)
    model.preview(outputTexFile='stokes_2d.tex')

def test_latex_2d_3():
    DIM = 2

    domain = Domain('Omega', dim=DIM)

    B1 = Boundary(r'\Gamma_1', domain)
    B2 = Boundary(r'\Gamma_2', domain)
    B3 = Boundary(r'\Gamma_3', domain)

    V = FunctionSpace('V', domain)

    x = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    # ...
    expr = dot(grad(v), grad(u))
    a_0 = BilinearForm((v,u), expr, name='a_0')

    expr = v*trace_0(u, B1)
    a_bnd = BilinearForm((v, u), expr, name='a_bnd')

    expr = a_0(v,u) + a_bnd(v,u)
    a = BilinearForm((v,u), expr, name='a')
    print(latex(a_0))
    print(latex(a_bnd))
    print(latex(a))
#    print(a)
    print('')
    # ...

def test_latex_2d_4():
    DIM = 2

    domain = Domain('Omega', dim=DIM)

    # ... abstract model
    V = VectorFunctionSpace('V', domain)
    W = FunctionSpace('W', domain)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')
    p = TestFunction(W, name='p')
    q = TestFunction(W, name='q')

    a = BilinearForm((v,u), inner(grad(v), grad(u)), name='a')
    b = BilinearForm((v,p), div(v)*p, name='b')
    A = BilinearForm(((v,q),(u,p)), a(v,u) - b(v,p) + b(u,q), name='A')
    # ...

    print(latex(A))
    print(latex(tensorize(A)))
    print('')
    # ...

def test_latex_2d_5():
    DIM = 2

    domain = Domain('Omega', dim=DIM)

    # ... abstract model
    W1 = VectorFunctionSpace('W1', domain)

    w1 = VectorTestFunction(W1, name='w1')

    F = VectorField(W1, 'F')

    # ...
    l1 = LinearForm(w1, dot(w1, F), name='l1')

    print(latex(l1))
    print('')
    # ...

    # ...
    l2 = LinearForm(w1, rot(w1)*rot(F) + div(w1)*div(F), name='l2')

    print(latex(l2))
    print('')
    # ...


####################
if __name__ == '__main__':
#    test_latex_1d()
#    test_latex_2d_1()
#    test_latex_3d_1()
#
#    test_latex_2d_2()
#    test_latex_3d_2()

#    test_latex_model_2d_1()

#    test_latex_2d_3()
#    test_latex_2d_4()
    test_latex_2d_5()

