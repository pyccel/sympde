# coding: utf-8

from sympy import Function

from sympde.core.expr import BilinearForm, LinearForm, Integral
from sympde.core.model import Model, Equation
from sympde.core import grad, dot, inner, cross, rot, curl, div, bracket
from sympde.core import FunctionSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import Field
from sympde.core import Domain
from sympde.core import Constant

DIM = 2

# ...
def test_model_2d_1():
    print('============ test_model_2d_1 ==============')

    domain = Domain('\Omega', dim=DIM)

    V = FunctionSpace('V', domain)
    x,y = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    psi = Field('psi', V)
    f = Function('f')

    expr = dot(grad(v), grad(u)) + bracket(u,psi) * v

    a = BilinearForm((v,u), expr, name='a')
    b = LinearForm(v, f(x,y)*v, name='b')

    forms = [a, b]
    equation = Equation(a(v,u), b(v))

    model = Model(domain, forms=forms, equation=equation)
    model.preview(outputTexFile='test_model_2d_1.tex')
# ...

# ...
def test_model_2d_2():
    print('============ test_model_2d_2 ==============')

    domain = Domain('\Omega', dim=DIM)

    V = FunctionSpace('V', domain)
    x,y = V.coordinates

    u = TestFunction(V, name='u')
    v = TestFunction(V, name='v')
    v_1 = TestFunction(V, name='v_1')
    v_2 = TestFunction(V, name='v_2')
    v_3 = TestFunction(V, name='v_3')
    v_4 = TestFunction(V, name='v_4')

    psi = TestFunction(V, name='psi')
    phi = TestFunction(V, name='phi')
    w   = TestFunction(V, name='omega')
    j   = TestFunction(V, name='J')

    psi_n = Field('psi_n',   V)
    phi_n = Field('phi_n',   V)
    w_n   = Field('omega_n', V)
    j_n   = Field('J_n',     V)
    j_c   = Field('J_c',     V)

    nu  = Constant('nu')
    eta = Constant('eta')
    dt  = Constant('dt')

    expr = psi*v - dt*bracket(psi,phi_n)*v #+ dt*nu*dot(grad(psi), grad(v))
    a_1 = BilinearForm((v,psi), expr, name='a_1')

    expr = (w*v -
            dt*bracket(w,phi_n)*v +
            dt*bracket(psi,j_n)*v +
            dt*eta*dot(grad(psi), grad(v)))
    a_2 = BilinearForm((v,w), expr, name='a_2')

    m = BilinearForm((v,u), u*v, name='m')

    expr = dot(grad(phi), grad(v))
    a_4 = BilinearForm((v,phi), expr, name='a_4')

    f = Function('f')
    b = LinearForm(v, f(x,y)*v, name='b')

    expr = a_1(v_1, psi) + m(v_1, j) + a_2(v_2, w) + m(v_3, j) + a_4(v_4, phi)
    a = BilinearForm(((v_1, v_2, v_3, v_4),(psi, w, j, phi)), expr, name='a')

    forms = [m, a_1, a_2, a_4, a, b]
    equation = Equation(a(v,psi), b(v))

    model = Model(domain, forms=forms, equation=equation)
    model.preview(outputTexFile='test_model_2d_2.tex')
# ...

# .....................................................
if __name__ == '__main__':

#    test_model_2d_1()
    test_model_2d_2()
