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

# .....................................................
if __name__ == '__main__':

    test_model_2d_1()
