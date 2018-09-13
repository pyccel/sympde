# coding: utf-8

from sympde.core.expr import BilinearForm, LinearForm, Integral
from sympde.core.model import Model, Equation
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import FunctionSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import Domain

DIM = 1

# ...
def test_model_1d_1():
    print('============ test_model_1d_1 ==============')

    domain = Domain('\Omega', dim=DIM)

    V = FunctionSpace('V', domain)
    x = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    expr = dot(grad(v), grad(u))

    a = BilinearForm((v,u), expr, name='a')
    b = LinearForm(v, x*(1.-x)*v, name='b')

    forms = [a, b]
    equation = Equation(a(v,u), b(v))

    model = Model(domain, forms=forms, equation=equation)
    model.preview(outputTexFile='test_model_1d_1.tex')
# ...

# .....................................................
if __name__ == '__main__':

    test_model_1d_1()
