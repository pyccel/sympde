# coding: utf-8

from sympde.core.expr import BilinearForm, LinearForm, Integral
from sympde.core.model import Model, Equation
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import FunctionSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction

DIM = 1

# ...
def test_model_1d_1():
    print('============ test_model_1d_1 ==============')

    V = FunctionSpace('V', ldim=DIM)
    x = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    expr = dot(grad(v), grad(u))

    a = BilinearForm((v,u), expr)
    b = LinearForm(v, x*(1.-x)*v)

    forms = {'a': a,
             'b': b}
    equations = (Equation(a, b),
                )

    model = Model(forms=forms, equations=equations)
    model.preview(outputTexFile='test_model_1d_1.tex')
# ...

# .....................................................
if __name__ == '__main__':

    test_model_1d_1()
