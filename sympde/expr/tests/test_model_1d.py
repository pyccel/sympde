# coding: utf-8

from sympy import Function

from sympde.expr import BilinearForm, LinearForm, Integral
from sympde.expr import Equation, DirichletBC
from sympde.expr import Model
from sympde.core import grad, dot
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import TestFunction
from sympde.topology import Domain

DIM = 1
domain = Domain('\Omega', dim=DIM)

#==============================================================================
def test_model_1d_1():

    V = FunctionSpace('V', domain)
    x = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    f = Function('f')

    expr = dot(grad(v), grad(u))

    a = BilinearForm((v,u), expr, name='a')
    b = LinearForm(v, f(x)*v, name='b')

    forms = [a, b]
    equation = Equation(a(v,u), b(v))

    model = Model(domain, forms=forms, equation=equation)

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
