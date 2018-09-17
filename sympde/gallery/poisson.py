# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Integer, Float

from sympde.core.expr import BilinearForm, LinearForm, Integral
from sympde.core.model import Model, Equation
from sympde.core import grad, dot
from sympde.core import FunctionSpace
from sympde.core import TestFunction

class Poisson(Model):
    """
    Represents a mathematical model for Poisson.

    Examples

    """
    def __new__(cls, domain, **kwargs):
        # ...
        dim = domain.dim
        if not(dim in [1, 2, 3]):
            raise ValueError('> only 1d, 2d and 3d models are possible')
        # ...

        # ... abstract model
        V = FunctionSpace('V', domain)

        v = TestFunction(V, name='v')
        u = TestFunction(V, name='u')

        a = BilinearForm((v,u), dot(grad(v), grad(u)), name='a')
        # ...

        # ... rhs as undefined function
        xyz = domain.coordinates
        f = Function('f')
        if domain.dim == 1: xyz = [xyz]
        l = LinearForm(v, f(*xyz)*v, name='l')
        # ...

        forms = [a, l]
        equation = Equation(a(v,u), l(v))

        obj = Model.__new__(cls, domain, forms=forms, equation=equation, **kwargs)

        obj._space = V
        obj._domain = domain

        return obj

    @property
    def space(self):
        return self._space
