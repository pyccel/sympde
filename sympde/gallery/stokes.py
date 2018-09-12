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
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import FunctionSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction

class Stokes(Model):
    """
    Represents a mathematical model for Stokes.

    Examples

    """
    def __new__(cls, **kwargs):
        # ...
        domain = kwargs.pop('domain', None)
        if domain is None:
            raise ValueError('> Expecting a domain entry')

        dim = domain.dim
        if not(dim in [2, 3]):
            raise ValueError('> only 2d and 3d models are possible')
        # ...

        # ... abstract model
        V = FunctionSpace('V', domain, is_block=True, shape=dim)
        W = FunctionSpace('W', domain)

        v = VectorTestFunction(V, name='v')
        u = VectorTestFunction(V, name='u')
        p = TestFunction(W, name='p')
        q = TestFunction(W, name='q')

        a = BilinearForm((v,u), inner(grad(v), grad(u)), name='a')
        b = BilinearForm((v,p), div(v)*p, name='b')
        A = BilinearForm(((v,q),(u,p)), a(v,u) - b(v,p) + b(u,q), name='A')
        # ...

        forms = [a, b, A]
        equation = Equation(A((v,u),(q,p)), None)

        obj = Model.__new__(cls, forms=forms, equation=equation, **kwargs)

        obj._spaces = [V, W]
        obj._domain = domain

        return obj

    @property
    def spaces(self):
        return self._spaces
