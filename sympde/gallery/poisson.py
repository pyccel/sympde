# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Integer, Float

from sympde.core.expr import BilinearForm, LinearForm, FunctionForm
from sympde.core.model import Model
from sympde.core import grad, dot
from sympde.core import FunctionSpace
from sympde.core import TestFunction

class Poisson(Model):
    """
    Represents a mathematical model for Poisson.

    Examples

    """
    def __new__(cls, dim, *args, **kwargs):

        # ... abstract model
        V = FunctionSpace('V', ldim=dim)

        v = TestFunction(V, name='v')
        u = TestFunction(V, name='u')

        expr = dot(grad(v), grad(u))

        a = BilinearForm((v,u), expr)
        # ...

        obj = Model.__new__(cls, a, *args, **kwargs)

        obj._space = V

        return obj

    @property
    def space(self):
        return self._space
