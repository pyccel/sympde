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
        dim = kwargs.pop('dim', None)
        if dim is None:
            raise ValueError('> Expecting a dim entry')

        if not(dim in [2, 3]):
            raise ValueError('> only 2d and 3d models are possible')
        # ...

        # ... abstract model
        V = FunctionSpace('V', ldim=dim, is_block=True, shape=dim)
        W = FunctionSpace('W', ldim=dim)

        v = VectorTestFunction(V, name='v')
        u = VectorTestFunction(V, name='u')
        p = TestFunction(W, name='p')
        q = TestFunction(W, name='q')

        a1 = BilinearForm((v,u), inner(grad(v), grad(u)))
        a2 = BilinearForm((v,p), div(v)*p)
        a  = BilinearForm(((v,q), (u,p)), a1(v,u) - a2(v,p) + a2(u,q))
        # ...

#        obj = Model.__new__(cls, a, a1, a2, **kwargs)
        obj = Model.__new__(cls, a1=a1, a2=a2, **kwargs)

        obj._spaces = [V, W]

        return obj

    @property
    def spaces(self):
        return self._spaces
