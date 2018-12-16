# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Integer, Float

from sympde.core.generic   import grad, dot, inner, cross, rot, curl, div
from sympde.expr.expr      import BilinearForm, LinearForm, Integral
from sympde.topology.space import FunctionSpace, VectorFunctionSpace
from sympde.topology.space import TestFunction
from sympde.topology.space import VectorTestFunction
from sympde.expr.equation  import Equation
from sympde.expr.model     import Model

class Stokes(Model):
    """
    Represents a mathematical model for Stokes.

    Examples

    """
    def __new__(cls, domain, **kwargs):
        # ...
        dim = domain.dim
        if not(dim in [2, 3]):
            raise ValueError('> only 2d and 3d models are possible')
        # ...

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

        # ... rhs as undefined function
        xyz = domain.coordinates
        f = Function('f')
        l = LinearForm((v,q), f(*xyz)*v, name='l')
        # ...

        forms = [a, b, A, l]
        equation = Equation(A((v,q),(u,p)), l(v,q))

        obj = Model.__new__(cls, domain, forms=forms, equation=equation, **kwargs)

        obj._spaces = [V, W]
        obj._domain = domain

        return obj

    @property
    def spaces(self):
        return self._spaces
