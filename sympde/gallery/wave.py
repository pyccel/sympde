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
from sympde.core import Constant
from sympde.core.derivatives import dx, dy, dz

class Wave_1d(Model):
    """
    Represents a mathematical model for the 1d Wave problem.

    Examples

    """
    def __new__(cls, **kwargs):
        # ...
        domain = kwargs.pop('domain', None)
        if domain is None:
            raise ValueError('> Expecting a domain entry')

        dim = domain.dim
        if not(dim in [1]):
            raise ValueError('> only 1d is possible')

        mapping = kwargs.pop('mapping', None)
        # ...

        # ... abstract model
        U = FunctionSpace('U', domain)
        V = FunctionSpace('V', domain)

        # trial functions
        u = TestFunction(U, name='u')
        f = TestFunction(V, name='f')

        # test functions
        v = TestFunction(U, name='v')
        w = TestFunction(V, name='w')

        rho = Constant('rho', real=True, label='mass density')
        dt = Constant('dt', real=True, label='time step')

        a = BilinearForm((v,u), v*u, mapping=mapping, name='a')
        b  = BilinearForm((v,u), dx(v)*u, mapping=mapping, name='b')

        expr = rho*a(v,u) + dt*b(v, f) + dt*b(w,u) + a(w,f)
        A = BilinearForm(((v,w), (u,f)), expr, mapping=mapping, name='A')
        # ...

        forms = [a, b, A]
        equation = Equation(A((v,w), (u,f)), None)

        obj = Model.__new__(cls, forms=forms, equation=equation, **kwargs)

        obj._space = V
        obj._domain = domain

        return obj

    @property
    def space(self):
        return self._space

class Wave(Model):
    """
    Represents a mathematical model for Wave.

    Examples

    """
    def __new__(cls, **kwargs):
        try:
            domain = kwargs['domain']
        except:
            raise ValueError('> domain must be provided')

        dim = domain.dim
        construct = eval('Wave_{}d'.format(dim))
        return construct(**kwargs)