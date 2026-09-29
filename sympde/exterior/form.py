# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Integer, Float
from sympy.core.singleton import Singleton
from sympy.core import Add, Mul
from sympy.core.singleton import S

from sympde.core.basic import _coeffs_registery
from sympde.core import LinearOperator

from .datatype import get_index_form

#==============================================================================
class DifferentialForm(Symbol):
    """
    Represents a differential form symbol.

    Examples

    """
    def __new__(cls, name, index, dim):
        if not isinstance(name, str):
            raise TypeError('> Expecting a string for name')

        assert(isinstance(dim, (int, Symbol)))

        index = get_index_form(index)

        obj = Symbol.__xnew__(cls, name)
        obj._index = index
        obj._dim = Integer(dim) if isinstance(dim, int) else dim
        return obj

    @property
    def index(self):
        return self._index

    @property
    def dim(self):
        return self._dim

    def _hashable_content(self):
        return Symbol._hashable_content(self) + (self.index, self.dim)
