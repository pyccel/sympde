# coding: utf-8

from sympy.core import Integer
from sympy.core import Symbol

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
