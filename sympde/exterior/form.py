# coding: utf-8

from sympy.core import Basic
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

        return Basic.__new__(cls, name, index, dim)

    @property
    def name(self):
        return self._args[0]

    @property
    def index(self):
        return self._args[1]

    @property
    def dim(self):
        return self._args[2]
