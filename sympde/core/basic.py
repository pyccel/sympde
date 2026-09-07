# coding: utf-8

from sympy        import Function
from sympy        import Number
from sympy        import NumberSymbol
from sympy.core   import Basic
from sympy.core   import Symbol
from sympy.tensor import IndexedBase

#==============================================================================
class Constant(Symbol):
    """
    Represents a constant symbol.

    Examples

    """
    _label = ''
    is_number = True
    def __new__(cls, *args, **kwargs):
        label = kwargs.pop('label', '')

        obj = Symbol.__new__(cls, *args, **kwargs)
        obj._label = label
        return obj

    @property
    def label(self):
        return self._label

#==============================================================================
class CalculusFunction(Function):
    """this class is needed to distinguish between functions and calculus
    functions when manipulating our expressions"""
    pass

#==============================================================================
class SymbolicMapping(IndexedBase):
    """
    Common root of the unified mapping hierarchy: a symbolic transformation of
    coordinates identified by a name and a pair of dimensions (logical ``ldim``
    to physical ``pdim``).

    A ``SymbolicMapping`` may be undefined (name and dimensions only) or carry
    more structure in a subclass. It stays callable on a *domain*, returning a
    symbolic mapped domain; point evaluation is the responsibility of
    ``DefinedMapping``.

    Lives here (rather than in ``sympde.topology.mapping``) so that leaf modules
    such as ``sympde.topology.derivatives`` can type-check against it without a
    circular import.
    """

# Deprecated alias for the pre-WP06 name; removed in WP06d-4.
BasicMapping = SymbolicMapping

#==============================================================================
class BasicDerivable(Basic):
    pass

#==============================================================================
_coeffs_registery = (int, float, complex, Number, NumberSymbol, Constant)
