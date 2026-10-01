# coding: utf-8

from sympy        import Function
from sympy        import Number
from sympy        import NumberSymbol
from sympy        import Integer
from sympy        import sympify
from sympy.core   import Basic
from sympy.core   import Symbol
from sympy.core.containers import Dict, Tuple
from sympy.core.symbol import Str
from sympy.tensor import IndexedBase


class _NoneArgument(Basic):
    """SymPy-compatible representation of ``None`` in expression arguments."""


_none = _NoneArgument()


def _sympify_argument(value):
    """Convert structural metadata to objects accepted in ``Basic.args``."""
    if isinstance(value, Basic):
        return value
    if value is None:
        return _none
    if isinstance(value, str):
        return Str(value)
    if isinstance(value, dict):
        return Dict(*((_sympify_argument(k), _sympify_argument(v))
                      for k, v in value.items()))
    if isinstance(value, (tuple, list)):
        return Tuple(*(_sympify_argument(v) for v in value))
    return sympify(value)


def _new_basic(cls, *args, **options):
    """Construct a ``Basic`` after normalizing all structural arguments."""
    args = tuple(_sympify_argument(arg) for arg in args)
    return Basic.__new__(cls, *args, **options)


def _is_none_argument(value):
    return value is _none


def _restore_optional_argument(value):
    """Restore an optional value stored in a symbolic argument list."""
    return None if _is_none_argument(value) else value


def _restore_string_argument(value):
    """Restore a Python string stored as a SymPy ``Str`` object."""
    return value.name if isinstance(value, Str) else value


def _restore_int_argument(value):
    """Restore a Python integer stored as a SymPy ``Integer`` object."""
    value = _restore_optional_argument(value)
    return int(value) if isinstance(value, Integer) else value


def _restore_container_argument(value):
    """Restore Python containers used as non-expression metadata."""
    value = _restore_optional_argument(value)
    if isinstance(value, Str):
        return value.name
    if isinstance(value, Dict):
        return {
            _restore_container_argument(k): _restore_container_argument(v)
            for k, v in value.items()
        }
    if isinstance(value, Tuple):
        return tuple(_restore_container_argument(v) for v in value)
    if isinstance(value, Integer):
        return int(value)
    return value

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
class BasicMapping(IndexedBase):
    """
    Represents a basic class for mapping.
    """
    pass

#==============================================================================
class BasicDerivable(Basic):
    pass

#==============================================================================
_coeffs_registery = (int, float, complex, Number, NumberSymbol, Constant)
