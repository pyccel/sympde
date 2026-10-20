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
    """SymPy-compatible representation of `None` in expression arguments."""


_none = _NoneArgument()


def _sympify_argument(value):
    """Convert structural metadata to objects accepted in `Basic.args`."""
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


def new_basic(cls, *args, **options):
    """Construct a `Basic` after normalizing all structural arguments."""
    args = tuple(_sympify_argument(arg) for arg in args)
    return Basic.__new__(cls, *args, **options)


def restore_optional_argument(value):
    """Restore an optional value stored in a symbolic argument list."""
    return None if value is _none else value


def restore_string_argument(value):
    """Restore a Python string stored as a SymPy `Str` object."""
    value = restore_optional_argument(value)
    return value.name if isinstance(value, Str) else value


def restore_int_argument(value):
    """Restore a Python integer stored as a SymPy `Integer` object."""
    value = restore_optional_argument(value)
    return int(value) if isinstance(value, Integer) else value


def restore_container_argument(value):
    """Restore Python containers used as non-expression metadata."""
    value = restore_optional_argument(value)
    if isinstance(value, Str):
        return value.name
    if isinstance(value, Dict):
        return {
            restore_container_argument(k): restore_container_argument(v)
            for k, v in value.items()
        }
    if isinstance(value, Tuple):
        return tuple(restore_container_argument(v) for v in value)
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
coefficient_types = (int, float, complex, Number, NumberSymbol, Constant)
