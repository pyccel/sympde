# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Integer, Float

# ...
class Constant(Symbol):
    """
    Represents a constant symbol.

    Examples

    """
    _label = ''
    def __new__(cls, *args, **kwargs):
        label = kwargs.pop('label', '')

        obj = Symbol.__new__(cls, *args, **kwargs)
        obj._label = label
        return obj

    @property
    def label(self):
        return self._label
# ...

# ...
class Field(Symbol):
    """
    Represents a Field variable.

    Examples

    """
    _space = None
    is_commutative = True
    def __new__(cls, name, space=None):
        obj =  Basic.__new__(cls, name)
        obj._space = space
        return obj

    @property
    def space(self):
        return self._space

    @property
    def name(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)
# ...

class CalculusFunction(Function):
    """this class is needed to distinguish between functions and calculus
    functions when manipulating our expressions"""
    pass

class DottedName(Basic):
    """
    Represents a dotted variable.

    Examples

    >>> from pyccel.ast.core import DottedName
    >>> DottedName('matrix', 'n_rows')
    matrix.n_rows
    >>> DottedName('pyccel', 'stdlib', 'parallel')
    pyccel.stdlib.parallel
    """
    def __new__(cls, *args):
        return Basic.__new__(cls, *args)

    @property
    def name(self):
        return self._args

    def __str__(self):
        return '.'.join(str(n) for n in self.name)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '.'.join(sstr(n) for n in self.name)


# ...
class Mapping(IndexedBase):
    """
    Represents a Mapping object.

    Examples

    """
    def __new__(cls, name, rdim, coordinates=None):
        if isinstance(rdim, (tuple, list, Tuple)):
            if not len(rdim) == 1:
                raise ValueError('> Expecting a tuple, list, Tuple of length 1')

            rdim = rdim[0]

        obj = IndexedBase.__new__(cls, name, shape=(rdim))

        if coordinates is None:
            _coordinates = [Symbol(name) for name in ['x', 'y', 'z'][:rdim]]
        else:
            if not isinstance(coordinates, (list, tuple, Tuple)):
                raise TypeError('> Expecting list, tuple, Tuple')

            for a in coordinates:
                if not isinstance(a, (str, Symbol)):
                    raise TypeError('> Expecting str or Symbol')

            _coordinates = [Symbol(name) for name in coordinates]

        obj._name = name
        obj._rdim = rdim
        obj._coordinates = _coordinates
        return obj

    @property
    def name(self):
        return self._name

    @property
    def rdim(self):
        return self._rdim

    @property
    def coordinates(self):
        if self.rdim == 1:
            return self._coordinates[0]
        else:
            return self._coordinates

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)
# ...

_coeffs_registery = (Integer, Float, Constant)
