# coding: utf-8

# TODO improve this file

from sympy.core import Basic
from sympy.core.containers import Tuple

class BasicMeasure(Basic):

    def _sympystr(self, printer):
        sstr = printer.doprint
        txt = ''.join(j for j in ['d{}'.format(sstr(i)) for i in self.args])
        return txt

class CanonicalMeasure(BasicMeasure):
    """
    Represents a canonical measure.

    Examples

    """
    def __new__(cls, ldim):
        if not isinstance(ldim, int):
            raise TypeError('> Expecting an integer')

        coordinates = ['x1', 'x2', 'x3'][:ldim]
        coordinates = Tuple(*coordinates)

        return Basic.__new__(cls, coordinates)

    @property
    def args(self):
        return self._args[0]

class CartesianMeasure(BasicMeasure):
    """
    Represents a cartesian measure.

    Examples

    """
    def __new__(cls, ldim):
        if not isinstance(ldim, int):
            raise TypeError('> Expecting an integer')

        coordinates = ['x', 'y', 'z'][:ldim]
        coordinates = Tuple(*coordinates)

        return Basic.__new__(cls, coordinates)

    @property
    def args(self):
        return self._args[0]

class Measure(BasicMeasure):
    """
    Represents a measure from coordinates.

    Examples

    """
    def __new__(cls, coordinates):
        if not isinstance(coordinates, (tuple, list, Tuple)):
            coordinates = [coordinates]

        coordinates = Tuple(*coordinates)

        return Basic.__new__(cls, coordinates)

    @property
    def args(self):
        return self._args[0]
