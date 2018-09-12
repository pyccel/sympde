# coding: utf-8

from sympy.core import Basic
from sympy.core.containers import Tuple

class BasicDomain(Basic):
    _name = None
    _dim = None

    @property
    def name(self):
        return self._name

    @property
    def dim(self):
        return self._dim

class Domain(BasicDomain):
    """
    Represents an undefined domain.

    Examples

    """
    _name = None
    _dim = None
    def __new__(cls, name, dim):
        obj = Basic.__new__(cls)
        obj._name = name
        obj._dim = dim
        return obj

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(name))

class Line(BasicDomain):
    """
    Represents a 1D line.

    Examples

    """
    _dim = 1
    _name = 'Line'
    def __new__(cls, xmin=0, xmax=1):
        return Basic.__new__(cls, Tuple(xmin, xmax))

    @property
    def bounds(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        xmin, xmax = self.bounds
        return '[{xmin}, {xmax}]'.format(xmin=sstr(xmin),
                                         xmax=sstr(xmax))

class Square(BasicDomain):
    """
    Represents a 2D square.

    Examples

    """
    _dim = 2
    _name = 'Square'

    def __new__(cls, xmin=[0,0], xmax=[1,1]):
        if not isinstance(xmin, (list, tuple, Tuple)):
            raise TypeError('> Expecting a list, tuple or Tuple')

        if not isinstance(xmax, (list, tuple, Tuple)):
            raise TypeError('> Expecting a list, tuple or Tuple')

        if not len(xmin) == 2:
            raise ValueError('> xmin must be of length 2')

        if not len(xmax) == 2:
            raise ValueError('> xmax must be of length 2')

        xmin = Tuple(*xmin)
        xmax = Tuple(*xmax)

        return Basic.__new__(cls, Tuple(xmin, xmax))

    @property
    def bounds(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        xmin, xmax = self.bounds
        xmins = [sstr(i) for i in xmin]
        xmaxs = [sstr(i) for i in xmax]

        intervals = ['[{a}, {b}]'.format(a=xmin,b=xmax) for xmin, xmax in zip(xmins, xmaxs)]
        return ' * '.join(i for i in intervals)

class Cube(BasicDomain):
    """
    Represents a 3D cube.

    Examples

    """
    _dim = 3
    _name = 'Cube'
    def __new__(cls, xmin=[0,0,0], xmax=[1,1,1]):
        if not isinstance(xmin, (list, tuple, Tuple)):
            raise TypeError('> Expecting a list, tuple or Tuple')

        if not isinstance(xmax, (list, tuple, Tuple)):
            raise TypeError('> Expecting a list, tuple or Tuple')

        if not len(xmin) == 3:
            raise ValueError('> xmin must be of length 3')

        if not len(xmax) == 3:
            raise ValueError('> xmax must be of length 3')

        xmin = Tuple(*xmin)
        xmax = Tuple(*xmax)

        return Basic.__new__(cls, Tuple(xmin, xmax))

    @property
    def bounds(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        xmin, xmax = self.bounds
        xmins = [sstr(i) for i in xmin]
        xmaxs = [sstr(i) for i in xmax]

        intervals = ['[{a}, {b}]'.format(a=xmin,b=xmax) for xmin, xmax in zip(xmins, xmaxs)]
        return ' * '.join(i for i in intervals)
