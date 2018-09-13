# coding: utf-8

from sympy.core import Basic
from sympy.core.containers import Tuple
from sympy.tensor import IndexedBase

class BasicDomain(Basic):
    _dim = None

    @property
    def dim(self):
        return self._dim

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

class Domain(BasicDomain):
    """
    Represents an undefined domain.

    Examples

    """
    def __new__(cls, name, dim):
        obj = Basic.__new__(cls, name)
        obj._dim = dim
        return obj

    @property
    def name(self):
        return self._args[0]

class Boundary(BasicDomain):
    """
    Represents an undefined boundary over a domain.

    Examples

    """
    def __new__(cls, name, domain):
        obj = Basic.__new__(cls, name)
        obj._domain = domain
        return obj

    @property
    def name(self):
        return self._args[0]

    @property
    def domain(self):
        return self._domain

    @property
    def dim(self):
        return self.domain.dim

class BoundaryVector(IndexedBase):

    def __new__(cls, boundary):
        if not isinstance(boundary, Boundary):
            raise TypeError('> Expecting a Boundary')

        return Basic.__new__(cls, boundary)

    @property
    def boundary(self):
        return self._args[0]

    @property
    def domain(self):
        return self.boundary.domain

    @property
    def shape(self):
        return [self.boundary.dim]

class NormalVector(BoundaryVector):
    pass

class TangentVector(BoundaryVector):
    pass


class Line(BasicDomain):
    """
    Represents a 1D line.

    Examples

    """
    _dim = 1
    def __new__(cls, xmin=0, xmax=1):
        return Basic.__new__(cls, Tuple(xmin, xmax))

    @property
    def bounds(self):
        return self._args[0]

    @property
    def name(self):
        return 'Line'

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

    @property
    def name(self):
        return 'Square'

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

    @property
    def name(self):
        return 'Cube'

    def _sympystr(self, printer):
        sstr = printer.doprint
        xmin, xmax = self.bounds
        xmins = [sstr(i) for i in xmin]
        xmaxs = [sstr(i) for i in xmax]

        intervals = ['[{a}, {b}]'.format(a=xmin,b=xmax) for xmin, xmax in zip(xmins, xmaxs)]
        return ' * '.join(i for i in intervals)
