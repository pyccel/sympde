# coding: utf-8


from collections import OrderedDict
from collections import abc

from sympy.core import Basic, Symbol
from sympy.core.containers import Tuple
from sympy.tensor import IndexedBase

from .basic import BasicDomain, InteriorDomain, Interval, Boundary
from .basic import ProductDomain
from .domain import Domain



#==============================================================================
class Line(Domain):
    def __new__(cls, name=None):
        if name is None:
            name = 'Ix'

        x  = Symbol('x')
        Ix = Interval(name, coordinate=x)

        Gamma_1 = Boundary('Gamma_1', Ix, axis=0, ext=-1)
        Gamma_2 = Boundary('Gamma_2', Ix, axis=0, ext=1)

        return Domain('Line', interiors=[Ix],
                      boundaries=[Gamma_1, Gamma_2])

#==============================================================================
class Square(Domain):
    def __new__(cls, name=None):
        if name is None:
            name = 'Square'

        x  = Symbol('x')
        Ix = Interval('Ix', coordinate=x)

        y  = Symbol('y')
        Iy = Interval('Iy', coordinate=y)

        interior = ProductDomain(Ix, Iy, name=name)

        boundaries = []
        i = 1
        for axis in range(interior.dim):
            for ext in [-1, 1]:
                name = 'Gamma_{}'.format(i)
                Gamma = Boundary(name, interior, axis=axis, ext=ext)
                boundaries += [Gamma]

                i += 1

        interior = InteriorDomain(interior)

        return Domain(name, interiors=[interior],
                      boundaries=boundaries)


#==============================================================================
class Cube(Domain):
    def __new__(cls, name=None):
        if name is None:
            name = 'Cube'

        x  = Symbol('x')
        Ix = Interval('Ix', coordinate=x)

        y  = Symbol('y')
        Iy = Interval('Iy', coordinate=y)

        z  = Symbol('z')
        Iz = Interval('Iz', coordinate=z)

        interior = ProductDomain(Ix, Iy, Iz, name=name)

        boundaries = []
        i = 1
        for axis in range(interior.dim):
            for ext in [-1, 1]:
                name = 'Gamma_{}'.format(i)
                I = interior.domains[axis]
                Gamma = Boundary(name, I, axis=axis, ext=ext)
                boundaries += [Gamma]

                i += 1

        interior = InteriorDomain(interior)

        return Domain(name, interiors=[interior],
                      boundaries=boundaries)

