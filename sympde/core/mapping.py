# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Integer, Float
from sympy import Indexed, IndexedBase, Matrix, ImmutableDenseMatrix

from .basic import BasicMapping
from .derivatives import dx, dy, dz
from .derivatives import (Grad_1d, Div_1d,
                          Grad_2d, Curl_2d, Rot_2d, Div_2d,
                          Grad_3d, Curl_3d, Div_3d)
from .generic import Dot, Inner, Cross
from .generic import Grad, Rot, Curl, Div
from .generic import _generic_ops

from .algebra import (Dot_1d,
                      Dot_2d, Inner_2d, Cross_2d,
                      Dot_3d, Inner_3d, Cross_3d)


# ...
class Mapping(BasicMapping):
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


class MappingApplication(Function):
    nargs = None

    def __new__(cls, *args, **options):

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

class Jacobian(MappingApplication):
    """

    Examples

    """

    @classmethod
    def eval(cls, F):

        if not isinstance(F, Mapping):
            raise TypeError('> Expecting a Mapping object')

        rdim = F.rdim

        F = [F[i] for i in range(0, F.rdim)]
        F = Tuple(*F)

        if rdim == 1:
            return Grad_1d(F)

        elif rdim == 2:
            return Grad_2d(F)

        elif rdim == 3:
            return Grad_3d(F)

class DetJacobian(MappingApplication):
    """

    Examples

    """

    @classmethod
    def eval(cls, F):

        if not isinstance(F, Mapping):
            raise TypeError('> Expecting a Mapping object')

        J = Matrix(Jacobian(F))
        dim = F.rdim

        if dim == 2:
            det = J[0,0]* J[1,1] - J[0,1]* J[1,0]
            return det

        elif dim == 3:
            det = (J[0, 0]*J[1, 1]*J[2, 2] -
                   J[0, 0]*J[1, 2]*J[2, 1] -
                   J[0, 1]*J[1, 0]*J[2, 2] +
                   J[0, 1]*J[1, 2]*J[2, 0] +
                   J[0, 2]*J[1, 0]*J[2, 1] -
                   J[0, 2]*J[1, 1]*J[2, 0])
            return det

        else:
            return J.det()

class InvDetJacobian(MappingApplication):
    """

    Examples

    """

    @classmethod
    def eval(cls, F):

        if not isinstance(F, Mapping):
            raise TypeError('> Expecting a Mapping object')

        det = DetJacobian(F)
        return 1/det

class Covariant(MappingApplication):
    """

    Examples

    """

    @classmethod
    def eval(cls, F, v):

        if not isinstance(F, Mapping):
            raise TypeError('> Expecting a Mapping')

        if not isinstance(v, (tuple, list, Tuple, Matrix)):
            raise TypeError('> Expecting a tuple, list, Tuple, Matrix')

        v = Matrix(v)
        J = Matrix(Jacobian(F))
        dim = F.rdim
        if dim == 1:
            j_inv = 1/J[0,0]
            v = v[0,0]
            v = [j_inv * v]
            return Tuple(*v)

        elif dim == 2:
            det = J[0,0]* J[1,1] - J[0,1]* J[1,0]
            J_inv = Matrix([[J[1,1], -J[0,1]], [-J[1,0], J[0,0]]])
            M = J_inv.transpose() / det
            v = M*v
            return Tuple(*v)

        elif dim == 3:
            det = (J[0, 0]*J[1, 1]*J[2, 2] -
                   J[0, 0]*J[1, 2]*J[2, 1] -
                   J[0, 1]*J[1, 0]*J[2, 2] +
                   J[0, 1]*J[1, 2]*J[2, 0] +
                   J[0, 2]*J[1, 0]*J[2, 1] -
                   J[0, 2]*J[1, 1]*J[2, 0])

            J_inv = Matrix([[(((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0]))*J[0, 0]*J[1, 1] - ((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*J[0, 2] - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*J[0, 1])*(-(J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*J[2, 0] + (J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0])*J[1, 0]))/(((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0]))*(J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*J[0, 0]), (-((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0]))*J[0, 0]*J[0, 1] + ((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*J[0, 2] - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*J[0, 1])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0])*J[0, 0])/(((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0]))*(J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*J[0, 0]), -((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*J[0, 2] - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*J[0, 1])/((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0]))],
[                                                                       (-((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0]))*J[1, 0] - (-(J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*J[2, 0] + (J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0])*J[1, 0])*(J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0]))/(((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0]))*(J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])),                                                                          (((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0]))*J[0, 0] + (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0])*J[0, 0])/(((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0]))*(J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])),                                                 -(J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*J[0, 0]/((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0]))],
[                                                                                                                                                                                                                                                                                                                  (-(J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*J[2, 0] + (J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0])*J[1, 0])/((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0])),                                                                                                                                                                                                                                                                                                                   -(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0])*J[0, 0]/((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0])),                                                  (J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*J[0, 0]/((J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0])*(J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0]) - (J[0, 0]*J[1, 2] - J[0, 2]*J[1, 0])*(J[0, 0]*J[2, 1] - J[0, 1]*J[2, 0]))]])

            M = J_inv.transpose() / det
            v = M*v
            return Tuple(*v)

        else:
            M = J.inv().transpose()
            v = M*v
            return Tuple(*v)

class Contravariant(MappingApplication):
    """

    Examples

    """

    @classmethod
    def eval(cls, F, v):

        if not isinstance(F, Mapping):
            raise TypeError('> Expecting a Mapping')

        if not isinstance(v, (tuple, list, Tuple, Matrix)):
            raise TypeError('> Expecting a tuple, list, Tuple, Matrix')

        v = Matrix(v)
        J = Matrix(Jacobian(F))
        j = J.det()
        v = [i/j for i in J*v]
        return Tuple(*v)

