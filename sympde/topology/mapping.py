# coding: utf-8

from sympy.core import Basic
from sympy.tensor import IndexedBase
from sympy.core import Symbol
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Matrix

from sympde.core.basic import BasicMapping
from sympde.core.generic import Dot, Inner, Cross
from sympde.core.generic import Grad, Rot, Curl, Div
from sympde.core.generic import _generic_ops
from sympde.core.algebra import (Dot_1d,
                                 Dot_2d, Inner_2d, Cross_2d,
                                 Dot_3d, Inner_3d, Cross_3d)

from .domain import Domain
from .derivatives import dx, dy, dz
from .derivatives import (Grad_1d, Div_1d,
                          Grad_2d, Curl_2d, Rot_2d, Div_2d,
                          Grad_3d, Curl_3d, Div_3d)


# ...
class Mapping(BasicMapping):
    """
    Represents a Mapping object.

    Examples

    """
    # TODO shall we keep rdim ?
    def __new__(cls, name, rdim, domain=None, coordinates=None):
        if isinstance(rdim, (tuple, list, Tuple)):
            if not len(rdim) == 1:
                raise ValueError('> Expecting a tuple, list, Tuple of length 1')

            rdim = rdim[0]

        if not isinstance(domain, Domain):
            raise TypeError('> Expecting a Domain object')

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
        obj._domain = domain
        obj._rdim = rdim
        obj._coordinates = _coordinates

        obj._jacobian = None
        obj._det_jacobian = None
        obj._covariant = None
        obj._contravariant = None
        obj._hessian = None

        return obj

    @property
    def name(self):
        return self._name

    @property
    def rdim(self):
        return self._rdim

    @property
    def domain(self):
        return self._domain

    @property
    def coordinates(self):
        if self.rdim == 1:
            return self._coordinates[0]
        else:
            return self._coordinates

    @property
    def jacobian(self):
        if self._jacobian is None:
            self._compute_jacobian()

        return self._jacobian

    @property
    def det_jacobian(self):
        if self._det_jacobian is None:
            self._compute_det_jacobian()

        return self._det_jacobian

    @property
    def covariant(self):
        if self._covariant is None:
            self._compute_covariant()

        return self._covariant

    @property
    def contravariant(self):
        if self._contravariant is None:
            self._compute_contravariant()

        return self._contravariant

    @property
    def hessian(self):
        if self._hessian is None:
            self._compute_hessian()

        return self._hessian

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

    def _compute_jacobian(self):
        M = Matrix(Jacobian(self))
        self._jacobian = M

    def _compute_det_jacobian(self):
        J = self.jacobian

        dim = self.rdim
        if dim == 2:
            det = J[0,0]* J[1,1] - J[0,1]* J[1,0]

        elif dim == 3:
            det = (J[0, 0]*J[1, 1]*J[2, 2] -
                   J[0, 0]*J[1, 2]*J[2, 1] -
                   J[0, 1]*J[1, 0]*J[2, 2] +
                   J[0, 1]*J[1, 2]*J[2, 0] +
                   J[0, 2]*J[1, 0]*J[2, 1] -
                   J[0, 2]*J[1, 1]*J[2, 0])

        else:
            det = J.det()

        self._det_jacobian = det

    def _compute_covariant(self):

        J = self.jacobian
        dim = self.rdim
        if dim == 1:
            M = 1/J[0,0]

        elif dim == 2:
            det = J[0,0]* J[1,1] - J[0,1]* J[1,0]
            J_inv = Matrix([[J[1,1], -J[0,1]], [-J[1,0], J[0,0]]])
            M = J_inv.transpose() / det

        elif dim == 3:
            det = (J[0, 0]*J[1, 1]*J[2, 2] -
                   J[0, 0]*J[1, 2]*J[2, 1] -
                   J[0, 1]*J[1, 0]*J[2, 2] +
                   J[0, 1]*J[1, 2]*J[2, 0] +
                   J[0, 2]*J[1, 0]*J[2, 1] -
                   J[0, 2]*J[1, 1]*J[2, 0])

            M = Matrix([[J[1, 1]*J[2, 2] - J[1, 2]*J[2, 1], J[0, 2]*J[2, 1] - J[0, 1]*J[2, 2], J[0, 1]*J[1, 2] - J[0, 2]*J[1, 1]],
                        [J[1, 2]*J[2, 0] - J[1, 0]*J[2, 2], J[0, 0]*J[2, 2] - J[0, 2]*J[2, 0], J[0, 2]*J[1, 0] - J[0, 0]*J[1, 2]],
                        [J[1, 0]*J[2, 1] - J[1, 1]*J[2, 0], J[0, 1]*J[2, 0] - J[0, 0]*J[2, 1], J[0, 0]*J[1, 1] - J[0, 1]*J[1, 0]]])

            M = 1/det * M.transpose()

        else:
            M = J.inv().transpose()

        self._covariant = M

    def _compute_contravariant(self):
        J = self.jacobian
        j = self.det_jacobian
        inv_j = 1/j
        n_rows, n_cols = J.shape
        for i_row in range(0, n_rows):
            for i_col in range(0, n_cols):
                J[i_row, i_col] *= inv_j

        self._contravariant = J

    def _compute_hessian(self):
        raise NotImplementedError('TODO')
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

        return F.det_jacobian

class Covariant(MappingApplication):
    """

    Examples

    """

    @classmethod
    def eval(cls, F, v):

        if not isinstance(v, (tuple, list, Tuple, Matrix)):
            raise TypeError('> Expecting a tuple, list, Tuple, Matrix')

        M = F.covariant
#        v = Matrix(v)
        dim = F.rdim
        if dim == 1:
            v = [M * v[0]]
            return Tuple(*v)

        else:
#            v = M*v
#            return Tuple(*v)

            n,m = M.shape
            w = []
            for i in range(0, n):
                w.append(0)

            for i in range(0, n):
                for j in range(0, m):
                    w[i] += M[i,j] * v[j]

            return Tuple(*w)

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

        M = F.contravariant
        v = Matrix(v)
        v = M*v
        return Tuple(*v)
