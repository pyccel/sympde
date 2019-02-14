# coding: utf-8

from collections  import OrderedDict

from sympy.core import Basic
from sympy.tensor import IndexedBase, Indexed
from sympy.core import Symbol
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Matrix
from sympy.core import Add, Mul, Pow
from sympy.core.singleton import S
from sympy.core.expr import AtomicExpr

from sympde.core.basic import BasicMapping
from sympde.core.algebra import (Dot_1d,
                                 Dot_2d, Inner_2d, Cross_2d,
                                 Dot_3d, Inner_3d, Cross_3d)

from sympde.core.basic import CalculusFunction
from sympde.core.basic import _coeffs_registery
from .basic import BasicDomain
from .domain import Domain
from .derivatives import dx, dy, dz
from .derivatives import dx1, dx2, dx3
from .derivatives import (Grad_1d, Div_1d,
                          Grad_2d, Curl_2d, Rot_2d, Div_2d,
                          Grad_3d, Curl_3d, Div_3d)
from .derivatives import LogicalGrad_1d, LogicalGrad_2d, LogicalGrad_3d
from .derivatives import _logical_partial_derivatives
from .derivatives import get_atom_logical_derivatives, get_index_logical_derivatives_atom

from .space import ScalarTestFunction, VectorTestFunction, IndexedTestTrial
from .space import ScalarField, VectorField, IndexedVectorField

#==============================================================================
class Mapping(BasicMapping):
    """
    Represents a Mapping object.

    Examples

    """
    # TODO shall we keep rdim ?
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
    def coordinates(self):
        if self.rdim == 1:
            return self._coordinates[0]
        else:
            return self._coordinates

    def __call__(self, domain):
        assert(isinstance(domain, BasicDomain))

        return MappedDomain(self, domain)

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


#==============================================================================
class MappedDomain(BasicDomain):
    def __new__(cls, mapping, domain):
        assert(isinstance(mapping, Mapping))
        assert(isinstance(domain, BasicDomain))

        return Basic.__new__(cls, mapping, domain)

    @property
    def mapping(self):
        return self._args[0]

    @property
    def domain(self):
        return self._args[1]

    @property
    def dim(self):
        return self.domain.dim

#==============================================================================
class SymbolicMappingExpr(AtomicExpr):

    def __new__(cls, mapping):
        assert(isinstance(mapping, Mapping))

        return Basic.__new__(cls, mapping)

    @property
    def mapping(self):
        return self._args[0]

class SymbolicDeterminant(SymbolicMappingExpr):

    def _sympystr(self, printer):
        sstr = printer.doprint
        mapping = sstr(self.mapping)
        return 'det({})'.format(mapping)

class SymbolicCovariant(SymbolicMappingExpr):

    def _sympystr(self, printer):
        sstr = printer.doprint
        mapping = sstr(self.mapping)
        return 'covariant({})'.format(mapping)

class SymbolicContravariant(SymbolicMappingExpr):

    def _sympystr(self, printer):
        sstr = printer.doprint
        mapping = sstr(self.mapping)
        return 'contravariant({})'.format(mapping)

#==============================================================================
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

#==============================================================================
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
            expr = LogicalGrad_1d(F)

        elif rdim == 2:
            expr = LogicalGrad_2d(F)

        elif rdim == 3:
            expr = LogicalGrad_3d(F)

        return expr

#==============================================================================
class DetJacobian(MappingApplication):
    """

    Examples

    """

    @classmethod
    def eval(cls, F):

        if not isinstance(F, Mapping):
            raise TypeError('> Expecting a Mapping object')

        return F.det_jacobian

#==============================================================================
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

#==============================================================================
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


#==============================================================================
class LogicalExpr(CalculusFunction):

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

    @classmethod
    def eval(cls, *_args):
        """."""

        if not _args:
            return

        if not len(_args) == 2:
            raise ValueError('Expecting two arguments')

        M    = _args[0]
        expr = _args[1]
        dim  = M.rdim # TODO this is not the dim of the domain

        if isinstance(expr, Add):
            args = [cls.eval(M, a) for a in expr.args]
            return Add(*args)

        elif isinstance(expr, Mul):
            args = [cls.eval(M, a) for a in expr.args]
            return Mul(*args)

        elif isinstance(expr, _coeffs_registery):
            return expr

        elif isinstance(expr, _logical_partial_derivatives):
            return expr

        elif isinstance(expr, (ScalarField, ScalarTestFunction, IndexedTestTrial, IndexedVectorField)):
            return expr

        elif isinstance(expr, (VectorField, VectorTestFunction)):
            raise NotImplementedError('')

        elif isinstance(expr, dx):
            arg = expr.args[0]
            arg = cls.eval(M, arg)

            # ...
            if dim == 1:
                lgrad_arg = LogicalGrad_1d(arg)

                if not isinstance(lgrad_arg, (list, tuple, Tuple, Matrix)):
                    lgrad_arg = Tuple(lgrad_arg)

            elif dim == 2:
                lgrad_arg = LogicalGrad_2d(arg)

            elif dim == 3:
                lgrad_arg = LogicalGrad_3d(arg)

            grad_arg = Covariant(M, lgrad_arg)
            # ...

            return grad_arg[0]

        elif isinstance(expr, dy):
            arg = expr.args[0]
            arg = cls.eval(M, arg)

            # ...
            if dim == 1:
                lgrad_arg = LogicalGrad_1d(arg)

            elif dim == 2:
                lgrad_arg = LogicalGrad_2d(arg)

            elif dim == 3:
                lgrad_arg = LogicalGrad_3d(arg)

            grad_arg = Covariant(M, lgrad_arg)
            # ...

            return grad_arg[1]

        elif isinstance(expr, dz):
            arg = expr.args[0]
            arg = cls.eval(M, arg)

            # ...
            if dim == 1:
                lgrad_arg = LogicalGrad_1d(arg)

            elif dim == 2:
                lgrad_arg = LogicalGrad_2d(arg)

            elif dim == 3:
                lgrad_arg = LogicalGrad_3d(arg)

            grad_arg = Covariant(M, lgrad_arg)
            # ...

            return grad_arg[2]

        return cls(M, expr, evaluate=False)

#==============================================================================
class SymbolicExpr(CalculusFunction):
    """returns a sympy expression where partial derivatives are converted into
    sympy Symbols."""

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

    @classmethod
    def eval(cls, *_args, **kwargs):
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        code = kwargs.pop('code', None)

        if isinstance(expr, Add):
            args = [cls.eval(a, code=code) for a in expr.args]
            return Add(*args)

        elif isinstance(expr, Mul):
            args = [cls.eval(a, code=code) for a in expr.args]
            return Mul(*args)

        elif isinstance(expr, Pow):
            b = expr.base
            e = expr.exp
            return Pow(cls(b, code=code), e)

        elif isinstance(expr, _coeffs_registery):
            return expr

        elif isinstance(expr, (list, tuple, Tuple)):
            expr = [cls.eval(a, code=code) for a in expr]
            return Tuple(*expr)

        elif isinstance(expr, Matrix):

            lines = []
            n_row,n_col = expr.shape
            for i_row in range(0,n_row):
                line = []
                for i_col in range(0,n_col):
                    line.append(cls.eval(expr[i_row, i_col], code=code))

                lines.append(line)

            return Matrix(lines)

        elif isinstance(expr, (ScalarField, ScalarTestFunction)):
            if code:
                name = '{name}_{code}'.format(name=expr.name, code=code)
            else:
                name = str(expr.name)

            return Symbol(name)

        elif isinstance(expr, (VectorField, VectorTestFunction)):
            raise NotImplementedError('')

        elif isinstance(expr, Indexed):
            base = expr.base
            if isinstance(base, Mapping):
                if expr.indices[0] == 0:
                    name = 'x'
                elif expr.indices[0] == 1:
                    name = 'y'
                elif expr.indices[0] == 2:
                    name = 'z'
                else:
                    raise ValueError('Wrong index')

            else:
                name =  '{base}_{i}'.format(base=base.name, i=expr.indices[0])

            if code:
                name = '{name}_{code}'.format(name=name, code=code)

            return Symbol(name)

        elif isinstance(expr, _logical_partial_derivatives):

            atom = get_atom_logical_derivatives(expr)
            indices = get_index_logical_derivatives_atom(expr, atom)

            code = None
            if indices:
                index = indices[0]
                code = ''
                index = OrderedDict(sorted(index.items()))

                for k,n in list(index.items()):
                    code += k*n

            return cls.eval(atom, code=code)

        return cls(expr, evaluate=False)

