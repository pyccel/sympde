# coding: utf-8
from abc import ABC, abstractmethod
from sympy                 import Indexed, IndexedBase, Idx
from sympy                 import Matrix, ImmutableDenseMatrix
from sympy                 import Function, Expr
from sympy                 import sympify
from sympy                 import cacheit
from sympy.core            import Basic
from sympy.core            import Symbol,Integer
from sympy.core            import Add, Mul, Pow
from sympy.core.numbers    import ImaginaryUnit
from sympy.core.containers import Tuple
from sympy                 import S
from sympy                 import sqrt, symbols
from sympy.core.exprtools  import factor_terms
from sympy.polys.polytools import parallel_poly_from_expr

from sympde.core              import constant
from sympde.core.basic        import BasicMapping
from sympde.core.basic        import CalculusFunction
from sympde.core.basic        import _coeffs_registery
from sympde.calculus.core     import PlusInterfaceOperator, MinusInterfaceOperator
from sympde.calculus.core     import grad, div, curl
from sympde.calculus.core     import dot, inner, outer, _diff_ops
from sympde.calculus.core     import has, DiffOperator
from sympde.calculus.matrices import MatrixSymbolicExpr, MatrixElement, SymbolicTrace, Inverse
from sympde.calculus.matrices import SymbolicDeterminant, Transpose

from .basic       import BasicDomain, Union, InteriorDomain
from .basic       import Boundary, Connectivity, Interface
from .domain      import Domain, NCubeInterior
from .domain      import NormalVector
from .space       import ScalarFunction, VectorFunction, IndexedVectorFunction
from .space       import Trace
from .datatype    import HcurlSpaceType, H1SpaceType, L2SpaceType, HdivSpaceType, UndefinedSpaceType
from .derivatives import dx, dy, dz, DifferentialOperator
from .derivatives import _partial_derivatives
from .derivatives import get_atom_derivatives, get_index_derivatives_atom
from .derivatives import _logical_partial_derivatives
from .derivatives import get_atom_logical_derivatives, get_index_logical_derivatives_atom
from .derivatives import LogicalGrad_1d, LogicalGrad_2d, LogicalGrad_3d

# TODO fix circular dependency between sympde.topology.domain and sympde.topology.mapping
# TODO fix circular dependency between sympde.expr.evaluation and sympde.topology.mapping

__all__ = (
    'BasicCallableMapping',
    'Contravariant',
    'Covariant',
    'InterfaceMapping',
    'InverseMapping',
    'Jacobian',
    'JacobianInverseSymbol',
    'JacobianSymbol',
    'MappedDomain',
    'MappingApplication',
    'MultiPatchMapping',
    'PullBack',
    'SymbolicWeightedVolume',
    'get_logical_test_function',
)

#==============================================================================
@cacheit
def cancel(f):
    try:
        f           = factor_terms(f, radical=True)
        p, q        = f.as_numer_denom()
        # TODO accelerate parallel_poly_from_expr
        (p, q), opt = parallel_poly_from_expr((p,q))
        c, P, Q     = p.cancel(q)
        return c*(P.as_expr()/Q.as_expr())
    except:
        return f

def get_logical_test_function(u):
    space           = u.space
    kind            = space.kind
    dim             = space.ldim
    logical_domain  = space.domain.logical_domain
    l_space         = type(space)(space.name, logical_domain, kind=kind)
    el              = l_space.element(u.name)
    return el

#==============================================================================
class BasicCallableMapping(ABC):
    """
    Transformation of coordinates, which can be evaluated.

    F: R^l -> R^p
    F(eta) = x

    with l <= p
    """
    @abstractmethod
    def __call__(self, *eta):
        """ Evaluate mapping at location eta. """

    @abstractmethod
    def jacobian(self, *eta):
        """ Compute Jacobian matrix at location eta. """

    @abstractmethod
    def jacobian_inv(self, *eta):
        """ Compute inverse Jacobian matrix at location eta.
            An exception should be raised if the matrix is singular.
        """

    @abstractmethod
    def metric(self, *eta):
        """ Compute components of metric tensor at location eta. """

    @abstractmethod
    def metric_det(self, *eta):
        """ Compute determinant of metric tensor at location eta. """

    @property
    @abstractmethod
    def ldim(self):
        """ Number of logical/parametric dimensions in mapping
            (= number of eta components).
        """

    @property
    @abstractmethod
    def pdim(self):
        """ Number of physical dimensions in mapping
            (= number of x components).
        """

#==============================================================================
class Mapping(BasicMapping):
    """
    Represents a Mapping object.

    Examples

    """
    _expressions  = None # used for analytical mapping
    _jac          = None
    _inv_jac      = None
    _constants    = None
    _callable_map = None
    _ldim         = None
    _pdim         = None

    def __new__(cls, name, dim=None, **kwargs):

        ldim        = kwargs.pop('ldim', cls._ldim)
        pdim        = kwargs.pop('pdim', cls._pdim)
        coordinates = kwargs.pop('coordinates', None)
        evaluate    = kwargs.pop('evaluate', True)

        dims = [dim, ldim, pdim]
        for i,d in enumerate(dims):
            if isinstance(d, (tuple, list, Tuple, Matrix, ImmutableDenseMatrix)):
                if not len(d) == 1:
                    raise ValueError('> Expecting a tuple, list, Tuple of length 1')
                dims[i] = d[0]

        dim, ldim, pdim = dims

        if dim is None:
            assert ldim is not None
            assert pdim is not None
            assert pdim >= ldim
        else:
            ldim = dim
            pdim = dim


        obj = IndexedBase.__new__(cls, name, shape=pdim)

        if not evaluate:
            return obj

        if coordinates is None:
            _coordinates = [Symbol(name) for name in ['x', 'y', 'z'][:pdim]]
        else:
            if not isinstance(coordinates, (list, tuple, Tuple)):
                raise TypeError('> Expecting list, tuple, Tuple')

            for a in coordinates:
                if not isinstance(a, (str, Symbol)):
                    raise TypeError('> Expecting str or Symbol')

            _coordinates = [Symbol(u) for u in coordinates]

        obj._name                = name
        obj._ldim                = ldim
        obj._pdim                = pdim
        obj._coordinates         = tuple(_coordinates)
        obj._jacobian            = kwargs.pop('jacobian', JacobianSymbol(obj))
        obj._is_minus            = None
        obj._is_plus             = None

        lcoords = ['x1', 'x2', 'x3'][:ldim]
        lcoords = [Symbol(i) for i in lcoords]
        obj._logical_coordinates = Tuple(*lcoords)
        # ...
        if not( obj._expressions is None ):
            coords = ['x', 'y', 'z'][:pdim]

            # ...
            args = []
            for i in coords:
                x = obj._expressions[i]
                x = sympify(x)
                args.append(x)

            args = Tuple(*args)
            # ...
            zero_coords = ['x1', 'x2', 'x3'][ldim:]

            for i in zero_coords:
                x = sympify(i)
                args = args.subs(x,0)
            # ...

            constants        = list(set(args.free_symbols) - set(lcoords))
            constants_values = {a.name:constant(a.name, dtype=float) for a in constants}
            # subs constants as sympde constant objects instead of Symbol
            constants_values.update( kwargs )
            d = {a:constants_values[a.name] for a in constants}
            args = args.subs(d)

            obj._expressions = args
            obj._constants   = tuple(a for a in constants if isinstance(constants_values[a.name], Symbol))

            args  = [obj[i] for i in range(pdim)]
            exprs = obj._expressions
            subs  = list(zip(_coordinates, exprs))

            if obj._jac is None and obj._inv_jac is None:
                obj._jac     = Jacobian(obj).subs(list(zip(args, exprs)))
                obj._inv_jac = obj._jac.inv() if pdim == ldim else None
            elif obj._inv_jac is None:
                obj._jac     = ImmutableDenseMatrix(sympify(obj._jac)).subs(subs)
                obj._inv_jac = obj._jac.inv() if pdim == ldim else None

            elif obj._jac is None:
                obj._inv_jac = ImmutableDenseMatrix(sympify(obj._inv_jac)).subs(subs)
                obj._jac     = obj._inv_jac.inv()
            else:
                obj._jac     = ImmutableDenseMatrix(sympify(obj._jac)).subs(subs)
                obj._inv_jac = ImmutableDenseMatrix(sympify(obj._inv_jac)).subs(subs)

        else:
            obj._jac     = Jacobian(obj)

        obj._metric     = obj._jac.T*obj._jac
        obj._metric_det = obj._metric.det()

        return obj

    #--------------------------------------------------------------------------
    # Callable mapping
    #--------------------------------------------------------------------------
    def get_callable_mapping(self):
        if self._callable_map is None:
            if self._expressions is None:
                msg = 'Cannot generate callable mapping without analytical expressions. '\
                      'A user-defined callable mapping of type `BasicCallableMapping` '\
                      'can be provided using the method `set_callable_mapping`.'
                raise ValueError(msg)

            from sympde.topology.callable_mapping import CallableMapping
            self._callable_map = CallableMapping(self)

        return self._callable_map

    def set_callable_mapping(self, F):

        if not isinstance(F, BasicCallableMapping):
            raise TypeError(
                f'F must be a BasicCallableMapping, got {type(F)} instead')

        self._callable_map = F

    #--------------------------------------------------------------------------
    @property
    def name( self ):
        return self._name

    @property
    def ldim( self ):
        return self._ldim

    @property
    def pdim( self ):
        return self._pdim

    @property
    def coordinates( self ):
        if self.pdim == 1:
            return self._coordinates[0]
        else:
            return self._coordinates

    @property
    def logical_coordinates( self ):
        if self.ldim == 1:
            return self._logical_coordinates[0]
        else:
            return self._logical_coordinates

    def __call__(self, domain):
        assert(isinstance(domain, BasicDomain))
        return MappedDomain(self, domain)

    @property
    def jacobian( self ):
        return self._jacobian

    @property
    def det_jacobian( self ):
        return self.jacobian.det()

    @property
    def is_analytical( self ):
        return not( self._expressions is None )

    @property
    def expressions( self ):
        return self._expressions

    @property
    def jacobian_expr( self ):
        return self._jac

    @property
    def jacobian_inv_expr( self ):
        if not self.is_analytical and self._inv_jac is None:
            self._inv_jac = self.jacobian_expr.inv()
        return self._inv_jac

    @property
    def metric_expr( self ):
        return self._metric

    @property
    def metric_det_expr( self ):
        return self._metric_det

    @property
    def constants( self ):
        return self._constants

    @property
    def is_minus( self ):
        return self._is_minus

    @property
    def is_plus( self ):
        return self._is_plus

    def set_plus_minus( self, **kwargs):
        minus = kwargs.pop('minus', False)
        plus  = kwargs.pop('plus', False)
        assert plus is not minus

        self._is_plus  = plus
        self._is_minus = minus

    def copy(self):
        obj = Mapping(self.name,
                     ldim=self.ldim,
                     pdim=self.pdim,
                     evaluate=False)

        obj._name                = self.name
        obj._ldim                = self.ldim
        obj._pdim                = self.pdim
        obj._coordinates         = self.coordinates
        obj._jacobian            = JacobianSymbol(obj)
        obj._logical_coordinates = self.logical_coordinates
        obj._expressions         = self._expressions
        obj._constants           = self._constants
        obj._jac                 = self._jac
        obj._inv_jac             = self._inv_jac
        obj._metric              = self._metric
        obj._metric_det          = self._metric_det
        obj.__callable_map       = self._callable_map
        obj._is_plus             = self._is_plus
        obj._is_minus            = self._is_minus
        return obj

    def _hashable_content(self):
        args = (self.name, self.ldim, self.pdim, self._coordinates, self._logical_coordinates,
                self._expressions, self._constants, self._is_plus, self._is_minus)
        return tuple([a for a in args if a is not None])

    def _eval_subs(self, old, new):
        return self

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

#==============================================================================
class InverseMapping(Mapping):
    def __new__(cls, mapping):
        assert isinstance(mapping, Mapping)
        name     = mapping.name
        ldim     = mapping.ldim
        pdim     = mapping.pdim
        coords   = mapping.logical_coordinates
        jacobian = mapping.jacobian.inv()
        return Mapping.__new__(cls, name, ldim=ldim, pdim=pdim, coordinates=coords, jacobian=jacobian)

#==============================================================================
class JacobianSymbol(MatrixSymbolicExpr):
    _axis = None
    def __new__(cls, mapping, axis=None):
        assert isinstance(mapping, Mapping)
        if axis is not None:
            assert isinstance(axis, (int, Integer))
        obj = MatrixSymbolicExpr.__new__(cls, mapping)
        obj._axis = axis
        return obj

    @property
    def mapping(self):
        return self._args[0]

    @property
    def axis(self):
        return self._axis

    def inv(self):
        return JacobianInverseSymbol(self.mapping, self.axis)

    def _hashable_content(self):
        if self.axis is not None:
            return (type(self).__name__, self.mapping, self.axis)
        else:
            return (type(self).__name__, self.mapping)

    def __hash__(self):
        return hash(self._hashable_content())

    def _eval_subs(self, old, new):
        if isinstance(new, Mapping):
            if self.axis is not None:
                obj = JacobianSymbol(new, self.axis)
            else:
                obj = JacobianSymbol(new)
            return obj
        return self
    def _sympystr(self, printer):
        sstr = printer.doprint
        if self.axis:
            return 'Jacobian({},{})'.format(sstr(self.mapping.name), self.axis)
        else:
            return 'Jacobian({})'.format(sstr(self.mapping.name))

#==============================================================================
class JacobianInverseSymbol(MatrixSymbolicExpr):
    _axis = None
    is_Matrix     = False
    def __new__(cls, mapping, axis=None):
        assert isinstance(mapping, Mapping)
        if axis is not None:
            assert isinstance(axis, int)
        obj = MatrixSymbolicExpr.__new__(cls, mapping)
        obj._axis = axis
        return obj

    @property
    def mapping(self):
        return self._args[0]

    @property
    def axis(self):
        return self._axis

    def _hashable_content(self):
        if self.axis is not None:
            return (type(self).__name__, self.mapping, self.axis)
        else:
            return (type(self).__name__, self.mapping)

    def __hash__(self):
        return hash(self._hashable_content())

    def _sympystr(self, printer):
        sstr = printer.doprint
        if self.axis:
            return 'Jacobian({},{})**(-1)'.format(sstr(self.mapping.name), self.axis)
        else:
            return 'Jacobian({})**(-1)'.format(sstr(self.mapping.name))

#==============================================================================
class InterfaceMapping(Mapping):
    """
    InterfaceMapping is used to represent a mapping in the interface.

    Attributes
    ----------
    minus : Mapping
        the mapping on the negative direction of the interface
    plus  : Mapping
        the mapping on the positive direction of the interface
    """

    def __new__(cls, minus, plus):
        assert isinstance(minus, Mapping)
        assert isinstance(plus,  Mapping)
        minus = minus.copy()
        plus  = plus.copy()

        minus.set_plus_minus(minus=True)
        plus.set_plus_minus(plus=True)

        name = '{}|{}'.format(str(minus.name), str(plus.name))
        obj  = Mapping.__new__(cls, name, ldim=minus.ldim, pdim=minus.pdim)
        obj._minus = minus
        obj._plus  = plus
        return obj

    @property
    def minus(self):
        return self._minus

    @property
    def plus(self):
        return self._plus

    @property
    def is_analytical(self):
        return self.minus.is_analytical and self.plus.is_analytical

    def _eval_subs(self, old, new):
        minus = self.minus.subs(old, new)
        plus  = self.plus.subs(old, new)
        return InterfaceMapping(minus, plus)

    def _eval_simplify(self, **kwargs):
        return self

#==============================================================================
class MultiPatchMapping(Mapping):

    def __new__(cls, dic):
        assert isinstance( dic, dict)
        return Basic.__new__(cls, dic)

    @property
    def mappings(self):
        return self.args[0]

    @property
    def is_analytical(self):
        return all(a.is_analytical for a in self.mappings.values())

    @property
    def ldim(self):
        return list(self.mappings.values())[0].ldim

    @property
    def pdim(self):
        return list(self.mappings.values())[0].pdim

    @property
    def is_analytical(self):
        return all(e.is_analytical for e in self.mappings.values())

    def _eval_subs(self, old, new):
        return self

    def _eval_simplify(self, **kwargs):
        return self

    def __hash__(self):
        return hash((*self.mappings.values(), *self.mappings.keys()))

    def _sympystr(self, printer):
        sstr = printer.doprint
        mappings = (sstr(i) for i in self.mappings.values())
        return 'MultiPatchMapping({})'.format(', '.join(mappings))

#==============================================================================
class MappedDomain(BasicDomain):
    """."""

    @cacheit
    def __new__(cls, mapping, logical_domain):
        assert(isinstance(mapping, Mapping))
        assert(isinstance(logical_domain, BasicDomain))
        if isinstance(logical_domain, Domain):
            kwargs = dict(
            dim            = logical_domain._dim,
            mapping        = mapping,
            logical_domain = logical_domain)
            boundaries     = logical_domain.boundary
            interiors      = logical_domain.interior

            if isinstance(interiors, Union):
                kwargs['interiors'] = Union(*[mapping(a) for a in interiors.args])
            else:
                kwargs['interiors'] = mapping(interiors)

            if isinstance(boundaries, Union):
                kwargs['boundaries'] = [mapping(a) for a in boundaries.args]
            elif boundaries:
                kwargs['boundaries'] = mapping(boundaries)

            interfaces =  logical_domain.connectivity.interfaces
            if interfaces:
                if isinstance(interfaces, Union):
                    interfaces = interfaces.args
                else:
                    interfaces = [interfaces]
                connectivity = {}
                for e in interfaces:
                    connectivity[e.name] = Interface(e.name, mapping(e.minus), mapping(e.plus))
                kwargs['connectivity'] = Connectivity(connectivity)

            name = '{}({})'.format(str(mapping.name), str(logical_domain.name))
            return Domain(name, **kwargs)

        elif isinstance(logical_domain, NCubeInterior):
            name  = logical_domain.name
            dim   = logical_domain.dim
            dtype = logical_domain.dtype
            min_coords = logical_domain.min_coords
            max_coords = logical_domain.max_coords
            name = '{}({})'.format(str(mapping.name), str(name))
            return NCubeInterior(name, dim, dtype, min_coords, max_coords, mapping, logical_domain)
        elif isinstance(logical_domain, InteriorDomain):
            name  = logical_domain.name
            dim   = logical_domain.dim
            dtype = logical_domain.dtype
            name = '{}({})'.format(str(mapping.name), str(name))
            return InteriorDomain(name, dim, dtype, mapping, logical_domain)
        elif isinstance(logical_domain, Boundary):
            name   = logical_domain.name
            axis   = logical_domain.axis
            ext    = logical_domain.ext
            domain = mapping(logical_domain.domain)
            return Boundary(name, domain, axis, ext, mapping, logical_domain)
        else:
            raise NotImplementedError('TODO')
#==============================================================================
class SymbolicWeightedVolume(Expr):
    """
    This class represents the symbolic weighted volume of a quadrature rule
    """
#TODO move this somewhere else
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

class PullBack(Expr):
    is_commutative = False

    def __new__(cls, u, mapping=None):
        if not isinstance(u, (VectorFunction, ScalarFunction)):
            raise TypeError('{} must be of type ScalarFunction or VectorFunction'.format(str(u)))

        if u.space.domain.mapping is None:
            raise ValueError('The pull-back can be performed only to mapped domains')

        space = u.space
        kind  = space.kind
        dim   = space.ldim
        el    = get_logical_test_function(u)

        if space.is_broken:
            assert mapping is not None
        else:
            mapping = space.domain.mapping

        J = mapping.jacobian
        if isinstance(kind, (UndefinedSpaceType, H1SpaceType)):
            expr = el

        elif isinstance(kind, HcurlSpaceType):
            expr = J.inv().T * el

        elif isinstance(kind, HdivSpaceType):
            expr = (J/J.det()) * el

        elif isinstance(kind, L2SpaceType):
            expr = el/J.det()

#        elif isinstance(kind, UndefinedSpaceType):
#            raise ValueError('kind must be specified in order to perform the pull-back transformation')
        else:
            raise ValueError("Unrecognized kind '{}' of space {}".format(kind, str(u.space)))

        obj       = Expr.__new__(cls, u)
        obj._expr = expr
        obj._kind = kind
        obj._test = el
        return obj

    @property
    def expr(self):
        return self._expr

    @property
    def kind(self):
        return self._kind

    @property
    def test(self):
        return self._test

#==============================================================================
class Jacobian(MappingApplication):
    r"""
    This class calculates the Jacobian of a mapping F
    where [J_{F}]_{i,j} =  \frac{\partial F_{i}}{\partial x_{j}}
    or simply J_{F} =  (\nabla F)^T

    """

    @classmethod
    def eval(cls, F):
        """
        this class methods computes the jacobian of a mapping

        Parameters:
        ----------
         F: Mapping
            mapping object

        Returns:
        ----------
         expr : ImmutableDenseMatrix
            the jacobian matrix
        """

        if not isinstance(F, Mapping):
            raise TypeError('> Expecting a Mapping object')

        if F.jacobian_expr is not None:
            return F.jacobian_expr

        pdim = F.pdim
        ldim = F.ldim

        F = [F[i] for i in range(0, F.pdim)]
        F = Tuple(*F)

        if ldim == 1:
            expr = LogicalGrad_1d(F)

        elif ldim == 2:
            expr = LogicalGrad_2d(F)

        elif ldim == 3:
            expr = LogicalGrad_3d(F)

        return expr.T

#==============================================================================
class Covariant(MappingApplication):
    """

    Examples

    """

    @classmethod
    def eval(cls, F, v):

        """
        This class methods computes the covariant transformation

        Parameters:
        ----------
         F: Mapping
            mapping object

         v: <tuple|list|Tuple|ImmutableDenseMatrix|Matrix>
            the basis function

        Returns:
        ----------
         expr : Tuple
            the covariant transformation
        """

        if not isinstance(v, (tuple, list, Tuple, ImmutableDenseMatrix, Matrix)):
            raise TypeError('> Expecting a tuple, list, Tuple, Matrix')

        assert F.pdim == F.ldim

        M   = Jacobian(F).inv().T
        dim = F.pdim

        if dim == 1:
            b = M[0,0] * v[0]
            return Tuple(b)
        else:
            n,m = M.shape
            w   = []
            for i in range(0, n):
                w.append(S.Zero)

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
        """
        This class methods computes the contravariant transformation

        Parameters:
        ----------
         F: Mapping
            mapping object

         v: <tuple|list|Tuple|ImmutableDenseMatrix|Matrix>
            the basis function

        Returns:
        ----------
         expr : Tuple
            the contravariant transformation
        """

        if not isinstance(F, Mapping):
            raise TypeError('> Expecting a Mapping')

        if not isinstance(v, (tuple, list, Tuple, ImmutableDenseMatrix, Matrix)):
            raise TypeError('> Expecting a tuple, list, Tuple, Matrix')

        M = Jacobian(F)
        M = M/M.det()
        v = Matrix(v)
        v = M*v
        return Tuple(*v)
