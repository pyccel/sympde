# coding: utf-8
import warnings
from abc import ABC, ABCMeta, abstractmethod
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

from sympde.core              import Constant
from sympde.core.basic        import SymbolicMapping
from sympde.core.basic        import CalculusFunction
from sympde.core.basic        import _coeffs_registery
from sympde.calculus.core     import PlusInterfaceOperator, MinusInterfaceOperator
from sympde.calculus.core     import grad, div, curl, laplace #, hessian
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
    'AnalyticMapping',
    'BasicCallableMapping',
    'Contravariant',
    'Covariant',
    'DefinedMapping',
    'InterfaceMapping',
    'InverseMapping',
    'Jacobian',
    'JacobianInverseSymbol',
    'JacobianSymbol',
    'LogicalExpr',
    'MappedDomain',
    'Mapping',
    'MappingApplication',
    'MultiPatchMapping',
    'PullBack',
    'StructuralMapping',
    'SymbolicExpr',
    'SymbolicMapping',
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

import numpy as np

def numpy_to_native_python(a):
    if isinstance(a, np.generic):
        return a.item()
    return a

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
class _MappingABCMeta(ABCMeta, type(SymbolicMapping)):
    """
    Metaclass merging ``abc.ABCMeta`` with sympy's metaclass
    (``ManagedProperties`` in sympy 1.9).

    A mapping class derived from ``IndexedBase`` (via ``SymbolicMapping``)
    already carries sympy's metaclass; declaring ``@abstractmethod`` members on
    it additionally requires ``ABCMeta``. Python rejects a class whose metaclass
    is not a subclass of every base's metaclass, so the two are merged here once
    and reused by ``DefinedMapping``, ``StructuralMapping``, and
    ``AnalyticMapping``. ``type(...)`` is used instead of importing the name so
    this keeps working if a future sympy renames its metaclass.
    """

#==============================================================================
# ``SymbolicMapping`` is defined in ``sympde.core.basic`` (imported above) so
# that leaf modules can type-check against it without importing this one.

#==============================================================================
class DefinedMapping(SymbolicMapping, BasicCallableMapping, metaclass=_MappingABCMeta):
    """
    Abstract base class for *point-evaluable* mappings.

    F: R^l -> R^p ,  F(eta) = x ,  with l <= p

    A concrete subclass (``AnalyticMapping`` in sympde, ``SplineMapping`` in
    psydac) can be evaluated on logical coordinates -- single points or arrays
    of points -- and returns physical coordinates. The point-evaluation
    interface (``__call__``, ``jacobian``, ``jacobian_inv``, ``metric``,
    ``metric_det``, ``ldim``, ``pdim``) is inherited verbatim from
    ``BasicCallableMapping``; every one must be implemented for a subclass to be
    instantiable, which is what guarantees the sympde and psydac concrete
    mappings are interchangeable.
    """

#==============================================================================
class StructuralMapping(SymbolicMapping, metaclass=_MappingABCMeta):
    """
    Abstract base class for *symbolic, non-point-evaluable* mappings.

    These objects (``InverseMapping``, ``InterfaceMapping``,
    ``MultiPatchMapping``) describe how patches are related or assembled and
    only make sense symbolically: they stay callable on a domain (returning a
    symbolic mapped domain, like any ``SymbolicMapping``) but reject point
    evaluation.

    Since WP06d-4a a ``StructuralMapping`` no longer inherits ``Mapping`` -- it
    owns the small symbolic surface it needs (``name``, ``is_plus``/``is_minus``,
    coordinates, ``__call__``, hashing) directly here. Each concrete subclass
    sets ``_name`` / ``_ldim`` / ``_pdim`` (and ``_coordinates`` /
    ``_logical_coordinates`` where meaningful) in its own ``__new__`` and
    implements ``ldim`` / ``pdim``.
    """

    _is_minus = None
    _is_plus  = None

    def __call__(self, *args, domain=None):
        """
        Call this structural mapping on a domain (positional or as the
        ``domain`` keyword, matching ``SymbolicMapping``'s signature) to get a
        symbolic mapped domain.

        Parameters
        ----------
        domain : BasicDomain
            The logical domain to map.

        Returns
        -------
        MappedDomain

        Raises
        ------
        TypeError
            If called with anything other than exactly one ``BasicDomain``
            argument -- a ``StructuralMapping`` is symbolic and not
            point-evaluable.
        """
        # Same domain-vs-point dispatch as AnalyticMapping.__call__; the
        # terminal action here is the domain call (Mapping.__call__'s body,
        # reproduced -- Mapping is no longer in the MRO to super() into).
        if domain is None and len(args) == 1 and isinstance(args[0], BasicDomain):
            domain = args[0]
        if domain is not None:
            assert domain.logical_domain is None
            assert domain.dim == self.ldim
            return MappedDomain(self, domain)
        raise TypeError(
            f"{type(self).__name__} is a StructuralMapping: it is symbolic "
            "and not point-evaluable.")

    @property
    def name(self):
        return self._name

    @property
    def is_minus(self):
        return self._is_minus

    @property
    def is_plus(self):
        return self._is_plus

    @property
    def coordinates(self):
        c = self._coordinates
        return c[0] if self.pdim == 1 else c

    @property
    def logical_coordinates(self):
        c = self._logical_coordinates
        return c[0] if self.ldim == 1 else c

    @property
    def jacobian_symbol(self):
        # lazy: multipatch pull-back (PullBack.__new__) asks for this on an
        # InterfaceMapping. InverseMapping overrides with a pre-inverted one.
        j = getattr(self, '_jacobian', None)
        if j is None:
            j = self._jacobian = JacobianSymbol(self)
        return j

    # A structural mapping carries no analytic `_expressions`; its symbolic
    # jacobian / metric are the same `Jacobian(self)` forms that Mapping.__new__
    # produced for it before WP06d-4a severed the base. Each concrete __new__
    # calls _init_symbolic_jacobian(obj) once (straight-line, like the old
    # Mapping.__new__ `else` branch -- a lazy property would recurse through
    # Jacobian.eval). psydac codegen (api/ast/fem.py) reads `jacobian_expr` on
    # multipatch mappings.
    _jac        = None
    _inv_jac    = None
    _metric     = None
    _metric_det = None

    @staticmethod
    def _init_symbolic_jacobian(obj):
        obj._jac        = Jacobian(obj)
        obj._metric     = obj._jac.T * obj._jac
        obj._metric_det = obj._metric.det()

    @property
    def jacobian_expr(self):
        return self._jac

    @property
    def jacobian_inv_expr(self):
        return self._inv_jac

    @property
    def metric_expr(self):
        return self._metric

    @property
    def metric_det_expr(self):
        return self._metric_det

    def _hashable_content(self):
        return (type(self).__name__, self._name, self.ldim, self.pdim,
                self._coordinates, self._logical_coordinates)

    @property
    @abstractmethod
    def ldim(self):
        """ Number of logical/parametric dimensions. """

    @property
    @abstractmethod
    def pdim(self):
        """ Number of physical dimensions. """

#==============================================================================
class Mapping(SymbolicMapping):
    """
    Represents a Mapping object.

    Now sits under ``SymbolicMapping`` (WP 06a) so that
    ``isinstance(_, SymbolicMapping)`` holds for every mapping in the hierarchy;
    ``SymbolicMapping`` is still a thin pass-through, so behaviour is unchanged.

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
            _coordinates = [Symbol(name, real=True) for name in ['x', 'y', 'z'][:pdim]]
        else:
            if not isinstance(coordinates, (list, tuple, Tuple)):
                raise TypeError('> Expecting list, tuple, Tuple')

            for a in coordinates:
                if not isinstance(a, (str, Symbol)):
                    raise TypeError('> Expecting str or Symbol')

            _coordinates = [Symbol(u, real=True) for u in coordinates]

        obj._name                = name
        obj._ldim                = ldim
        obj._pdim                = pdim
        obj._coordinates         = tuple(_coordinates)
        obj._jacobian            = kwargs.pop('jacobian', JacobianSymbol(obj))
        obj._is_minus            = None
        obj._is_plus             = None

        lcoords_names = ['x1', 'x2', 'x3'][:ldim]
        lcoords_symbols_real = [Symbol(i, real=True) for i in lcoords_names]
        lcoords_symbols_real_dict = {key: symbol for key,symbol in zip(lcoords_names, lcoords_symbols_real)}
        obj._logical_coordinates = Tuple(*lcoords_symbols_real) # coordinates must be symbols with real=True      
        lcoords_general_symbols = [Symbol(i) for i in lcoords_names] #list of symbols without real=True       

        # ...
        if not( obj._expressions is None ):
            coords_names = ['x', 'y', 'z'][:pdim]

            # ...
            args = []
            for i in coords_names:
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
            # get constants by subtracting coordinates from list of free symbols
            constants        = list(set(args.free_symbols) - set(lcoords_general_symbols))
            constants_values = {a.name:Constant(a.name) for a in constants}
            # subs constants as Constant objects instead of Symbol
            constants_values.update( kwargs )
            d = {a:numpy_to_native_python(constants_values[a.name]) for a in constants}
            args = args.subs(d)
            # subs coordinate symbols without real=True to symbols with real=True
            args = args.subs(lcoords_symbols_real_dict)

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
        # An explicitly attached callable (set_callable_mapping) always wins.
        if self._callable_map is not None:
            return self._callable_map

        if self._expressions is None:
            raise ValueError(
                'Cannot generate a callable mapping without analytical '
                'expressions. Attach a user-defined callable of type '
                '`BasicCallableMapping` with the method `set_callable_mapping`.')

        # Reachable only for a bare `Mapping` subclass that carries
        # `_expressions` but is not an `AnalyticMapping` (whose override returns
        # `self`). Since WP06c such a class is a mistake -- there is no longer a
        # `CallableMapping` wrapper to build.
        raise TypeError(
            f'{type(self).__name__} carries analytical expressions but is not '
            'an AnalyticMapping. Subclass AnalyticMapping (not Mapping) for a '
            'point-evaluable analytic mapping; an AnalyticMapping instance is '
            'its own callable mapping.')

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

    # Applying the mapping to a logical domain returns a mapped domain
    def __call__(self, domain):
        assert isinstance(domain, BasicDomain)
        assert domain.logical_domain is None
        assert domain.dim == self.ldim
        return MappedDomain(self, domain)

    @property
    def jacobian_symbol( self ):
        """Symbolic JacobianSymbol of this mapping."""
        return self._jacobian

    @property
    def jacobian( self ):
        # WP 02a: the name `jacobian` is being freed so that in 02b it can become
        # the numeric point-evaluation method of DefinedMapping. The symbolic
        # JacobianSymbol now lives on `jacobian_symbol`.
        warnings.warn(
            "Mapping.jacobian (symbolic) is deprecated; use "
            "Mapping.jacobian_symbol.",
            DeprecationWarning, stacklevel=2)
        return self.jacobian_symbol

    @property
    def det_jacobian( self ):
        return self.jacobian_symbol.det()

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
        # Use type(self), not Mapping, so a copy of an AnalyticMapping (or any
        # other Mapping subclass) keeps its concrete class -- and hence its
        # point-evaluation capability -- rather than being downgraded to a
        # plain, non-point-evaluable Mapping. Safe because evaluate=False
        # short-circuits Mapping.__new__ right after IndexedBase.__new__,
        # before any subclass-specific (_expressions) branching runs.
        obj = type(self)(self.name,
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
        obj._callable_map        = self._callable_map
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
class AnalyticMapping(Mapping, DefinedMapping, metaclass=_MappingABCMeta):
    """
    Analytic mapping: symbolic like ``Mapping`` (it carries ``_expressions``),
    and *directly point-evaluable* through the ``DefinedMapping`` interface.

    Point evaluation is done by the mapping itself: on first use it lambdifies
    its own analytic expressions into numpy callables (cached, one quantity at a
    time), so an ``AnalyticMapping`` and a psydac ``SplineMapping`` are
    interchangeable wherever the point-evaluation interface is expected, and
    ``get_callable_mapping()`` returns ``self``. A callable mapping attached
    explicitly with ``set_callable_mapping`` (e.g. a spline approximation) still
    wins: it is then used both by ``get_callable_mapping()`` and by every point
    call (``F(eta)``, ``F.jacobian(eta)``, ...).

    Examples
    --------
    >>> from sympde.topology import IdentityMapping
    >>> F = IdentityMapping('F', dim=2)
    >>> F(0.3, 0.4)                       # point evaluation -> physical coords
    (0.3, 0.4)
    >>> F.jacobian_symbol                 # symbolic JacobianSymbol (unchanged)
    Jacobian(F)
    """

    def _ensure_lambdified(self):
        """ Check this mapping can be point-evaluated and return the per-quantity
        lambdified-callable cache (created empty on first call; entries are
        filled in on demand by ``_lambdify``). """
        cache = getattr(self, '_lambdified', None)
        if cache is None:
            if self._expressions is None:
                raise ValueError('Cannot point-evaluate a mapping without '
                                 'analytical expressions.')
            if self._constants:
                raise ValueError(
                    f'{self.name} has unresolved symbolic constants '
                    f'{self._constants}; give them numeric values at '
                    'construction to point-evaluate it.')
            self._lambdified = cache = {}
        return cache

    def _lambdify(self, key):
        """ Lambdify one quantity -- ``'call'``, ``'jacobian'``,
        ``'jacobian_inv'``, ``'metric'`` or ``'metric_det'`` -- into a numpy
        callable, caching it on the instance. Done one quantity at a time so
        that e.g. ``F(eta)`` never builds ``jacobian_inv`` (which needs a
        square Jacobian -- surface mappings have ``ldim != pdim``). """
        cache = self._ensure_lambdified()
        if key not in cache:
            # Lazy import: sympde.utilities.utils imports sympde.topology,
            # so a module-level import here would be circular.
            from sympde.utilities.utils import lambdify_sympde
            v = self.logical_coordinates
            if key == 'call':
                cache[key] = tuple(lambdify_sympde(v, e) for e in self.expressions)
            else:
                expr = {'jacobian':     self.jacobian_expr,
                        'jacobian_inv': self.jacobian_inv_expr,
                        'metric':       self.metric_expr,
                        'metric_det':   self.metric_det_expr}[key]
                cache[key] = lambdify_sympde(v, expr)
        return cache[key]

    def _delegate_point_eval(self, name, *eta):
        """ Point-evaluate quantity ``name`` (``'call'`` | ``'jacobian'`` |
        ``'jacobian_inv'`` | ``'metric'`` | ``'metric_det'``) at ``eta``.

        A callable mapping attached with ``set_callable_mapping`` wins; otherwise
        this mapping's own lambdified analytic expressions are used. Keeps
        ``F(eta)`` / ``F.jacobian(eta)`` / ... consistent with
        ``get_callable_mapping()``. """
        cm = self.get_callable_mapping()
        if cm is not self:
            # cm is a SplineMapping or user-supplied BasicCallableMapping -- its
            # own methods, no recursion. It has no `.call`, so 'call' -> cm(*eta).
            return cm(*eta) if name == 'call' else getattr(cm, name)(*eta)
        if name == 'call':
            return tuple(f(*eta) for f in self._lambdify('call'))
        return self._lambdify(name)(*eta)

    def __call__(self, *args, domain=None):
        # A single BasicDomain (positional or `domain=`) -> symbolic
        # MappedDomain (unchanged Mapping behaviour). Anything else -> point
        # evaluation on logical coordinates. An unexpected keyword raises a
        # natural TypeError. Same dispatch shape as StructuralMapping.__call__.
        if domain is None and len(args) == 1 and isinstance(args[0], BasicDomain):
            domain, args = args[0], ()
        if domain is not None:
            return super().__call__(domain)
        return self._delegate_point_eval('call', *args)

    def jacobian(self, *eta):
        """ Jacobian matrix evaluated at the logical point(s) ``eta``. """
        return self._delegate_point_eval('jacobian', *eta)

    def jacobian_inv(self, *eta):
        """ Inverse Jacobian matrix evaluated at the logical point(s) ``eta``. """
        return self._delegate_point_eval('jacobian_inv', *eta)

    def metric(self, *eta):
        """ Metric tensor evaluated at the logical point(s) ``eta``. """
        return self._delegate_point_eval('metric', *eta)

    def metric_det(self, *eta):
        """ Determinant of the metric tensor at the logical point(s) ``eta``. """
        return self._delegate_point_eval('metric_det', *eta)

    def get_callable_mapping(self):
        # An AnalyticMapping *is* its own callable mapping. An explicitly
        # attached callable (set_callable_mapping) still wins.
        return self._callable_map if self._callable_map is not None else self

    # ldim / pdim: the concrete properties inherited from Mapping satisfy the
    # DefinedMapping / BasicCallableMapping abstract members.

#==============================================================================
class InverseMapping(StructuralMapping, metaclass=_MappingABCMeta):
    """ Symbolic inverse of a mapping: F^{-1}. Not point-evaluable. """

    def __new__(cls, mapping):
        assert isinstance(mapping, SymbolicMapping)
        # Build the IndexedBase directly (WP06d-4a: no longer routing through
        # Mapping.__new__, which wrapped `coordinates` Symbols in Symbol(...)
        # again -- the pre-existing `Symbol(Symbol(...))` TypeError).
        obj = IndexedBase.__new__(cls, mapping.name, shape=mapping.pdim)
        lcoords                  = mapping._logical_coordinates  # raw Tuple
        obj._name                = mapping.name
        obj._ldim                = mapping.ldim
        obj._pdim                = mapping.pdim
        obj._coordinates         = lcoords
        obj._logical_coordinates = lcoords
        obj._jacobian            = mapping.jacobian_symbol.inv()
        obj._is_minus            = None
        obj._is_plus             = None
        obj._base_mapping        = mapping
        cls._init_symbolic_jacobian(obj)
        return obj

    @property
    def ldim(self):
        return self._ldim

    @property
    def pdim(self):
        return self._pdim

    @property
    def jacobian_symbol(self):
        return self._jacobian

    @property
    def is_analytical(self):
        return self._base_mapping.is_analytical

    def copy(self):
        return InverseMapping(self._base_mapping)

#==============================================================================
class JacobianSymbol(MatrixSymbolicExpr):
    _axis = None
    def __new__(cls, mapping, axis=None):
        assert isinstance(mapping, SymbolicMapping)   # incl. structural mappings (WP06d-4a)
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
        if isinstance(new, SymbolicMapping):
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
        assert isinstance(mapping, SymbolicMapping)   # incl. structural mappings (WP06d-4a)
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
class InterfaceMapping(StructuralMapping, metaclass=_MappingABCMeta):
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
        # Mapping, not SymbolicMapping: interface legs are always concrete patch
        # mappings -- the body below needs Mapping's set_plus_minus() / copy().
        assert isinstance(minus, Mapping)
        assert isinstance(plus,  Mapping)
        minus = minus.copy()
        plus  = plus.copy()

        minus.set_plus_minus(minus=True)
        plus.set_plus_minus(plus=True)

        name       = '{}|{}'.format(str(minus.name), str(plus.name))
        ldim, pdim = minus.ldim, minus.pdim

        # WP06d-4a: build the IndexedBase directly instead of via Mapping.__new__.
        obj = IndexedBase.__new__(cls, name, shape=pdim)
        obj._name                = name
        obj._ldim                = ldim
        obj._pdim                = pdim
        obj._coordinates         = tuple(Symbol(u, real=True) for u in ['x', 'y', 'z'][:pdim])
        obj._logical_coordinates = Tuple(*(Symbol(u, real=True) for u in ['x1', 'x2', 'x3'][:ldim]))
        obj._is_minus            = None
        obj._is_plus             = None
        obj._minus               = minus
        obj._plus                = plus
        cls._init_symbolic_jacobian(obj)
        return obj

    @property
    def minus(self):
        return self._minus

    @property
    def plus(self):
        return self._plus

    @property
    def ldim(self):
        return self._ldim

    @property
    def pdim(self):
        return self._pdim

    @property
    def is_analytical(self):
        return self.minus.is_analytical and self.plus.is_analytical

    def copy(self):
        return InterfaceMapping(self._minus, self._plus)

    def _eval_subs(self, old, new):
        minus = self.minus.subs(old, new)
        plus  = self.plus.subs(old, new)
        return InterfaceMapping(minus, plus)

    def _eval_simplify(self, **kwargs):
        return self

#==============================================================================
class MultiPatchMapping(StructuralMapping, metaclass=_MappingABCMeta):

    def __new__(cls, dic):
        assert isinstance( dic, dict)
        obj = Basic.__new__(cls, dic)
        # WP06d-4a bug fix: MultiPatchMapping.__new__ never set _name, so any
        # domain call or `==` that reached `name` raised AttributeError.
        # An empty dict is degenerate but constructed fine pre-WP06d-4a.
        first = next(iter(dic.values()), None)
        obj._name                = '|'.join(str(m.name) for m in dic.values())
        obj._coordinates         = first._coordinates         if first is not None else None
        obj._logical_coordinates = first._logical_coordinates if first is not None else None
        obj._is_minus            = None
        obj._is_plus             = None
        # NB: no _init_symbolic_jacobian -- MultiPatchMapping uses Basic.__new__
        # (not IndexedBase), so it can't be indexed; pre-WP06d-4a it never went
        # through Mapping.__new__ either, so jacobian_expr stays None.
        return obj

    @property
    def mappings(self):
        return self.args[0]

    @property
    def ldim(self):
        return list(self.mappings.values())[0].ldim

    @property
    def pdim(self):
        return list(self.mappings.values())[0].pdim

    @property
    def is_analytical(self):
        return all(e.is_analytical for e in self.mappings.values())

    def copy(self):
        return MultiPatchMapping(dict(self.mappings))

    def _eval_subs(self, old, new):
        return self

    def _eval_simplify(self, **kwargs):
        return self

    def _hashable_content(self):
        return (type(self).__name__, *self.mappings.keys(), *self.mappings.values())

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
        # SymbolicMapping, not Mapping: since WP06d-4a the structural mappings
        # (InterfaceMapping / MultiPatchMapping / InverseMapping) are callable
        # on a domain but no longer subclass Mapping.
        assert(isinstance(mapping, SymbolicMapping))
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

        J = mapping.jacobian_symbol
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

        # SymbolicMapping, not Mapping: WP06d-4a's structural mappings ask for
        # their own Jacobian(self) at construction (via
        # StructuralMapping._init_symbolic_jacobian), and they are no longer
        # Mapping subclasses.
        if not isinstance(F, SymbolicMapping):
            raise TypeError('> Expecting a SymbolicMapping object')

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

        # SymbolicMapping, not Mapping: consistent with Covariant.eval (no
        # guard) and Jacobian.eval -- WP06d-4a's structural mappings are valid
        # here (Hdiv push-forward over a multipatch/broken domain).
        if not isinstance(F, SymbolicMapping):
            raise TypeError('> Expecting a SymbolicMapping')

        if not isinstance(v, (tuple, list, Tuple, ImmutableDenseMatrix, Matrix)):
            raise TypeError('> Expecting a tuple, list, Tuple, Matrix')

        M = Jacobian(F)
        M = M/M.det()
        v = Matrix(v)
        v = M*v
        return Tuple(*v)

#==============================================================================
class LogicalExpr(CalculusFunction):

    def __new__(cls, expr, domain, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(expr, domain, **options)
        else:
            r = None

        if r is None:
            obj = Basic.__new__(cls, expr, domain)
            return obj
        else:
            return r

    @property
    def expr(self):
        return self._args[0]

    @property
    def domain(self):
        return self._args[1]

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

    @classmethod
    def eval(cls, expr, domain, **options):
        """."""

        from sympde.expr.evaluation import TerminalExpr, DomainExpression
        from sympde.expr.expr import BilinearForm, LinearForm, BasicForm, Norm
        from sympde.expr.expr import Integral

        types = (ScalarFunction, VectorFunction, DifferentialOperator, Trace, Integral)

        mapping   = domain.mapping
        dim       = domain.dim
        assert mapping

        # TODO this is not the dim of the domain
        l_coords  = ['x1', 'x2', 'x3'][:dim]
        ph_coords = ['x', 'y', 'z']

        if not has(expr, types):
            if has(expr, DiffOperator):
                return cls( expr, domain, evaluate=False)
            else:
                syms = symbols(ph_coords[:dim], real=True)
                if isinstance(mapping, InterfaceMapping):
                    mapping = mapping.minus
                    # here we assume that the two mapped domains
                    # are identical in the interface so we choose one of them
                Ms   = [mapping[i] for i in range(dim)]
                expr = expr.subs(list(zip(syms, Ms)))

                if mapping.is_analytical:
                    expr = expr.subs(list(zip(Ms, mapping.expressions)))
                return expr

        if isinstance(expr, Symbol) and expr.name in l_coords:
            return expr

        if isinstance(expr, Symbol) and expr.name in ph_coords:
            return mapping[ph_coords.index(expr.name)]

        elif isinstance(expr, Add):
            args = [cls.eval(a, domain) for a in expr.args]
            v    =  S.Zero
            for i in args:
                v += i
            n,d = v.as_numer_denom()
            return n/d

        elif isinstance(expr, Mul):
            args = [cls.eval(a, domain) for a in expr.args]
            v    =  S.One
            for i in args:
                v *= i
            return v

        elif isinstance(expr, _logical_partial_derivatives):
            if mapping.is_analytical:
                Ms   = [mapping[i] for i in range(dim)]
                expr = expr.subs(list(zip(Ms, mapping.expressions)))
            return expr

        elif isinstance(expr, IndexedVectorFunction):
            el = cls.eval(expr.base, domain)
            el = TerminalExpr(el, domain=domain.logical_domain)
            return el[expr.indices[0]]

        elif isinstance(expr, MinusInterfaceOperator):
            mapping = mapping.minus
            newexpr = PullBack(expr.args[0], mapping)
            test    = newexpr.test
            newexpr = newexpr.expr.subs(test, MinusInterfaceOperator(test))
            return newexpr

        elif isinstance(expr, PlusInterfaceOperator):
            mapping = mapping.plus
            newexpr = PullBack(expr.args[0], mapping)
            test    = newexpr.test
            newexpr = newexpr.expr.subs(test, PlusInterfaceOperator(test))
            return newexpr

        elif isinstance(expr, (VectorFunction, ScalarFunction)):
            return PullBack(expr, mapping).expr

        elif isinstance(expr, Transpose):
            arg = cls(expr.arg, domain)
            return Transpose(arg)
            
        elif isinstance(expr, grad):
            arg = expr.args[0]
            if isinstance(mapping, InterfaceMapping):
                if isinstance(arg, MinusInterfaceOperator):
                    a     = arg.args[0]
                    mapping = mapping.minus
                elif isinstance(arg, PlusInterfaceOperator):
                    a = arg.args[0]
                    mapping = mapping.plus
                else:
                    raise TypeError(arg)

                arg = type(arg)(cls.eval(a, domain))
            else:
                arg = cls.eval(arg, domain)

            return mapping.jacobian_symbol.inv().T*grad(arg)

        elif isinstance(expr, curl):
            arg = expr.args[0]
            if isinstance(mapping, InterfaceMapping):
                if isinstance(arg, MinusInterfaceOperator):
                    arg     = arg.args[0]
                    mapping = mapping.minus
                elif isinstance(arg, PlusInterfaceOperator):
                    arg = arg.args[0]
                    mapping = mapping.plus
                else:
                    raise TypeError(arg)

            if isinstance(arg, VectorFunction):
                arg = PullBack(arg, mapping)
            else:
                arg = cls.eval(arg, domain)

            if isinstance(arg, PullBack) and isinstance(arg.kind, HcurlSpaceType):
                J   = mapping.jacobian_symbol
                arg = arg.test
                if isinstance(expr.args[0], (MinusInterfaceOperator, PlusInterfaceOperator)):
                    arg = type(expr.args[0])(arg)
                if expr.is_scalar:
                    return (1/J.det())*curl(arg)

                return (J/J.det())*curl(arg)
            else:
                raise NotImplementedError('TODO')

        elif isinstance(expr, div):
            arg = expr.args[0]
            if isinstance(mapping, InterfaceMapping):
                if isinstance(arg, MinusInterfaceOperator):
                    arg     = arg.args[0]
                    mapping = mapping.minus
                elif isinstance(arg, PlusInterfaceOperator):
                    arg = arg.args[0]
                    mapping = mapping.plus
                else:
                    raise TypeError(arg)

            if isinstance(arg, (ScalarFunction, VectorFunction)):
                arg = PullBack(arg, mapping)
            else:

                arg = cls.eval(arg, domain)

            if isinstance(arg, PullBack) and isinstance(arg.kind, HdivSpaceType):
                J   = mapping.jacobian_symbol
                arg = arg.test
                if isinstance(expr.args[0], (MinusInterfaceOperator, PlusInterfaceOperator)):
                    arg = type(expr.args[0])(arg)
                return (1/J.det())*div(arg)
            elif isinstance(arg, PullBack):
                return SymbolicTrace(mapping.jacobian_symbol.inv().T*grad(arg.test))
            else:
                raise NotImplementedError('TODO')

        elif isinstance(expr, laplace):
            arg = expr.args[0]
            v   = cls.eval(grad(arg), domain)
            v   = mapping.jacobian_symbol.inv().T*grad(v)
            return SymbolicTrace(v)

#        elif isinstance(expr, hessian):
#           arg = expr.args[0]
#            if isinstance(mapping, InterfaceMapping):
#                if isinstance(arg, MinusInterfaceOperator):
#                    arg     = arg.args[0]
#                    mapping = mapping.minus
#                elif isinstance(arg, PlusInterfaceOperator):
#                    arg = arg.args[0]
#                    mapping = mapping.plus
#                else:
#                    raise TypeError(arg)
#            v   = cls.eval(grad(expr.args[0]), domain)
#            v   = mapping.jacobian.inv().T*grad(v)
#            return v

        elif isinstance(expr, (dot, inner, outer)):
            args = [cls.eval(arg, domain) for arg in expr.args]
            return type(expr)(*args)

        elif isinstance(expr, _diff_ops):
            raise NotImplementedError('TODO')

        # TODO MUST BE MOVED AFTER TREATING THE CASES OF GRAD, CURL, DIV IN FEEC
        elif isinstance(expr, (Matrix, ImmutableDenseMatrix)):
            n_rows, n_cols = expr.shape
            lines          = []
            for i_row in range(0, n_rows):
                line = []
                for i_col in range(0, n_cols):
                    line.append(cls.eval(expr[i_row,i_col], domain))
                lines.append(line)
            return type(expr)(lines)

        elif isinstance(expr, dx):
            if expr.atoms(PlusInterfaceOperator):
                mapping = mapping.plus
            elif expr.atoms(MinusInterfaceOperator):
                mapping = mapping.minus

            arg = expr.args[0]
            arg = cls(arg, domain, evaluate=True)

            if isinstance(arg, PullBack):
                arg = TerminalExpr(arg, domain=domain.logical_domain)
            elif isinstance(arg, MatrixElement):
                arg = TerminalExpr(arg, domain=domain.logical_domain)
            # ...
            if dim == 1:
                lgrad_arg = LogicalGrad_1d(arg)

                if not isinstance(lgrad_arg, (list, tuple, Tuple, Matrix)):
                    lgrad_arg = Tuple(lgrad_arg)

            elif dim == 2:
                lgrad_arg = LogicalGrad_2d(arg)

            elif dim == 3:
                lgrad_arg = LogicalGrad_3d(arg)
            
            grad_arg = Covariant(mapping, lgrad_arg)
            expr = grad_arg[0]
            return expr

        elif isinstance(expr, dy):
            if expr.atoms(PlusInterfaceOperator):
                mapping = mapping.plus
            elif expr.atoms(MinusInterfaceOperator):
                mapping = mapping.minus

            arg = expr.args[0]
            arg = cls(arg, domain, evaluate=True)
            if isinstance(arg, PullBack):
                arg = TerminalExpr(arg, domain=domain.logical_domain)
            elif isinstance(arg, MatrixElement):
                arg = TerminalExpr(arg, domain=domain.logical_domain)

            # ..p
            if dim == 1:
                lgrad_arg = LogicalGrad_1d(arg)

            elif dim == 2:
                lgrad_arg = LogicalGrad_2d(arg)

            elif dim == 3:
                lgrad_arg = LogicalGrad_3d(arg)

            grad_arg = Covariant(mapping, lgrad_arg)

            expr = grad_arg[1]
            return expr

        elif isinstance(expr, dz):
            if expr.atoms(PlusInterfaceOperator):
                mapping = mapping.plus
            elif expr.atoms(MinusInterfaceOperator):
                mapping = mapping.minus

            arg = expr.args[0]
            arg = cls(arg, domain, evaluate=True)
            if isinstance(arg, PullBack):
                arg = TerminalExpr(arg, domain=domain.logical_domain)
            elif isinstance(arg, MatrixElement):
                arg = TerminalExpr(arg, domain=domain.logical_domain)
            # ...
            if dim == 1:
                lgrad_arg = LogicalGrad_1d(arg)

            elif dim == 2:
                lgrad_arg = LogicalGrad_2d(arg)

            elif dim == 3:
                lgrad_arg = LogicalGrad_3d(arg)

            grad_arg = Covariant(mapping, lgrad_arg)

            expr = grad_arg[2]

            return expr

        elif isinstance(expr, (Symbol, Indexed)):
            return expr

        elif isinstance(expr, NormalVector):
            return expr

        elif isinstance(expr, Pow):
            b = expr.base
            e = expr.exp
            expr =  Pow(cls(b, domain), cls(e, domain))
            return expr

        elif isinstance(expr, Trace):
            e = cls.eval(expr.expr, domain)
            bd = expr.boundary.logical_domain
            order = expr.order
            return Trace(e, bd, order)

        elif isinstance(expr, Integral):
            domain  = expr.domain
            mapping = domain.mapping


            assert domain is not None

            if expr.is_domain_integral:
                J   = mapping.jacobian_symbol
                det = sqrt((J.T*J).det())
            else:
                axis = domain.axis
                J    = JacobianSymbol(mapping, axis=axis)
                det  = sqrt((J.T*J).det())

            body   = cls.eval(expr.expr, domain)*det
            domain  = domain.logical_domain
            return Integral(body, domain)

        elif isinstance(expr, BilinearForm):
            tests   = [get_logical_test_function(a) for a in expr.test_functions]
            trials  = [get_logical_test_function(a) for a in expr.trial_functions]
            body    = cls.eval(expr.expr, domain)
            return BilinearForm((trials, tests), body)

        elif isinstance(expr, LinearForm):
            tests   = [get_logical_test_function(a) for a in expr.test_functions]
            body    = cls.eval(expr.expr, domain)
            return LinearForm(tests, body)

        elif isinstance(expr, Norm):
            kind           = expr.kind
            exponent       = expr.exponent
            e              = cls.eval(expr.expr, domain)
            domain         = domain.logical_domain
            norm           = Norm(e, domain, kind, evaluate=False)
            norm._exponent = exponent
            return norm

        elif isinstance(expr, DomainExpression):
            domain  = expr.target
            J       = domain.mapping.jacobian_symbol
            newexpr = cls.eval(expr.expr, domain)
            newexpr = TerminalExpr(newexpr, domain=domain)
            domain  = domain.logical_domain
            det     = TerminalExpr(sqrt((J.T*J).det()), domain=domain)
            return DomainExpression(domain, ImmutableDenseMatrix([[newexpr*det]]))
            
        elif isinstance(expr, Function):
            args = [cls.eval(a, domain) for a in expr.args]
            return type(expr)(*args)

        return cls(expr, domain, evaluate=False)

#==============================================================================
class SymbolicExpr(CalculusFunction):
    """returns a sympy expression where partial derivatives are converted into
    sympy Symbols."""

    @cacheit
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
    @cacheit
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
            v = Add(*args)
            return v

        elif isinstance(expr, Mul):
            args = [cls.eval(a, code=code) for a in expr.args]
            v    = Mul(*args)
            return v

        elif isinstance(expr, Pow):
            b = expr.base
            e = expr.exp
            v = Pow(cls.eval(b, code=code), e)
            return v

        elif isinstance(expr, _coeffs_registery):
            return expr

        elif isinstance(expr, (list, tuple, Tuple)):
            expr = [cls.eval(a, code=code) for a in expr]
            return Tuple(*expr)

        elif isinstance(expr, (Matrix, ImmutableDenseMatrix)):

            lines = []
            n_row,n_col = expr.shape
            for i_row in range(0,n_row):
                line = []
                for i_col in range(0,n_col):
                    line.append(cls.eval(expr[i_row, i_col], code=code))

                lines.append(line)

            return type(expr)(lines)

        elif isinstance(expr, (ScalarFunction, VectorFunction)):
            if code:
                name = '{name}_{code}'.format(name=expr.name, code=code)
            else:
                name = str(expr.name)

            return Symbol(name)

        elif isinstance(expr, ( PlusInterfaceOperator, MinusInterfaceOperator)):
            return cls.eval(expr.args[0], code=code)

        elif isinstance(expr, Indexed):
            base = expr.base
            # SymbolicMapping, not Mapping: WP06d-4a's structural mappings are
            # IndexedBase-backed but no longer Mapping (mirrors the L1800 check).
            if isinstance(base, SymbolicMapping):
                if expr.indices[0] == 0:
                    name = 'x'
                elif expr.indices[0] == 1:
                    name = 'y'
                elif expr.indices[0] == 2:
                    name = 'z'
                else:
                    raise ValueError('Wrong index')

                if base.is_plus:
                    name = name + '_plus'
            else:
                name =  '{base}_{i}'.format(base=base.name, i=expr.indices[0])

            if code:
                name = '{name}_{code}'.format(name=name, code=code)

            return Symbol(name)

        elif isinstance(expr, _partial_derivatives):
            atom = get_atom_derivatives(expr)
            indices = get_index_derivatives_atom(expr, atom)
            code = None
            if indices:
                index = indices[0]
                code = ''
                index =dict(sorted(index.items()))

                for k,n in list(index.items()):
                    code += k*n
            return cls.eval(atom, code=code)

        elif isinstance(expr, _logical_partial_derivatives):
            atom = get_atom_logical_derivatives(expr)
            indices = get_index_logical_derivatives_atom(expr, atom)
            code = None
            if indices:
                index = indices[0]
                code = ''
                index = dict(sorted(index.items()))
                for k,n in list(index.items()):
                    code += k*n
            return cls.eval(atom, code=code)

        elif isinstance(expr, SymbolicMapping):
            return Symbol(expr.name)

        # ... this must be done here, otherwise codegen for FEM will not work
        elif isinstance(expr, Symbol):
            return expr

        elif isinstance(expr, IndexedBase):
            return expr

        elif isinstance(expr, Indexed):
            return expr

        elif isinstance(expr, Idx):
            return expr

        elif isinstance(expr, Function):
            args = [cls.eval(a, code=code) for a in expr.args]
            return type(expr)(*args)

        elif isinstance(expr, ImaginaryUnit):
            return expr


        elif isinstance(expr, SymbolicWeightedVolume):
            mapping = expr.args[0]
            if isinstance(mapping, InterfaceMapping):
                mapping = mapping.minus
            name = 'wvol_{mapping}'.format(mapping=mapping)

            return Symbol(name)

        elif isinstance(expr, SymbolicDeterminant):
            name = 'det_{}'.format(str(expr.args[0]))
            return Symbol(name)

        elif isinstance(expr, PullBack):
            return cls.eval(expr.expr, code=code)

        # Expression must always be translated to Sympy!
        # TODO: check if we should use 'sympy.sympify(expr)' instead
        else:
            raise NotImplementedError('Cannot translate to Sympy: {}'.format(expr))
