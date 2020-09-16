# coding: utf-8

from collections  import OrderedDict

from sympy import Indexed, IndexedBase, Idx
from sympy import Matrix, ImmutableDenseMatrix
from sympy import Function, Expr
from sympy import sympify
from sympy import cacheit

from sympy.core import Basic
from sympy.core import Symbol
from sympy.core import Add, Mul, Pow

from sympy.core.expr       import AtomicExpr
from sympy.core.numbers    import ImaginaryUnit
from sympy.core.containers import Tuple
from sympy                 import S

from sympde.calculus.core  import PlusInterfaceOperator, MinusInterfaceOperator
from sympde.calculus.core  import grad, div, rot, curl, laplace, hessian
from sympde.calculus.core  import dot, inner, outer, _diff_ops
from sympde.calculus.core  import has, DiffOperator

from sympde.calculus.matrices import MatrixSymbolicExpr, MatrixElement, SymbolicTrace
from sympde.calculus.matrices import SymbolicDeterminant

from sympde.core       import Constant
from sympde.core.basic import BasicMapping
from sympde.core.basic import CalculusFunction
from sympde.core.basic import _coeffs_registery


from .basic       import BasicDomain, Union, InteriorDomain
from .basic       import Boundary, Connectivity, Interface
from .domain      import Domain, NCubeInterior
from .domain      import NormalVector
from .space       import ScalarTestFunction, VectorTestFunction, IndexedTestTrial
from .space       import ScalarField, VectorField, IndexedVectorField
from .space       import Trace, trace_0, trace_1
from .datatype    import HcurlSpaceType, H1SpaceType, L2SpaceType, HdivSpaceType, UndefinedSpaceType
from .derivatives import dx, dy, dz, DifferentialOperator
from .derivatives import _partial_derivatives
from .derivatives import get_atom_derivatives, get_index_derivatives_atom
from .derivatives import _logical_partial_derivatives
from .derivatives import get_atom_logical_derivatives, get_index_logical_derivatives_atom
from .derivatives import LogicalGrad_1d, LogicalGrad_2d, LogicalGrad_3d

# TODO fix circular dependency between sympde.topology.domain and sympde.topology.mapping
# TODO fix circular dependency between sympde.expr.evaluation and sympde.topology.mapping
#==============================================================================
from sympy                   import sqrt, symbols
from sympy.core.exprtools    import factor_terms
from sympy.polys.polytools   import parallel_poly_from_expr
from sympy                   import cacheit

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
class Mapping(BasicMapping):
    """
    Represents a Mapping object.

    Examples

    """
    _expressions = None # used for analytical mapping
    _rdim        = None
    # TODO shall we keep rdim ?
    def __new__(cls, name, rdim=None, coordinates=None, **kwargs):
        if isinstance(rdim, (tuple, list, Tuple)):
            if not len(rdim) == 1:
                raise ValueError('> Expecting a tuple, list, Tuple of length 1')

            rdim = rdim[0]
        elif rdim is None:
            rdim = cls._rdim

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

        obj._name                = name
        obj._rdim                = rdim
        obj._coordinates         = _coordinates
        obj._jacobian            = kwargs.pop('jacobian', JacobianSymbol(obj))

        lcoords = ['x1', 'x2', 'x3'][:rdim]
        lcoords = [Symbol(i) for i in lcoords]
        obj._logical_coordinates = Tuple(*lcoords)
        # ...
        if not( obj._expressions is None ):
            coords = ['x', 'y', 'z'][:rdim]

            # ...
            args = []
            for i in coords:
                x = obj._expressions[i]
                x = sympify(x)
                args.append(x)

            args = Tuple(*args)
            # ...
            zero_coords = ['x1', 'x2', 'x3'][rdim:]

            for i in zero_coords:
                x = sympify(i)
                args = args.subs(x,0)
            # ...

            constants = list(set(args.free_symbols) - set(lcoords))
            # subs constants as Constant objects instead of Symbol
            d = {}
            for i in constants:
                # TODO shall we add the type?
                # by default it is real
                if i.name in kwargs:
                    d[i] = kwargs[i.name]
                else:
                    d[i] = Constant(i.name)

            args = args.subs(d)
            # ...

            obj._expressions = args
        # ...

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

    @property
    def logical_coordinates(self):
        if self.rdim == 1:
            return self._logical_coordinates[0]
        else:
            return self._logical_coordinates

    def __call__(self, domain):
        assert(isinstance(domain, BasicDomain))
        return MappedDomain(self, domain)

    @property
    def jacobian(self):
        return self._jacobian

    @property
    def det_jacobian(self):
        return self.jacobian.det()

    @property
    def is_analytical(self):
        return not( self._expressions is None )

    @property
    def expressions(self):
        return self._expressions

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

class InverseMapping(Mapping):
    def __new__(cls, mapping):
        assert isinstance(mapping, Mapping)
        name     = mapping.name
        rdim     = mapping.rdim
        coords   = mapping.logical_coordinates
        jacobian = mapping.jacobian.inv()
        return Mapping.__new__(cls, name, rdim , coordinates=coords, jacobian=jacobian)

class JacobianSymbol(MatrixSymbolicExpr):

    def __new__(cls, mapping, axis=None):
        assert isinstance(mapping, Mapping)
        if axis is not None:
            assert isinstance(axis, int)
        return MatrixSymbolicExpr.__new__(cls, mapping, axis)

    @property
    def mapping(self):
        return self._args[0]

    @property
    def axis(self):
        return self._args[1]

    def _hashable_content(self):
        return (self.mapping, self.axis)

    def _sympystr(self, printer):
        sstr = printer.doprint
        if self.axis:
            return 'Jacobian({},{})'.format(sstr(self.mapping.name), self.axis)
        else:
            return 'Jacobian({})'.format(sstr(self.mapping.name))

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
        name = '{}|{}'.format(str(minus), str(plus))
        obj  = Mapping.__new__(cls, name, rdim=minus.rdim)
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

class MultiPatchMapping(Mapping):

    def __new__(cls, dic):
        assert isinstance( dic, dict)
        return Basic.__new__(cls, dict)

    def mappings(self):
        return self.args[0]

    @property
    def is_analytical(self):
        return all(a.is_analytical for a in self.mappings.values())

    @property
    def rdim(self):
        return self.args[0].rdim
#==============================================================================
class IdentityMapping(Mapping):
    """
    Represents an identity 1D/2D/3D Mapping object.

    Examples

    """
    _expressions = {'x': 'x1',
                    'y': 'x2',
                    'z': 'x3'}

#==============================================================================
class AffineMapping(Mapping):
    """
    Represents a 1D/2D/3D Affine Mapping object.

    Examples

    """
    _expressions = {'x': 'c1 + a11*x1 + a12*x2 + a13*x3',
                    'y': 'c2 + a21*x1 + a22*x2 + a23*x3',
                    'z': 'c3 + a31*x1 + a32*x2 + a33*x3'}

#==============================================================================
class PolarMapping(Mapping):
    """
    Represents a Polar 2D Mapping object (Annulus).

    Examples

    """
    _expressions = {'x': 'c1 + (rmin*(1-x1)+rmax*x1)*cos(x2)',
                    'y': 'c2 + (rmin*(1-x1)+rmax*x1)*sin(x2)'}

    _rdim        = 2
#==============================================================================
class TargetMapping(Mapping):
    """
    Represents a Target 2D Mapping object.

    Examples

    """
    _expressions = {'x': 'c1 + (1-k)*x1*cos(x2) - D*x1**2',
                    'y': 'c2 + (1+k)*x1*sin(x2)'}

    _rdim        = 2
#==============================================================================
class CzarnyMapping(Mapping):
    """
    Represents a Czarny 2D Mapping object.

    Examples

    """
    _expressions = {'x': '(1 - sqrt( 1 + eps*(eps + 2*x1*cos(x2)) )) / eps',
                    'y': 'c2 + (b / sqrt(1-eps**2/4) * x1 * sin(x2)) /'
                        '(2 - sqrt( 1 + eps*(eps + 2*x1*cos(x2)) ))'}

    _rdim        = 2
#==============================================================================
class CollelaMapping(Mapping):
    """
    Represents a Collela 2D Mapping object.

    Examples

    """
    _expressions = {'x': '2.*(x1 + eps*sin(2.*pi*k1*x1)*sin(2.*pi*k2*x2)) - 1.',
                    'y': '2.*(x2 + eps*sin(2.*pi*k1*x1)*sin(2.*pi*k2*x2)) - 1.'}

    _rdim        = 2
#==============================================================================
class TorusMapping(Mapping):
    """
    Represents a Torus 3D Mapping object.

    Examples

    """
    _expressions = {'x': '(R0+x1*cos(x2))*cos(x3)',
                    'y': '(R0+x1*cos(x2))*sin(x3)',
                    'z': 'x1*sin(x2)'}

    _rdim        = 3
#==============================================================================
class TwistedTargetMapping(Mapping):
    """
    Represents a Twisted Target 3D Mapping object.

    Examples

    """
    _expressions = {'x': 'c1 + (1-k)*x1*cos(x2) - D*x1**2',
                    'y': 'c2 + (1+k)*x1*sin(x2)',
                    'z': 'c3 + x3*x1**2*sin(2*x2)'}

    _rdim        = 3

#==============================================================================
class MappedDomain(BasicDomain):
    """."""

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
            return Domain(logical_domain.name, **kwargs)

        elif isinstance(logical_domain, NCubeInterior):
            name  = logical_domain.name
            dim   = logical_domain.dim
            dtype = logical_domain.dtype
            min_coords = logical_domain.min_coords
            max_coords = logical_domain.max_coords
            return NCubeInterior(name, dim, dtype, min_coords, max_coords, mapping, logical_domain)
        elif isinstance(logical_domain, InteriorDomain):
            name  = logical_domain.name
            dim   = logical_domain.dim
            dtype = logical_domain.dtype
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

    @cacheit
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
        if not isinstance(u, (VectorTestFunction, ScalarTestFunction, ScalarField, VectorField)):
            raise TypeError('{} must be of type ScalarTestFunction or VectorTestFunction'.format(str(u)))

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

        J               = mapping.jacobian
        if isinstance(kind, (UndefinedSpaceType, H1SpaceType)):
            expr =  el
        elif isinstance(kind, HdivSpaceType):
            expr  =  (J/J.det())*ImmutableDenseMatrix(tuple(el[i] for i in range(dim)))
        elif isinstance(kind , HcurlSpaceType):
            expr  = J.inv().T*ImmutableDenseMatrix(tuple(el[i] for i in range(dim)))
        elif isinstance(kind, L2SpaceType):
            expr = J.det()*el
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
    """
    This class calculates the Jacobian of a mapping

    """

    @classmethod
    @cacheit
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
class Covariant(MappingApplication):
    """

    Examples

    """

    @classmethod
    @cacheit
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

        M   = Jacobian(F).inv().T
        dim = F.rdim
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
    @cacheit
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

#==============================================================================
class LogicalExpr(CalculusFunction):

    @cacheit
    def __new__(cls, expr, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(expr, **options)
            if options.pop('subs', False):
                mapping = options['mapping']
                if mapping.is_analytical and not isinstance(mapping, InterfaceMapping):
                    for i in range(mapping.rdim):
                        r = r.subs(mapping[i], mapping.expressions[i])

                elif mapping.is_analytical and isinstance(mapping, InterfaceMapping):
                    M1 = mapping.minus
                    M2 = mapping.plus
                    for i in range(M1.rdim):
                        r = r.subs(M1[i], M1.expressions[i])
                        r = r.subs(M2[i], M2.expressions[i])
        else:
            r = None

        if r is None:
            obj          = Basic.__new__(cls, expr)
            obj._mapping = options.pop('mapping', None)
            obj._dim     = options.pop('dim', None)
            return obj
        else:
            return r

    @property
    def expr(self):
        return self._args[0]

    @property
    def mapping(self):
        return self._mapping

    @property
    def dim(self):
        return self._dim

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

    @classmethod
    @cacheit
    def eval(cls, expr, mapping=None, dim=None, **options):
        """."""

        from sympde.expr.evaluation import TerminalExpr, DomainExpression
        from sympde.expr.expr import BilinearForm, LinearForm, BasicForm, Norm
        from sympde.expr.expr import Integral

        types = (ScalarTestFunction, VectorTestFunction, VectorField, ScalarField, DifferentialOperator, Trace, Integral)

        # TODO this is not the dim of the domain
        l_coords  = ['x1', 'x2', 'x3'][:dim]
        ph_coords = ['x', 'y', 'z']

        if not has(expr, types):
            if has(expr, DiffOperator):
                return cls( expr, mapping=mapping, dim=dim, evaluate=False)
            else:
                syms = symbols(ph_coords[:dim])
                if isinstance(mapping, InterfaceMapping):
                    mapping = mapping.minus
                    # here we assume that the two mapped domains
                    # are identical in the interface so we choose one of them
                Ms   = [mapping[i] for i in range(dim)]
                return expr.subs(zip(syms, Ms))

        if isinstance(expr, Symbol) and expr.name in l_coords:
            return expr

        if isinstance(expr, Symbol) and expr.name in ph_coords:
            return mapping[ph_coords.index(expr.name)]

        elif isinstance(expr, Add):
            args = [cls.eval(a, mapping=mapping, dim=dim) for a in expr.args]
            v    =  S.Zero
            for i in args:
                v += i
            n,d = v.as_numer_denom()
            return n/d

        elif isinstance(expr, Mul):
            args = [cls.eval(a, mapping=mapping, dim=dim) for a in expr.args]
            v    =  S.One
            for i in args:
                v *= i
            return v

        elif isinstance(expr, _logical_partial_derivatives):
            return expr

        elif isinstance(expr, (IndexedTestTrial, IndexedVectorField)):
            el = cls.eval(expr.base, mapping=mapping, dim=dim)
            return el[expr.indices[0]]

        elif isinstance(expr, MinusInterfaceOperator):
            mapping = mapping.minus
            return MinusInterfaceOperator(PullBack(expr.args[0], mapping).expr)

        elif isinstance(expr, PlusInterfaceOperator):
            mapping = mapping.plus
            return PlusInterfaceOperator(PullBack(expr.args[0], mapping).expr)

        elif isinstance(expr, (ScalarField, VectorField, 
                               VectorTestFunction, 
                               ScalarTestFunction)):
            return PullBack(expr, mapping).expr

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

                arg = type(arg)(cls.eval(a, mapping=mapping, dim=dim))
            else:
                arg = cls.eval(arg, mapping=mapping, dim=dim)

            return mapping.jacobian.inv().T*grad(arg)

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

            if isinstance(arg,(VectorField, VectorTestFunction)):
                arg = PullBack(arg, mapping)
            else:
                arg = cls.eval(arg, mapping=mapping, dim=dim)
                
            if isinstance(arg, PullBack) and isinstance(arg.kind, HcurlSpaceType):
                J   = mapping.jacobian
                arg = arg.test
                if isinstance(expr.args[0], (MinusInterfaceOperator, PlusInterfaceOperator)):
                    arg = type(expr.args[0])(arg)
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

            if isinstance(arg,(VectorField, VectorTestFunction)):
                arg = PullBack(arg, mapping)
            else:
                arg = cls.eval(arg, mapping=mapping, dim=dim)
            if isinstance(arg, PullBack) and isinstance(arg.kind, HdivSpaceType):
                J   = mapping.jacobian
                arg = arg.test
                if isinstance(expr.args[0], (MinusInterfaceOperator, PlusInterfaceOperator)):
                    arg = type(expr.args[0])(arg)
                return (1/J.det())*div(arg)
            else:
                raise NotImplementedError('TODO')

        elif isinstance(expr, laplace):
            arg = expr.args[0]
            v   = cls.eval(grad(arg), mapping=mapping, dim=dim)
            v   = mapping.jacobian.inv()*grad(v)
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
#            v   = cls.eval(grad(expr.args[0]), mapping=mapping, dim=dim)
#            v   = mapping.jacobian.inv().T*grad(v)
#            return v
          
        elif isinstance(expr, (dot, inner, outer)):
            args = [cls.eval(arg, mapping=mapping, dim=dim) for arg in expr.args]
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
                    line.append(cls.eval(expr[i_row,i_col], mapping=mapping, dim=dim))

                lines.append(line)

            return type(expr)(lines)

        elif isinstance(expr, dx):
            if expr.atoms(PlusInterfaceOperator):
                mapping = mapping.plus
            elif expr.atoms(MinusInterfaceOperator):
                mapping = mapping.minus

            arg = expr.args[0]
            arg = cls(arg, mapping=mapping, dim=dim, evaluate=True)

            if isinstance(arg, PullBack):
                arg = TerminalExpr(arg, dim=dim)
            elif isinstance(arg, MatrixElement):
                arg = TerminalExpr(arg, dim=dim, logical=True)
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
            arg = cls(arg, mapping=mapping, dim=dim, evaluate=True)
            if isinstance(arg, PullBack):
                arg = TerminalExpr(arg, dim=dim)
            elif isinstance(arg, MatrixElement):
                arg = TerminalExpr(arg, dim=dim, logical=True)

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
            arg = cls(arg, mapping=mapping, dim=dim, evaluate=True)
            if isinstance(arg, PullBack):
                arg = TerminalExpr(arg, dim=dim)
            elif isinstance(arg, MatrixElement):
                arg = TerminalExpr(arg, dim=dim, logical=True)
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
            expr =  Pow(cls(b, mapping=mapping, dim=dim), cls(e, mapping=mapping, dim=dim))
            return expr

        elif isinstance(expr, Trace):
            e = cls.eval(expr.expr, mapping=mapping, dim=dim)
            bd = expr.boundary.logical_domain
            order = expr.order
            return Trace(e, bd, order)

        elif isinstance(expr, Integral):
            domain  = expr.domain
            mapping = domain.mapping
            domain  = domain.logical_domain
            dim     = domain.dim
            assert domain is not None

            if expr.is_domain_integral:
                J   = mapping.jacobian
                det = sqrt((J.T*J).det())
            else:
                axis = expr.domain.axis
                J    = JacobianSymbol(mapping, axis)
                det  = sqrt((J.T*J).det())
            body   = cls.eval(expr.expr, mapping=mapping, dim=dim)*det
            return Integral(body, domain)

        elif isinstance(expr, BilinearForm):
            tests  = [get_logical_test_function(a) for a in expr.test_functions]
            trials = [get_logical_test_function(a) for a in expr.trial_functions]
            body   = cls.eval(expr.expr, mapping=mapping, dim=dim)
            return BilinearForm((trials, tests), body)

        elif isinstance(expr, LinearForm):
            tests  = [get_logical_test_function(a) for a in expr.test_functions]
            body   = cls.eval(expr.expr, mapping=mapping, dim=dim)
            return LinearForm(tests, body)
        elif isinstance(expr, Norm):
            kind           = expr.kind
            domain         = expr.domain.logical_domain
            exponent       = expr.exponent
            e              = cls.eval(expr.expr, mapping=mapping, dim=dim)
            norm           = Norm(e, domain, kind, evaluate=False)
            norm._exponent = exponent
            return norm
        elif isinstance(expr, DomainExpression):
            domain  = expr.target.logical_domain
            mapping = expr.target.mapping
            dim     = domain.dim
            J       = mapping.jacobian
            newexpr = cls.eval(expr.expr, mapping=mapping, dim=dim)
            det     = TerminalExpr(sqrt((J.T*J).det()), dim=dim, logical=True)
            
            return DomainExpression(domain, ImmutableDenseMatrix([[newexpr*det]]))
            
        elif isinstance(expr, Function):
            args = [cls.eval(a, mapping=mapping, dim=dim) for a in expr.args]
            return type(expr)(*args)
        return cls(expr, mapping=mapping, dim=dim, evaluate=False)

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

        elif isinstance(expr, (ScalarField, ScalarTestFunction, VectorField, VectorTestFunction)):
            if code:
                name = '{name}_{code}'.format(name=expr.name, code=code)
            else:
                name = str(expr.name)

            return Symbol(name)

        elif isinstance(expr, ( PlusInterfaceOperator, MinusInterfaceOperator)):
            return cls.eval(expr.args[0], code=code)

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

        elif isinstance(expr, _partial_derivatives):
            atom = get_atom_derivatives(expr)
            indices = get_index_derivatives_atom(expr, atom)
            code = None
            if indices:
                index = indices[0]
                code = ''
                index = OrderedDict(sorted(index.items()))

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
                index = OrderedDict(sorted(index.items()))
                for k,n in list(index.items()):
                    code += k*n
            return cls.eval(atom, code=code)

        elif isinstance(expr, Mapping):
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

