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
from sympde.calculus.core  import grad, div, rot, curl, dot, inner, outer, _diff_ops
from sympde.calculus.matrices import MatrixSymbolicExpr, MatrixElement
from sympy                 import Determinant
from sympde.core       import Constant
from sympde.core.basic import BasicMapping
from sympde.core.basic import CalculusFunction
from sympde.core.basic import _coeffs_registery
from sympde.topology   import NormalVector

from .basic       import BasicDomain, Union, InteriorDomain
from .domain      import Domain
from .space       import ScalarTestFunction, VectorTestFunction, IndexedTestTrial
from .space       import ScalarField, VectorField, IndexedVectorField
from .datatype    import HcurlSpaceType, H1SpaceType, L2SpaceType, HdivSpaceType, UndefinedSpaceType
from .derivatives import dx, dy, dz
from .derivatives import _partial_derivatives
from .derivatives import get_atom_derivatives, get_index_derivatives_atom
from .derivatives import _logical_partial_derivatives
from .derivatives import get_atom_logical_derivatives, get_index_logical_derivatives_atom
from .derivatives import LogicalGrad_1d, LogicalGrad_2d, LogicalGrad_3d

#==============================================================================
from sympy.core.exprtools    import factor_terms
from sympy.polys.polytools   import parallel_poly_from_expr
from sympy                   import cacheit

@cacheit
def cancel(f):
    f           = factor_terms(f, radical=True)
    p, q        = f.as_numer_denom()
    # TODO accelerate parallel_poly_from_expr
    (p, q), opt = parallel_poly_from_expr((p,q))
    c, P, Q     = p.cancel(q)
    return c*(P.as_expr()/Q.as_expr())

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

    def __new__(cls, mapping):
        assert isinstance(mapping, Mapping)
        return Expr.__new__(cls, mapping)

    @property
    def mapping(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
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
        return Basic.__new__(cls, minus, plus)

    @property
    def minus(self):
        return self._args[0]

    @property
    def plus(self):
        return self._args[1]

    @property
    def is_analytical(self):
        return self.minus.is_analytical and self.plus.is_analytical

    @property
    def rdim(self):
        return self.minus.rdim

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
            boundaries     = logical_domain.boundary,
            interiors      = logical_domain.interior
            if isinstance(interiors, Union):
                kwargs['interiors'] = Union(*[mapping(a) for a in interiors.args])
            else:
                kwargs['interiors'] = mapping(interiors)

            if isinstance(boundaries, Union):
                kwargs['boundaries'] = Union(*[mapping(a) for a in boundaries])

            interfaces =  logical_domain.connectivity.interfaces
            if interfaces:
                if isinstance(interfaces, Union):
                    interfaces = interfaces.args
                else:
                    interfaces = [interfaces]
                connectivity = {}
                for e in interfaces:
                    connectivity[e.name] = (mapping(e.minus), mapping(e.plus))
                kwargs['connectivity'] = Connectivity(connectivity)
            return Domain(logical_domain.name, **kwargs)

        elif isinstance(logical_domain, InteriorDomain):
            name  = logical_domain.name
            dim   = logical_domain.dim
            dtype = logical_domain.dtype
            return InteriorDomain(name, dim, dtype, mapping, logical_domain)

#==============================================================================
class SymbolicWeightedVolume(Expr):
    pass
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

    def __new__(cls, u):
        if not isinstance(u, (VectorTestFunction, ScalarTestFunction)):
            raise TypeError('{} must be of type ScalarTestFunction or VectorTestFunction'.format(str(u)))

        space           = u.space
        kind            = space.kind
        dim             = space.ldim
        logical_domain  = space.domain.logical_domain
        l_space         = type(space)(space.name, logical_domain, kind=kind)
        el              = l_space.element(u.name)
        J               = space.domain.mapping.jacobian
        if isinstance(kind, H1SpaceType):
            expr =  el
        elif isinstance(kind, HdivSpaceType):
            A     = J/J.det()
            expr  =  A*ImmutableDenseMatrix(tuple(el[i] for i in range(dim)))
        elif isinstance(kind , HcurlSpaceType):
            expr  = J.inv().T*ImmutableDenseMatrix(tuple(el[i] for i in range(dim)))
        elif isinstance(kind, L2SpaceType):
            expr = J.det()*el
        elif isinstance(kind, UndefinedSpaceType):
            raise ValueError('kind must be specified in order to perform the pull-back transformation')
        else:
            raise NotImplementedError('TODO')

        return Expr.__new__(cls, expr, kind, el)

    @property
    def expr(self):
        return self._args[0]

    @property
    def kind(self):
        return self._args[1]

    @property
    def test(self):
        return self._args[2]
#==============================================================================
class Jacobian(MappingApplication):
    """
    Examples

    """

    @classmethod
    @cacheit
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
class Covariant(MappingApplication):
    """

    Examples

    """

    @classmethod
    @cacheit
    def eval(cls, F, v):

        if not isinstance(v, (tuple, list, Tuple, Matrix)):
            raise TypeError('> Expecting a tuple, list, Tuple, Matrix')

        M = Jacobian(F).inv().T
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

        if not isinstance(F, Mapping):
            raise TypeError('> Expecting a Mapping')

        if not isinstance(v, (tuple, list, Tuple, Matrix)):
            raise TypeError('> Expecting a tuple, list, Tuple, Matrix')

        M = Jacobian(F)
        M = M/M.det()
        v = Matrix(v)
        v = M*v
        return Tuple(*v)

#==============================================================================
class LogicalExpr(CalculusFunction):

    @cacheit
    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            #expr = args[0]
            #M    = list(expr.atoms(MappedDomain))
            #assert len(M) == 1
            #M    = M[0].mapping
            # TODO remove mapping from args
            r = cls.eval(*args)
            M = args[0]
            if M.is_analytical and not isinstance(M, InterfaceMapping):
                for i in range(M.rdim):
                    r = r.subs(M[i], M.expressions[i])

            elif isinstance(M, InterfaceMapping):
                M1 = M.minus
                M2 = M.plus
                for i in range(M1.rdim):
                    r = r.subs(M1[i], M1.expressions[i])
                    r = r.subs(M2[i], M2.expressions[i])
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
    def eval(cls, *_args):
        """."""

        if not _args:
            return

        if not len(_args) == 2:
            raise ValueError('Expecting two arguments')

        from sympde.expr.evaluation import TerminalExpr

        M         = _args[0]
        expr      = _args[1]
        dim       = M.rdim # TODO this is not the dim of the domain
        l_coords  = ['x1', 'x2', 'x3'][:dim]
        ph_coords = ['x', 'y', 'z']

        if isinstance(expr, Symbol) and expr.name in l_coords:
            return expr
        elif isinstance(expr, Symbol) and expr.name in ph_coords:
            if M.is_analytical:
                return M[ph_coords.index(expr.name)]
            return expr

        elif isinstance(expr, Add):
            args = [cls.eval(M, a) for a in expr.args]
            v    =  Add(*args)
            n,d = v.as_numer_denom()
            return cancel(n/d)

        elif isinstance(expr, Mul):
            args = [cls.eval(M, a) for a in expr.args]
            v    =  Mul(*args)
            return v

        elif isinstance(expr, _coeffs_registery):
            return expr

        elif isinstance(expr, _logical_partial_derivatives):
            if M.is_analytical:
                arg = cls.eval(M, expr._args[0])
                op  = expr
                return op.eval(arg)
            else:
                return expr
        elif isinstance(expr, (IndexedTestTrial, IndexedVectorField)):
            el = cls.eval(M, expr.base)
            if isinstance(el, PullBack):
                el = el.expr
            return el[expr.indices[0]]
        elif isinstance(expr, (VectorField, VectorTestFunction, ScalarField, ScalarTestFunction)):
            return PullBack(expr)
        elif isinstance(expr, grad):
            arg = cls.eval(M, expr.args[0])
            if isinstance(arg, PullBack):
                arg = arg.expr
            return M.jacobian.inv().T*grad(arg)

        elif isinstance(expr, curl):
            arg = cls.eval(M, expr.args[0])
            if isinstance(arg, PullBack) and isinstance(arg.kind, HcurlSpaceType):
                J   = M.jacobian
                return (J/J.det())*curl(arg.test)
            else:
                raise NotImplementedError('TODO')

        elif isinstance(expr, div):
            arg = expr.args[0]
            arg = cls.eval(M, expr.args[0])
            if isinstance(arg, PullBack) and isinstance(arg.kind, HdivSpaceType):
                J   = M.jacobian
                return (1/J.det())*div(arg.test)
            else:
                raise NotImplementedError('TODO')
        elif isinstance(expr, (dot, inner, outer)):
            args = [cls.eval(M, arg) for arg in expr.args]
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
                    line.append(cls.eval(M, expr[i_row,i_col]))

                lines.append(line)

            return type(expr)(lines)

        elif isinstance(expr, dx):
            if expr.atoms(PlusInterfaceOperator):
                M = M.plus
            elif expr.atoms(MinusInterfaceOperator):
                M = M.minus

            arg = expr.args[0]
            arg = cls(M, arg, evaluate=True)
            if isinstance(arg, MatrixSymbolicExpr) or isinstance(arg, MatrixElement):
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

            grad_arg = Covariant(M, lgrad_arg)

            expr = grad_arg[0]
            return expr

        elif isinstance(expr, dy):
            if expr.atoms(PlusInterfaceOperator):
                M = M.plus
            elif expr.atoms(MinusInterfaceOperator):
                M = M.minus

            arg = expr.args[0]
            arg = cls(M, arg, evaluate=True)
            if isinstance(arg, MatrixSymbolicExpr) or isinstance(arg, MatrixElement):
                arg = TerminalExpr(arg, dim=dim, logical=True)

            # ..p
            if dim == 1:
                lgrad_arg = LogicalGrad_1d(arg)

            elif dim == 2:
                lgrad_arg = LogicalGrad_2d(arg)

            elif dim == 3:
                lgrad_arg = LogicalGrad_3d(arg)

            grad_arg = Covariant(M, lgrad_arg)

            expr = grad_arg[1]
            return expr

        elif isinstance(expr, dz):
            if expr.atoms(PlusInterfaceOperator):
                M = M.plus
            elif expr.atoms(MinusInterfaceOperator):
                M = M.minus

            arg = expr.args[0]
            arg = cls(M, arg, evaluate=True)
            if isinstance(arg, MatrixSymbolicExpr) or isinstance(arg, MatrixElement):
                arg = TerminalExpr(arg, dim=dim, logical=True)
            # ...
            if dim == 1:
                lgrad_arg = LogicalGrad_1d(arg)

            elif dim == 2:
                lgrad_arg = LogicalGrad_2d(arg)

            elif dim == 3:
                lgrad_arg = LogicalGrad_3d(arg)

            grad_arg = Covariant(M, lgrad_arg)

            expr = grad_arg[2]

            return expr

        elif isinstance(expr, (Symbol, Indexed)):
            return expr
        elif isinstance(expr, NormalVector):
            return expr
        elif isinstance(expr, Function):
            func = expr.func
            args = expr.args
            args = [cls(M, a) for a in args]
            return func(*args)

        elif isinstance(expr, Pow):
            b = expr.base
            e = expr.exp
            expr =  Pow(cls(M, b), cls(M, e))
            return expr

        return cls(M, expr, evaluate=False)

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
            return expr

        elif isinstance(expr, ImaginaryUnit):
            return expr
        # ...

        # Expression must always be translated to Sympy!
        # TODO: check if we should use 'sympy.sympify(expr)' instead
        else:
            raise NotImplementedError('Cannot translate to Sympy: {}'.format(expr))

