# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr,Mul, Add
from sympy.core.expr import AtomicExpr
from sympy.core.containers import Tuple
from sympy import Integer, Float, expand
from sympy import Function
from sympy.core.numbers import Zero as sy_Zero
from sympy.core.singleton import S

from sympde.core.utils import expand_name_patterns
from sympde.core.utils import random_string
from sympde.core.basic import CalculusFunction
from sympde.core.basic import _coeffs_registery
from .basic import BasicDomain, Union, Interval
from .datatype import SpaceType, dtype_space_registry
from .datatype import RegularityType, dtype_regularity_registry

#==============================================================================
def element_of(space, name):
    """ Create a single element of a given space (possibly a ProductSpace).

    Parameters
    ----------
    space : ScalarFunctionSpace | VectorFunctionSpace | ProductSpace
        Function space from which a single element should be created.

    name : str | iterable
        If space is ProductSpace, 'name' must either be an explicit list of
        function names, or a single pattern string that will be expanded into
        such a list. Otherwise, 'name' must be a simple string.

    Results
    -------
    res : ScalarTestFunction | VectorTestFunction | iterable
        Single element taken from the given space. If space is ProductSpace,
        an element is a list of functions; otherwise, it is a single function.

    """
    if not isinstance(space, BasicFunctionSpace):
        raise TypeError(space)
    names = expand_name_patterns(name)
    return _recursive_element_of(space, names)

#----------------------------------------
def _recursive_element_of(space, names):

    if isinstance(names, str):
        name = names
        return space.element(name)

    elif isinstance(space, ProductSpace):
        spaces = space.spaces
        result = [_recursive_element_of(s, n) for s, n in zip(spaces, names)]
        return type(names)(result)

    else:
        msg = "To create multiple elements of same space, use 'elements_of'."
        raise ValueError(msg)

#==============================================================================
def elements_of(space, names):
    """ Create multiple elements of same space (possibly a ProductSpace).

    Parameters
    ----------
    space : ScalarFunctionSpace | VectorFunctionSpace | ProductSpace

    names : str | iterable
        Pattern or list of patterns from which a list of function names is
        produced.

    Results
    -------
    res : iterable
        Multiple elements taken from the given space. If space is ProductSpace,
        each element is a list of functions; otherwise, each element is a
        single function.

    """
    if not isinstance(space, BasicFunctionSpace):
        raise TypeError(space)
    names = expand_name_patterns(names, seq=True)
    return _recursive_elements_of(space, names)

#----------------------------------------
def _recursive_elements_of(space, names):

    if isinstance(names, str):
        name = names
        return space.element(name)

    elif isinstance(space, ProductSpace):
        spaces = space.spaces
        result = [_recursive_elements_of(s, n) for s, n in zip(spaces, names)]
        return type(names)(result)

    else:
        result = [_recursive_elements_of(space, n) for n in names]
        return type(names)(result)

#==============================================================================
class BasicFunctionSpace(Basic):
    """
    Represents a basic continuous Function space.

    Examples

    """
    _domain     = None
    _shape      = None
    _kind       = None
    _regularity = None # TODO pass it as arg to __new__
    _is_broken  = None
    def __new__(cls, name, domain, shape, kind):

        if not isinstance(domain, BasicDomain):
            raise TypeError('> Expecting a BasicDomain object for domain')

        obj = Basic.__new__(cls)
        obj._name   = name
        obj._domain = domain
        obj._shape  = shape

        # ...
        if kind is None:
            kind = 'undefined'

        if isinstance(kind, str):
            kind_str = kind.lower()
            assert(kind_str in ['h1', 'hcurl', 'hdiv', 'l2', 'undefined'])

            kind = dtype_space_registry[kind_str]
        elif not isinstance(kind, SpaceType):
            raise TypeError('Expecting kind to be of SpaceType')

        kind_str  = kind.name
        obj._kind = kind

        # ...
        if not(kind_str == 'undefined'):
            obj._regularity = dtype_regularity_registry[kind_str]
        # ...

        # ...
        # TODO remove this if => bug in tensor form
        if isinstance(domain, Interval):
            is_broken = False

        else:
            is_broken = len(domain) > 1

        obj._is_broken = is_broken
        # ...

        return obj

    @property
    def name(self):
        return self._name

    @property
    def domain(self):
        return self._domain

    @property
    def ldim(self):
        return self.domain.dim

    @property
    def shape(self):
        return self._shape

    @property
    def kind(self):
        return self._kind

    @property
    def is_broken(self):
        return self._is_broken

    @property
    def regularity(self):
        return self._regularity

    @property
    def coordinates(self):
        return self.domain.coordinates

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

    def __mul__(self, other):
        return ProductSpace(self, other)

    def __hash__(self):
        return hash((self.name, self.domain, self.shape, self.kind))
#==============================================================================
class ScalarFunctionSpace(BasicFunctionSpace):
    """
    Represents a basic continuous scalar Function space.
    """
    def __new__(cls, name, domain, kind=None):
        shape = 1
        return BasicFunctionSpace.__new__(cls, name, domain, shape, kind)

    def element(self, name):
        return ScalarTestFunction(self, name)

    def field(self, name):
        return ScalarField(self, name)

#==============================================================================
class VectorFunctionSpace(BasicFunctionSpace):
    """
    Represents a basic continuous vector Function space.
    """
    def __new__(cls, name, domain, kind=None):
        shape = domain.dim
        return BasicFunctionSpace.__new__(cls, name, domain, shape, kind)

    def element(self, name):
        return VectorTestFunction(self, name)

    def field(self, name):
        return VectorField(self, name)

#=============================================================================
class Derham:
    """."""
    def __init__(self, domain, sequence=None):
        shape = domain.dim
        self._V0  = None
        self._V1  = None
        self._V2  = None
        self._V3  = None

        if shape == 1:
            spaces = [ScalarFunctionSpace('H1', domain, kind='H1'),
                      ScalarFunctionSpace('L2', domain, kind='L2')]

            self._V0  = spaces[0]
            self._V1  = spaces[1]

        elif shape == 2:
            assert sequence is not None

            space = sequence[1]
            spaces = [ScalarFunctionSpace('H1', domain, kind='H1'),
                      VectorFunctionSpace(space, domain, kind=space),
                      ScalarFunctionSpace('L2', domain, kind='L2')]

            self._V0  = spaces[0]
            self._V1  = spaces[1]
            self._V2  = spaces[2]

        elif shape == 3:
            spaces = [ScalarFunctionSpace('H1', domain, kind='H1'),
                      VectorFunctionSpace('Hcurl', domain, kind='Hcurl'),
                      VectorFunctionSpace('Hdiv', domain, kind='Hdiv'),
                      ScalarFunctionSpace('L2', domain, kind='L2')]

            self._V0  = spaces[0]
            self._V1  = spaces[1]
            self._V2  = spaces[2]
            self._V3  = spaces[3]


        self._spaces = spaces
        self._domain = domain
        self._shape  = shape


    # ...
    @property
    def spaces(self):
        return self._spaces

    @property
    def domain(self):
        return self._domain

    @property
    def shape(self):
        return self._shape

    @property
    def V0(self):
        return self._V0

    @property
    def V1(self):
        return self._V1

    @property
    def V2(self):
        return self._V2

    @property
    def V3(self):
        return self._V3

    @property
    def spaces(self):
        if self.shape == 1:
            return (self._V0, self._V1)

        elif self.shape == 2:
            return (self._V0, self._V1, self._V2)

        elif self.shape == 3:
            return (self._V0, self._V1, self._V2, self._V3)

#==============================================================================
# TODO must check that all spaces have the same domain
#     for the moment this class is not used
class ProductSpace(BasicFunctionSpace):
    """
    Represents a product of continuous Sobolev spaces.

    Examples

    """
    def __new__(cls, *spaces):

        # ...
        if not (isinstance(spaces, (tuple, list, Tuple))):
            raise TypeError('> Expecting a tuple, list or Tuple')
        # ...

        # ...
        args = []
        for V in spaces:
            if isinstance(V, ProductSpace):
                args += [W for W in V.spaces]

            else:
                args += [V]

        spaces = Tuple(*args)
        # ...

        # ... all spaces must have the same domain
        domain = spaces[0].domain

        #TODO uncomment
        #for space in spaces:
        #    if not(space.domain is domain):
        #        raise ValueError('> all spaces must have the same domain')

        ldim = domain.dim
        # ...

        # ...
        shape = sum([i.shape for i in spaces])
        # ...

        # ...
        name = ''.join(i.name for i in spaces)
        # ...

        # ...
        def _get_name(V):
            if V.ldim == 1:
                return V.coordinates.name
            else:
                return [i.name for i in V.coordinates]

        coordinates = unique([_get_name(i) for i in spaces])
        for i in coordinates:
            if not isinstance(i, str):
                raise TypeError('> Expecting a string')

        coordinates = spaces[0].coordinates
        if isinstance(coordinates, Symbol):
            coordinates = [coordinates]
        # ...

        # ...
        obj = Basic.__new__(cls, spaces)

        obj._shape = shape
        obj._coordinates = coordinates
        obj._name = name
        # ...

        return obj

    @property
    def spaces(self):
        return self._args[0]

    @property
    def name(self):
        return self._name

    @property
    def domain(self):
        return self.spaces[0].domain

    @property
    def ldim(self):
        return self.spaces[0].ldim

    @property
    def shape(self):
        return self._shape

    def element(self, name):
        raise NotImplementedError('TODO')

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

#==============================================================================
class ScalarTestFunction(Symbol):
    """
    Represents a test function as an element of a fem space.

    Examples

    >>> from sympde.codegen.core import SplineFemSpace
    >>> from sympde.codegen.core import ScalarTestFunction
    >>> V = SplineFemSpace('V')
    >>> phi = ScalarTestFunction(V, 'phi')
    """
    is_commutative = True
    _space         = None
    _projection_of = None

    def __new__(cls, space, name):
        if not isinstance(space, ScalarFunctionSpace):
            raise ValueError('Expecting a ScalarFunctionSpace')
        obj = Expr.__new__(cls)
        obj._space = space
        obj._name  = name
        return obj

    @property
    def space(self):
        return self._space

    @property
    def name(self):
        return self._name

    @property
    def ldim(self):
        return self.space.ldim

    @property
    def projection_of(self):
        return self._projection_of

    def duplicate(self, name):
        return ScalarTestFunction(self.space, name)

    def set_as_projection(self, expr):
        self._projection_of = expr

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

    def __hash__(self):
        return hash((self.name, self.space))
#==============================================================================
# this class is needed, otherwise sympy will convert VectorTestFunction to
# IndexedBase
class IndexedTestTrial(Indexed):
    """Represents a mathematical object with indices.

    """
    is_commutative = True
    is_Indexed = True
    is_symbol = True
    is_Atom = True

    def __new__(cls, base, *args, **kw_args):

        if isinstance(base, VectorTestFunction):
            pass

        elif isinstance(base, Add):
            return Add(*[cls(b, *args, **kw_args) for b in base.args])

        elif isinstance(base, Mul):
            scalar_types = (*_coeffs_registery, ScalarTestFunction)
            scalars = [s for s in base.args if isinstance(s, scalar_types)]
            others  = [s for s in base.args if s not in scalars]
            return Mul(*scalars) * Mul(*[cls(b, *args, **kw_args) for b in others])

        elif isinstance(base, Trace):
            expr     = cls(base.expr, *args, **kw_args)
            boundary = base.boundary
            order    = base.order
            return Trace(expr, boundary, order)

        else:
            raise ValueError('Expecting VectorTestFunction, Trace, or Add/Mul object')

        if not args:
            raise IndexException("Indexed needs at least one index.")

        return Expr.__new__(cls, base, *args, **kw_args)

    # free_symbols is redefined otherwise an expression u[0].free_symbols will
    # give the error:  AttributeError: 'int' object has no attribute 'free_symbols'
    @property
    def free_symbols(self):
        base_free_symbols = self.base.free_symbols
        symbolic_indices = [i for i in self.indices if isinstance(i, Basic)]
        if len(symbolic_indices) > 0:
            raise ValueError('symbolic indices not yet available')

        return base_free_symbols

    @property
    def ldim(self):
        return self.base.space.ldim

    def __hash__(self):
        return hash(self._args)

#==============================================================================
class VectorTestFunction(Symbol, IndexedBase):
    """
    Represents a vector test function as an element of a fem space.

    Examples

    """
    is_commutative = False
    _space         = None
    _projection_of = None

    def __new__(cls, space, name):
        if not isinstance(space, VectorFunctionSpace):
            raise ValueError('Expecting a VectorFunctionSpace')
        obj        = Expr.__new__(cls)
        obj._space = space
        obj._name = name
        return obj

    @property
    def space(self):
        return self._space

    @property
    def name(self):
        return self._name

    @property
    def shape(self):
        # we return a list to make it compatible with IndexedBase sympy object
        return [self.space.shape]

    @property
    def ldim(self):
        return self.space.ldim

    @property
    def projection_of(self):
        return self._projection_of

    def __getitem__(self, *args):

        if self.shape and len(self.shape) != len(args):
            raise IndexException("Rank mismatch.")

        if not(len(args) == 1):
            raise ValueError('expecting exactly one argument')

        args = list(args)
        for i in range(len(args)):
            if isinstance(args[i], int):
                args[i] = Integer(args[i])

        obj = IndexedTestTrial(self, *args)
        return obj

    def duplicate(self, name):
        return VectorTestFunction(self.space, name)

    def set_as_projection(self, expr):
        self._projection_of = expr

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

    def __hash__(self):
        return hash((self.name, self.space))

#==============================================================================
# this is implemented as a function, it would be better to have it as a class
def TestFunction(space, name=None):

    if isinstance(space,ScalarFunctionSpace):
        return ScalarTestFunction(space, name=name)

    elif isinstance(space, VectorFunctionSpace):
        return VectorTestFunction(space, name=name)

    elif isinstance(space, ProductSpace):

        if not(name is None):
            assert(isinstance(name, (tuple, list, Tuple)))
            assert(len(name) == len(space.spaces))

        else:
            name = [None for V in space.spaces]

        args = []
        for V, n in zip(space.spaces, name):
            args += [TestFunction(V, name=n)]

        return Tuple(*args)

    else:
        raise TypeError('Wrong space type. given {}'.format(type(space)))

#==============================================================================
class ScalarField(Symbol):
    """
    Represents a ScalarField variable.

    Examples

    """
    _space         = None
    is_commutative = True
    _projection_of = None
    def __new__(cls, space, name=None):
        if not isinstance(space, ScalarFunctionSpace):
            raise ValueError('Expecting a ScalarFunctionSpace')

        obj = Symbol.__new__(cls, name)
        obj._space = space
        return obj

    @property
    def space(self):
        return self._space

    @property
    def ldim(self):
        return self.space.ldim

    @property
    def projection_of(self):
        return self._projection_of

    def set_as_projection(self, expr):
        self._projection_of = expr

    def __eq__(self, a):
        if isinstance(a, ScalarField):
            return a.space == self.space and a.name == self.name
        return False

    def __hash__(self):
        return hash((self.name, self.space))

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

#==============================================================================
class VectorField(Symbol, IndexedBase):
    """
    Represents a vector field as an element of a fem space.

    Examples

    """
    is_commutative = False
    _space = None
    _projection_of = None
    def __new__(cls, space, name=None):
        if not isinstance(space, VectorFunctionSpace):
            raise ValueError('Expecting a VectorFunctionSpace')

        obj = Basic.__new__(cls, name)
        obj._space = space
        return obj

    @property
    def space(self):
        return self._space

    @property
    def name(self):
        return self._args[0]

    @property
    def shape(self):
        # we return a list to make it compatible with IndexedBase sympy object
        return [self.space.shape]

    @property
    def ldim(self):
        return self.space.ldim

    def __getitem__(self, *args):

        if self.shape and len(self.shape) != len(args):
            raise IndexException("Rank mismatch.")

        if not(len(args) == 1):
            raise ValueError('expecting exactly one argument')

        args = list(args)
        for i in range(len(args)):
            if isinstance(args[i], int):
                args[i] = Integer(args[i])

        assumptions ={}
        obj = IndexedVectorField(self, *args)
        return obj

    def duplicate(self, name):
        return VectorField(self.space, name)

    @property
    def projection_of(self):
        return self._projection_of

    def set_as_projection(self, expr):
        self._projection_of = expr

#==============================================================================
# this class is needed, otherwise sympy will convert VectorTestFunction to
# IndexedBase
class IndexedVectorField(Indexed):
    """Represents a mathematical object with indices.

    """
    is_commutative = True
    is_Indexed = True
    is_symbol = True
    is_Atom = True

    def __new__(cls, base, *args, **kw_args):
        assert(isinstance(base, VectorField))

        if not args:
            raise IndexException("Indexed needs at least one index.")

        return Expr.__new__(cls, base, *args, **kw_args)

    # free_symbols is redefined otherwise an expression u[0].free_symbols will
    # give the error:  AttributeError: 'int' object has no attribute 'free_symbols'
    @property
    def free_symbols(self):
        base_free_symbols = self.base.free_symbols
        symbolic_indices = [i for i in self.indices if isinstance(i, Basic)]
        if len(symbolic_indices) > 0:
            raise ValueError('symbolic indices not yet available')

        return base_free_symbols

    @property
    def ldim(self):
        return self.base.space.ldim

#==============================================================================
# TODO Expr or AtomicExpr?
class Trace(AtomicExpr):
    """
    Represents the trace over a boundary and a space function

    """
    is_commutative = True
    def __new__(cls, expr, boundary, order=0, **options):
#        # TODO these tests are not working for the moment for Grad(u)
#        if not expr.atoms((ScalarTestFunction, VectorTestFunction, ScalarField)):
#            raise TypeError('> Wrong type for expr')
#
#        if not(expr.space.domain is boundary.domain):
#            raise ValueError('> Space and boundary domains must be the same')
        evaluate = options.pop('evaluate',True)
        if evaluate:
            return cls.eval(expr, boundary, order)
        return Basic.__new__(cls, expr, boundary, order)

    @property
    def expr(self):
        return self._args[0]

    @property
    def boundary(self):
        return self._args[1]

    @property
    def order(self):
        return self._args[2]

    @classmethod
    def eval(cls, expr, boundary, order):

        if not isinstance(expr, Tuple):
            expr = expand(expr)

        if not isinstance(expr, (Expr, Tuple)):
            raise TypeError('only Expr are accepted')

        if isinstance(expr, sy_Zero):
            return sy_Zero

        if isinstance(expr, Add):
            args = [cls.eval(a, boundary, order) for a in expr.args]
            return expr._new_rawargs(*args)

        if isinstance(expr, Mul):
            args = expr.args
            coeffs = [a for a in args if isinstance(a, _coeffs_registery)]
            a      = Mul(*coeffs)
            args   = [a for a in args if not(a in coeffs)]

            funcs  = [a for a in args if isinstance(a, Function) and not isinstance(a, CalculusFunction)]
            a      = a*Mul(*funcs)
            args   = [a for a in args if not(a in funcs)]

            b      = cls(Mul(*args), boundary, order, evaluate=False)

            return a*b

        if isinstance(boundary, Union):
            expr = [Integral.eval(expr, d, order) for d in boundary.args]
            return Add(*expr)

        return cls(expr, boundary, order, evaluate=False)

#==============================================================================
# ... user friendly functions
trace_0 = lambda x, B: Trace(x, B, order=0)
trace_1 = lambda x, B: Trace(x, B, order=1)

_is_sympde_atom = lambda a: isinstance(a, (ScalarTestFunction, VectorTestFunction,
                                                   ScalarField, VectorField))

_is_test_function = lambda a: isinstance(a, (ScalarTestFunction, VectorTestFunction))
_is_field         = lambda a: isinstance(a, (ScalarField, VectorField))


#==============================================================================
class Projector(Basic):
    """
    Represents a Projector over a function space.

    Examples

    """
    _kind = None
    def __new__(cls, space, kind=None):

        if not isinstance(space, BasicFunctionSpace):
            raise TypeError('> Expecting a BasicFunctionSpace object for space')

        obj = Basic.__new__(cls, space)
        obj._kind = kind

        return obj

    @property
    def space(self):
        return self._args[0]

    @property
    def kind(self):
        return self._kind

    def __call__(self, expr):
        V = self.space

        if isinstance(expr, (ScalarField, VectorField, ScalarTestFunction, VectorTestFunction)):
            if expr.space is V:
                return expr

            elif isinstance(expr.projection_of, Projection):
                if expr.projection_of.projector.space is V:
                    return expr

        if isinstance(V, (ScalarFunctionSpace, VectorFunctionSpace)):
            name = 'Proj_' + random_string( 4 )
            F = element_of(V, name)
        else:
            raise TypeError('Only scalar and vector space are handled')

        F.set_as_projection(Projection(self, expr))
        return F

#==============================================================================
class Projection(AtomicExpr):
    """
    Represents a projection

    Examples

    """
    _kind = None
    def __new__(cls, projector, expr):

        return Basic.__new__(cls, projector, expr)

    @property
    def projector(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]
