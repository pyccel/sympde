# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.expr import AtomicExpr
from sympy.core.containers import Tuple

from sympde.core.utils import random_string
from .basic import BasicDomain
from .datatype import SpaceType, dtype_space_registry
from .datatype import RegularityType, dtype_regularity_registry


def element_of_space(space, name):
    assert isinstance(space, BasicFunctionSpace)
    if isinstance(name, (list,tuple,Tuple)):
        return [Element(space,nm) for nm in name]
    return Element(space, name)

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
    def __new__(cls, name, domain, shape, kind):

        if not isinstance(domain, BasicDomain):
            raise TypeError('> Expecting a BasicDomain object for domain')

        obj = Basic.__new__(cls, name, domain)
        obj._shape = shape

        # ...
        if kind is None:
            kind = 'undefined'

        assert(isinstance(kind, (str, SpaceType)))

        kind_str = kind.lower()
        assert(kind_str in ['h1', 'hcurl', 'hdiv', 'l2', 'undefined'])

        kind = dtype_space_registry[kind_str]

        if not isinstance(kind, SpaceType):
            raise TypeError('Expecting kind to be of SpaceType')

        obj._kind = kind
        # ...

        # ...
        if not(kind_str == 'undefined'):
            obj._regularity = dtype_regularity_registry[kind_str]
        # ...

        return obj

    @property
    def name(self):
        return self._args[0]

    @property
    def domain(self):
        return self._args[1]

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

#==============================================================================
class FunctionSpace(BasicFunctionSpace):
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
            spaces = [FunctionSpace('H1', domain, kind='H1'),
                        FunctionSpace('L2', domain, kind='L2')]
            
            self._V0  = spaces[0]
            self._V1  = spaces[1]

        elif shape == 2:
            assert sequence is not None
            
            space = sequence[1]
            spaces = [FunctionSpace('H1', domain, kind='H1'),
                        VectorFunctionSpace(space, domain, kind=space),
                        FunctionSpace('L2', domain, kind='L2')]

            self._V0  = spaces[0]
            self._V1  = spaces[1]
            self._V2  = spaces[2]
                        
        elif shape == 3:
            spaces = [FunctionSpace('H1', domain, kind='H1'),
                        VectorFunctionSpace('Hcurl', domain, kind='Hcurl'),
                        VectorFunctionSpace('Hdiv', domain, kind='Hcurl'),
                        FunctionSpace('L2', domain, kind='L2')]
                   
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
        for space in spaces:
            if not(space.domain is domain):
                raise ValueError('> all spaces must have the same domain')

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
    _space = None
    is_commutative = True
    def __new__(cls, space, name=None):
        obj = Symbol.__new__(cls, name)
        obj._space = space
        return obj

    @property
    def space(self):
        return self._space

    @property
    def ldim(self):
        return self.space.ldim

    def duplicate(self, name):
        return ScalarTestFunction(self.space, name)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)


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
        assert(isinstance(base, VectorTestFunction))

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
class VectorTestFunction(Symbol, IndexedBase):
    """
    Represents a vector test function as an element of a fem space.

    Examples

    """
    is_commutative = True
    _space = None
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

        assumptions ={}
        obj = IndexedTestTrial(self, *args)
        return obj

    def duplicate(self, name):
        return VectorTestFunction(self.space, name)


#==============================================================================
# this is implemented as a function, it would be better to have it as a class
def TestFunction(space, name=None):

    if isinstance(space, FunctionSpace):
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
# TODO not tested yet
# this is implemented as a function, it would be better to have it as a class
def Field(space, name=None):

    if isinstance(space, FunctionSpace):
        return ScalarField(space, name=name)

    elif isinstance(space, VectorFunctionSpace):
        return VectorField(space, name=name)

    elif isinstance(space, ProductSpace):

        if not(name is None):
            assert(isinstance(name, (tuple, list, Tuple)))
            assert(len(name) == len(space.spaces))

        else:
            name = [None for V in space.spaces]

        args = []
        for V, n in zip(space.spaces, name):
            args += [Field(V, name=n)]

        return Tuple(*args)

    else:
        raise TypeError('Wrong space type. given {}'.format(type(space)))


class Element(Symbol):
    """
    Represents an element of a fem space.

    Examples

    >>> from sympde.codegen.core import SplineFemSpace
    >>> from sympde.codegen.core import ScalarTestFunction
    >>> V = SplineFemSpace('V')
    >>> phi = Element(V, 'phi')
    """
    _space = None
    is_commutative = True
    def __new__(cls, space, name=None):
        obj = Symbol.__new__(cls, name)
        obj._space = space
        obj._iterable = False
        return obj

    @property
    def space(self):
        return self._space

    @property
    def ldim(self):
        return self.space.ldim
        
    @property
    def shape(self):
        # we return a list to make it compatible with IndexedBase sympy object
        return [self.space.shape]


    def duplicate(self, name):
        return Element(self.space, name)
        
    def __getitem__(self, *args):
        if self.shape and len(self.shape) != len(args):
            raise IndexException("Rank mismatch.")

        if not(len(args) == 1):
            raise ValueError('expecting exactly one argument')

        assumptions ={}
        obj = IndexedElement(self, *args)
        return obj

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)
    
class IndexedElement(Indexed):
    """Represents a mathematical object with indices.

    """
    is_commutative = True
    is_Indexed = True
    is_symbol = True
    is_Atom = True

    def __new__(cls, base, *args, **kw_args):
        if isinstance(base, (VectorTestFunction, VectorField)):
            # this is a hack when sympy do the substitution 
            # we return the right object
            return base.__getitem__(*args)
            
        assert isinstance(base, Element)

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
# TODO to be removed
class Unknown(ScalarTestFunction):
    """
    Represents an unknown function

    """
    def __new__(cls, name, domain):
        space_name = 'space_{}'.format(abs(hash(name)))
        V = FunctionSpace(space_name, domain)
        return ScalarTestFunction.__new__(cls, V, name)


#==============================================================================
# TODO to be removed
class VectorUnknown(VectorTestFunction):
    """
    Represents an unknown function

    """
    def __new__(cls, name, domain, shape):
        space_name = 'space_{}'.format(abs(hash(name)))
        V = VectorFunctionSpace(space_name, domain)
        return VectorTestFunction.__new__(cls, V, name)

#==============================================================================
class ScalarField(Symbol):
    """
    Represents a ScalarField variable.

    Examples

    """
    _space = None
    is_commutative = True
    _projection_of = None
    def __new__(cls, space, name=None):
        if not isinstance(space, FunctionSpace):
            raise ValueError('Expecting a FunctionSpace')

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

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

#==============================================================================
class VectorField(Symbol, IndexedBase):
    """
    Represents a vector field as an element of a fem space.

    Examples

    """
    is_commutative = True
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
class Trace(Expr):
    """
    Represents the trace over a boundary and a space function

    """
    def __new__(cls, expr, boundary, order=0):
#        # TODO these tests are not working for the moment for Grad(u)
#        if not expr.atoms((ScalarTestFunction, VectorTestFunction, ScalarField)):
#            raise TypeError('> Wrong type for expr')
#
#        if not(expr.space.domain is boundary.domain):
#            raise ValueError('> Space and boundary domains must be the same')

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

        name = 'Proj_' + random_string( 4 )
        if isinstance(V, FunctionSpace):
            F = ScalarField(V, name)

        elif isinstance(V, VectorFunctionSpace):
            F = VectorField(V, name)

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
