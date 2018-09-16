# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple

from .basic import Field
from .geometry import BasicDomain

# ...
class BasicFunctionSpace(Basic):
    """
    Represents a basic continuous Function space.

    Examples

    """
    _domain = None
    _shape = None
    def __new__(cls, name, domain, shape):

        if not isinstance(domain, BasicDomain):
            raise TypeError('> Expecting a BasicDomain object for domain')

        obj = Basic.__new__(cls, name, domain)
        obj._shape = shape

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
    def coordinates(self):
        return self.domain.coordinates

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)
# ...

# ...
class FunctionSpace(BasicFunctionSpace):
    """
    Represents a basic continuous scalar Function space.
    """
    def __new__(cls, name, domain):
        shape = 1
        return BasicFunctionSpace.__new__(cls, name, domain, shape)

class VectorFunctionSpace(BasicFunctionSpace):
    """
    Represents a basic continuous vector Function space.
    """
    def __new__(cls, name, domain):
        shape = domain.dim
        return BasicFunctionSpace.__new__(cls, name, domain, shape)
# ...

# ...
# TODO must check that all spaces have the same domain
#     for the moment this class is not used
class ProductSpace(FunctionSpace):
    """
    Represents a product of continuous Sobolev spaces.

    Examples

    """
    def __new__(cls, *spaces):

        # ...
        if not (isinstance(spaces, (tuple, list, Tuple))):
            raise TypeError('> Expecting a tuple, list or Tuple')

        spaces = Tuple(*spaces)
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

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)
# ...

class TestFunction(Symbol):
    """
    Represents a test function as an element of a fem space.

    Examples

    >>> from sympde.codegen.core import SplineFemSpace
    >>> from sympde.codegen.core import TestFunction
    >>> V = SplineFemSpace('V')
    >>> phi = TestFunction(V, 'phi')
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
        return TestFunction(self.space, name)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)


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


class Unknown(TestFunction):
    """
    Represents an unknown function

    """
    def __new__(cls, name, domain):
        space_name = 'space_{}'.format(abs(hash(name)))
        V = FunctionSpace(space_name, domain)
        return TestFunction.__new__(cls, V, name)


class VectorUnknown(VectorTestFunction):
    """
    Represents an unknown function

    """
    def __new__(cls, name, domain, shape):
        space_name = 'space_{}'.format(abs(hash(name)))
        V = VectorFunctionSpace(space_name, domain)
        return VectorTestFunction.__new__(cls, V, name)


class Trace(Basic):
    """
    Represents the trace over a boundary and a space function

    """
    def __new__(cls, expr, boundary, order=0):
#        # TODO these tests are not working for the moment for Grad(u)
#        if not expr.atoms((TestFunction, VectorTestFunction, Field)):
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

# ... user friendly functions
trace_0 = lambda x, B: Trace(x, B, order=0)
trace_1 = lambda x, B: Trace(x, B, order=1)
