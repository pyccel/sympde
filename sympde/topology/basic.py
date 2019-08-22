# coding: utf-8


from collections import OrderedDict
from collections import abc

from sympy.core import Basic, Symbol
from sympy.core.containers import Tuple
from sympy.tensor import IndexedBase

#==============================================================================
class BasicDomain(Basic):
    _dim         = None
    _name        = None
    _coordinates = None

    @property
    def name(self):
        return self._name

    @property
    def dim(self):
        return self._dim

    @property
    def coordinates(self):
        dim = self.dim
        if self._coordinates is None:
            xyz = ['x', 'y', 'z'][:dim]
            xyz = [Symbol(i) for i in xyz]
            self._coordinates = xyz

        if dim == 1:
            return self._coordinates[0]
        else:
            return self._coordinates

#==============================================================================
class InteriorDomain(BasicDomain):
    """
    Represents an undefined interior domain.

    Examples

    """
    def __new__(cls, name, dim=None):
        target = None
        if not isinstance(name, str):
            target = name
            name   = name.name

        if not( target is None ):
            if isinstance(target, ProductDomain):
                dim = target.dim

        obj = Basic.__new__(cls, name)

        obj._dim    = dim
        obj._target = target

        return obj

    @property
    def name(self):
        return self._args[0]

    @property
    def target(self):
        return self._target

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

    def todict(self):
        name   = str(self.name)
        d = {'name': name}

        return OrderedDict(sorted(d.items()))


#==============================================================================
# TODO remove redundancy
class Union(BasicDomain):
    def __new__(cls, *args):
        args = Tuple(*args)
        if not all( [isinstance(i, BasicDomain) for i in args] ):
            raise TypeError('arguments must be of BasicDomain type')

        assert(len(args) > 1)

        dim = args[0].dim
        dims = [a.dim for a in args[1:]]
        #TODO
     #   if not all( [d == dim for d in dims]):
     #       raise ValueError('arguments must have the same dimension')

        # ...
        unions = [a for a in args if isinstance(a, Union)]
        args   = [a for a in args if not isinstance(a, Union)]
        for union in unions:
            args += list(union.as_tuple())
        # ...

        # sort domains by name
        args = sorted(args, key=lambda x: x.name)

        return Basic.__new__(cls, *args)

    @property
    def dim(self):
        return self._args[0].dim

    def __len__(self):
        return len(self._args)

    def complement(self, arg):
        if isinstance(arg, Union):
            arg = arg._args

        elif isinstance(arg, BasicDomain):
            arg = [arg]

        ls = [i for i in self._args if not(i in arg)]
        if len(ls) > 1:
            return Union(*ls)
        else:
            return ls[0]

    def __sub__(self, other):
        return self.complement(other)

    def todict(self):
        return [i.todict() for i in self._args]

    def as_tuple(self):
        ls = [i for i in self._args]
        return tuple(ls)

#==============================================================================
class ProductDomain(BasicDomain):
    def __new__(cls, *args, name=None):
        args = Tuple(*args)
        if not all( [isinstance(i, BasicDomain) for i in args] ):
            raise TypeError('arguments must be of BasicDomain type')

        assert(len(args) > 1)

        obj = Basic.__new__(cls, *args)
        obj._dim = sum(i.dim for i in args)
        obj._name = name
        return obj

    @property
    def domains(self):
        return self._args

#==============================================================================
class Interval(InteriorDomain):
    """
    Represents a 1D interval.

    Examples

    """
    _dim = 1
    def __new__(cls, name=None, coordinate=None, bounds=None):
        if name is None:
            name = 'Interval'

        if bounds is None:
            bounds = (0, 1)

        obj = Basic.__new__(cls, name)
        if coordinate:
            obj._coordinates = [coordinate]

        obj._bounds = bounds

        return obj

    @property
    def name(self):
        return self._args[0]

    @property
    def bounds(self):
        return self._bounds

#==============================================================================
class Boundary(BasicDomain):
    """
    Represents an undefined boundary over a domain.

    Examples

    """
    def __new__(cls, name, domain, axis=None, ext=None):

        if not( axis is None ):
            assert(isinstance(axis, int))

        if not( ext is None ):
            assert(isinstance(ext, int))

        obj = Basic.__new__(cls, name)
        obj._domain = domain
        obj._axis = axis
        obj._ext = ext
        return obj

    @property
    def name(self):
        return self._args[0]

    @property
    def domain(self):
        return self._domain

    @property
    def dim(self):
        return self.domain.dim

    @property
    def axis(self):
        return self._axis

    @property
    def ext(self):
        return self._ext

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

    def __add__(self, other):
        if isinstance(other, ComplementBoundary):
            raise TypeError('> Cannot add complement of boundary')

        return Union(self, other)

    def todict(self):
        domain = str(self.domain.name)
        name   = str(self.name)
        axis   = str(self.axis)
        ext    = str(self.ext)

        d = {'patch': domain,
             'name':  name,
             'axis':  axis,
             'ext':   ext}

        return OrderedDict(sorted(d.items()))

    def __hash__(self):
        return hash(self.name)

    def __lt__(self, other):
        #add this method to avoid sympy error in Basic.compare
        return 0

    def __eq__(self, other):
        if not isinstance(other, Boundary):
            return False

        return ( self.name == other.name ) and ( self.domain is other.domain )

#==============================================================================
class Interface(BasicDomain):
    """
    Represents an interface between two subdomains through two boundaries.

    Examples

    """
    def __new__(cls, edge, bnd_minus, bnd_plus):
        assert(isinstance(edge, Edge))
        assert(isinstance(bnd_minus, Boundary))
        assert(isinstance(bnd_plus, Boundary))

        return Basic.__new__(cls, edge.name, bnd_minus, bnd_plus)

    @property
    def name(self):
        return self._args[0]

    @property
    def minus(self):
        return self._args[1]

    @property
    def plus(self):
        return self._args[2]

    def _sympystr(self, printer):
        sstr = printer.doprint

#        name  = self.name
#        minus = self.minus
#        plus  = self.plus
#        minus = '{domain}.{bnd}'.format( domain = sstr(minus.domain),
#                                         bnd    = sstr(minus) )
#        plus = '{domain}.{bnd}'.format( domain = sstr(plus.domain),
#                                         bnd    = sstr(plus) )
#        pattern = 'Interface( {name}; {minus}, {plus} )'
#        return pattern.format( name  = sstr(self.name),
#                               minus = minus,
#                               plus = plus )

        return '{}'.format(sstr(self.name))


#==============================================================================
class Edge(object):
    def __init__(self, name):
        self._name = name

    @property
    def name(self):
        return self._name

    def __lt__(self, other):
        return self.name.__lt__(other.name)


#==============================================================================
class Connectivity(abc.Mapping):
    _patches  = []

    def __init__(self, data=None):
        if data is None:
            data = {}

        else:
            assert( isinstance( data, (dict, OrderedDict)) )

        self._data = data

    @property
    def patches(self):
        return self._patches

    @property
    def interfaces(self):
        ls = []
        data = OrderedDict(sorted(self._data.items()))
        for k, (minus, plus) in data.items():
            ls.append(Interface(k, minus, plus))

        if len(ls) == 1:
            return ls[0]
        else:
            return Union(*ls)

    def todict(self):
        # ... create the connectivity
        connectivity = {}
        data = OrderedDict(sorted(self._data.items()))
        for edge, pair in data.items():
            connectivity[edge.name] = [bnd.todict() for bnd in pair]

        connectivity = OrderedDict(sorted(connectivity.items()))
        # ...

        return connectivity

    def __setitem__(self, key, value):
        if isinstance(key, str):
            key = Edge(key)

        assert( isinstance( key, Edge ) )
        assert( isinstance( value, (tuple, list)  ) )
        assert( len(value) in [1, 2] )
        assert( all( [isinstance( P, Boundary ) for P in value ] ) )

        self._data[key] = value

    # ==========================================
    #  abstract methods
    # ==========================================
    def __getitem__(self, key):
        return self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __hash__(self):
        return hash(tuple(self._data.values()))

    def __lt__(self, other):
        #add this method to avoid sympy error in Basic.compare
        return 0

    # ==========================================
