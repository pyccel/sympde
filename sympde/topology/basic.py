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
            if self.mapping is None:
                xyz = ['x1', 'x2', 'x3'][:dim]
            else:
                xyz = ['x', 'y', 'z'][:dim]

            xyz = [Symbol(i) for i in xyz]
            self._coordinates = xyz

        if dim == 1:
            return self._coordinates[0]
        else:
            return self._coordinates

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

#==============================================================================
class InteriorDomain(BasicDomain):
    """
    Represents an undefined interior domain.

    Examples

    """
    def __new__(cls, name, dim=None, dtype=None, mapping=None, logical_domain=None):
        target = None
        if not isinstance(name, str):
            target = name
            name   = name.name

        if not( target is None ):
            dim = target.dim

        assert mapping is None and logical_domain is None or \
        mapping is not None and logical_domain  is not None

        obj = Basic.__new__(cls, name)

        obj._dim            = dim
        obj._target         = target
        obj._dtype          = dtype
        obj._mapping        = mapping
        obj._logical_domain = logical_domain

        return obj

    @property
    def name(self):
        return self.args[0]

    @property
    def target(self):
        return self._target

    @property
    def mapping(self):
        return self._mapping

    @property
    def logical_domain(self):
        return self._logical_domain

    @property
    def dtype(self):
        return self._dtype

    @property
    def dim(self):
        return self._dim

    def _sympystr(self, printer):
        sstr = printer.doprint
        if self.mapping:
            return '{}({})'.format(sstr(self.mapping.name), sstr(self.name))
        else:
            return '{}'.format(sstr(self.name))

    def todict(self):
        name   = str(self.name)
        d = {'name': name}

        return OrderedDict(sorted(d.items()))


#==============================================================================
# TODO remove redundancy
class Union(BasicDomain):

    def __new__(cls, *args):

        # Discard empty Unions (represented as None) from args
        args = Tuple(*[a for a in args if a is not None])

        # Verify types
        if not all(isinstance(a, BasicDomain) for a in args):
            raise TypeError('arguments must be of BasicDomain type')

        # Verify dimensionality
        if len({a.dim for a in args}) > 1:
            dims = ', '.join(str(a.dim) for a in args)
            msg  = 'arguments must have the same dimension, '\
                   'given [{}] instead'.format(dims)
            raise ValueError(msg)

        # Flatten arguments into a single list of domains
        unions = [a for a in args if     isinstance(a, Union)]
        args   = [a for a in args if not isinstance(a, Union)]
        for union in unions:
            args += list(union.as_tuple())

        # remove duplicates and sort domains by their string representation
        args = sorted(set(args), key=str)

        # a. If the required Union contains no domains, return None;
        # b. If it contains a single domain, return the domain itself;
        # c. If it contains multiple domains, create a Union object.
        if not args:
            obj = None
        elif len(args) == 1:
            obj = args[0]
        else:
            obj       = Basic.__new__(cls, *args)
            obj.index = 0
        return obj

    @property
    def dim(self):
        return self.args[0].dim

    def __len__(self):
        return len(self.args)

    @property
    def coordinates(self):
        coords = self.args[0].coordinates
        assert all(e.coordinates == coords for e in self)
        return coords

    def complement(self, arg):
        if isinstance(arg, Union):
            arg = arg.args
        elif isinstance(arg, BasicDomain):
            arg = [arg]
        elif arg is None:
            return self
        else:
            TypeError('Invalid argument {}'.format(arg))

        return Union(*[i for i in self.args if (i not in arg)])

    def __sub__(self, other):
        return self.complement(other)

    def todict(self):
        return [i.todict() for i in self.args]

    def as_tuple(self):
        ls = [i for i in self.args]
        return tuple(ls)

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        try:
            result = self.args[self.index]
        except IndexError:
            raise StopIteration
        self.index += 1
        return result

    def _sympystr(self, printer):
        sstr = printer.doprint
        args = ', '.join(sstr(a) for a in self.args)
        return 'Union({})'.format(args)

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
        return self.args

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
        return self.args[0]

    @property
    def bounds(self):
        return self._bounds


#==============================================================================
class Boundary(BasicDomain):
    """
    Represents an undefined boundary over a domain.

    Examples

    """
    def __new__(cls, name, domain, axis=None, ext=None, mapping=None, logical_domain=None):

        if axis is not None:
            assert isinstance(axis, int)

        if ext is not None:
            assert isinstance(ext, int)

        obj                 = Basic.__new__(cls, name, domain, axis, ext)
        obj._mapping        = mapping
        obj._logical_domain = logical_domain

        return obj

    @property
    def name(self):
        return self.args[0]

    @property
    def domain(self):
        return self.args[1]

    @property
    def axis(self):
        return self.args[2]

    @property
    def ext(self):
        return self.args[3]

    @property
    def mapping(self):
        return self._mapping

    @property
    def logical_domain(self):
        return self._logical_domain

    @property
    def dim(self):
        return self.domain.dim

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}_{}'.format(sstr(self.domain),sstr(self.name))

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

#==============================================================================
class Interface(BasicDomain):
    """
    Represents an interface between two subdomains through two boundaries.

    Examples

    """
    def __new__(cls, name, bnd_minus, bnd_plus, mapping=None, logical_domain=None):

        if not isinstance(name     , str    ): raise TypeError(name)
        if not isinstance(bnd_minus, Boundary): raise TypeError(bnd_minus)
        if not isinstance(bnd_plus , Boundary): raise TypeError(bnd_plus)

        if bnd_minus.dim != bnd_plus.dim:
            raise TypeError('Dimension mismatch: {} != {}'.format(
                bnd_minus.dim, bnd_plus.dim))

        assert mapping is None and logical_domain is None or\
               mapping is not None and logical_domain is not None

        assert bnd_minus.axis == bnd_plus.axis
        obj = Basic.__new__(cls, name, bnd_minus, bnd_plus)
        obj._mapping        = mapping
        obj._logical_domain = logical_domain
        return obj

    @property
    def dim(self):
        return self.minus.dim

    @property
    def name(self):
        return self.args[0]

    @property
    def minus(self):
        return self.args[1]

    @property
    def plus(self):
        return self.args[2]

    @property
    def axis(self):
        return self.plus.axis

    @property
    def mapping(self):
        return self._mapping

    @property
    def logical_domain(self):
        return self._logical_domain

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
            for k,v in data.items():
                assert( isinstance( k, str ) )
                assert( isinstance(v, Interface) )
        self._data = data

    @property
    def patches(self):
        return self._patches

    @property
    def interfaces(self):
        ls = []
        data = OrderedDict(sorted(self._data.items()))
        for _,v in data.items():
            ls.append(v)
        return Union(*ls)

    def todict(self):
        # ... create the connectivity
        connectivity = {}
        data = OrderedDict(sorted(self._data.items()))
        for name, v in data.items():
            connectivity[name] = [v.minus.todict(), v.plus.todict()]
        connectivity = OrderedDict(sorted(connectivity.items()))
        # ...

        return connectivity

    def __setitem__(self, key, value):

        assert( isinstance( key, str ) )
        assert( isinstance(value, Interface) )
        assert( str(value.name) == key )

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


