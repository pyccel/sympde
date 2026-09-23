# coding: utf-8



from collections import abc
from numbers import Integral

from sympy.core import Basic, Symbol, Expr, Integer
from sympy.core.containers import Tuple
from sympy.tensor import IndexedBase


def _as_integer(value, name):
    """Return an integer input without silently truncating other scalars."""
    if isinstance(value, bool) or not isinstance(value, (Integral, Integer)):
        raise TypeError(f'{name} must be an integer')
    return int(value)

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

            xyz = [Symbol(i, real=True) for i in xyz]
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
        return '{}'.format(sstr(self.name))

    def todict(self):
        return {'name': str(self.logical_domain.name if self.logical_domain else self.name ),
                'mapping':str(self.mapping.name if self.mapping else None)}


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
            axis = _as_integer(axis, 'boundary axis')
            if not 0 <= axis < domain.dim:
                raise ValueError(
                    f'boundary axis must be between 0 and {domain.dim - 1}')

        if ext is not None:
            ext = _as_integer(ext, 'boundary extremity')
            if ext not in (-1, 1):
                raise ValueError('boundary extremity must be either -1 or 1')

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
    def normal_axis(self):
        """Logical axis normal to this patch side.

        This is the descriptive multipatch alias of :attr:`axis`.
        """
        return self.axis

    @property
    def ext(self):
        return self.args[3]

    @property
    def patch(self):
        """Patch interior owning this side; an alias of :attr:`domain`."""
        return self.domain

    @property
    def patch_dim(self):
        """Logical dimension of the patch owning this side."""
        return self.dim

    @property
    def intrinsic_dim(self):
        """Intrinsic dimension of this codimension-one side."""
        return self.dim - 1

    @property
    def mapping(self):
        return self._mapping

    @property
    def logical_domain(self):
        return self._logical_domain

    @property
    def dim(self):
        return self.domain.dim

    @property
    def adjacent_boundaries(self):
        boundaries = [a for a in self.domain.boundary if a.axis !=self.axis]
        return Union(*boundaries)

    def rotate(self, *directions):
        assert len(directions) == self.dim-1

        if self.dim == 2:
            if directions[0] == 1:
                return self
            elif directions[0] == -1:
                return self.domain.get_boundary(axis=self.axis, ext=-self.ext)

            else:
                raise TypeError('must be int')

        raise NotImplementedError('only 2d case is available')

    def join(self, boundary, orientation, *, name=None, logical_name=None):
        if not isinstance(boundary, Boundary):
            raise TypeError('boundary must be a Boundary')

        if self.mapping and boundary.mapping:
            # Imported lazily to keep the basic topology layer independent of
            # the mapping/domain import cycle.
            from sympde.topology.mapping import InterfaceMapping
            int_map            = InterfaceMapping(self.mapping , boundary.mapping)
            a,b                = self.logical_domain, boundary.logical_domain
            l_name             = logical_name or '{l}|{r}'.format(l=a.domain.name, r=b.domain.name)
            int_logical_domain = Interface(l_name, a, b, orientation=orientation)
        else:
            int_map            = None
            int_logical_domain = None

        name = name or '{l}|{r}'.format(l=self.domain.name, r=boundary.domain.name)
        interface = Interface(name, self, boundary,
                              mapping=int_map,
                              logical_domain=int_logical_domain,
                              orientation=orientation)
        return interface

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}_{}'.format(sstr(self.domain),sstr(self.name))

    def __add__(self, other):
        if isinstance(other, ComplementBoundary):
            raise TypeError('> Cannot add complement of boundary')

        return Union(self, other)

    def todict(self):
        name = self.domain.logical_domain.name if self.logical_domain else self.domain.name
        mapping = self.domain.mapping.name if self.domain.mapping else 'None'
        d = {'axis'  : str(self.axis),
                'ext'  : str(self.ext),
                'name' : str(self.name),
                'patch': str(name),
                'mapping':str(mapping)}
        return d
#==============================================================================
class CornerBoundary(BasicDomain):
    """
    Represents a vertex as the intersection of a patch's boundary faces.

    """
    def __new__(cls, *boundaries):
        assert all(isinstance(i, Boundary) for i in boundaries)
        assert all(i.domain==boundaries[0].domain for i in boundaries)

        boundaries = sorted(boundaries, key=lambda x:x.axis)
        obj = Basic.__new__(cls, *boundaries)
        obj._domain = boundaries[0].domain
        return obj

    @property
    def boundaries(self):
        return self._args

    @property
    def domain(self):
        return self._domain

    @property
    def coordinates(self):
        coords = [None]*self.domain.dim
        for b in self.boundaries:
            coords[b.axis] = (b.ext + 1)//2
        return tuple(coords)

    @property
    def logical_domain(self):
        boundaries = [a.logical_domain for a in self.boundaries]
        if boundaries[0]:
            return CornerBoundary(*boundaries)
        else:
            return None

    def _sympystr(self, printer):
        sstr = printer.doprint
        boundaries = ', '.join(sstr(b) for b in self.boundaries)
        return 'CornerBoundary({})'.format(boundaries)

#==============================================================================
class CornerInterface(BasicDomain):
    """
    Represents a vertex shared by multiple patches.

    """
    def __new__(cls, *corners):
        assert all(isinstance(i, CornerBoundary) for i in corners)
        corners = sorted(corners, key=lambda x:x.domain.name)
        return Basic.__new__(cls, *corners)

    @property
    def corners(self):
        return self._args

    @property
    def logical_domain(self):
        corners = [a.logical_domain for a in self.corners]
        if corners[0]:
            return CornerInterface(*corners)
        else:
            return None

    def __len__(self):
        return len(self.corners)

    def _sympystr(self, printer):
        sstr = printer.doprint
        corners = ', '.join(sstr(b) for b in self.corners)
        return 'CornerInterface({})'.format(corners)


# Descriptive multipatch aliases.  The original class names remain canonical
# SymPy and serialization names for backwards compatibility.
PatchVertex = CornerBoundary
SharedVertex = CornerInterface


#==============================================================================
def _interface_axis_map(dim, minus_axis, plus_axis, orientation):
    """Decode a validated compact orientation into a signed axis map."""
    dim = int(dim)
    minus_tangents = tuple(axis for axis in range(dim) if axis != minus_axis)
    plus_tangents = tuple(axis for axis in range(dim) if axis != plus_axis)

    if dim == 1:
        return ()
    if dim == 2:
        return ((minus_tangents[0], plus_tangents[0], orientation),)

    flag, sign1, sign2 = orientation
    plus_positions = (0, 1) if flag == 1 else (1, 0)
    return tuple(
        (minus_axis, plus_tangents[plus_position], direction)
        for minus_axis, plus_position, direction
        in zip(minus_tangents, plus_positions, (sign1, sign2))
    )


#==============================================================================
class Interface(BasicDomain):
    """
    Represents an interface between two subdomains through two boundaries.

    Parameters
    ----------
    name : str
        Name of the interface.

    bnd_minus : Boundary
        Boundary on the "minus" side of the interface.

    bnd_plus : Boundary
        Boundary on the "plus" side of the interface.

    orientation : None, int, or tuple[int, int, int]
        Compact orientation from the minus face to the plus face. Let ``M``
        and ``P`` be the ordered tangential axes obtained by removing the
        normal axis of the minus and plus face, respectively, from
        ``range(dim)``.

        In 1D the orientation is ``None``. In 2D it is one sign, and maps
        ``M[0]`` to ``P[0]`` with that direction. In 3D it is
        ``(flag, sign1, sign2)``. If ``flag`` is ``+1``, ``M[0]`` and ``M[1]``
        map to ``P[0]`` and ``P[1]``. If ``flag`` is ``-1``, they map to
        ``P[1]`` and ``P[0]``. ``sign1`` and ``sign2`` give the respective
        directions. Every flag and sign is either ``-1`` or ``+1``.

    mapping : Mapping, optional
        Mapping from the logical domain to the physical domain, if available.

    logical_domain : BasicDomain, optional
        Logical domain associated with the interface, if available. It should
        be consistent with the mapping if provided.

    Notes
    -----
    ``minus`` and ``plus`` are ordered labels: they define the direction in
    which ``orientation`` is read. They do not describe geometric position or
    mesh size. The two faces may have different normal axes; those axes are
    part of ``bnd_minus`` and ``bnd_plus`` and are deliberately not repeated
    in the orientation.

    This compact orientation convention follows the multipatch convention
    used by GeoPDEs; see its `multipatch geometry specification
    <https://github.com/rafavzqz/geopdes/blob/master/geopdes/doc/geo_specs_mp_v21.txt#L193-L237>`_
    and T. Dokken, E. Quak, V. Skytt, *Requirements from Isogeometric
    Analysis for Changes in Product Design Ontologies* (2010).

    :attr:`axis_map` exposes the decoded convention as
    ``(minus_axis, plus_axis, direction)`` triples.

    Examples
    --------
    A cube interface may exchange its two tangential axes and reverse one of
    them, even when the faces have different normal axes:

    >>> from sympde.topology import Cube, Interface
    >>> A = Cube('A')
    >>> B = Cube('B')
    >>> interface = Interface(
    ...     'A|B',
    ...     A.get_boundary(axis=0, ext=+1),
    ...     B.get_boundary(axis=1, ext=-1),
    ...     orientation=(-1, +1, -1),
    ... )
    >>> interface.orientation
    (-1, 1, -1)
    >>> interface.axis_map
    ((1, 2, 1), (2, 0, -1))

    """
    def __new__(cls, name, bnd_minus, bnd_plus, orientation, *, mapping=None,
                logical_domain=None):

        if not isinstance(name     , str     ): raise TypeError(name)
        if not isinstance(bnd_minus, Boundary): raise TypeError(bnd_minus)
        if not isinstance(bnd_plus , Boundary): raise TypeError(bnd_plus)

        # Check that the dimensions of the two boundaries are the same
        if bnd_minus.dim != bnd_plus.dim:
            raise ValueError(f'Dimension mismatch between boundaries: {bnd_minus.dim} != {bnd_plus.dim}')
        else:
            # Number of logical dimensions of the interface
            ldim = bnd_minus.dim

        # Mapping and logical domain must be provided together or not at all
        assert mapping is None and logical_domain is None or\
               mapping is not None and logical_domain is not None

        # If provided, check that mapping is consistent with the boundaries
        if mapping is not None:
            from sympde.topology.mapping import Mapping
            if not isinstance(mapping, Mapping):
                raise TypeError(f'mapping must be of type Mapping, got {type(mapping)} instead')
            if mapping.ldim != ldim:
                raise ValueError(f'mapping should have logical dimension = {ldim}, got {mapping.ldim} instead')

        # If provided, check that logical domain is consistent with the boundaries
        if logical_domain is not None:
            if not isinstance(logical_domain, BasicDomain):
                raise TypeError(f'logical_domain must be of type BasicDomain, got {type(logical_domain)} instead')
            if logical_domain.dim != ldim:
                raise ValueError(f'Logical domain should have dimension = {ldim}, got {logical_domain.dim} instead')

        # Validate the public input and convert it to the canonical SymPy
        # representation stored in ``args``.  ``axis_map`` can consequently
        # decode ``self.orientation`` without validating it a second time.
        if ldim == 1:
            if orientation is not None:
                raise ValueError('a 1D interface orientation must be None')
            orientation_arg = None
        elif ldim == 2:
            try:
                orientation = _as_integer(
                    orientation, 'a 2D interface orientation')
            except TypeError as error:
                raise TypeError(
                    'a 2D interface orientation must be an integer') \
                    from error
            if orientation not in (-1, 1):
                raise ValueError('a 2D interface orientation must be +1 or -1')
            orientation_arg = Integer(orientation)
        elif ldim == 3:
            if (not isinstance(orientation, (tuple, list, Tuple)) or
                    len(orientation) != 3):
                raise TypeError(
                    'a 3D interface orientation must be '
                    '(flag, sign1, sign2)')
            try:
                orientation = tuple(
                    _as_integer(value, 'a 3D interface orientation value')
                    for value in orientation)
            except TypeError as error:
                message = 'a 3D interface orientation must contain integers'
                raise TypeError(message) from error
            if any(value not in (-1, 1) for value in orientation):
                raise ValueError(
                    'each 3D interface orientation value must be +1 or -1')
            orientation_arg = Tuple(*orientation)
        else:
            raise ValueError(f'unsupported interface dimension: {ldim}')

        obj = Basic.__new__(cls, name, bnd_minus, bnd_plus, orientation_arg)
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
    def minus_side(self):
        """Descriptive alias of :attr:`minus`."""
        return self.minus

    @property
    def plus(self):
        return self.args[2]

    @property
    def plus_side(self):
        """Descriptive alias of :attr:`plus`."""
        return self.plus

    @property
    def orientation(self):
        """Compact dimension-specific orientation of this interface."""
        if self.dim == 1:
            return None
        if self.dim == 2:
            return int(self.args[3])
        return tuple(int(value) for value in self.args[3])

    @property
    def axis_map(self):
        """Tangential-axis correspondence from the minus to the plus side."""
        return _interface_axis_map(
            self.dim, self.minus.axis, self.plus.axis, self.orientation)

    @property
    def patch_dim(self):
        """Logical dimension of each adjacent patch."""
        return self.dim

    @property
    def intrinsic_dim(self):
        """Intrinsic dimension of the interface."""
        return self.dim - 1

    @property
    def mapping(self):
        return self._mapping

    @property
    def logical_domain(self):
        return self._logical_domain

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))


#==============================================================================
class Connectivity(abc.Mapping):
    _patches  = []

    def __init__(self, data=None):
        if data is None:
            data = {}
        else:
            assert isinstance(data, dict)
            for k, v in data.items():
                assert isinstance(k, str)
                assert isinstance(v, Interface)
        self._data = data

    @property
    def patches(self):
        return self._patches

    @property
    def interfaces(self):
        ls = []
        data = dict(sorted(self._data.items()))
        for _,v in data.items():
            ls.append(v)
        return Union(*ls)

    def todict(self):
        # ... create the connectivity
        connectivity = {}
        data = dict(sorted(self._data.items()))
        for name, v in data.items():
            orientation = v.orientation
            if isinstance(orientation, tuple):
                orientation = list(orientation)
            connectivity[name] = [
                v.minus.todict(), v.plus.todict(), orientation]
        connectivity = dict(sorted(connectivity.items()))
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
