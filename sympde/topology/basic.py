# coding: utf-8



from collections import abc
from math import prod
from numbers import Integral

from sympy.core import Basic, Symbol, Integer
from sympy.core.containers import Tuple

from sympde.core.basic import (
    new_basic,
    restore_container_argument,
    restore_int_argument,
    restore_optional_argument,
    restore_string_argument,
)


def _as_integer(value, name):
    """Return an integer input without silently truncating other scalars."""
    if isinstance(value, bool) or not isinstance(value, (Integral, Integer)):
        raise TypeError(f'{name} must be an integer')
    return int(value)

#==============================================================================
class BasicDomain(Basic):
    _ldim        = None
    _name        = None
    _coordinates = None

    @property
    def name(self):
        return self._name

    @property
    def ldim(self):
        """Logical dimension of the domain."""
        return self._ldim

    @property
    def pdim(self):
        """Physical dimension of the domain: the physical dimension of its
        mapping, or its logical dimension if it has no mapping."""
        mapping = getattr(self, 'mapping', None)
        return self.ldim if mapping is None else int(mapping.pdim)

    @property
    def coordinates(self):
        """The logical coordinates of a domain without a mapping, and the
        physical coordinates of a mapped domain."""
        if self._coordinates is None:
            if self.mapping is None:
                xyz = ['x1', 'x2', 'x3'][:self.ldim]
            else:
                xyz = ['x', 'y', 'z'][:self.pdim]

            xyz = [Symbol(i, real=True) for i in xyz]
            self._coordinates = xyz

        if len(self._coordinates) == 1:
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
    def __new__(cls, name, ldim=None, dtype=None, mapping=None,
                logical_domain=None, target=None):
        name = restore_string_argument(name)
        ldim = restore_int_argument(ldim)
        dtype = restore_container_argument(dtype)
        mapping = restore_optional_argument(mapping)
        logical_domain = restore_optional_argument(logical_domain)
        target = restore_optional_argument(target)

        if not isinstance(name, str):
            target = name
            name   = name.name

        if not( target is None ):
            ldim = target.ldim

        assert mapping is None and logical_domain is None or \
        mapping is not None and logical_domain  is not None

        obj = new_basic(cls, name, ldim, dtype, mapping,
                         logical_domain, target)

        obj._ldim           = ldim
        obj._target         = target
        obj._dtype          = dtype
        obj._mapping        = mapping
        obj._logical_domain = logical_domain

        return obj

    @property
    def name(self):
        return self.args[0].name

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
        if len({a.ldim for a in args}) > 1:
            dims = ', '.join(str(a.ldim) for a in args)
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
    def ldim(self):
        return self.args[0].ldim

    @property
    def pdim(self):
        return self.args[0].pdim

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
        obj._ldim = sum(i.ldim for i in args)
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

    _ldim = 1

    def __new__(cls, name=None, coordinate=None, bounds=None):
        name = restore_string_argument(name)
        coordinate = restore_optional_argument(coordinate)
        bounds = restore_container_argument(bounds)

        if name is None:
            name = 'Interval'

        if bounds is None:
            bounds = (0, 1)

        obj = new_basic(cls, name, coordinate, bounds)
        if coordinate:
            obj._coordinates = [coordinate]

        obj._bounds = bounds

        return obj

    @property
    def name(self):
        return self.args[0].name

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

        name = restore_string_argument(name)
        axis = restore_int_argument(axis)
        ext = restore_int_argument(ext)
        mapping = restore_optional_argument(mapping)
        logical_domain = restore_optional_argument(logical_domain)

        if axis is not None:
            axis = _as_integer(axis, 'boundary axis')
            if not 0 <= axis < domain.ldim:
                raise ValueError(
                    f'boundary axis must be between 0 and {domain.ldim - 1}')

        if ext is not None:
            ext = _as_integer(ext, 'boundary extremity')
            if ext not in (-1, 1):
                raise ValueError('boundary extremity must be either -1 or 1')

        obj                 = new_basic(cls, name, domain, axis, ext,
                                         mapping, logical_domain)
        obj._mapping        = mapping
        obj._logical_domain = logical_domain

        return obj

    @property
    def name(self):
        return self.args[0].name

    @property
    def domain(self):
        """The patch (interior domain) of this side."""
        return self.args[1]

    @property
    def axis(self):
        return restore_int_argument(self.args[2])

    @property
    def ext(self):
        return restore_int_argument(self.args[3])

    @property
    def intrinsic_dim(self):
        """Intrinsic dimension of this codimension-one side."""
        return self.ldim - 1

    @property
    def mapping(self):
        return self._mapping

    @property
    def logical_domain(self):
        return self._logical_domain

    @property
    def ldim(self):
        """Logical dimension of the patch owning this side.

        The dimension of the side itself is :attr:`intrinsic_dim`.
        """
        return self.domain.ldim

    @property
    def pdim(self):
        """Physical dimension of the patch owning this side."""
        return self.domain.pdim

    @property
    def adjacent_boundaries(self):
        boundaries = [a for a in self.domain.boundary if a.axis !=self.axis]
        return Union(*boundaries)

    def rotate(self, *directions):
        assert len(directions) == self.ldim-1

        if self.ldim == 2:
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
class InterfaceTransform:
    """
    Full logical coordinate transform from the minus to the plus side of an
    interface.

    :attr:`axis_map` only relates the tangential axes of the two faces. This
    mixin completes it with the normal axes, so that every logical axis of the
    minus patch is mapped to an axis of the plus patch with a direction. It
    is used by :class:`Interface` and by
    :class:`~sympde.topology.incidence.IndexedInterface`, which must provide
    ``minus`` and ``plus`` with ``axis`` and ``ext`` attributes, and
    ``axis_map``.

    Going out of the minus patch through its face means going into the plus
    patch through its face, hence the normal axes have the direction
    ``-minus.ext * plus.ext``.
    """

    @property
    def normal_direction(self):
        """Direction of the plus normal axis along the minus normal axis,
        ``-minus.ext * plus.ext``."""
        return -int(self.minus.ext) * int(self.plus.ext)

    @property
    def signed_permutation(self):
        """``(plus_axis, direction)`` for every axis of the minus patch,
        normal axis included, indexed by the minus axis."""
        ldim = len(self.axis_map) + 1
        permutation = [None] * ldim
        permutation[self.minus.axis] = (int(self.plus.axis), self.normal_direction)
        for minus_axis, plus_axis, direction in self.axis_map:
            permutation[minus_axis] = (int(plus_axis), int(direction))
        return tuple(permutation)

    @property
    def inverse_signed_permutation(self):
        """``(minus_axis, direction)`` for every axis of the plus patch,
        normal axis included, indexed by the plus axis."""
        inverse = [None] * len(self.signed_permutation)
        for minus_axis, (plus_axis, direction) in enumerate(self.signed_permutation):
            inverse[plus_axis] = (minus_axis, direction)
        return tuple(inverse)

    @property
    def permutation_matrix(self):
        """Signed permutation matrix ``P`` as a tuple of rows, such that a
        direction ``d`` of the minus patch is ``P @ d`` in the plus patch."""
        permutation = self.signed_permutation
        ldim = len(permutation)
        rows = [[0] * ldim for _ in range(ldim)]
        for minus_axis, (plus_axis, direction) in enumerate(permutation):
            rows[plus_axis][minus_axis] = direction
        return tuple(tuple(row) for row in rows)

    @property
    def handedness(self):
        """Determinant of :attr:`permutation_matrix`: +1 if the transform
        preserves the orientation of the logical axes, -1 otherwise."""
        permutation = self.signed_permutation
        sign = 1
        for _, direction in permutation:
            sign *= direction
        # Parity of the permutation, from the lengths of its cycles
        visited = [False] * len(permutation)
        for start in range(len(permutation)):
            length = 0
            axis = start
            while not visited[axis]:
                visited[axis] = True
                axis = permutation[axis][0]
                length += 1
            if length and length % 2 == 0:
                sign = -sign
        return sign

    def map_point(self, eta, minus_bounds=None, plus_bounds=None, *, inverse=False):
        """
        Map logical coordinates of the minus patch to the plus patch.

        The map is affine: tangential coordinates are rescaled from the minus
        to the plus interval, reversed if their direction is -1, and the
        minus face is mapped to the plus face. Points inside the minus patch
        are mapped outside the plus patch, across its face.

        Parameters
        ----------
        eta : sequence
            Logical coordinates of the minus patch, one per axis. They may be
            numbers, NumPy arrays or SymPy expressions.

        minus_bounds, plus_bounds : sequence of (lower, upper), optional
            Bounds of the logical coordinates of each patch, one pair per
            axis. Default: ``(0, 1)`` for every axis.

        inverse : bool, default=False
            If true, map coordinates of the plus patch to the minus patch.

        Returns
        -------
        tuple
            Logical coordinates in the other patch, one per axis.
        """
        if inverse:
            source, target = self.plus, self.minus
            source_bounds, target_bounds = plus_bounds, minus_bounds
            permutation = self.inverse_signed_permutation
        else:
            source, target = self.minus, self.plus
            source_bounds, target_bounds = minus_bounds, plus_bounds
            permutation = self.signed_permutation

        ldim = len(permutation)
        if len(eta) != ldim:
            raise ValueError(f'expected {ldim} coordinates, got {len(eta)}')
        source_bounds = ((0, 1),) * ldim if source_bounds is None else tuple(source_bounds)
        target_bounds = ((0, 1),) * ldim if target_bounds is None else tuple(target_bounds)

        image = [None] * ldim
        for source_axis, (target_axis, direction) in enumerate(permutation):
            source_lower, source_upper = source_bounds[source_axis]
            target_lower, target_upper = target_bounds[target_axis]
            # A point of the source with coordinate source_anchor along this
            # axis is mapped to target_anchor along target_axis
            if source_axis == source.axis:
                source_anchor = source_lower if source.ext == -1 else source_upper
                target_anchor = target_lower if target.ext == -1 else target_upper
            else:
                source_anchor = source_lower if direction == 1 else source_upper
                target_anchor = target_lower
            image[target_axis] = target_anchor + direction * (
                (eta[source_axis] - source_anchor) * (target_upper - target_lower)
                / (source_upper - source_lower))
        return tuple(image)

    def map_component(self, k, kind, minus_bounds=None, plus_bounds=None, *, inverse=False):
        """
        Relate a logical component of a vector field on the minus side to
        the plus side.

        The logical components are those of the pull-back of a physical
        field, as defined by :class:`~sympde.topology.mapping.PullBack`: the
        covariant components for H(curl), and the contravariant Piola
        components, with the signed Jacobian determinant, for H(div). Only
        the components that a conforming field keeps continuous across the
        interface are related: the tangential ones for H(curl), the normal
        one for H(div). At a point of the minus face and its image by
        :meth:`map_point`, with ``(j, factor)`` the returned value,
        ``u_minus[k] = factor * u_plus[j]``.

        Parameters
        ----------
        k : int
            Component of the minus side, or of the plus side if ``inverse``
            is true.

        kind : str
            ``'hcurl'`` or ``'hdiv'``.

        minus_bounds, plus_bounds : sequence of (lower, upper), optional
            Bounds of the logical coordinates of each patch, one pair per
            axis. Default: ``(0, 1)`` for every axis.

        inverse : bool, default=False
            If true, relate component ``k`` of the plus side to the minus
            side: ``u_plus[k] = factor * u_minus[j]``.

        Returns
        -------
        tuple[int, float]
            The component ``j`` of the other side and the factor.

        Raises
        ------
        TypeError
            If ``k`` is not an integer.

        ValueError
            If ``k`` is not a component, if ``kind`` is unknown or ``'l2'``,
            or if component ``k`` is not continuous for this kind: the
            normal component for H(curl), a tangential component for H(div).
            These components, and L2 densities, have no relation across the
            interface, unless the geometry is smooth across it.
        """
        if inverse:
            source, permutation = self.plus, self.inverse_signed_permutation
            source_bounds, target_bounds = plus_bounds, minus_bounds
        else:
            source, permutation = self.minus, self.signed_permutation
            source_bounds, target_bounds = minus_bounds, plus_bounds

        ldim = len(permutation)
        k = _as_integer(k, 'the component k')
        if not 0 <= k < ldim:
            raise ValueError(f'the component k must be in range({ldim}), got {k}')

        if kind == 'hcurl':
            if k == source.axis:
                raise ValueError('the normal component of an H(curl) field is '
                                 'not continuous across an interface')
        elif kind == 'hdiv':
            if k != source.axis:
                raise ValueError('a tangential component of an H(div) field is '
                                 'not continuous across an interface')
        elif kind == 'l2':
            raise ValueError('the densities of L2 fields have no relation across an '
                             'interface, unless the geometry is smooth across it')
        else:
            raise ValueError(f"kind must be 'hcurl' or 'hdiv', got {kind!r}")

        source_bounds = ((0, 1),) * ldim if source_bounds is None else tuple(source_bounds)
        target_bounds = ((0, 1),) * ldim if target_bounds is None else tuple(target_bounds)
        source_lengths = [upper - lower for lower, upper in source_bounds]
        target_lengths = [upper - lower for lower, upper in target_bounds]

        # The Jacobian of map_point is A = diag(target_lengths) @ P @ diag(1/source_lengths),
        # whose tangential columns relate the tangential columns of the
        # Jacobians of the two mappings on the interface. Hence the covariant
        # components satisfy u_source = A^T u_target, and the normal
        # contravariant Piola component u_source = det(A) A^{-1} u_target.
        j, direction = permutation[k]
        if kind == 'hcurl':
            factor = direction * target_lengths[j] / source_lengths[k]
        else:
            factor = (self.handedness * direction
                      * prod(target_lengths) / prod(source_lengths)
                      * source_lengths[k] / target_lengths[j])
        return j, factor


#==============================================================================
class Interface(BasicDomain, InterfaceTransform):
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
        ``range(ldim)``.

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
    ``(minus_axis, plus_axis, direction)`` triples. The properties and
    methods inherited from :class:`InterfaceTransform` complete it with the
    normal axes, e.g. :attr:`~InterfaceTransform.signed_permutation`,
    :meth:`~InterfaceTransform.map_point` and
    :meth:`~InterfaceTransform.map_component`.

    Examples
    --------
    A cube interface may exchange its two tangential axes and reverse one of
    them, even when the faces have different normal axes:

    >>> from sympde.api import Cube, Interface
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
    >>> interface.signed_permutation
    ((1, 1), (2, 1), (0, -1))
    >>> interface.map_point((1, 0.25, 0.5))
    (0.5, 0.0, 0.25)
    >>> interface.map_component(0, 'hdiv')
    (1, -1.0)

    """
    def __new__(cls, name, bnd_minus, bnd_plus, orientation, mapping=None,
                logical_domain=None):

        name = restore_string_argument(name)
        orientation = restore_container_argument(orientation)
        mapping = restore_optional_argument(mapping)
        logical_domain = restore_optional_argument(logical_domain)

        if not isinstance(name     , str     ): raise TypeError(name)
        if not isinstance(bnd_minus, Boundary): raise TypeError(bnd_minus)
        if not isinstance(bnd_plus , Boundary): raise TypeError(bnd_plus)

        # Check that the dimensions of the two boundaries are the same
        if bnd_minus.ldim != bnd_plus.ldim:
            raise ValueError(f'Dimension mismatch between boundaries: {bnd_minus.ldim} != {bnd_plus.ldim}')
        else:
            # Number of logical dimensions of the interface
            ldim = bnd_minus.ldim

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
            if logical_domain.ldim != ldim:
                raise ValueError(f'Logical domain should have dimension = {ldim}, got {logical_domain.ldim} instead')

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

        obj = new_basic(cls, name, bnd_minus, bnd_plus, orientation_arg,
                         mapping, logical_domain)
        obj._mapping        = mapping
        obj._logical_domain = logical_domain
        return obj

    @property
    def ldim(self):
        """Logical dimension of the patches joined by this interface.

        The dimension of the interface itself is :attr:`intrinsic_dim`.
        """
        return self.minus.ldim

    @property
    def pdim(self):
        """Physical dimension of the patches joined by this interface."""
        return self.minus.pdim

    @property
    def name(self):
        return self.args[0].name

    @property
    def minus(self):
        return self.args[1]

    @property
    def plus(self):
        return self.args[2]

    @property
    def orientation(self):
        """Compact dimension-specific orientation of this interface."""
        if self.ldim == 1:
            return None
        if self.ldim == 2:
            return int(self.args[3])
        return restore_container_argument(self.args[3])

    @property
    def axis_map(self):
        """Tangential-axis correspondence from the minus to the plus side."""
        return _interface_axis_map(
            self.ldim, self.minus.axis, self.plus.axis, self.orientation)

    @property
    def intrinsic_dim(self):
        """Intrinsic dimension of the interface."""
        return self.ldim - 1

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
    """Read-only mapping from interface names to interfaces."""

    def __init__(self, data=None):
        if data is None:
            data = {}
        else:
            assert isinstance(data, dict)
            for k, v in data.items():
                assert isinstance(k, str)
                assert isinstance(v, Interface)
        self._data = dict(data)

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

        return connectivity

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
