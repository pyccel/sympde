#---------------------------------------------------------------------------#
# This file is part of SymPDE.                                              #
#---------------------------------------------------------------------------#
"""Index-based views of the topology of a multipatch domain.

The objects of this module are plain Python data, meant for discretization
codes such as PSYDAC: patches are given by their position in
``domain.patches``, and the vertices and edges of a patch by their
position in the parametric box of the patch, with 0 for the lower and 1 for
the upper end of an axis.
"""

from dataclasses import dataclass
from itertools import combinations, product
from types import MappingProxyType

__all__ = (
    'IndexedInterface',
    'PatchSide',
    'TopologyEdge',
    'TopologyEdgeIncidence',
    'TopologyVertex',
    'TopologyVertexIncidence',
    'collect_edges',
    'collect_vertices',
    'index_exterior_sides',
    'index_interfaces',
)

#==============================================================================
# Interfaces and exterior sides
#==============================================================================
@dataclass(frozen=True)
class PatchSide:
    """A face of the patch `patch_index`: a side of an interface, or an
    exterior side of the domain."""

    patch_index: int  # position of the patch in domain.patches
    axis: int         # logical axis normal to the face
    ext: int          # -1 for the lower face, +1 for the upper face


@dataclass(frozen=True)
class IndexedInterface:
    """An interface of a multipatch domain, with patch indices.

    `orientation` is the compact orientation of the interface (None in 1D,
    +1 or -1 in 2D, (flag, sign1, sign2) in 3D), and `axis_map` its decoded
    form: one entry (minus_axis, plus_axis, direction) per tangential axis of
    the minus side.
    """

    name: str
    minus: PatchSide
    plus: PatchSide
    orientation: object
    axis_map: tuple


def _patch_indices(domain):
    return {patch: index for index, patch in enumerate(domain.patches)}


def _side(boundary, patch_indices):
    if boundary.axis is None or boundary.ext is None:
        raise ValueError(
            f'the side {boundary} has no axis or extremity')
    return PatchSide(
        patch_indices[boundary.domain], int(boundary.axis), int(boundary.ext))


def index_interfaces(domain):
    """Return the interfaces of a domain with patch indices, keyed by name.

    This is the value of :attr:`Domain.indexed_interfaces`.

    Parameters
    ----------
    domain : sympde.topology.domain.Domain
        A logical or mapped domain.

    Returns
    -------
    Mapping[str, IndexedInterface]
        Read-only mapping, ordered by interface name like
        ``domain.interface_map``.
    """
    patch_indices = _patch_indices(domain)
    return MappingProxyType({
        str(name): IndexedInterface(
            name=str(name),
            minus=_side(interface.minus, patch_indices),
            plus=_side(interface.plus, patch_indices),
            orientation=interface.orientation,
            axis_map=interface.axis_map,
        )
        for name, interface in domain.interface_map.items()
    })


def index_exterior_sides(domain):
    """Return the exterior sides of a domain with patch indices.

    This is the value of :attr:`Domain.indexed_exterior_sides`.

    Parameters
    ----------
    domain : sympde.topology.domain.Domain
        A logical or mapped domain.

    Returns
    -------
    tuple[PatchSide, ...]
        The faces of the patches that are not part of an interface, in the
        order of ``domain.exterior_sides``.
    """
    patch_indices = _patch_indices(domain)
    return tuple(_side(side, patch_indices) for side in domain.exterior_sides)

#==============================================================================
# Vertices and edges
#==============================================================================
@dataclass(frozen=True)
class TopologyVertexIncidence:
    """One patch-local representative of a topological vertex."""

    patch_index: int
    patch_name: str
    logical_corner: tuple  # 0 or 1 for every axis


@dataclass(frozen=True)
class TopologyVertex:
    """A topological vertex and all patch-local corners identified with it."""

    index: int
    incidences: tuple
    is_boundary: bool

    @property
    def patch_indices(self):
        """Indices of incident patches, without duplicates."""
        return tuple(dict.fromkeys(i.patch_index for i in self.incidences))

    @property
    def patches(self):
        """Names of incident patches, without duplicates."""
        return tuple(dict.fromkeys(i.patch_name for i in self.incidences))


@dataclass(frozen=True)
class TopologyEdgeIncidence:
    """One patch-local representative of a topological edge.

    `logical_edge` has None for the axis along the edge, and 0 or 1 for the
    other axes. `direction` is +1 if the increasing coordinate along the edge
    has the same direction as in the first incidence of the edge, and -1
    otherwise.
    """

    patch_index: int
    patch_name: str
    logical_edge: tuple
    direction: int

    @property
    def axis(self):
        """Logical axis along the edge."""
        return self.logical_edge.index(None)


@dataclass(frozen=True)
class TopologyEdge:
    """A topological edge and all patch-local edges identified with it."""

    index: int
    incidences: tuple
    is_boundary: bool

    @property
    def patch_indices(self):
        """Indices of incident patches, without duplicates."""
        return tuple(dict.fromkeys(i.patch_index for i in self.incidences))

    @property
    def patches(self):
        """Names of incident patches, without duplicates."""
        return tuple(dict.fromkeys(i.patch_name for i in self.incidences))


def _patch_name(patch):
    logical_patch = patch.logical_domain
    return str(logical_patch.name if logical_patch is not None else patch.name)


def _box_corners_and_edges(ndim, dim):
    """All corners (`dim` 0) or edges (`dim` 1) of an `ndim`-dimensional box.

    Each one is given by its position: a tuple with None for the axis along
    an edge, and 0 (lower end) or 1 (upper end) for the other axes.
    """
    for free_axes in combinations(range(ndim), dim):
        fixed_axes = [axis for axis in range(ndim) if axis not in free_axes]
        for values in product((0, 1), repeat=len(fixed_axes)):
            position = [None] * ndim
            for axis, value in zip(fixed_axes, values):
                position[axis] = value
            yield tuple(position)


def _map_to_plus_side(position, interface):
    """Map a corner or edge of the minus face of an interface to the plus face.

    Return None if `position` is not on the minus face, and otherwise its
    position in the plus patch and the sign relating the directions of the
    two edges (+1 for corners).
    """
    minus = interface.minus
    plus = interface.plus
    if position[minus.axis] != (0 if minus.ext == -1 else 1):
        return None
    image = [None] * len(position)
    image[plus.axis] = 0 if plus.ext == -1 else 1
    sign = 1
    for minus_axis, plus_axis, direction in interface.axis_map:
        value = position[minus_axis]
        if value is None:
            sign *= direction
        else:
            image[plus_axis] = value if direction == 1 else 1 - value
    return tuple(image), sign


def _sort_key(key):
    patch_index, position = key
    return patch_index, tuple(
        -1 if value is None else value for value in position)


def _identify_through_interfaces(domain, dim):
    """Group the corners (`dim` 0) or edges (`dim` 1) of all patches that the
    interfaces identify, directly or through other patches.

    Return a list of (keys, signs, is_boundary) for every group, in
    deterministic order. `keys` are the (patch_index, position) of the
    members of the group, sorted, and `signs` their directions relative to
    the first one.
    """
    ndim = int(domain.ldim)
    patch_indices = _patch_indices(domain)
    interfaces = tuple(domain.interface_map.values())
    for interface in interfaces:
        for side in (interface.minus, interface.plus):
            if side.axis is None or side.ext is None:
                raise ValueError(
                    f'the interface side {side} has no axis or extremity, so '
                    'its vertices and edges are not defined')

    keys = [
        (patch_index, position)
        for patch_index in range(len(domain.patches))
        for position in _box_corners_and_edges(ndim, dim)
    ]
    # Union-find in which every key also stores the sign of its direction
    # relative to its parent
    parent = {key: (key, 1) for key in keys}

    def find(key):
        root, sign = parent[key]
        if root != key:
            root, root_sign = find(root)
            sign *= root_sign
            parent[key] = (root, sign)
        return root, sign

    def union(first, second, sign):
        # direction(second) = sign * direction(first)
        first_root, first_sign = find(first)
        second_root, second_sign = find(second)
        relative_sign = first_sign * sign * second_sign
        if first_root == second_root:
            if relative_sign != 1:
                raise ValueError(
                    'the interfaces identify an edge with itself in the '
                    'opposite direction')
            return
        if _sort_key(second_root) < _sort_key(first_root):
            first_root, second_root = second_root, first_root
        parent[second_root] = (first_root, relative_sign)

    glued_faces = set()
    for interface in interfaces:
        minus_index = patch_indices[interface.minus.domain]
        plus_index = patch_indices[interface.plus.domain]
        glued_faces.add(
            (minus_index, int(interface.minus.axis), int(interface.minus.ext)))
        glued_faces.add(
            (plus_index, int(interface.plus.axis), int(interface.plus.ext)))
        for position in _box_corners_and_edges(ndim, dim):
            image = _map_to_plus_side(position, interface)
            if image is not None:
                union((minus_index, position), (plus_index, image[0]),
                      image[1])

    groups = {}
    for key in keys:
        groups.setdefault(find(key)[0], []).append(key)

    def faces(key):
        patch_index, position = key
        return ((patch_index, axis, -1 if value == 0 else 1)
                for axis, value in enumerate(position) if value is not None)

    def group_key(group):
        return _sort_key(min(group, key=_sort_key))

    result = []
    for group in sorted(groups.values(), key=group_key):
        group = sorted(group, key=_sort_key)
        signs = tuple(find(key)[1] * find(group[0])[1] for key in group)
        is_boundary = any(
            face not in glued_faces for key in group for face in faces(key))
        result.append((tuple(group), signs, is_boundary))
    return result


def collect_vertices(domain):
    """Return all vertices of a domain and their patch incidences.

    This is the value of :attr:`Domain.vertices`. Each returned vertex records
    all of its patch-local representatives, including unshared patch corners,
    and whether it lies on the exterior boundary.

    Parameters
    ----------
    domain : sympde.topology.domain.Domain
        A logical or mapped domain, whose interface sides have an axis and an
        extremity.

    Returns
    -------
    tuple[TopologyVertex, ...]
        Topological vertices in deterministic patch/corner order.
    """
    names = tuple(_patch_name(patch) for patch in domain.patches)
    return tuple(
        TopologyVertex(
            index=index,
            incidences=tuple(
                TopologyVertexIncidence(
                    patch_index, names[patch_index], position)
                for patch_index, position in keys),
            is_boundary=is_boundary,
        )
        for index, (keys, _, is_boundary)
        in enumerate(_identify_through_interfaces(domain, 0)))


def collect_edges(domain):
    """Return all edges of a 3D domain and their patch incidences.

    This is the value of :attr:`Domain.edges`. An edge where several patches
    meet has one incidence per patch edge identified with it, through any
    chain of interfaces.

    Parameters
    ----------
    domain : sympde.topology.domain.Domain
        A three-dimensional logical or mapped domain, whose interface sides
        have an axis and an extremity.

    Returns
    -------
    tuple[TopologyEdge, ...]
        Topological edges in deterministic patch/edge order. The incidences
        of an edge record their direction relative to the first one, which
        fixes the sign of the tangential components shared along the edge.

    Raises
    ------
    ValueError
        If the domain is not three-dimensional, or if the interfaces identify
        an edge with itself in the opposite direction.
    """
    if int(domain.ldim) != 3:
        raise ValueError(
            'edges are only defined for 3D domains; in 2D, the edges of the '
            'patches are their sides, given by domain.indexed_interfaces and '
            'domain.indexed_exterior_sides')
    names = tuple(_patch_name(patch) for patch in domain.patches)
    return tuple(
        TopologyEdge(
            index=index,
            incidences=tuple(
                TopologyEdgeIncidence(
                    patch_index, names[patch_index], position, sign)
                for (patch_index, position), sign in zip(keys, signs)),
            is_boundary=is_boundary,
        )
        for index, (keys, signs, is_boundary)
        in enumerate(_identify_through_interfaces(domain, 1)))
