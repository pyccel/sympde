#---------------------------------------------------------------------------#
# This file is part of SymPDE.                                              #
#---------------------------------------------------------------------------#
"""Inspection and plotting helpers for symbolic multipatch topology."""

from dataclasses import dataclass
from itertools import product
import numbers
import sys

import numpy as np

__all__ = (
    'TopologyVertexIncidence',
    'TopologyVertex',
    'collect_topology_vertices',
    'format_topology',
    'print_topology',
    'plot_domain',
)


@dataclass(frozen=True)
class TopologyVertexIncidence:
    """One patch-local representative of a topological vertex."""

    patch_index: int
    patch_name: str
    logical_corner: tuple


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


def _patch_name(patch):
    logical_patch = patch.logical_domain
    return str(logical_patch.name if logical_patch is not None else patch.name)


def _interfaces(domain):
    return tuple(domain.interface_map.values())


def _shared_vertices(domain):
    return domain.shared_vertices


def _check_topology_dimension(domain):
    if int(domain.dim) not in (2, 3):
        raise ValueError(
            'topology vertex collection and annotated topology plotting '
            'support only two- and three-dimensional domains')


def collect_topology_vertices(domain):
    """Return all vertices of a 2D or 3D domain and their patch incidences.

    Unlike ``domain.corners``, which stores only equivalence classes created
    by interface identifications, this function also includes unshared patch
    corners. Each returned vertex records all of its patch-local
    representatives and whether it lies on the exterior boundary.

    Parameters
    ----------
    domain : sympde.topology.domain.Domain
        A two- or three-dimensional logical or mapped domain.

    Returns
    -------
    tuple[TopologyVertex, ...]
        Topological vertices in deterministic patch/corner order.
    """
    _check_topology_dimension(domain)

    ndim = int(domain.dim)
    patches = domain.patches
    patch_indices = {patch: index for index, patch in enumerate(patches)}
    patch_names = tuple(_patch_name(patch) for patch in patches)
    keys = [
        (patch_index, tuple(corner))
        for patch_index in range(len(patches))
        for corner in product((0, 1), repeat=ndim)
    ]
    parent = {key: key for key in keys}

    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(first, second):
        first_root = find(first)
        second_root = find(second)
        if first_root != second_root:
            parent[second_root] = first_root

    for shared_vertex in _shared_vertices(domain):
        members = [
            (patch_indices[corner.domain], tuple(map(int, corner.coordinates)))
            for corner in shared_vertex.corners
        ]
        for member in members[1:]:
            union(members[0], member)

    groups = {}
    for key in keys:
        groups.setdefault(find(key), []).append(key)

    glued_faces = set()
    for interface in _interfaces(domain):
        for face in (interface.minus, interface.plus):
            glued_faces.add(
                (patch_indices[face.domain], int(face.axis), int(face.ext)))

    ordered_groups = sorted(
        (sorted(group) for group in groups.values()),
        key=lambda group: group[0],
    )

    vertices = []
    for index, group in enumerate(ordered_groups):
        incidences = tuple(
            TopologyVertexIncidence(
                patch_index=patch_index,
                patch_name=patch_names[patch_index],
                logical_corner=corner,
            )
            for patch_index, corner in group
        )
        local_faces = (
            (
                incidence.patch_index,
                axis,
                -1 if incidence.logical_corner[axis] == 0 else 1,
            )
            for incidence in incidences
            for axis in range(ndim)
        )
        is_boundary = any(face not in glued_faces for face in local_faces)
        vertices.append(TopologyVertex(index, incidences, is_boundary))

    return tuple(vertices)


def format_topology(domain):
    """Return a human-readable summary of patches, interfaces, and vertices."""
    _check_topology_dimension(domain)

    patches = domain.patches
    lines = [
        f'=== {domain.name} ===',
        'patches: ' + ', '.join(_patch_name(patch) for patch in patches),
        'interfaces:',
    ]

    interfaces = _interfaces(domain)
    if not interfaces:
        lines.append('  (none)')
    for index, interface in enumerate(interfaces):
        minus = interface.minus
        plus = interface.plus
        lines.extend((
            (
                f'  I{index} ({interface.name}): '
                f'{_patch_name(minus.domain)}'
                f'(axis={minus.axis}, ext={int(minus.ext):+d}) <-> '
                f'{_patch_name(plus.domain)}'
                f'(axis={plus.axis}, ext={int(plus.ext):+d})'
            ),
            (
                f'      orientation={interface.orientation}, '
                f'axis_map={interface.axis_map}'
            ),
        ))

    lines.append('vertices:')
    for vertex in collect_topology_vertices(domain):
        incidences = ', '.join(
            f'{incidence.patch_name}{incidence.logical_corner}'
            for incidence in vertex.incidences
        )
        location = 'boundary' if vertex.is_boundary else 'interior'
        lines.append(
            f'  V{vertex.index}: {location}; patches={vertex.patches}; '
            f'local corners=[{incidences}]')

    return '\n'.join(lines)


def print_topology(domain, file=None):
    """Print :func:`format_topology` to ``file`` or standard output."""
    print(format_topology(domain), file=sys.stdout if file is None else file)


#==============================================================================
# Plotting style
#==============================================================================
# Okabe-Ito palette, which is safe for colour-blind readers
_PATCH_COLORS = (
    '#E69F00', '#56B4E9', '#009E73', '#F0E442',
    '#0072B2', '#D55E00', '#CC79A7', '#999999',
)
_INTERFACE_COLORMAP = 'tab10'
_PLUS_LINESTYLE = (0, (3, 2))

_STYLE = {
    # Patches
    'patch_alpha'         : 0.3,   # 2D patch fill
    'surface_alpha'       : 0.6,   # 3D patch faces, geometry plot
    'face_alpha'          : 0.08,  # 3D exterior faces, topology plot
    'interface_alpha'     : 0.3,   # 3D minus faces of the interfaces
    'boundary_width'      : 1.8,   # exterior boundary
    'edge_width'          : 0.8,   # other patch edges
    'edge_color'          : '0.45',
    'isoline_width'       : 0.5,
    'isoline_color'       : '0.55',
    'isoline_count'       : 7,     # interior isolines per direction
    # Interface sides
    'side_width'          : 2.2,
    'side_inset'          : 0.06,  # fraction of the patch width, in 2D
    'arrow_size'          : 14,
    'connector_width'     : 1.2,
    # Labels and vertices
    'patch_fontsize'      : 11,
    'label_fontsize'      : 9,
    'vertex_fontsize'     : 7,
    'vertex_size'         : 5,
    'halo_width'          : 2.5,
    # Figures
    'figsize_2d'          : (7.5, 6.5),
    'figsize_3d'          : (8.0, 7.0),
    'view'                : (24, -58),  # elevation and azimuth in 3D
    'min_box_aspect'      : 0.3,   # smallest 3D box side, relative
}


def _halo():
    import matplotlib.patheffects as path_effects
    return [path_effects.withStroke(
        linewidth=_STYLE['halo_width'], foreground='white')]


#==============================================================================
# Sampling of the geometry
#==============================================================================
def _logical_limits(patch):
    logical_patch = patch.logical_domain
    parameter_domain = logical_patch if logical_patch is not None else patch
    return (
        np.asarray(parameter_domain.min_coords, dtype=float),
        np.asarray(parameter_domain.max_coords, dtype=float),
    )


def _map_points(patch, logical_points):
    """Map an array of logical points, of shape (n, dim), to the physical
    space."""
    logical_points = np.asarray(logical_points, dtype=float)
    if patch.mapping is None:
        return logical_points
    mapping = patch.mapping.get_callable_mapping()
    components = mapping(*logical_points.T)
    return np.column_stack([
        np.broadcast_to(np.asarray(c, dtype=float), len(logical_points))
        for c in components])


def _physical_dimension(domain):
    dimensions = {
        int(patch.mapping.pdim) if patch.mapping is not None
        else int(patch.dim)
        for patch in domain.patches
    }
    if len(dimensions) != 1:
        raise ValueError('all plotted patches must have one physical dimension')
    return dimensions.pop()


def _face_points(patch, axis, ext, samples, inset=0.0):
    """Logical points on a face of a patch, of shape (samples, ..., dim).

    The face is moved into the patch by the fraction `inset` of the patch
    width along its normal axis.
    """
    lower, upper = _logical_limits(patch)
    dim = len(lower)
    tangents = [a for a in range(dim) if a != axis]
    grids = np.meshgrid(
        *[np.linspace(lower[a], upper[a], samples) for a in tangents],
        indexing='ij')
    width = upper[axis] - lower[axis]
    position = lower[axis] + inset * width if ext == -1 else \
               upper[axis] - inset * width
    points = np.empty((*grids[0].shape, dim))
    points[..., axis] = position
    for a, grid in zip(tangents, grids):
        points[..., a] = grid
    return points


def _map_face(patch, axis, ext, samples, inset=0.0):
    points = _face_points(patch, axis, ext, samples, inset)
    shape = points.shape
    return _map_points(patch, points.reshape(-1, shape[-1])).reshape(shape)


def _patch_center(patch):
    lower, upper = _logical_limits(patch)
    return _map_points(patch, [(lower + upper) / 2])[0]


def _face_center(face):
    lower, upper = _logical_limits(face.domain)
    point = (lower + upper) / 2
    point[face.axis] = lower[face.axis] if face.ext == -1 else upper[face.axis]
    return _map_points(face.domain, [point])[0]


def _corner_point(patch, logical_corner):
    lower, upper = _logical_limits(patch)
    point = np.where(np.asarray(logical_corner) == 0, lower, upper)
    return _map_points(patch, [point])[0]


def _patch_polygon(patch, samples):
    """Closed boundary of a 2D patch, counter-clockwise in logical space."""
    return np.vstack((
        _map_face(patch, 1, -1, samples),
        _map_face(patch, 0, +1, samples),
        _map_face(patch, 1, +1, samples)[::-1],
        _map_face(patch, 0, -1, samples)[::-1],
    ))


#==============================================================================
# Plot model: what is drawn, independently of Matplotlib
#==============================================================================
@dataclass(frozen=True)
class _PlotInterface:
    """An interface with the points needed to draw it."""

    index: int
    interface: object
    color: tuple
    minus_center: np.ndarray
    plus_center: np.ndarray
    separated: bool   # the two sides are drawn at different places


@dataclass(frozen=True)
class _PlotModel:
    """The sampled geometry and topology of a domain."""

    domain: object
    dim: int
    refinement: int
    patches: tuple
    patch_colors: tuple
    faces: tuple           # per patch, {(axis, ext): sampled face}
    exterior_faces: frozenset  # (patch index, axis, ext)
    interfaces: tuple      # _PlotInterface
    vertices: tuple        # (TopologyVertex, tuple of physical points)
    extent: np.ndarray     # size of the bounding box

    @property
    def scale(self):
        return max(float(np.max(self.extent)), 1.0)


def _patch_neighbours(domain, patches):
    """Pairs of distinct patches that share an interface or a vertex."""
    index = {patch: i for i, patch in enumerate(patches)}
    neighbours = {i: set() for i in range(len(patches))}

    def connect(group):
        for i in group:
            neighbours[i].update(j for j in group if j != i)

    for interface in _interfaces(domain):
        connect({index[interface.minus.domain], index[interface.plus.domain]})
    for vertex in collect_topology_vertices(domain):
        connect(set(vertex.patch_indices))
    return neighbours


def _color_patches(domain, patches):
    """Give every patch its own colour if the palette is large enough, and
    otherwise give neighbouring patches different colours, when possible."""
    if len(patches) <= len(_PATCH_COLORS):
        return _PATCH_COLORS[:len(patches)]
    neighbours = _patch_neighbours(domain, patches)
    colors = {}
    for i in range(len(patches)):
        used = {colors[j] for j in neighbours[i] if j in colors}
        free = [c for c in range(len(_PATCH_COLORS)) if c not in used]
        colors[i] = free[0] if free else i % len(_PATCH_COLORS)
    return tuple(_PATCH_COLORS[colors[i]] for i in range(len(patches)))


def _build_model(domain, refinement, topology):
    import matplotlib.pyplot as plt

    dim = int(domain.dim)
    patches = domain.patches
    faces = tuple(
        {(axis, ext): _map_face(patch, axis, ext, refinement)
         for axis in range(dim) for ext in (-1, 1)}
        for patch in patches)
    points = np.vstack([
        face.reshape(-1, dim) for patch_faces in faces
        for face in patch_faces.values()])
    extent = np.ptp(points, axis=0)

    patch_index = {patch: i for i, patch in enumerate(patches)}
    glued_faces = {
        (patch_index[face.domain], int(face.axis), int(face.ext))
        for interface in _interfaces(domain)
        for face in (interface.minus, interface.plus)}
    exterior_faces = frozenset(
        (i, axis, ext) for i in range(len(patches))
        for axis in range(dim) for ext in (-1, 1)
        if (i, axis, ext) not in glued_faces)

    if not topology:
        return _PlotModel(
            domain, dim, refinement, patches,
            _color_patches(domain, patches),
            faces, exterior_faces, (), (), extent)

    scale = max(float(np.max(extent)), 1.0)
    interface_colors = plt.get_cmap(_INTERFACE_COLORMAP).colors
    interfaces = []
    for index, interface in enumerate(_interfaces(domain)):
        minus_center = _face_center(interface.minus)
        plus_center = _face_center(interface.plus)
        separated = (
            np.linalg.norm(minus_center - plus_center) > 1e-10 * scale)
        interfaces.append(_PlotInterface(
            index, interface,
            interface_colors[index % len(interface_colors)],
            minus_center, plus_center, separated))

    vertices = []
    for vertex in collect_topology_vertices(domain):
        locations = {}
        for incidence in vertex.incidences:
            point = _corner_point(
                patches[incidence.patch_index], incidence.logical_corner)
            locations.setdefault(tuple(np.round(point, decimals=12)), point)
        vertices.append((vertex, tuple(locations.values())))

    return _PlotModel(
        domain, dim, refinement, patches, _color_patches(domain, patches),
        faces, exterior_faces, tuple(interfaces), tuple(vertices), extent)


#==============================================================================
# Drawing helpers shared by 2D and 3D
#==============================================================================
def _new_axes(dim, ax, topology):
    import matplotlib.pyplot as plt

    if ax is None:
        figsize = _STYLE['figsize_2d'] if dim == 2 else _STYLE['figsize_3d']
        figure = plt.figure(
            figsize=figsize if topology else None, layout='constrained')
        ax = figure.add_subplot(projection='3d' if dim == 3 else None)
        return figure, ax

    if hasattr(ax, 'zaxis') != (dim == 3):
        raise ValueError(
            f'a {dim}D Matplotlib axes is required for a {dim}D domain')
    return ax.figure, ax


def _text(ax, point, text, **options):
    """Write a label with a white halo, in 2D or 3D."""
    options.setdefault('ha', 'center')
    options.setdefault('va', 'center')
    ax.text(*point, text, path_effects=_halo(), **options)


def _draw_interface_labels(ax, model):
    for item in model.interfaces:
        style = dict(
            color=item.color, fontsize=_STYLE['label_fontsize'],
            fontweight='bold', zorder=7,
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                      edgecolor=item.color, linewidth=0.6))
        if not item.separated:
            _text(ax, item.minus_center, f'I{item.index}', **style)
            continue
        _text(ax, item.minus_center, f'I{item.index}-', **style)
        _text(ax, item.plus_center, f'I{item.index}+', **style)


def _draw_vertices(ax, model, labels):
    for vertex, points in model.vertices:
        face_color = 'white' if vertex.is_boundary else 'black'
        for point in points:
            ax.plot(*([coordinate] for coordinate in point), marker='o',
                    linestyle='', markersize=_STYLE['vertex_size'],
                    markerfacecolor=face_color, markeredgecolor='black',
                    zorder=8)
            if labels:
                _text(ax, point, f'V{vertex.index}', ha='left', va='bottom',
                      fontsize=_STYLE['vertex_fontsize'], color='0.2',
                      zorder=9)


def _draw_legend(ax, dim):
    from matplotlib.lines import Line2D

    side = 'side' if dim == 2 else 'face arrows'
    handles = [
        Line2D([], [], color='0.3', linewidth=_STYLE['side_width'],
               label=f'minus {side}'),
        Line2D([], [], color='0.3', linewidth=_STYLE['side_width'],
               linestyle=_PLUS_LINESTYLE, label=f'plus {side}'),
    ]
    if dim == 2:
        handles.append(Line2D(
            [], [], color='black', linewidth=_STYLE['boundary_width'],
            label='exterior boundary'))
    handles += [
        Line2D([], [], marker='o', linestyle='', markerfacecolor='white',
               markeredgecolor='black', label='boundary vertex'),
        Line2D([], [], marker='o', linestyle='', markerfacecolor='black',
               markeredgecolor='black', label='interior vertex'),
    ]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.0),
              ncol=3, fontsize=8, frameon=False)


#==============================================================================
# 2D plots
#==============================================================================
def _draw_patches_2d(ax, model, isolines, topology, labels):
    for i, patch in enumerate(model.patches):
        polygon = _patch_polygon(patch, model.refinement)
        ax.fill(*polygon.T, facecolor=model.patch_colors[i],
                alpha=_STYLE['patch_alpha'], linewidth=0, zorder=1)
        if isolines:
            _draw_isolines_2d(ax, patch, model.refinement)
        for (axis, ext), curve in model.faces[i].items():
            exterior = (i, axis, ext) in model.exterior_faces
            ax.plot(*curve.T, zorder=3, **(
                dict(color='black', linewidth=_STYLE['boundary_width'])
                if exterior or not topology else
                dict(color=_STYLE['edge_color'],
                     linewidth=_STYLE['edge_width'])))
        if labels:
            _text(ax, _patch_center(patch), _patch_name(patch),
                  fontsize=_STYLE['patch_fontsize'], fontweight='bold',
                  zorder=6)


def _draw_isolines_2d(ax, patch, samples):
    lower, upper = _logical_limits(patch)
    fractions = np.linspace(0, 1, _STYLE['isoline_count'] + 2)[1:-1]
    for axis in range(2):
        tangent = 1 - axis
        points = np.empty((samples, 2))
        points[:, tangent] = np.linspace(lower[tangent], upper[tangent], samples)
        for fraction in fractions:
            points[:, axis] = lower[axis] + fraction * (upper[axis] - lower[axis])
            ax.plot(*_map_points(patch, points).T,
                    color=_STYLE['isoline_color'],
                    linewidth=_STYLE['isoline_width'], zorder=2)


def _draw_side_2d(ax, face, color, linestyle, samples, direction=1):
    """Draw an interface side inside its patch, with an arrow along its
    tangential coordinate, reversed if `direction` is -1."""
    from matplotlib.patches import FancyArrowPatch

    curve = _map_face(face.domain, int(face.axis), int(face.ext), samples,
                      inset=_STYLE['side_inset'])
    ax.plot(*curve.T, color=color, linewidth=_STYLE['side_width'],
            linestyle=linestyle, solid_capstyle='butt', zorder=4)
    start, end = curve[int(0.45 * samples)], curve[int(0.62 * samples)]
    if direction == -1:
        start, end = end, start
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle='-|>', mutation_scale=_STYLE['arrow_size'],
        color=color, linewidth=0, zorder=5))


def _draw_interfaces_2d(ax, model):
    from matplotlib.patches import FancyArrowPatch

    for item in model.interfaces:
        # The arrow of the plus side shows where the coordinate of the minus
        # side goes, so that it is reversed for a reversing interface
        (_, _, direction), = item.interface.axis_map
        _draw_side_2d(ax, item.interface.minus, item.color, '-',
                      model.refinement)
        _draw_side_2d(ax, item.interface.plus, item.color, _PLUS_LINESTYLE,
                      model.refinement, direction)
        if item.separated:
            ax.add_patch(FancyArrowPatch(
                item.minus_center, item.plus_center, arrowstyle='-',
                connectionstyle='arc3,rad=0.25', linestyle=':',
                linewidth=_STYLE['connector_width'], color=item.color,
                zorder=2))


def _plot_domain_2d(model, ax, isolines, topology, patch_labels,
                    interface_labels, vertex_labels, legend):
    _draw_patches_2d(ax, model, isolines, topology, topology and patch_labels)
    if topology:
        _draw_interfaces_2d(ax, model)
        if interface_labels:
            _draw_interface_labels(ax, model)
        _draw_vertices(ax, model, vertex_labels)
        if legend:
            _draw_legend(ax, 2)
        ax.set_axis_off()
        ax.margins(0.03)
    else:
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    ax.set_aspect('equal', adjustable='datalim')


#==============================================================================
# 3D plots
#==============================================================================
def _surface_edges(surface):
    return (surface[0], surface[-1], surface[:, 0], surface[:, -1])


def _draw_patches_3d(ax, model, isolines, topology, labels):
    for i, patch in enumerate(model.patches):
        color = model.patch_colors[i]
        for (axis, ext), surface in model.faces[i].items():
            if not topology:
                ax.plot_surface(*np.moveaxis(surface, -1, 0), color=color,
                                alpha=_STYLE['surface_alpha'], linewidth=0,
                                shade=True)
            elif (i, axis, ext) in model.exterior_faces:
                ax.plot_surface(*np.moveaxis(surface, -1, 0), color=color,
                                alpha=_STYLE['face_alpha'], linewidth=0,
                                shade=False)
            for edge in _surface_edges(surface):
                ax.plot(*edge.T, color=_STYLE['edge_color'],
                        linewidth=_STYLE['edge_width'])
            if isolines:
                stride = max(1, model.refinement // 5)
                ax.plot_wireframe(*np.moveaxis(surface, -1, 0),
                                  rstride=stride, cstride=stride,
                                  color=_STYLE['isoline_color'],
                                  linewidth=_STYLE['isoline_width'],
                                  alpha=0.5)
        if labels:
            _text(ax, _patch_center(patch), _patch_name(patch),
                  fontsize=_STYLE['patch_fontsize'], fontweight='bold')


def _draw_side_arrows_3d(ax, face, color, linestyle, alpha, arrows):
    """Draw arrows on a face along its tangential axes.

    Each arrow is given by `(axis, direction, label)`: it follows the logical
    `axis` of the face, reversed if `direction` is -1, and is labelled with
    the coordinate eta of index `label`.
    """
    lower, upper = _logical_limits(face.domain)
    outward = 0.10 * (_face_center(face) - _patch_center(face.domain))
    for tangent, direction, label in arrows:
        start = (lower + upper) / 2
        start[face.axis] = (lower if face.ext == -1 else upper)[face.axis]
        end = start.copy()
        width = upper[tangent] - lower[tangent]
        start[tangent] = lower[tangent] + 0.25 * width
        end[tangent] = lower[tangent] + 0.72 * width
        start, end = _map_points(face.domain, [start, end])
        if direction == -1:
            start, end = end, start
        delta = end - start
        ax.quiver(*start, *delta, color=color, linestyles=linestyle,
                  linewidth=1.5, arrow_length_ratio=0.18, alpha=alpha)
        _text(ax, end + 0.08 * delta + outward, rf'$\eta_{label + 1}$',
              color=color, fontsize=_STYLE['label_fontsize'])


def _draw_interfaces_3d(ax, model):
    for item in model.interfaces:
        interface = item.interface
        surface = model.faces[model.patches.index(interface.minus.domain)][
            (int(interface.minus.axis), int(interface.minus.ext))]
        ax.plot_surface(*np.moveaxis(surface, -1, 0), color=item.color,
                        alpha=_STYLE['interface_alpha'], linewidth=0,
                        shade=False)
        # The arrows of the plus face show where the coordinates of the
        # minus face go, with the labels of the minus coordinates
        axis_map = interface.axis_map
        _draw_side_arrows_3d(
            ax, interface.minus, item.color, '-', 1.0,
            [(minus_axis, 1, minus_axis) for minus_axis, _, _ in axis_map])
        _draw_side_arrows_3d(
            ax, interface.plus, item.color, '--', 0.7,
            [(plus_axis, direction, minus_axis)
             for minus_axis, plus_axis, direction in axis_map])
        if item.separated:
            ax.plot(*np.column_stack((item.minus_center, item.plus_center)),
                    color=item.color, linestyle=':',
                    linewidth=_STYLE['connector_width'])


def _plot_domain_3d(model, ax, isolines, topology, patch_labels,
                    interface_labels, vertex_labels, legend):
    _draw_patches_3d(ax, model, isolines, topology, topology and patch_labels)
    if topology:
        _draw_interfaces_3d(ax, model)
        if interface_labels:
            _draw_interface_labels(ax, model)
        _draw_vertices(ax, model, vertex_labels)
        if legend:
            _draw_legend(ax, 3)
        ax.view_init(*_STYLE['view'])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    extent = np.maximum(model.extent, 1e-12)
    ax.set_box_aspect(np.maximum(extent,
                                 _STYLE['min_box_aspect'] * extent.max()))


#==============================================================================
# Entry point
#==============================================================================
def plot_domain(
        domain, draw=True, refinement=None, isolines=False, ax=None,
        topology=False,
        patch_labels=True, interface_labels=True, vertex_labels=True,
        legend=True):
    """Plot a 2D or 3D domain, optionally with topology annotations.

    This is the common plotting entry point for ordinary geometry plots and
    annotated multipatch topology plots. Both modes use the same sampling of
    the mapped patch edges and faces, and colour neighbouring patches
    differently.

    With ``topology=True``, each side of an interface is drawn inside its own
    patch, the minus side solid and the plus side dashed, in a colour shared
    by both sides. The arrows of the minus side show its increasing
    tangential coordinates, and those of the plus side where these go on the
    plus side: they point the other way for a reversing interface. In 3D, the
    arrows are labelled with the coordinates of the minus side.
    A dotted connector joins the two sides of an interface when they are drawn
    at different places, as for periodic or self-interfaces. The exterior
    boundary is drawn in black, and the vertices are marked: white for
    boundary vertices, black for interior ones. :func:`print_topology` gives
    the same information as text, including the axis map of every interface.

    Parameters
    ----------
    domain : sympde.topology.domain.Domain
        Two- or three-dimensional logical or mapped domain with equal logical
        and physical dimensions.
    draw : bool, default=True
        Call ``matplotlib.pyplot.show`` after constructing the figure.
    refinement : int or None
        Samples per edge. Defaults to 40 in 2D and 15 in 3D.
    isolines : bool, default=False
        Draw logical isolines or face wireframes.
    ax : matplotlib.axes.Axes, optional
        Existing 2D or 3D axes on which to draw.
    topology : bool, default=False
        Add symbolic multipatch topology annotations.
    patch_labels, interface_labels, vertex_labels, legend : bool
        Control annotations when ``topology=True``. They are ignored in the
        ordinary geometry mode.

    Returns
    -------
    matplotlib.figure.Figure
        The containing figure.
    """
    import matplotlib.pyplot as plt

    dim = int(domain.dim)
    if dim not in (2, 3):
        raise ValueError('domain plotting supports only dimensions 2 and 3')
    if _physical_dimension(domain) != dim:
        raise ValueError(
            'domain plotting currently requires equal logical and physical '
            'dimensions')
    if refinement is None:
        refinement = 40 if dim == 2 else 15
    if isinstance(refinement, bool) or not isinstance(
            refinement, numbers.Integral) or refinement < 2:
        raise ValueError('refinement must be an integer greater than one')

    model = _build_model(domain, int(refinement), topology)
    figure, ax = _new_axes(dim, ax, topology)
    plot = _plot_domain_2d if dim == 2 else _plot_domain_3d
    plot(model, ax, isolines, topology, patch_labels, interface_labels,
         vertex_labels, legend)
    ax.set_title(str(domain.name).replace('_', ' '))

    if draw:
        plt.show()
    return figure
