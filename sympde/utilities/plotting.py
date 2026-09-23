#---------------------------------------------------------------------------#
# This file is part of SymPDE.                                              #
#---------------------------------------------------------------------------#
"""Inspection and plotting helpers for symbolic multipatch topology."""

from dataclasses import dataclass
from itertools import product
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Patch
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
    domain : sympde.topology.Domain
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


def _logical_limits(patch):
    logical_patch = patch.logical_domain
    parameter_domain = logical_patch if logical_patch is not None else patch
    return (
        np.asarray(parameter_domain.min_coords, dtype=float),
        np.asarray(parameter_domain.max_coords, dtype=float),
    )


def _map_points(patch, logical_points):
    logical_points = np.asarray(logical_points, dtype=float)
    if patch.mapping is None:
        return logical_points
    mapping = patch.mapping.get_callable_mapping()
    return np.asarray(
        [mapping(*point) for point in logical_points], dtype=float)


def _physical_dimension(domain):
    dimensions = {
        int(patch.mapping.pdim) if patch.mapping is not None
        else int(patch.dim)
        for patch in domain.patches
    }
    if len(dimensions) != 1:
        raise ValueError('all plotted patches must have one physical dimension')
    return dimensions.pop()


def _patch_polygon(patch, samples):
    lower, upper = _logical_limits(patch)
    x_values = np.linspace(lower[0], upper[0], samples)
    y_values = np.linspace(lower[1], upper[1], samples)
    boundary = np.vstack((
        np.column_stack((x_values, np.full(samples, lower[1]))),
        np.column_stack((np.full(samples, upper[0]), y_values)),
        np.column_stack((x_values[::-1], np.full(samples, upper[1]))),
        np.column_stack((np.full(samples, lower[0]), y_values[::-1])),
    ))
    return _map_points(patch, boundary)


def _face_curve(face, samples):
    patch = face.domain
    lower, upper = _logical_limits(patch)
    tangent = 1 - int(face.axis)
    values = np.linspace(lower[tangent], upper[tangent], samples)
    points = np.empty((samples, 2))
    points[:, tangent] = values
    points[:, face.axis] = (
        lower[face.axis] if face.ext == -1 else upper[face.axis])
    return _map_points(patch, points)


def _corner_point(patch, logical_corner):
    lower, upper = _logical_limits(patch)
    point = np.where(np.asarray(logical_corner) == 0, lower, upper)
    return _map_points(patch, [point])[0]


def _draw_parameter_arrow(axis, curve, color, linestyle):
    axis.plot(
        curve[:, 0], curve[:, 1], color=color, linewidth=5,
        linestyle=linestyle, solid_capstyle='round', zorder=4,
    )
    start = curve[len(curve) * 3 // 10]
    end = curve[len(curve) * 7 // 10]
    axis.annotate(
        '', xy=end, xytext=start,
        arrowprops={
            'arrowstyle': '-|>', 'color': color, 'linewidth': 1.8,
        },
        zorder=6,
    )


def _draw_isolines(axis, patch, samples, count=9):
    lower, upper = _logical_limits(patch)
    fractions = np.linspace(0.0, 1.0, count)[1:-1]
    for fixed_axis in range(2):
        tangent = 1 - fixed_axis
        tangent_values = np.linspace(
            lower[tangent], upper[tangent], samples)
        for fraction in fractions:
            points = np.empty((samples, 2))
            points[:, tangent] = tangent_values
            points[:, fixed_axis] = (
                lower[fixed_axis]
                + fraction * (upper[fixed_axis] - lower[fixed_axis]))
            curve = _map_points(patch, points)
            axis.plot(
                curve[:, 0], curve[:, 1], color='darkgrey',
                linewidth=0.6, alpha=0.7, zorder=2,
            )


def _plot_domain_2d(
        domain, draw=True, refinement=40, isolines=False, ax=None,
        topology=False,
        patch_labels=True, interface_labels=True, vertex_labels=True,
        legend=True):
    """Plot a 2D domain, optionally with symbolic topology annotations.

    Interface sides use a common color. The minus side is solid, the plus
    side is dashed, and arrows show each side's increasing native tangential
    coordinate. A dotted connector is drawn when identified sides occupy
    different physical locations, which is useful for cut-open quotient
    topologies and self-identifications.

    Parameters
    ----------
    domain : sympde.topology.Domain
        Two-dimensional logical or mapped domain.
    draw : bool, default=True
        Call ``matplotlib.pyplot.show`` after constructing the figure.
    refinement : int, default=40
        Number of sample points used for each curved edge.
    isolines : bool, default=False
        Draw logical-coordinate isolines inside each patch.
    ax : matplotlib.axes.Axes, optional
        Existing axes on which to draw. New axes are created when omitted.
    patch_labels, interface_labels, vertex_labels, legend : bool
        Control topology annotations.

    Returns
    -------
    matplotlib.figure.Figure
        The containing figure.
    """
    _check_topology_dimension(domain)
    physical_dim = _physical_dimension(domain)
    if physical_dim != 2:
        raise ValueError('2D domain plotting requires physical dimension 2')
    if not isinstance(refinement, int) or refinement < 2:
        raise ValueError('refinement must be an integer greater than one')

    created_axes = ax is None
    if created_axes:
        if topology:
            figure, axis = plt.subplots(figsize=(10, 6.5))
        else:
            figure, axis = plt.subplots()
    else:
        axis = ax
        figure = axis.figure
        if hasattr(axis, 'zaxis'):
            raise ValueError('a 2D Matplotlib axes is required for a 2D domain')

    patches = domain.patches
    patch_by_index = dict(enumerate(patches))
    patch_colors = plt.get_cmap('Pastel1').colors
    interface_colors = plt.get_cmap('Dark2').colors
    legend_handles = []
    all_points = []

    for index, patch in enumerate(patches):
        polygon = _patch_polygon(patch, refinement)
        all_points.append(polygon)
        color = patch_colors[index % len(patch_colors)]
        if topology:
            axis.fill(
                polygon[:, 0], polygon[:, 1], facecolor=color,
                edgecolor='black', linewidth=1.2, alpha=0.8, zorder=1,
            )
        else:
            for normal_axis in range(2):
                for ext in (-1, 1):
                    face = patch.get_boundary(axis=normal_axis, ext=ext)
                    curve = _face_curve(face, refinement)
                    axis.plot(
                        curve[:, 0], curve[:, 1], color='black',
                        linewidth=1.0)
        if isolines:
            _draw_isolines(axis, patch, refinement)
        if topology and patch_labels:
            lower, upper = _logical_limits(patch)
            center = _map_points(patch, [(lower + upper) / 2.0])[0]
            axis.text(
                *center, _patch_name(patch), ha='center', va='center',
                fontsize=14, fontweight='bold', zorder=3,
            )
        if topology and legend:
            legend_handles.append(Patch(
                facecolor=color, edgecolor='black', label=_patch_name(patch)))

    point_extent = np.ptp(np.vstack(all_points), axis=0)
    scale = max(float(np.max(point_extent)), 1.0)
    interface_text = []
    interfaces = _interfaces(domain) if topology else ()
    for index, interface in enumerate(interfaces):
        color = interface_colors[index % len(interface_colors)]
        minus_curve = _face_curve(interface.minus, refinement)
        plus_curve = _face_curve(interface.plus, refinement)
        _draw_parameter_arrow(axis, minus_curve, color, '-')
        _draw_parameter_arrow(axis, plus_curve, color, '--')

        minus_midpoint = minus_curve[len(minus_curve) // 2]
        plus_midpoint = plus_curve[len(plus_curve) // 2]
        if interface_labels:
            axis.text(
                *minus_midpoint, f'I{index}-', color=color,
                fontsize=9, zorder=7)
            axis.text(
                *plus_midpoint, f'I{index}+', color=color,
                fontsize=9, zorder=7)

        if np.linalg.norm(minus_midpoint - plus_midpoint) > 1e-10 * scale:
            axis.add_patch(FancyArrowPatch(
                minus_midpoint, plus_midpoint,
                arrowstyle='-', connectionstyle='arc3,rad=0.25',
                linestyle=':', linewidth=1.8, color=color, zorder=2,
            ))

        interface_text.append(_format_interface_axis_map(index, interface))

    vertex_colors = plt.get_cmap('tab20').colors
    vertices = collect_topology_vertices(domain) if topology else ()
    for vertex in vertices:
        color = vertex_colors[vertex.index % len(vertex_colors)]
        plotted_locations = set()
        for incidence in vertex.incidences:
            point = _corner_point(
                patch_by_index[incidence.patch_index],
                incidence.logical_corner)
            location = tuple(np.round(point, decimals=12))
            if location in plotted_locations:
                continue
            plotted_locations.add(location)
            axis.scatter(
                *point, s=42, color=color, edgecolor='black', zorder=8)
            if vertex_labels:
                axis.annotate(
                    f'V{vertex.index}', point, xytext=(5, 5),
                    textcoords='offset points', fontsize=8,
                    color='black', zorder=9,
                )

    if topology and legend:
        axis.legend(
            handles=legend_handles, loc='upper left',
            bbox_to_anchor=(1.02, 1.0))
    if topology and interface_labels:
        axis.text(
            1.02, 0.72,
            'solid: minus side\n'
            'dashed: plus side\n'
            'arrows: increasing local coordinate\n\n'
            + '\n'.join(interface_text),
            transform=axis.transAxes, va='top', fontsize=9,
        )

    if topology:
        axis.set_title(str(domain.name).replace('_', ' '))
    axis.set_aspect('equal', adjustable='box')
    axis.set_xlabel('physical x' if topology else 'X')
    axis.set_ylabel('physical y' if topology else 'Y')
    if topology:
        axis.margins(0.12)
    if topology and created_axes and (legend or interface_labels):
        figure.subplots_adjust(right=0.72)

    if draw:
        plt.show()
    return figure


def _face_surface(face, samples):
    """Return mapped points on a 3D patch face and its tangent axes."""
    patch = face.domain
    lower, upper = _logical_limits(patch)
    normal_axis = int(face.axis)
    tangent_axes = tuple(axis for axis in range(3) if axis != normal_axis)
    first_values = np.linspace(
        lower[tangent_axes[0]], upper[tangent_axes[0]], samples)
    second_values = np.linspace(
        lower[tangent_axes[1]], upper[tangent_axes[1]], samples)
    first, second = np.meshgrid(
        first_values, second_values, indexing='ij')

    points = np.empty((samples * samples, 3))
    points[:, normal_axis] = (
        lower[normal_axis] if face.ext == -1 else upper[normal_axis])
    points[:, tangent_axes[0]] = first.ravel()
    points[:, tangent_axes[1]] = second.ravel()
    surface = _map_points(patch, points).reshape(samples, samples, 3)
    return surface, tangent_axes


def _surface_boundary_curves(surface):
    return (
        surface[0, :, :],
        surface[-1, :, :],
        surface[:, 0, :],
        surface[:, -1, :],
    )


def _face_center(face):
    lower, upper = _logical_limits(face.domain)
    point = (lower + upper) / 2.0
    point[face.axis] = (
        lower[face.axis] if face.ext == -1 else upper[face.axis])
    return _map_points(face.domain, [point])[0]


def _face_parameter_arrow(face, tangent_axis):
    lower, upper = _logical_limits(face.domain)
    start = (lower + upper) / 2.0
    end = start.copy()
    start[face.axis] = (
        lower[face.axis] if face.ext == -1 else upper[face.axis])
    end[face.axis] = start[face.axis]
    width = upper[tangent_axis] - lower[tangent_axis]
    start[tangent_axis] = lower[tangent_axis] + 0.25 * width
    end[tangent_axis] = lower[tangent_axis] + 0.72 * width
    return _map_points(face.domain, [start, end])


def _draw_3d_parameter_arrows(
        axis, face, color, linestyle, alpha, label_at_start_axes=()):
    lower, upper = _logical_limits(face.domain)
    patch_center = _map_points(face.domain, [(lower + upper) / 2.0])[0]
    outward_offset = 0.10 * (_face_center(face) - patch_center)

    for tangent_axis in range(3):
        if tangent_axis == face.axis:
            continue
        start, end = _face_parameter_arrow(face, tangent_axis)
        delta = end - start
        axis.plot(
            [start[0], end[0]], [start[1], end[1]],
            [start[2], end[2]],
            color=color, linestyle=linestyle, linewidth=1.8,
            alpha=alpha, zorder=7,
        )
        axis.quiver(
            start[0], start[1], start[2],
            delta[0], delta[1], delta[2],
            color=color, linewidth=1.3, arrow_length_ratio=0.18,
            alpha=alpha, zorder=8,
        )
        if tangent_axis in label_at_start_axes:
            label_position = start - 0.08 * delta + outward_offset
        else:
            label_position = end + 0.08 * delta + outward_offset
        axis.text(
            label_position[0], label_position[1], label_position[2],
            f'eta{tangent_axis + 1}',
            color=color, fontsize=7, zorder=9,
        )


def _format_interface_axis_map(index, interface):
    mappings = '; '.join(
        f'eta{minus_axis + 1} -> '
        f'{"+" if direction == 1 else "-"}eta{plus_axis + 1}'
        for minus_axis, plus_axis, direction in interface.axis_map
    )
    return (
        f'I{index}: {_patch_name(interface.minus.domain)} '
        f'(axis={interface.minus.axis}, ext={int(interface.minus.ext):+d}) '
        f'<-> {_patch_name(interface.plus.domain)} '
        f'(axis={interface.plus.axis}, ext={int(interface.plus.ext):+d})\n'
        f'    axis map: {mappings}'
    )


def _plot_domain_3d(
        domain, draw=True, refinement=15, isolines=False, ax=None,
        topology=False,
        patch_labels=True, interface_labels=True, vertex_labels=True,
        legend=True):
    """Plot a 3D domain, optionally with symbolic topology annotations."""
    if not isinstance(refinement, int) or refinement < 2:
        raise ValueError('refinement must be an integer greater than one')

    created_axes = ax is None
    if created_axes:
        figure = plt.figure(figsize=(12, 8) if topology else None)
        axis = figure.add_subplot(111, projection='3d')
    else:
        axis = ax
        figure = axis.figure
        if not hasattr(axis, 'zaxis'):
            raise ValueError('a 3D Matplotlib axes is required for a 3D domain')

    patches = domain.patches
    patch_by_index = dict(enumerate(patches))
    patch_colors = plt.get_cmap('Pastel1').colors
    interface_colors = plt.get_cmap('Dark2').colors
    legend_handles = []
    all_points = []

    for patch_index, patch in enumerate(patches):
        color = patch_colors[patch_index % len(patch_colors)]
        for normal_axis in range(3):
            for ext in (-1, 1):
                face = patch.get_boundary(axis=normal_axis, ext=ext)
                surface, _ = _face_surface(face, refinement)
                all_points.append(surface.reshape(-1, 3))
                if topology:
                    axis.plot_surface(
                        surface[:, :, 0], surface[:, :, 1], surface[:, :, 2],
                        color=color, alpha=0.10, linewidth=0,
                        shade=False, zorder=1,
                    )
                    for curve in _surface_boundary_curves(surface):
                        axis.plot(
                            curve[:, 0], curve[:, 1], curve[:, 2],
                            color='0.35', linewidth=0.7, alpha=0.65,
                            zorder=2,
                        )
                else:
                    axis.plot_surface(
                        surface[:, :, 0], surface[:, :, 1], surface[:, :, 2],
                        color='c', alpha=0.7,
                    )
                if isolines:
                    stride = max(1, refinement // 5)
                    axis.plot_wireframe(
                        surface[:, :, 0], surface[:, :, 1],
                        surface[:, :, 2],
                        rstride=stride, cstride=stride,
                        color='0.45', linewidth=0.35, alpha=0.35,
                        zorder=2,
                    )

        if topology and patch_labels:
            lower, upper = _logical_limits(patch)
            center = _map_points(patch, [(lower + upper) / 2.0])[0]
            axis.text(
                center[0], center[1], center[2], _patch_name(patch),
                ha='center', va='center', fontsize=12,
                fontweight='bold', zorder=5,
            )
        if topology and legend:
            legend_handles.append(Patch(
                facecolor=color, edgecolor='0.35',
                label=_patch_name(patch)))

    point_cloud = np.vstack(all_points)
    point_extent = np.ptp(point_cloud, axis=0)
    scale = max(float(np.max(point_extent)), 1.0)
    interface_text = []
    interfaces = _interfaces(domain) if topology else ()
    for index, interface in enumerate(interfaces):
        color = interface_colors[index % len(interface_colors)]
        minus_surface, _ = _face_surface(interface.minus, refinement)
        plus_surface, _ = _face_surface(interface.plus, refinement)

        for surface, alpha, linestyle in (
                (minus_surface, 0.34, '-'), (plus_surface, 0.20, '--')):
            axis.plot_surface(
                surface[:, :, 0], surface[:, :, 1], surface[:, :, 2],
                color=color, alpha=alpha, linewidth=0,
                shade=False, zorder=4,
            )
            for curve in _surface_boundary_curves(surface):
                axis.plot(
                    curve[:, 0], curve[:, 1], curve[:, 2],
                    color=color, linestyle=linestyle, linewidth=2.0,
                    alpha=0.9, zorder=5,
                )

        _draw_3d_parameter_arrows(
            axis, interface.minus, color, '-', 1.0)
        plus_labels_at_start = {
            plus_axis
            for _, plus_axis, direction in interface.axis_map
            if direction == 1
        }
        _draw_3d_parameter_arrows(
            axis, interface.plus, color, '--', 0.72,
            label_at_start_axes=plus_labels_at_start)

        minus_center = _face_center(interface.minus)
        plus_center = _face_center(interface.plus)
        centers_coincide = (
            np.linalg.norm(minus_center - plus_center) <= 1e-10 * scale)
        if interface_labels:
            if centers_coincide:
                axis.text(
                    minus_center[0], minus_center[1], minus_center[2],
                    f'I{index} -/+', color=color, fontsize=9,
                    fontweight='bold', zorder=10)
            else:
                axis.text(
                    minus_center[0], minus_center[1], minus_center[2],
                    f'I{index}-', color=color, fontsize=9,
                    fontweight='bold', zorder=10)
                axis.text(
                    plus_center[0], plus_center[1], plus_center[2],
                    f'I{index}+', color=color, fontsize=9,
                    fontweight='bold', zorder=10)
                axis.plot(
                    [minus_center[0], plus_center[0]],
                    [minus_center[1], plus_center[1]],
                    [minus_center[2], plus_center[2]],
                    color=color, linestyle=':', linewidth=1.7,
                    alpha=0.8, zorder=3,
                )
        interface_text.append(_format_interface_axis_map(index, interface))

    vertex_colors = plt.get_cmap('tab20').colors
    vertices = collect_topology_vertices(domain) if topology else ()
    for vertex in vertices:
        color = vertex_colors[vertex.index % len(vertex_colors)]
        plotted_locations = set()
        for incidence in vertex.incidences:
            point = _corner_point(
                patch_by_index[incidence.patch_index],
                incidence.logical_corner)
            location = tuple(np.round(point, decimals=12))
            if location in plotted_locations:
                continue
            plotted_locations.add(location)
            axis.scatter(
                point[0], point[1], point[2], s=28, color=color,
                edgecolor='black', depthshade=False, zorder=11)
            if vertex_labels:
                axis.text(
                    point[0], point[1], point[2], f'V{vertex.index}',
                    fontsize=7, color='black', zorder=12)

    if topology and legend:
        axis.legend(
            handles=legend_handles, loc='upper left',
            bbox_to_anchor=(1.02, 1.0))
    if topology and interface_labels:
        axis.text2D(
            1.02, 0.72,
            'solid: minus face\n'
            'dashed: plus face\n'
            'arrows: native tangential coordinates\n\n'
            + '\n'.join(interface_text),
            transform=axis.transAxes, va='top', fontsize=8,
        )

    if topology:
        axis.set_title(str(domain.name).replace('_', ' '))
    axis.set_xlabel('physical x' if topology else 'X')
    axis.set_ylabel('physical y' if topology else 'Y')
    axis.set_zlabel('physical z' if topology else 'Z')
    axis.set_box_aspect(np.maximum(point_extent, 1e-12))
    if topology:
        axis.view_init(elev=24, azim=-58)
    if topology and created_axes and (legend or interface_labels):
        figure.subplots_adjust(right=0.72)

    if draw:
        plt.show()
    return figure


def plot_domain(
        domain, draw=True, refinement=None, isolines=False, ax=None,
        topology=False,
        patch_labels=True, interface_labels=True, vertex_labels=True,
        legend=True):
    """Plot a 2D or 3D domain, optionally with topology annotations.

    This is the common plotting entry point for ordinary geometry plots and
    annotated multipatch topology plots. Both modes use the same mapped-edge
    and mapped-face sampling. Setting ``topology=True`` adds patch labels,
    oriented interface sides, parameter arrows, interface descriptions, and
    topological vertices.

    In topology mode, the 2D view draws oriented interface edges and cut-open
    seam connectors. The lightweight 3D view draws translucent patch
    boundaries, highlights both sides of every interface, displays native
    tangential-coordinate arrows and the signed ``axis_map``, and labels
    topological vertices.

    Parameters
    ----------
    domain : sympde.topology.Domain
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
    logical_dim = int(domain.dim)
    if logical_dim not in (2, 3):
        raise ValueError('domain plotting supports only dimensions 2 and 3')
    physical_dim = _physical_dimension(domain)
    if logical_dim != physical_dim:
        raise ValueError(
            'domain plotting currently requires equal logical and physical '
            'dimensions')

    if not topology:
        patch_labels = False
        interface_labels = False
        vertex_labels = False
        legend = False

    options = dict(
        domain=domain,
        draw=draw,
        isolines=isolines,
        ax=ax,
        patch_labels=patch_labels,
        interface_labels=interface_labels,
        vertex_labels=vertex_labels,
        legend=legend,
        topology=topology,
    )
    if logical_dim == 2:
        options['refinement'] = 40 if refinement is None else refinement
        return _plot_domain_2d(**options)

    options['refinement'] = 15 if refinement is None else refinement
    return _plot_domain_3d(**options)
