import numpy as np
import pytest

from sympde.topology.analytical_mapping import AffineMapping, IdentityMapping
from sympde.topology.basic import Boundary, InteriorDomain, Interface
from sympde.topology.domain import Cube, Domain, Square
from sympde.topology.incidence import IndexedInterface, PatchSide
from sympde.topology.multipatch_gallery import (
    build_multipatch_domain_2d,
    build_multipatch_domain_3d,
)


def _make_oriented_cubes():
    patch_a = IdentityMapping('F_A3', dim=3)(Cube('A3'))
    patch_b = AffineMapping(
        'F_B3', dim=3,
        c1=1, c2=0, c3=1,
        a11=0, a12=1, a13=0,
        a21=0, a22=0, a23=1,
        a31=-1, a32=0, a33=0,
    )(Cube('B3'))
    return Domain.join(
        [patch_a, patch_b],
        [((0, 0, 1), (1, 1, -1), (-1, +1, -1))],
        'oriented_cubes',
    )


def _physical_point(patch, logical_point):
    """Map a point of the unit box to the physical space of a patch."""
    lower = np.asarray(patch.logical_domain.min_coords, dtype=float)
    upper = np.asarray(patch.logical_domain.max_coords, dtype=float)
    point = lower + np.asarray(logical_point, dtype=float) * (upper - lower)
    mapping = patch.mapping.get_callable_mapping()
    return np.asarray(mapping(*point), dtype=float)

#==============================================================================
def test_vertices_shared_and_self_interfaces():
    vertices = build_multipatch_domain_2d('three_patch_shared_vertex').vertices
    interior_vertices = [vertex for vertex in vertices if not vertex.is_boundary]

    assert len(interior_vertices) == 1
    center = interior_vertices[0]
    assert center.patch_indices == (0, 1, 2)
    assert all(
        incidence.logical_corner == (0, 0)
        for incidence in center.incidences)

    patch = Square('P')
    self_connected = Domain.join(
        [patch],
        [((0, 0, -1), (0, 0, 1), +1)],
        'self_connected',
    )
    self_vertices = self_connected.vertices
    assert len(self_vertices) == 2
    assert all(len(vertex.incidences) == 2 for vertex in self_vertices)
    assert all(vertex.patches == ('P',) for vertex in self_vertices)
    assert all(vertex.is_boundary for vertex in self_vertices)


def test_vertices_3d():
    vertices = _make_oriented_cubes().vertices
    shared_vertices = [
        vertex for vertex in vertices if len(vertex.incidences) == 2]

    assert len(vertices) == 12
    assert len(shared_vertices) == 4
    assert all(vertex.patches == ('A3', 'B3') for vertex in shared_vertices)

#==============================================================================
def test_edges_around_an_interior_edge():
    domain = build_multipatch_domain_3d('three_patch_shared_edge')
    edges = domain.edges

    # 3 cubes have 36 edges, 3 interfaces identify 4 pairs of them each, and
    # the central edge is part of all three interfaces
    assert len(edges) == 25
    central = [edge for edge in edges if len(edge.incidences) == 3]
    assert len(central) == 1
    assert central[0].patch_indices == (0, 1, 2)
    assert not central[0].is_boundary
    assert all(edge.is_boundary for edge in edges if edge is not central[0])


def test_edges_directions():
    edges = _make_oriented_cubes().edges
    shared = [edge for edge in edges if len(edge.incidences) == 2]

    assert len(edges) == 20
    assert len(shared) == 4
    # The first incidence of an edge is the reference for the directions
    assert all(edge.incidences[0].direction == 1 for edge in edges)
    # The orientation (-1, +1, -1) reverses the second tangential axis of A3,
    # which maps to the first axis of B3
    reversed_edges = {
        (edge.incidences[0].axis, edge.incidences[1].axis)
        for edge in shared if edge.incidences[1].direction == -1}
    assert reversed_edges == {(2, 0)}

    # In 2D, the edges of the patches are their sides, given by the interfaces
    strip = build_multipatch_domain_2d('two_patch_twisted_strip')
    with pytest.raises(ValueError, match='indexed_interfaces'):
        _ = strip.edges


@pytest.mark.parametrize('name', [
    'two_cube', 'oriented_two_cube', 'three_patch_shared_edge', 'torus'])
def test_vertices_and_edges_coincide_geometrically(name):
    domain = build_multipatch_domain_3d(name)

    for vertex in domain.vertices:
        points = [
            _physical_point(domain.patches[i.patch_index], i.logical_corner)
            for i in vertex.incidences]
        assert all(np.allclose(p, points[0], atol=1e-10) for p in points)

    for edge in domain.edges:
        segments = []
        for incidence in edge.incidences:
            patch = domain.patches[incidence.patch_index]
            ends = [
                [value if value is not None else end
                 for value in incidence.logical_edge]
                for end in (0, 1)]
            start, stop = (_physical_point(patch, end) for end in ends)
            if incidence.direction == -1:
                start, stop = stop, start
            segments.append((start, stop))
        first = segments[0]
        assert all(np.allclose(s, first[0], atol=1e-10) and
                   np.allclose(t, first[1], atol=1e-10)
                   for s, t in segments)


def test_topology_requires_interface_axes():
    A = InteriorDomain('A', ldim=2)
    B = InteriorDomain('B', ldim=2)
    domain = Domain('Omega', interiors=[A, B], connectivity=[
        Interface('I', Boundary('G1', A), Boundary('G2', B), 1)])

    with pytest.raises(ValueError, match='no axis or extremity'):
        _ = domain.indexed_interfaces
    with pytest.raises(ValueError, match='no axis or extremity'):
        _ = domain.vertices

#==============================================================================
def test_patch_index():
    A, B = Square('A'), Square('B')
    F = IdentityMapping('F', dim=2)
    G = IdentityMapping('G', dim=2)
    domain = Domain.join([F(A), G(B)], [((0, 0, 1), (1, 0, -1), 1)], 'Omega')
    interface = domain.interface_map['A|B']

    assert domain.patch_index(domain.patches[1]) == 1
    assert domain.patch_index(F(A)) == 0
    assert domain.patch_index(B.interior) == 1
    assert domain.patch_index(interface.minus) == 0
    assert domain.patch_index(interface.plus) == 1

    with pytest.raises(TypeError, match='two patches'):
        domain.patch_index(interface)
    with pytest.raises(ValueError, match='not a patch'):
        domain.patch_index(Square('C'))


def test_indexed_interfaces():
    A, B = Cube('A'), Cube('B')
    domain = Domain.join(
        [A, B],
        [((A, 0, +1), (B, 1, -1), (-1, +1, -1)),
         ((A, 2, +1), (B, 2, -1), (+1, +1, +1), 'top')],
        'Omega')
    interfaces = domain.indexed_interfaces

    assert tuple(interfaces) == tuple(domain.interface_map) == ('A|B', 'top')
    assert interfaces['A|B'] == IndexedInterface(
        name='A|B',
        minus=PatchSide(patch_index=0, axis=0, ext=+1),
        plus=PatchSide(patch_index=1, axis=1, ext=-1),
        orientation=(-1, +1, -1),
        axis_map=((1, 2, +1), (2, 0, -1)),
    )
    assert interfaces['top'].minus == PatchSide(0, 2, +1)

    # Every face of the two cubes is an exterior side or a side of one
    # interface
    exterior = domain.indexed_exterior_sides
    assert len(exterior) == 2 * 6 - 2 * len(interfaces)
    assert len(set(exterior)) == len(exterior)
    glued = {side for i in interfaces.values() for side in (i.minus, i.plus)}
    assert set(exterior) | glued == {
        PatchSide(patch, axis, ext)
        for patch in range(2) for axis in range(3) for ext in (-1, 1)}
    assert exterior == tuple(
        PatchSide(domain.patch_index(side), side.axis, side.ext)
        for side in domain.exterior_sides)

    with pytest.raises(TypeError):
        interfaces['other'] = interfaces['top']
