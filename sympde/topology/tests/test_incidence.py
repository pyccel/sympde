import itertools

import numpy as np
from sympy import lambdify
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
    with pytest.raises(ValueError, match='not a patch'):
        domain.patch_index([A])


def test_cached_topology():
    A, B = Square('A'), Square('B')
    domain = Domain.join([A, B], [((0, 0, 1), (1, 0, -1), 1)], 'Omega')

    # Computed once, then shared
    assert domain.indexed_interfaces is domain.indexed_interfaces
    assert domain.patch_handedness is domain.patch_handedness
    assert domain.patch_handedness == (1, 1)

    # Rebuilding the domain gives the same topology
    rebuilt = domain.func(*domain.args)
    assert rebuilt.indexed_interfaces == domain.indexed_interfaces
    assert rebuilt.patch_index(B) == domain.patch_index(B)


def test_patch_index_of_logical_patch():
    # The mapped patches sort as (F1(C), F2(B), F3(A)), in the reverse order
    # of the logical patches (A, B, C)
    A, B, C = Square('A'), Square('B'), Square('C')
    F1, F2, F3 = (IdentityMapping(f'F{i}', dim=2) for i in (1, 2, 3))
    domain = Domain.join(
        [F3(A), F2(B), F1(C)],
        [((0, 0, 1), (1, 0, -1), 1),
         ((1, 0, 1), (2, 0, -1), 1)],
        'Omega')

    for logical, mapped in [(A, F3(A)), (B, F2(B)), (C, F1(C))]:
        index = domain.patch_index(mapped)
        assert domain.patch_index(logical) == index
        assert domain.patch_index(logical.interior) == index
        assert domain.patches[index].logical_domain == logical.interior

    # A has a boundary side at (0, -1), C has an interface there
    assert domain.get_boundary(0, -1, patch=A) == domain.get_boundary(
        0, -1, patch=F3(A))
    assert domain.get_boundary(0, -1, patch=A).domain == F3(A).interior
    with pytest.raises(ValueError, match='could not find boundary'):
        domain.get_boundary(0, -1, patch=C)

    assert domain.get_interface(A, B) == (domain.interface_map['A|B'],)
    assert domain.get_interface(A, C) == ()


def test_get_interface():
    A, B, C = Square('A'), Square('B'), Square('C')
    domain = Domain.join(
        [A, B, C],
        [((0, 0, 1), (1, 0, -1), 1),
         ((0, 1, 1), (1, 1, -1), 1),
         ((1, 1, +1), (2, 0, -1), 1),
         ((2, 1, -1), (2, 1, 1), 1)],
        'Omega')

    # Two interfaces between A and B, in either order and any form
    assert domain.get_interface(A, B) == (
        domain.interface_map['A|B'], domain.interface_map['A|B#2'])
    assert domain.get_interface(B.interior, A) == domain.get_interface(A, B)
    assert domain.get_interface(B, C) == (domain.interface_map['B|C'],)

    # A self-interface, and no interface
    assert domain.get_interface(C, C) == (domain.interface_map['C|C'],)
    assert domain.get_interface(A, C) == ()


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

    # Every face of the two cubes is a boundary side or a side of one
    # interface
    exterior = domain.indexed_boundary_sides
    assert len(exterior) == 2 * 6 - 2 * len(interfaces)
    assert len(set(exterior)) == len(exterior)
    glued = {side for i in interfaces.values() for side in (i.minus, i.plus)}
    assert set(exterior) | glued == {
        PatchSide(patch, axis, ext)
        for patch in range(2) for axis in range(3) for ext in (-1, 1)}
    assert exterior == tuple(
        PatchSide(domain.patch_index(side), side.axis, side.ext)
        for side in domain.boundary_sides)

    with pytest.raises(TypeError):
        interfaces['other'] = interfaces['top']

#==============================================================================
@pytest.mark.parametrize('ldim, name', [
    (2, 'square_9'), (2, 'annulus_4'), (2, 'curved_L_shape'), (2, 'pretzel'),
    (2, 'pretzel_f'), (2, 'pretzel_annulus'), (2, 'self_connected_patch'),
    (2, 'three_patch_shared_vertex'),
    (3, 'two_cube'), (3, 'torus'), (3, 'oriented_two_cube'),
    (3, 'three_patch_shared_edge'),
])
def test_patch_handedness_follows_jacobian_determinants(ldim, name):
    build = {2: build_multipatch_domain_2d, 3: build_multipatch_domain_3d}[ldim]
    domain = build(name)

    signs = []
    for patch in domain.patches:
        logical_patch = patch.logical_domain
        center = [(lower + upper) / 2 for lower, upper
                  in zip(logical_patch.min_coords, logical_patch.max_coords)]
        mapping = patch.mapping
        jacobian = lambdify(mapping.logical_coordinates, mapping.jacobian_expr)(*center)
        signs.append(int(np.sign(np.linalg.det(np.array(jacobian, dtype=float)))))

    assert domain.is_orientable
    assert domain.patch_handedness == domain.logical_domain.patch_handedness
    assert domain.patch_handedness == tuple(signs[0] * sign for sign in signs)

#==============================================================================
def test_patch_handedness_non_orientable():
    # The twisted strip is a Möbius strip made of two patches
    strip = build_multipatch_domain_2d('two_patch_twisted_strip')
    assert not strip.is_orientable
    assert strip.patch_handedness is None
    # A non-orientable domain stays non-orientable once cached
    assert strip.patch_handedness is None

    # A square whose opposite sides are identified with the same or the
    # reversed direction is a cylinder or a Möbius strip
    cylinder = Domain.join([Square('P')], [((0, 0, -1), (0, 0, +1), +1)], 'cylinder')
    moebius  = Domain.join([Square('P')], [((0, 0, -1), (0, 0, +1), -1)], 'moebius')
    assert cylinder.patch_handedness == (1,)
    assert moebius.patch_handedness is None

    # Identifying a face with its opposite through a reflection
    cube = Domain.join([Cube('P')], [((0, 0, -1), (0, 0, +1), (+1, +1, -1))], 'cube')
    assert not cube.is_orientable

    # Each connected component starts with the sign +1. Joining two upper
    # faces, or a face normal to x2 with one normal to x1 without reversing
    # the tangent, reflects the logical axes.
    patches = [Square('A'), Square('B'), Square('C'), Square('D')]
    domain = Domain.join(patches, [((0, 0, +1), (1, 0, +1), +1),
                                   ((2, 1, +1), (3, 0, -1), +1)], 'two_components')
    assert [interface.handedness for interface in domain.interface_map.values()] == [-1, -1]
    assert domain.patch_handedness == (1, -1, 1, -1)

#==============================================================================
@pytest.mark.parametrize('ldim, name', [
    (2, 'square_9'), (2, 'pretzel_f'), (2, 'two_patch_twisted_strip'),
    (2, 'self_connected_patch'), (2, 'three_patch_shared_vertex'),
    (3, 'two_cube'), (3, 'torus'), (3, 'oriented_two_cube'),
    (3, 'three_patch_shared_edge'),
])
def test_reverse_lookups(ldim, name):
    build = {2: build_multipatch_domain_2d, 3: build_multipatch_domain_3d}[ldim]
    domain = build(name)
    npatches = len(domain.patches)

    # Every corner of every patch belongs to exactly one vertex
    corners = list(itertools.product((0, 1), repeat=ldim))
    for patch_index in range(npatches):
        for corner in corners:
            vertex = domain.vertex_of(patch_index, corner)
            assert (patch_index, corner) in {(i.patch_index, i.logical_corner)
                                             for i in vertex.incidences}
    assert sum(len(v.incidences) for v in domain.vertices) == npatches * len(corners)

    # Every side of every patch is a side of one interface or a boundary side
    for patch_index in range(npatches):
        for axis, ext in itertools.product(range(ldim), (-1, +1)):
            side = PatchSide(patch_index, axis, ext)
            entity = domain.side_of(patch_index, axis, ext)
            if isinstance(entity, IndexedInterface):
                assert side in (entity.minus, entity.plus)
                assert entity in domain.indexed_interfaces.values()
            else:
                assert entity == side
                assert entity in domain.indexed_boundary_sides

    if ldim == 3:
        for edge in domain.edges:
            for incidence in edge.incidences:
                assert domain.edge_of(incidence.patch_index, incidence.logical_edge) \
                    == (edge, incidence.direction)
        assert sum(len(e.incidences) for e in domain.edges) == npatches * 12
    else:
        with pytest.raises(ValueError, match='3D'):
            domain.edge_of(0, (None, 0))

#==============================================================================
def test_reverse_lookups_with_repeated_patch():
    # A cube whose opposite faces are all identified is a 3-torus made of a
    # single patch: all its corners are one vertex, and the four parallel
    # edges along every axis are one edge.
    cube = Domain.join([Cube('P')], [((0, axis, -1), (0, axis, +1), (+1, +1, +1))
                                     for axis in range(3)], 'three_torus')

    vertex, = cube.vertices
    assert vertex.patch_indices == (0,)
    assert len(vertex.incidences) == 8
    assert all(cube.vertex_of(0, corner) is vertex
               for corner in itertools.product((0, 1), repeat=3))

    assert len(cube.edges) == 3
    for axis in range(3):
        logical_edges = [tuple(None if i == axis else value for i, value in
                               enumerate(position)) for position in
                         itertools.product((0, 1), repeat=3) if position[axis] == 0]
        edges = {cube.edge_of(0, logical_edge) for logical_edge in logical_edges}
        assert len(edges) == 1
        (edge, direction), = edges
        assert direction == 1
        assert edge.patch_indices == (0,) and len(edge.incidences) == 4

    # Both sides of a self-interface belong to the same interface
    assert cube.side_of(0, 1, -1) is cube.side_of(0, 1, +1)

#==============================================================================
def test_reverse_lookups_reject_unknown_entities():
    domain = build_multipatch_domain_3d('two_cube')

    with pytest.raises(ValueError, match='corner'):
        domain.vertex_of(0, (0, 0))
    with pytest.raises(ValueError, match='corner'):
        domain.vertex_of(2, (0, 0, 0))
    with pytest.raises(ValueError, match='edge'):
        domain.edge_of(0, (0, 0, 0))
    with pytest.raises(ValueError, match='side'):
        domain.side_of(0, 3, -1)
