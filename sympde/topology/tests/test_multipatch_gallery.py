import numpy as np
import pytest

from sympde.topology import TransposedPolarMapping
from sympde.topology.multipatch_gallery import (
    build_annulus_3,
    build_annulus_4,
    build_cartesian_multipatch_domain_2d,
    build_curved_l_shape,
    build_multipatch_domain_2d,
    build_multipatch_domain_3d,
    build_oriented_two_cube,
    build_pretzel,
    build_pretzel_annulus,
    build_pretzel_debug,
    build_pretzel_f,
    build_self_connected_patch,
    build_square_2,
    build_square_4,
    build_square_6,
    build_square_8,
    build_square_9,
    build_three_patch_shared_edge,
    build_three_patch_shared_vertex,
    build_torus_3d,
    build_two_cube_3d,
    build_two_patch_twisted_strip,
    main,
    plot_multipatch_domain,
)


def _physical_point(side, parameter):
    patch = side.domain
    logical_patch = patch.logical_domain
    logical_point = [
        minimum + parameter * (maximum - minimum)
        for minimum, maximum in zip(
            logical_patch.min_coords, logical_patch.max_coords)
    ]
    logical_point[side.axis] = (
        logical_patch.min_coords[side.axis]
        if side.ext == -1 else logical_patch.max_coords[side.axis]
    )
    substitutions = dict(zip(
        patch.mapping.logical_coordinates, logical_point))
    return np.array([
        float(expression.subs(substitutions))
        for expression in patch.mapping.expressions
    ])


@pytest.mark.parametrize(
    'builder, name, npatches, ninterfaces, nboundaries',
    [
        (build_square_2,         'square_2',          2,  1,  6),
        (build_square_4,         'square_4',          4,  4,  8),
        (build_square_6,         'square_6',          6,  7, 10),
        (build_square_8,         'square_8',          8,  8, 16),
        (build_square_9,         'square_9',          9, 12, 12),
        (build_annulus_3,        'annulus_3',         3,  3,  6),
        (build_annulus_4,        'annulus_4',         4,  4,  8),
        (build_curved_l_shape,   'curved_L_shape',    3,  2,  8),
        (build_pretzel,          'pretzel',          11, 13, 18),
        (build_pretzel_f,        'pretzel_f',        18, 20, 32),
        (build_pretzel_annulus,  'pretzel_annulus',   9,  9, 18),
        (build_pretzel_debug,    'pretzel_debug',     2,  1,  6),
    ],
)
def test_multipatch_domain_builder(
        builder, name, npatches, ninterfaces, nboundaries):
    domain = builder()

    assert str(domain.name) == name
    assert len(domain.patches) == npatches
    assert len(domain.interface_map) == ninterfaces
    assert len(domain.exterior_sides) == nboundaries
    assert all(interface.orientation == +1
               for interface in domain.interface_map.values())
    for interface in domain.interface_map.values():
        for parameter in (0, 0.5, 1):
            np.testing.assert_allclose(
                _physical_point(interface.minus, parameter),
                _physical_point(interface.plus, parameter),
                atol=1e-14,
            )


@pytest.mark.parametrize('builder', [
    build_square_2,
    build_square_6,
    build_square_8,
    build_square_9,
    build_annulus_3,
    build_annulus_4,
    build_curved_l_shape,
    build_pretzel,
    build_pretzel_f,
    build_pretzel_annulus,
    build_pretzel_debug,
])
def test_numbered_gallery_domains_use_patch_names(builder):
    names = tuple(str(patch.logical_domain.name) for patch in builder().patches)

    assert len(names) == len(set(names))
    assert all(name.startswith('P') for name in names)


@pytest.mark.parametrize('mapping', ['identity', 'polar'])
def test_build_rectangular_cartesian_multipatch_domain(mapping):
    layout = np.array([
        [1, None, 5],
        [2,    3, 4],
    ], dtype=object)

    domain = build_cartesian_multipatch_domain_2d(
        layout, (0, 3), (0, 2), mapping=mapping)

    assert tuple(str(p.logical_domain.name) for p in domain.patches) == (
        'Log_0_0', 'Log_0_2', 'Log_1_0', 'Log_1_1', 'Log_1_2')
    assert len(domain.interface_map) == 4
    assert len(domain.exterior_sides) == 12

    patches = {str(p.logical_domain.name): p for p in domain.patches}
    top_left = patches['Log_0_0'].logical_domain
    bottom_right = patches['Log_1_2'].logical_domain
    assert top_left.min_coords == (0, 1)
    assert top_left.max_coords == (1, 2)
    assert bottom_right.min_coords == (2, 0)
    assert bottom_right.max_coords == (3, 1)


@pytest.mark.parametrize('layout, message', [
    ([1, 2], 'two-dimensional'),
    (np.empty((0, 2)), 'at least one layout position'),
    ([[None]], 'at least one non-None entry'),
])
def test_build_cartesian_multipatch_domain_2d_rejects_invalid_layout(
        layout, message):
    with pytest.raises(ValueError, match=message):
        build_cartesian_multipatch_domain_2d(layout, (0, 1), (0, 1))


def test_build_cartesian_multipatch_domain_2d_rejects_unknown_mapping():
    with pytest.raises(ValueError, match="'identity' or 'polar'"):
        build_cartesian_multipatch_domain_2d(
            [[1]], (0, 1), (0, 1), mapping='unknown')


def test_transposed_polar_mapping_is_public():
    mapping = TransposedPolarMapping(
        'F', dim=2, c1=0, c2=0, rmin=0, rmax=1)

    angle, radius = mapping.logical_coordinates
    x, y = mapping.expressions
    assert x.subs({angle: 0, radius: 2}) == 2
    assert y.subs({angle: 0, radius: 2}) == 0


@pytest.mark.parametrize('name', ['unknown', None])
def test_build_multipatch_domain_2d_rejects_unknown_name(name):
    with pytest.raises(ValueError, match='unknown 2D multipatch domain'):
        build_multipatch_domain_2d(name)


@pytest.mark.parametrize('name', ['annulus_4', 'pretzel'])
def test_build_multipatch_domain_2d_rejects_invalid_radii(name):
    with pytest.raises(ValueError, match='0 < r_min < r_max'):
        build_multipatch_domain_2d(name, r_min=2, r_max=1)


def test_build_multipatch_domain_2d_dispatches_names_and_radii():
    square = build_multipatch_domain_2d('square_2')
    annulus = build_multipatch_domain_2d(
        'annulus_3', r_min=1, r_max=3)

    assert str(square.name) == 'square_2'
    assert str(annulus.name) == 'annulus_3'
    assert annulus.patches[0].logical_domain.min_coords[0] == 1
    assert annulus.patches[0].logical_domain.max_coords[0] == 3


def test_build_two_cube_3d():
    domain = build_two_cube_3d()

    assert str(domain.name) == 'mapped_two_cube_3d'
    assert len(domain.patches) == 2
    assert len(domain.interface_map) == 1
    assert len(domain.exterior_sides) == 10

    interface, = domain.interface_map.values()
    assert interface.orientation == (-1, +1, -1)
    assert interface.axis_map == ((1, 2, +1), (2, 0, -1))


@pytest.mark.parametrize(
    'arguments, name, ninterfaces, nboundaries',
    [
        ({}, 'hollow_torus_3d', 8, 8),
        (
            {'toroidal_angle': 1.5 * np.pi},
            'open_hollow_torus_3d', 6, 12,
        ),
        (
            {'hollow': False, 'close_torus': False},
            'open_solid_torus_3d', 6, 12,
        ),
    ],
)
def test_build_torus_3d(
        arguments, name, ninterfaces, nboundaries):
    domain = build_torus_3d(**arguments)

    assert str(domain.name) == name
    assert len(domain.patches) == 4
    assert len(domain.interface_map) == ninterfaces
    assert len(domain.exterior_sides) == nboundaries
    assert all(interface.orientation == (+1, +1, +1)
               for interface in domain.interface_map.values())


@pytest.mark.parametrize('arguments, error, message', [
    ({'minor_bounds': (0.8, 0.35)}, ValueError, 'minor radius'),
    ({'hollow': 1}, TypeError, 'hollow must be a bool'),
    ({'toroidal_angle': 0}, ValueError, '0 < angle <= 2\\*pi'),
    (
        {'toroidal_angle': np.pi, 'close_torus': True},
        ValueError,
        'only a full 2\\*pi sweep can be closed',
    ),
])
def test_build_torus_3d_rejects_invalid_parameters(
        arguments, error, message):
    with pytest.raises(error, match=message):
        build_torus_3d(**arguments)


@pytest.mark.parametrize(
    'builder, name, npatches, ninterfaces, nboundaries',
    [
        (
            build_two_patch_twisted_strip,
            'two_patch_twisted_strip', 2, 2, 4,
        ),
        (build_self_connected_patch, 'self_connected_patch', 1, 1, 2),
        (
            build_three_patch_shared_vertex,
            'three_patch_shared_vertex', 3, 3, 6,
        ),
        (build_oriented_two_cube, 'oriented_two_cube', 2, 1, 10),
        (
            build_three_patch_shared_edge,
            'three_patch_shared_edge', 3, 3, 12,
        ),
    ],
)
def test_diagnostic_domain_builder(
        builder, name, npatches, ninterfaces, nboundaries):
    domain = builder()

    assert str(domain.name) == name
    assert len(domain.patches) == npatches
    assert len(domain.interface_map) == ninterfaces
    assert len(domain.exterior_sides) == nboundaries


@pytest.mark.parametrize('name', ['unknown', None])
def test_build_multipatch_domain_3d_rejects_unknown_name(name):
    with pytest.raises(ValueError, match='unknown 3D multipatch domain'):
        build_multipatch_domain_3d(name)


def test_build_multipatch_domain_3d_dispatches_names_and_arguments():
    two_cube = build_multipatch_domain_3d('two_cube')
    torus = build_multipatch_domain_3d(
        'torus', hollow=False, close_torus=False)

    assert str(two_cube.name) == 'mapped_two_cube_3d'
    assert str(torus.name) == 'open_solid_torus_3d'


@pytest.mark.parametrize('domain_name, dimension', [
    ('square_2', 2),
    ('two_cube', 3),
])
def test_plot_multipatch_domain(domain_name, dimension, tmp_path):
    output = tmp_path / f'{domain_name}.png'
    domain, figure = plot_multipatch_domain(
        domain_name, output=output, show=False, topology=True)

    assert int(domain.dim) == dimension
    assert figure.axes
    assert output.is_file()


def test_multipatch_gallery_main(tmp_path):
    output = tmp_path / 'square_2.png'

    main([
        'square_2',
        '--topology',
        '--no-show',
        '--output',
        str(output),
    ])

    assert output.is_file()
