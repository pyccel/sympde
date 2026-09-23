from io import StringIO
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import pytest

from sympy import Matrix, symbols, Array
from sympy import S
import sympde.utilities as utilities
import sympde.utilities.utils as utilities_utils
from sympde.topology import (
    AffineMapping,
    Cube,
    Domain,
    IdentityMapping,
    Square,
)
from sympde.utilities import (
    collect_topology_vertices,
    format_topology,
    plot_domain,
    print_topology,
)
from sympde.utilities.plotting import plot_domain as plotting_plot_domain
from sympde.utilities.utils import lambdify_sympde


def test_import_sympde_does_not_import_matplotlib():
    code = """
import sys
import sympde
assert 'matplotlib' not in sys.modules
import sympde.utilities
assert 'matplotlib' not in sys.modules
from sympde.utilities import plot_domain
assert 'matplotlib' in sys.modules
"""
    subprocess.run([sys.executable, '-c', code], check=True)


def test_lambdify_sympde_1d():
    x = symbols("x")

    independent_expr = S.One
    dependent_expr = x

    array_i_expr = Matrix([[1, 0], [0, 1]])
    array_d_expr = Matrix([[x, 0], [x, x]])

    f_i = lambdify_sympde(x, independent_expr)
    f_d = lambdify_sympde(x, dependent_expr)

    f_ia = lambdify_sympde(x, array_i_expr)
    f_da = lambdify_sympde(x, array_d_expr)

    scalar_input_i = f_i(0)
    scalar_input_d = f_d(0)

    scalar_input_ia = f_ia(0)
    scalar_input_da = f_da(0)

    assert np.array_equal(scalar_input_i, 1)
    assert np.array_equal(scalar_input_d, 0)

    assert np.array_equal(scalar_input_ia, np.eye(2))
    assert np.array_equal(scalar_input_da, np.zeros((2, 2)))

    array_input_i = f_i(np.linspace(0, 1, 10))
    array_input_d = f_d(np.linspace(0, 1, 10))

    array_input_ia = f_ia(np.linspace(0, 1, 10))
    array_input_da = f_da(np.linspace(0, 1, 10))

    assert np.array_equal(array_input_i, np.ones(10))
    assert np.array_equal(array_input_d, np.linspace(0, 1, 10))

    expected_a_ia = np.zeros((2, 2, 10)) + np.eye(2)[:, :, None]
    assert np.array_equal(array_input_ia, expected_a_ia)

    expected_a_da = np.zeros((2, 2, 10))
    expected_a_da[0, 0, ...] = np.linspace(0, 1, 10)[None, None, Ellipsis]
    expected_a_da[1, ...] = np.linspace(0, 1, 10)[None, None, Ellipsis]
    assert np.array_equal(array_input_da, expected_a_da)


def test_lambdify_sympde_2d():
    x,y = symbols("x, y")

    independent_expr = S.One * 2
    dependent_expr = x + y
    semi_dependent_expr = x + 3

    array_i_expr = Matrix([[0, 3, 5], [5, 3, 4]])
    array_d_expr = Matrix([[x + 3, 3], [y, x]])
    array_s_expr = Array([[x, x], [x, 5]])

    f_i = lambdify_sympde([x, y], independent_expr)
    f_d = lambdify_sympde([x, y], dependent_expr)
    f_s = lambdify_sympde([x, y], semi_dependent_expr)

    f_ia = lambdify_sympde([x, y], array_i_expr)
    f_da = lambdify_sympde([x, y], array_d_expr)
    f_sa = lambdify_sympde([x, y], array_s_expr)

    scalar_input_i = f_i(0, 0)
    scalar_input_d = f_d(0, 0)
    scalar_input_s = f_s(0, 0)

    scalar_input_ia = f_ia(0, 0)
    scalar_input_da = f_da(0, 0)
    scalar_input_sa = f_sa(0, 0)

    assert np.array_equal(scalar_input_i, 2)
    assert np.array_equal(scalar_input_d, 0)
    assert np.array_equal(scalar_input_s, 3)

    assert np.array_equal(scalar_input_ia, np.array([[0, 3, 5], [5, 3, 4]]))
    assert np.array_equal(scalar_input_da, np.array([[3, 3], [0, 0]]))
    assert np.array_equal(scalar_input_sa, np.array([[0, 0], [0, 5]]))

    dense_input = np.meshgrid(np.linspace(0, 1, 10),  np.linspace(0, 1, 5), sparse=False)
    sparse_input = np.meshgrid(np.linspace(0, 1, 10),  np.linspace(0, 1, 5), sparse=True)

    dense_input_i = f_i(*dense_input)
    dense_input_d = f_d(*dense_input)
    dense_input_s = f_s(*dense_input)

    dense_input_ia = f_ia(*dense_input)
    dense_input_da = f_da(*dense_input)
    dense_input_sa = f_sa(*dense_input)

    assert np.array_equal(dense_input_i, np.full((5, 10), 2))
    assert np.array_equal(dense_input_d, dense_input[0] + dense_input[1])
    assert np.array_equal(dense_input_s, np.zeros((5, 10)) + (dense_input[0] + 3))

    expected_ia = np.zeros((2, 3, 5, 10)) + np.array([[0, 3, 5], [5, 3, 4]])[..., None, None]
    expected_da = np.zeros((2, 2, 5, 10))
    expected_da[0, 0, ...] = dense_input[0] + 3
    expected_da[0, 1, ...] = 3
    expected_da[1, 0, ...] = dense_input[1]
    expected_da[1, 1, ...] = dense_input[0]

    expected_sa = np.zeros((2, 2, 5, 10))
    expected_sa[:, :, ...] = dense_input[0]
    expected_sa[1, 1, ...] = 5

    assert np.array_equal(dense_input_ia, expected_ia)
    assert np.array_equal(dense_input_da, expected_da)
    assert np.array_equal(dense_input_sa, expected_sa)

    sparse_input_i = f_i(*sparse_input)
    sparse_input_d = f_d(*sparse_input)
    sparse_input_s = f_s(*sparse_input)

    sparse_input_ia = f_ia(*sparse_input)
    sparse_input_da = f_da(*sparse_input)
    sparse_input_sa = f_sa(*sparse_input)

    assert np.array_equal(sparse_input_i, np.full((5, 10), 2))
    assert np.array_equal(sparse_input_d, dense_input[0] + dense_input[1])
    assert np.array_equal(sparse_input_s, np.zeros((5, 10)) + (dense_input[0] + 3))

    assert np.array_equal(sparse_input_ia, expected_ia)
    assert np.array_equal(sparse_input_da, expected_da)
    assert np.array_equal(sparse_input_sa, expected_sa)


def _make_twisted_strip():
    patch_a = IdentityMapping('F_A', dim=2)(Square('A'))
    patch_b = AffineMapping(
        'F_B', dim=2,
        c1=1, c2=0,
        a11=1, a12=0,
        a21=0, a22=1,
    )(Square('B'))
    return Domain.join(
        [patch_a, patch_b],
        [
            ((0, 0, 1), (1, 0, -1), +1),
            ((0, 0, -1), (1, 0, 1), -1),
        ],
        'twisted_strip',
    )


def _make_three_patch_vertex():
    patches = [Square(name) for name in ('A', 'B', 'C')]
    return Domain.join(
        patches,
        [
            ((0, 0, -1), (1, 1, -1), +1),
            ((1, 0, -1), (2, 1, -1), +1),
            ((2, 0, -1), (0, 1, -1), +1),
        ],
        'three_patch_vertex',
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


def test_collect_topology_vertices_shared_and_self_interfaces():
    vertices = collect_topology_vertices(_make_three_patch_vertex())
    interior_vertices = [vertex for vertex in vertices if not vertex.is_boundary]

    assert len(interior_vertices) == 1
    center = interior_vertices[0]
    assert center.patch_indices == (0, 1, 2)
    assert center.patches == ('A', 'B', 'C')
    assert all(
        incidence.logical_corner == (0, 0)
        for incidence in center.incidences)

    patch = Square('P')
    self_connected = Domain.join(
        [patch],
        [((0, 0, -1), (0, 0, 1), +1)],
        'self_connected',
    )
    self_vertices = collect_topology_vertices(self_connected)
    assert len(self_vertices) == 2
    assert all(len(vertex.incidences) == 2 for vertex in self_vertices)
    assert all(vertex.patches == ('P',) for vertex in self_vertices)
    assert all(vertex.is_boundary for vertex in self_vertices)


def test_collect_and_format_3d_topology():
    domain = _make_oriented_cubes()
    vertices = collect_topology_vertices(domain)
    shared_vertices = [
        vertex for vertex in vertices if len(vertex.incidences) == 2]

    assert len(vertices) == 12
    assert len(shared_vertices) == 4
    assert all(vertex.patches == ('A3', 'B3') for vertex in shared_vertices)

    text = format_topology(domain)
    assert 'patches: A3, B3' in text
    assert 'orientation=(-1, 1, -1)' in text
    assert 'axis_map=((1, 2, 1), (2, 0, -1))' in text


def test_format_and_print_topology():
    domain = _make_twisted_strip()
    text = format_topology(domain)

    assert 'patches: A, B' in text
    assert 'I0' in text and 'I1' in text
    assert 'orientation=1' in text
    assert 'orientation=-1' in text
    assert 'axis_map=((1, 1, -1),)' in text
    assert 'vertices:' in text

    stream = StringIO()
    print_topology(domain, file=stream)
    assert stream.getvalue() == text + '\n'


def test_plot_domain_2d_integration():
    assert plot_domain is plotting_plot_domain
    assert not hasattr(utilities, 'plot_topology')
    assert not hasattr(utilities_utils, 'plot_domain')

    domain = _make_twisted_strip()
    figure = plot_domain(
        domain, draw=False, refinement=8, isolines=True, topology=True)
    labels = {text.get_text() for text in figure.axes[0].texts}
    assert {'A', 'B', 'I0-', 'I0+', 'I1-', 'I1+'}.issubset(labels)
    assert any(label.startswith('V') for label in labels)
    assert any('\n    axis map:' in label for label in labels)
    assert figure.axes[0].get_title() == 'twisted strip'
    plt.close(figure)

    figure = plot_domain(Square('single'), draw=False, refinement=8)
    assert len(figure.axes[0].lines) == 4
    plt.close(figure)


def test_plot_domain_3d_integration():
    domain = _make_oriented_cubes()
    figure = plot_domain(
        domain, draw=False, refinement=4, isolines=True, topology=True)
    axis = figure.axes[0]
    labels = {text.get_text() for text in axis.texts}
    assert {'A3', 'B3', 'I0 -/+'}.issubset(labels)
    assert any(label.startswith('V') for label in labels)
    assert any('\n    axis map:' in label for label in labels)
    assert axis.get_zlabel() == 'physical z'
    assert figure.axes[0].get_title() == 'oriented cubes'
    plt.close(figure)

    figure = plot_domain(Cube('single3d'), draw=False, refinement=4)
    axis = figure.axes[0]
    assert axis.get_zlabel() == 'Z'
    assert len(axis.collections) == 6
    plt.close(figure)
