import pytest
import numpy as np

from sympy import Matrix, symbols, Array
from sympy import S
from sympde.utilities.utils import lambdify_sympde




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
