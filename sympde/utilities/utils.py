import numpy as np
import itertools as it
from sympy import lambdify

from mpl_toolkits.mplot3d import *
import matplotlib.pyplot as plt

from sympde.topology import IdentityMapping, InteriorDomain, MultiPatchMapping
from sympde.topology.analytical_mapping import TorusMapping

def lambdify_sympde(variables, expr):
    """
    Custom lambify function that covers the
    shortcomings of sympy's lambdify. Most notably,
    this function uses numpy broadcasting rules to
    compute the shape of the output.

    Parameters
    ----------
    variables : sympy.core.symbol.Symbol or list of sympy.core.symbol.Symbol
        variables that appear in the expression
    expr :
        Sympy expression

    Returns
    -------
    lambda_f : callable
        Lambdified function built using numpy.

    Notes
    -----
    Compared to Sympy's lambdify, this function
    is capable of properly handling constant values,
    and array_like structures where not all components
    depend on all variables. See below.

    Examples
    --------
    >>> import numpy as np
    >>> from sympy import symbols,  Matrix
    >>> from sympde.utilities.utils import lambdify_sympde
    >>> x, y = symbols("x,y")
    >>> expr = Matrix([[x, x + y], [0, y]])
    >>> f = lambdify_sympde([x,y], expr)
    >>> f(np.array([[0, 1]]), np.array([[2], [3]]))
    array([[[[0., 1.],
             [0., 1.]],

            [[2., 3.],
             [3., 4.]]],


           [[[0., 0.],
             [0., 0.]],

            [[2., 2.],
             [3., 3.]]]])
    """
    array_expr = np.asarray(expr)
    scalar_shape = array_expr.shape
    if scalar_shape == ():
        f = lambdify(variables, expr, 'numpy')
        def f_vec_sc(*XYZ):
            b = np.broadcast(*XYZ)
            if b.ndim == 0:
                return f(*XYZ)
            temp = np.asarray(f(*XYZ))
            if b.shape == temp.shape:
                return temp

            result = np.zeros(b.shape)
            result[...] = temp
            return result
        return f_vec_sc

    else:
        scalar_functions = {}
        for multi_index in it.product(*tuple(range(s) for s in scalar_shape)):
            scalar_functions[multi_index] = lambdify(variables, array_expr[multi_index], 'numpy')

        def f_vec_v(*XYZ):
            b = np.broadcast(*XYZ)
            result = np.zeros(scalar_shape + b.shape)
            for multi_index in it.product(*tuple(range(s) for s in scalar_shape)):
                result[multi_index] = scalar_functions[multi_index](*XYZ)
            return result
        return f_vec_v


def plot_domain(domain, draw=True):
    """
    Plots a 2D  or 3D domain using matplotlib

    Parameters
    ----------
    domain : sympde.topology.Domain
        Domain to plot

    draw : bool
        if true, plt.show() will be called.
    """
    pdim = domain.dim if domain.mapping is None else domain.mapping.pdim
    if pdim == 2:
        plot_2d(domain, draw=draw)
    elif pdim ==3:
        plot_3d(domain, draw=draw)


def plot_2d(domain, draw=True):
    """
    Plot a 2D domain

    Parameters
    ----------
    domain : sympde.topology.Domain
        Domain to plot

    draw : bool
        if true, plt.show() will be called.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if isinstance(domain.interior, InteriorDomain):
        plot_2d_single_patch(domain.interior, domain.mapping, ax)
    else:
        if isinstance(domain.mapping, MultiPatchMapping):
            for patch, mapping in domain.mapping.mappings.items():
                plot_2d_single_patch(patch, mapping, ax)
        else:
            for interior in domain.interior.as_tuple():
                plot_2d_single_patch(interior, interior.mapping, ax)

    ax.set_aspect('equal', adjustable='box')
    if draw:
        plt.show()

def plot_3d(domain, draw=True):
    """
    Plot a 3D domain

    Parameters
    ----------
    domain : sympde.topology.Domain
        Domain to plot

    draw : bool
        if true, plt.show() will be called.
    """
    mapping = domain.mapping

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    if isinstance(domain.interior, InteriorDomain):
        plot_3d_single_patch(domain.interior, domain.mapping, ax)
    else:
        if isinstance(domain.mapping, MultiPatchMapping):
            for patch, mapping in domain.mapping.mappings.items():
                plot_3d_single_patch(patch, mapping, ax)
        else:
            for interior in domain.interior.as_tuple():
                plot_3d_single_patch(interior, interior.mapping, ax)

    if draw:
        plt.show()

def plot_3d_single_patch(patch, mapping, ax):
    """
    Plot a singe patch in a 3D domain

    Parameters
    ----------
    patch : sympde.topology.InteriorDomain

    mapping : sympde.topology.mapping

    ax : mpl_toolkits.mplot3d.axes3d.Axes3D
        Axes object on which the patch is drawn.
    """
    if mapping is None:
        mapping = IdentityMapping('Id', dim=3)

    map_call = mapping.get_callable_mapping()
    refinement = 21
    mesh_grid = np.meshgrid(
        *[np.linspace(patch.min_coords[i],
                      patch.max_coords[i],
                      num=refinement,
                      endpoint=True) for i in range(3)],
        indexing='ij',
        sparse=True
    )

    XX, YY, ZZ = map_call(*mesh_grid)

    for i in range(0, XX.shape[-1], 2):
        ax.plot_wireframe(XX[:, :, i], YY[:, :, i], ZZ[:, :, i], color='k', cstride=20, rstride=20)
    for j in range(0, XX.shape[-2], 2):
        ax.plot_wireframe(XX[:, j, :], YY[:, j, :], ZZ[:, j, :], color='k', cstride=20, rstride=20)
    for k in range(0, XX.shape[-3], 2):
        ax.plot_wireframe(XX[k, :, :], YY[k, :, :], ZZ[k, :, :], color='k',  cstride=20, rstride=20)


def plot_2d_single_patch(patch, mapping, ax):
    """
    Plot a singe patch in a 2D domain

    Parameters
    ----------
    patch : sympde.topology.InteriorDomain

    mapping : sympde.topology.mapping

    ax : matplotlib.axes.Axes
        Axes object on which the patch is drawn.
    """
    if mapping is None:
        mapping = IdentityMapping('Id', dim=3)

    map_call = mapping.get_callable_mapping()
    refinement = 41
    linspace_0 = np.linspace(patch.min_coords[0], patch.max_coords[0], refinement, endpoint=True)
    linspace_1 = np.linspace(patch.min_coords[1], patch.max_coords[1], refinement, endpoint=True)

    X_00, Y_00 = map_call(linspace_0, np.full(refinement, linspace_1[0]))
    X_01, Y_01 = map_call(linspace_0, np.full(refinement, linspace_1[-1]))
    X_10, Y_10 = map_call(np.full(refinement, linspace_0[0]), linspace_1)
    X_11, Y_11 = map_call(np.full(refinement, linspace_0[-1]), linspace_1)


    ax.plot(X_00, Y_00, 'k')
    ax.plot(X_01, Y_01, 'k')
    ax.plot(X_10, Y_10, 'k')
    ax.plot(X_11, Y_11, 'k')

if __name__ == '__main__':

    from sympde.topology import Square, PolarMapping, Cube
    A = Square('A', bounds1=(0, 1), bounds2=(0, np.pi/2))
    F = PolarMapping('F', c1=0, c2=0, rmin=0.5, rmax=1)
    Omega = F(A)

    plot_domain(Omega)
