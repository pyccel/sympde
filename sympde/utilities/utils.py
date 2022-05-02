import numpy as np
import itertools as it
from sympy import lambdify

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
