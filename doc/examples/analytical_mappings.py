"""Work with built-in and user-defined analytical mappings."""

import numpy as np
from sympy import pi, sin, symbols

from sympde.api import CollelaMapping2D, Constant, Domain, LogicalExpr, Mapping


# Inspect a built-in mapping symbolically.
mapping = CollelaMapping2D('M', dim=2)
mapped_domain = mapping(Domain('Omega', dim=2))
x1, x2 = symbols('x1, x2', real=True)
eps, k1, k2 = [Constant(name) for name in ('eps', 'k1', 'k2')]

expected_x = 2.0 * (
    x1 + eps * sin(2.0 * pi * k1 * x1) * sin(2.0 * pi * k2 * x2)
) - 1.0
assert LogicalExpr(mapping[0], mapped_domain) == expected_x


# Fix its parameters and evaluate it numerically.
numerical_mapping = CollelaMapping2D(
    'M_numeric',
    dim=2,
    eps=0.1,
    k1=1.0,
    k2=1.0,
).get_callable_mapping()
np.testing.assert_allclose(numerical_mapping(0.25, 0.5), (-0.5, 0.0))


# Define a custom mapping by specifying its coordinate expressions.
class ComplexSquareMapping(Mapping):
    """Map ``(x1, x2)`` through the square of a complex number."""

    _expressions = {
        'x': 'A * (x1**2 - x2**2)',
        'y': '2 * A * x1 * x2',
    }


custom_mapping = ComplexSquareMapping('F', dim=2, A=1.0).get_callable_mapping()
np.testing.assert_allclose(custom_mapping(0.5, 0.25), (0.1875, 0.25))
