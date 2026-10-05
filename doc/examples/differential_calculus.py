"""Verify representative symbolic vector-calculus identities."""

from sympy import expand

from sympde.api import Constant, Domain, ScalarFunctionSpace, VectorFunctionSpace
from sympde.api import curl, div, elements_of, grad, laplace


domain = Domain('Omega', dim=3)
scalar_space = ScalarFunctionSpace('V', domain)
vector_space = VectorFunctionSpace('W', domain)

alpha, beta = [Constant(name) for name in ('alpha', 'beta')]
f, g, h = elements_of(scalar_space, names='f, g, h')
field_f, field_g, field_h = elements_of(vector_space, names='F, G, H')

# Gradient linearity and product rules.
assert grad(f + g) == grad(f) + grad(g)
assert grad(alpha * f + beta * g) == alpha * grad(f) + beta * grad(g)
assert grad(f * g) == f * grad(g) + g * grad(f)
assert grad(f / g) == -f * grad(g) / g**2 + grad(f) / g
assert expand(grad(f * g * h)) == (
    f * g * grad(h) + f * h * grad(g) + g * h * grad(f)
)

# Linearity of curl, Laplace, and divergence operators.
assert curl(field_f + field_g) == curl(field_f) + curl(field_g)
assert curl(alpha * field_f + beta * field_g) == (
    alpha * curl(field_f) + beta * curl(field_g)
)
assert laplace(f + g) == laplace(f) + laplace(g)
assert laplace(alpha * f + beta * g) == alpha * laplace(f) + beta * laplace(g)
assert div(field_f + field_g) == div(field_f) + div(field_g)
assert div(alpha * field_f + beta * field_g) == (
    alpha * div(field_f) + beta * div(field_g)
)

# Exact-sequence identities.
assert curl(grad(h)) == 0
assert div(curl(field_h)) == 0
