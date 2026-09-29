"""Linearize the residual of a nonlinear Poisson problem."""

from sympy import exp

from sympde.api import BilinearForm, Domain, LinearForm, ScalarFunctionSpace
from sympde.api import dot, elements_of, grad, integral, linearize


domain = Domain('Omega', dim=2)
space = ScalarFunctionSpace('V', domain)
field, increment, test = elements_of(space, names='u, delta_u, v')


def integrate(expression):
    """Integrate an expression over the problem domain."""

    return integral(domain, expression)


residual = LinearForm(
    test,
    integrate(dot(grad(test), grad(field)) - 4.0 * exp(-field) * test),
)

jacobian = linearize(residual, field, trials=increment)
assert isinstance(jacobian, BilinearForm)
assert jacobian(increment, test) == integrate(
    dot(grad(test), grad(increment))
    + 4.0 * exp(-field) * increment * test
)
