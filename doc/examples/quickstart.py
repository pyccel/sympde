"""Build a scalar Poisson problem with the current public API."""

from sympy import cos

from sympde.api import BilinearForm, LinearForm, ScalarFunctionSpace, Square
from sympde.api import TerminalExpr, div, dot, elements_of, find, grad
from sympde.api import integral, latex


# [domains-start]
domain = Square('Omega')
space = ScalarFunctionSpace('V', domain)
trial, test = elements_of(space, names='u, v')

strong_expression = -div(grad(trial)) + trial
assert strong_expression.has(trial)
# [domains-end]


# [forms-start]
bilinear_form = BilinearForm(
    (trial, test),
    integral(domain, dot(grad(trial), grad(test)) + trial * test),
)

x, y = domain.coordinates
linear_form = LinearForm(test, integral(domain, cos(x + y) * test))

equation = find(
    trial,
    forall=test,
    lhs=bilinear_form(trial, test),
    rhs=linear_form(test),
)
# [forms-end]


# [lowering-start]
coordinate_expression = TerminalExpr(
    dot(grad(trial), grad(test)),
    domain,
)
assert coordinate_expression.has(trial, test)
# [lowering-end]


# [printing-start]
latex_equation = latex(equation)
print(latex_equation)
# [printing-end]
