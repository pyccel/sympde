"""Apply exterior-calculus operators to symbolic differential forms."""

from sympy import Symbol

from sympde.api import Constant, DifferentialForm, Domain, VectorFunctionSpace
from sympde.api import d, delta, element_of, hodge, ip, jp, latex, ld, wedge


dimension = Symbol('n')
alpha = Constant('alpha')

u_0, v_0 = [
    DifferentialForm(name, index=0, dim=dimension)
    for name in ('u_0', 'v_0')
]
u_1, v_1 = [
    DifferentialForm(name, index=1, dim=dimension)
    for name in ('u_1', 'v_1')
]
u_2 = DifferentialForm('u_2', index=2, dim=dimension)
u_n = DifferentialForm('u_n', index=dimension, dim=dimension)

# The exterior derivative is linear, nilpotent, and zero on top forms.
assert d(d(u_0)) == 0
assert d(u_0 + v_0) == d(u_0) + d(v_0)
assert d(alpha * u_0 + v_0) == alpha * d(u_0) + d(v_0)
assert d(u_n) == 0

# Its adjoint is linear and zero on zero-forms.
assert delta(u_0) == 0
assert delta(u_1 + v_1) == delta(u_1) + delta(v_1)
assert delta(alpha * u_1 + v_1) == alpha * delta(u_1) + delta(v_1)

domain = Domain('Omega', dim=3)
vector_space = VectorFunctionSpace('W', domain)
beta = element_of(vector_space, name='beta')

expressions = (
    d(u_1) + u_2,
    wedge(d(u_1), u_2),
    hodge(u_1),
    ip(beta, u_1),
    jp(beta, u_1),
    ld(beta, u_1),
)

for expression in expressions:
    print(latex(expression))
