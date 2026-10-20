"""Write the expression trees of symbolic expressions in the DOT language."""

from sympy.printing.dot import dotprint

from sympde.api import Domain, ScalarFunctionSpace, div, dx, dy, elements_of, grad


domain = Domain('Omega', dim=2)
space = ScalarFunctionSpace('V', domain)
u, v = elements_of(space, names='u, v')


# [operator-start]
def L(w):
    return -div(grad(w)) + w


operator_graph = dotprint(L(u))
# [operator-end]


# [derivative-start]
expr = dx(dy(u * v))
assert expr == u * dx(dy(v)) + v * dx(dy(u)) + dx(u) * dy(v) + dx(v) * dy(u)

derivative_graph = dotprint(expr)
# [derivative-end]


# Graphs drawn in the documentation, by file name
graphs = {
    'operator': operator_graph,
    'derivative': derivative_graph,
}
