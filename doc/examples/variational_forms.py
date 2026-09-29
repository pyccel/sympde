"""Construct scalar variational forms on single- and multi-patch domains."""

from sympy import cos

from sympde.api import BilinearForm, Constant, Domain, LinearForm, ScalarFunctionSpace
from sympde.api import Square, Dn, dot, elements_of, find, grad, integral
from sympde.api import jump, latex, minus, plus


# A Poisson problem on a single square.
domain = Square('Omega')
space = ScalarFunctionSpace('V', domain)
trial, test = elements_of(space, names='u, v')

bilinear_form = BilinearForm(
    (trial, test),
    integral(domain, dot(grad(trial), grad(test))),
)

x, y = domain.coordinates
linear_form = LinearForm(test, integral(domain, cos(x + y) * test))
equation = find(
    trial,
    forall=test,
    lhs=bilinear_form(trial, test),
    rhs=linear_form(test),
)

print(latex(equation))


# A Nitsche-type interface contribution on two joined patches.
patch_a = Square('A')
patch_b = Square('B')
domain = Domain.join(
    [patch_a, patch_b],
    [((0, 0, 1), (1, 0, -1), 1)],
    'Omega',
)

space = ScalarFunctionSpace('V', domain)
trial, test = elements_of(space, names='u, v')
interface = domain.interfaces
kappa = Constant('kappa')

interface_integrand = (
    -jump(trial) * jump(Dn(test))
    + kappa * jump(trial) * jump(test)
    + plus(Dn(trial)) * minus(test)
    + minus(Dn(trial)) * plus(test)
)

multi_patch_form = BilinearForm(
    (trial, test),
    integral(domain, dot(grad(trial), grad(test)))
    + integral(interface, interface_integrand),
)

print(latex(multi_patch_form))
