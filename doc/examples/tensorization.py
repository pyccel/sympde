"""Expose the tensor-product structure of a scalar bilinear form."""

from sympde.api import BilinearForm, Constant, Domain, ScalarFunctionSpace
from sympde.api import TensorExpr, dot, elements_of, grad, integral


domain = Domain('Omega', dim=2)
space = ScalarFunctionSpace('V', domain)
trial, test = elements_of(space, names='u, v')
coefficient = Constant('mu', real=True)

form = BilinearForm(
    (trial, test),
    integral(
        domain,
        coefficient * trial * test + dot(grad(trial), grad(test)),
    ),
)

tensor_expression = TensorExpr(form, domain=domain)
print(tensor_expression)
