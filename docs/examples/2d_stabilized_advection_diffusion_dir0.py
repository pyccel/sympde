from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import ScalarFunctionSpace, Square, element_of
from sympde.calculus import grad, dot
from sympde.core     import Constant

from sympy import pi, sin, Tuple

domain = Square()

V = ScalarFunctionSpace('V', domain)

x,y = domain.coordinates

u,v = [element_of(V, name=i) for i in ['u', 'v']]

kappa = Constant('kappa', is_real=True)
b1 = 1.
b2 = 0.
b = Tuple(b1, b2)

# bilinear form
expr = kappa * dot(grad(v), grad(u)) + dot(b, grad(u)) * v
a = BilinearForm((u,v), integral(domain, expr))

# linear form
f = x*y
l = LinearForm(v, integral(domain, f*v))

# Dirichlet boundary conditions
bc = [EssentialBC(u, 0, domain.boundary)]

# Variational problem
equation   = find(u, forall=v, lhs=a(u, v), rhs=l(v), bc=bc)
