from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import ScalarFunctionSpace, Square, element_of
from sympde.calculus import grad, dot, laplace
from sympde.core     import Constant

from sympy import Tuple

kappa = Constant('kappa', is_real=True)
tau   = Constant('tau', is_real=True)

b1 = 1.
b2 = 0.
b = Tuple(b1, b2)

# Differential operator in the strong form
L = lambda w: -kappa*laplace(w) + dot(b, grad(w))

domain = Square()

V = ScalarFunctionSpace('V', domain)

x,y = domain.coordinates

u,v = [element_of(V, name=i) for i in ['u', 'v']]

# bilinear form
expr = kappa * dot(grad(v), grad(u)) + dot(b, grad(u)) * v
a = BilinearForm((u,v), integral(domain, expr))

# linear form
f = x*y
l = LinearForm(v, integral(domain, f*v))

# ...
s_supg = BilinearForm((v,u), integral(domain, tau * L(u) * dot(b, grad(v))))

l_supg = LinearForm(v, integral(domain, tau * f * dot(b, grad(v))))
# ...

# ...
s_gls = BilinearForm((v,u), integral(domain, tau * L(u) * L(v)))

l_gls = LinearForm(v, integral(domain, tau * f * L(v)))
# ...

# Dirichlet boundary conditions
bc = [EssentialBC(u, 0, domain.boundary)]

# Variational problem
equation_supg = find(u, forall=v, lhs=a(u, v) + s_supg(u,v), rhs=l(v) + l_supg(v), bc=bc)
equation_gls  = find(u, forall=v, lhs=a(u, v) + s_gls(u,v), rhs=l(v) + l_gls(v), bc=bc)
