from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import ScalarFunctionSpace, Square, element_of
from sympde.calculus import laplace, dot

domain = Square()

V = ScalarFunctionSpace('V', domain)

x,y = domain.coordinates

u,v = [element_of(V, name=i) for i in ['u', 'v']]

# bilinear form
a = BilinearForm((u,v), integral(domain , laplace(v) * laplace(u)))

# linear form
f = 4
l = LinearForm(v, integral(domain, f*v))

# Dirichlet boundary conditions
bc = [EssentialBC(u, 0, domain.boundary)]

# Variational problem
equation   = find(u, forall=v, lhs=a(u, v), rhs=l(v), bc=bc)
