from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import ScalarFunctionSpace, Square, element_of
from sympde.calculus import grad, dot

domain = Square()

V = ScalarFunctionSpace('V', domain)

x,y = domain.coordinates

u,v = [element_of(V, name=i) for i in ['u', 'v']]

# bilinear form
a = BilinearForm((u,v), integral(domain , dot(grad(v), grad(u))))

# linear form
l = LinearForm(v, integral(domain, 4*v))

# boundary condition
B_dirichlet_0 = domain.boundary

# Dirichlet boundary conditions
bc = [EssentialBC(u, 0, B_dirichlet_0)]

# Variational problem
equation   = find(u, forall=v, lhs=a(u, v), rhs=l(v), bc=bc)
