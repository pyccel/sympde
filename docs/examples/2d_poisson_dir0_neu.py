from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import ScalarFunctionSpace, Square, element_of
from sympde.topology import NormalVector
from sympde.calculus import grad, dot

domain = Square()
x,y = domain.coordinates

V = ScalarFunctionSpace('V', domain)

B_neumann   = domain.get_boundary(axis=0, ext=-1)
B_dirichlet = domain.boundary.complement(B_neumann)

u,v = [element_of(V, name=i) for i in ['u', 'v']]
nn = NormalVector('nn')

# bilinear form
a = BilinearForm((u,v), integral(domain , dot(grad(v), grad(u))))

# linear form
f = x*y
l = LinearForm(v, integral(domain, f*v))

# Boundary term for the Neumann BC
g = x**2 + y**2
ln = LinearForm(v, integral(B_neumann, v * dot(grad(g), nn)))

# Dirichlet boundary conditions
bc = [EssentialBC(u, 0, B_dirichlet)]

# Variational problem
equation   = find(u, forall=v, lhs=a(u, v), rhs=l(v) + ln(v), bc=bc)
