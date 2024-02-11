from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import ScalarFunctionSpace, Square, element_of
from sympde.topology import NormalVector, Union
from sympde.calculus import grad, dot

from sympy import pi, cos, exp

domain = Square()
x,y = domain.coordinates

V = ScalarFunctionSpace('V', domain)

B_neumann     = domain.get_boundary(axis=0, ext=-1)
B_dirichlet_0 = domain.get_boundary(axis=0, ext=1)
B_dirichlet_i = domain.boundary.complement(Union(B_neumann, B_dirichlet_0))

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

# Homogeneous Dirichlet boundary conditions
bc  = [EssentialBC(u, 0, B_dirichlet_0)]

# Inhomogeneous Dirichlet boundary conditions
g2 = cos(pi*x) * exp(-2*y)
bc += [EssentialBC(u, g2, B_dirichlet_i)]

# Variational problem
equation   = find(u, forall=v, lhs=a(u, v), rhs=l(v) + ln(v), bc=bc)
