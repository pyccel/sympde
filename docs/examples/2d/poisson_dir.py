from sympde.calculus import grad, dot
from sympde.topology import ScalarFunctionSpace
from sympde.topology import element_of
from sympde.topology import NormalVector
from sympde.topology import Square
from sympde.topology import Union
from sympde.expr     import BilinearForm, LinearForm, integral
from sympde.expr     import Norm, SemiNorm
from sympde.expr     import find, EssentialBC

domain = Square()

B_dirichlet_0 = Union(*[domain.get_boundary(**kw) for kw in dir_zero_boundary])
B_dirichlet_i = Union(*[domain.get_boundary(**kw) for kw in dir_nonzero_boundary])
B_dirichlet   = Union(B_dirichlet_0, B_dirichlet_i)
B_neumann = domain.boundary.complement(B_dirichlet)

V  = ScalarFunctionSpace('V', domain)
u  = element_of(V, name='u')
v  = element_of(V, name='v')
nn = NormalVector('nn')

# Bilinear form a: V x V --> R
a = BilinearForm((u, v), integral(domain, dot(grad(u), grad(v))))

# Linear form l: V --> R
l0 = LinearForm(v, integral(domain, f * v))
if B_neumann:
    l1 = LinearForm(v, integral(B_neumann, v * dot(grad(solution), nn)))
    l  = LinearForm(v, l0(v) + l1(v))
else:
    l = l0

# Dirichlet boundary conditions
bc = []
if B_dirichlet_0:  bc += [EssentialBC(u,        0, B_dirichlet_0)]
if B_dirichlet_i:  bc += [EssentialBC(u, solution, B_dirichlet_i)]

# Variational model
equation = find(u, forall=v, lhs=a(u, v), rhs=l(v), bc=bc)

# Error norms
error  = u - solution
l2norm =     Norm(error, domain, kind='l2')
h1norm = SemiNorm(error, domain, kind='h1')

