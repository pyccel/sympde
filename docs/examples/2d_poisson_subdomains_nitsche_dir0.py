from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import ScalarFunctionSpace, Square, Domain, element_of
from sympde.calculus import grad, dot
from sympde.calculus import jump, avg, minus, plus, Dn
from sympde.core import Constant

# ... create a domain as the union of two subdomains
A = Square('A',bounds1=(0, 0.5), bounds2=(0, 1))
B = Square('B',bounds1=(0.5, 1.), bounds2=(0, 1))

connectivity = [((0,0,1),(1,0,-1))]
subdomains = [A,B]
domain = Domain.join(subdomains, connectivity, 'domain')
# ...

# internal interafaces of the domain
I = domain.interfaces

kappa = Constant('kappa', is_real=True)

V = ScalarFunctionSpace('V', domain)

x,y = domain.coordinates

u,v = [element_of(V, name=i) for i in ['u', 'v']]

# bilinear form
a = BilinearForm((u,v),
                   integral(domain, dot(grad(u),grad(v)))
                 + integral(I, kappa * jump(u)*jump(v) - Dn(u)*jump(v) - jump(u)*Dn(v))
                )

# linear form
f = 4
l = LinearForm(v, integral(domain, f*v))

# Dirichlet boundary conditions
bc = [EssentialBC(u, 0, domain.boundary)]

# Variational problem
equation   = find(u, forall=v, lhs=a(u, v), rhs=l(v), bc=bc)
