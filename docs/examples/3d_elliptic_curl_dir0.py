from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import VectorFunctionSpace, Cube, element_of
from sympde.calculus import curl, dot
from sympde.core     import Constant

from sympy import pi, sin, Tuple

mu = Constant('mu', is_real=True)

domain = Cube()

V = VectorFunctionSpace('V', domain, kind='Hcurl')

x,y,z = domain.coordinates

u,v = [element_of(V, name=i) for i in ['u', 'v']]

# bilinear form
a = BilinearForm((u,v), integral(domain, dot(curl(v), curl(u)) + mu * dot(u,v)))

# linear form
f1 = sin(pi*x)*sin(pi*y)*sin(pi*z)
f2 = sin(pi*x)*sin(pi*y)*sin(pi*z)
f3 = sin(pi*x)*sin(pi*y)*sin(pi*z)
f = Tuple(f1, f2, f3)

l = LinearForm(v, integral(domain, dot(f,v)))

# Dirichlet boundary conditions
bc = [EssentialBC(u, 0, domain.boundary)]

# Variational problem
equation   = find(u, forall=v, lhs=a(u, v), rhs=l(v), bc=bc)
