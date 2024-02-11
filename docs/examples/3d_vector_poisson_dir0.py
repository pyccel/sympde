from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import VectorFunctionSpace, Cube, element_of
from sympde.calculus import grad, inner, dot

from sympy import pi, sin, Tuple

domain = Cube()

V = VectorFunctionSpace('V', domain)

x,y,z = domain.coordinates

u,v = [element_of(V, name=i) for i in ['u', 'v']]

# bilinear form
a = BilinearForm((u,v), integral(domain , inner(grad(v), grad(u))))

# linear form
f1 = 3*pi**2*sin(pi*x)*sin(pi*y)*sin(pi*z)
f2 = 3*pi**2*sin(pi*x)*sin(pi*y)*sin(pi*z)
f3 = 3*pi**2*sin(pi*x)*sin(pi*y)*sin(pi*z)
f = Tuple(f1, f2, f3)

l = LinearForm(v, integral(domain, dot(f,v)))

# Dirichlet boundary conditions
bc = [EssentialBC(u, 0, domain.boundary)]

# Variational problem
equation   = find(u, forall=v, lhs=a(u, v), rhs=l(v), bc=bc)
