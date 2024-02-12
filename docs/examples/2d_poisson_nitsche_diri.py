from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import ScalarFunctionSpace, Square, element_of
from sympde.calculus import grad, dot, Dn
from sympde.core import Constant

from sympy import pi, sin

kappa = Constant('kappa', is_real=True)

domain = Square()

V = ScalarFunctionSpace('V', domain)

x,y = domain.coordinates

u,v = [element_of(V, name=i) for i in ['u', 'v']]

# bilinear form
a = BilinearForm((u,v), integral(domain, dot(grad(v), grad(u)) + kappa*u*v - u*Dn(v) - v*Dn(u)))

# linear form
g = sin(pi*x)*sin(pi*y)
f = 2*pi**2*sin(pi*x)*sin(pi*y)

l = LinearForm(v, integral(domain, f*v + kappa*g*v - g*Dn(v)))

# Variational problem
equation   = find(u, forall=v, lhs=a(u, v), rhs=l(v))
