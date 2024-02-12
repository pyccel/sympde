from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import (ScalarFunctionSpace, VectorFunctionSpace, Square,
                             elements_of)
from sympde.calculus import grad, dot, div, inner
from sympde.core     import Constant

from sympy import pi, cos, Tuple

domain = Square()
x, y   = domain.coordinates

V1 = VectorFunctionSpace('V1', domain, kind='H1')
V2 = ScalarFunctionSpace('V2', domain, kind='L2')

# rhs
fx = -x**2*(x - 1)**2*(24*y - 12) - 4*y*(x**2 + 4*x*(x - 1) + (x - 1)**2)*(2*y**2 - 3*y + 1) - 2*pi*cos(2*pi*x)
fy = 4*x*(2*x**2 - 3*x + 1)*(y**2 + 4*y*(y - 1) + (y - 1)**2) + y**2*(24*x - 12)*(y - 1)**2 + 2*pi*cos(2*pi*y)
f  = Tuple(fx, fy)

u, v = elements_of(V1, names='u, v')
p, q = elements_of(V2, names='p, q')

# bilinear form
a  = BilinearForm(((u, p), (v, q)), integral(domain, inner(grad(u), grad(v)) - div(u)*q - p*div(v)) )

# linear form
l  = LinearForm((v, q), integral(domain, dot(f, v)))

# Dirichlet boundary conditions
bc = EssentialBC(u, 0, domain.boundary)

equation = find((u, p), forall=(v, q), lhs=a((u, p), (v, q)), rhs=l(v, q), bc=bc)
