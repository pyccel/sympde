from sympde.expr import BilinearForm, LinearForm, integral
from sympde.expr     import find, EssentialBC
from sympde.topology import (ScalarFunctionSpace, VectorFunctionSpace, Square,
                             element_of)
from sympde.calculus import grad, dot, div
from sympde.core     import Constant

domain = Square()

V1 = VectorFunctionSpace('V1', domain, kind='Hdiv')
V2 = ScalarFunctionSpace('V2', domain, kind='L2')

x,y = domain.coordinates

# rhs
f = -2*x*(1-x) -2*y*(1-y)

u,v = [element_of(V1, name=i) for i in ['u', 'v']]
p,q = [element_of(V2, name=i) for i in ['p', 'q']]

# bilinear form
a  = BilinearForm(((u,p),(v,q)), integral(domain, dot(u,v) - p*div(v) + div(u)*q) )

# linear form
l  = LinearForm((v,q), integral(domain, f*q))

# Variational problem
equation = find([u,p], forall=[v,q], lhs=a((u,p),(v,q)), rhs=l(v,q))
