# coding: utf-8

from sympy import cos
#from sympy import exp

from sympde.core     import Constant
from sympde.calculus import grad, dot
#from sympde.topology import dx, dy
from sympde.topology import ScalarFunctionSpace
from sympde.topology import element_of
from sympde.topology import Domain, Boundary
from sympde.expr     import BilinearForm, LinearForm, integral
from sympde.expr     import TerminalExpr
from sympde.expr     import find
from sympde.expr     import EssentialBC
#from sympde.expr     import NewtonIteration

DIM = 2
domain = Domain('Omega', dim=DIM)

#==============================================================================
def test_find_2d_1():

    V = ScalarFunctionSpace('V', domain)
    U = ScalarFunctionSpace('U', domain)

    v = element_of(V, name='v')
    u = element_of(U, name='u')

    x,y = domain.coordinates

    alpha = Constant('alpha')
    kappa = Constant('kappa', real=True)
    eps   = Constant('eps', real=True)

    B1 = Boundary(r'\Gamma_1', domain)
    B2 = Boundary(r'\Gamma_2', domain)
    B3 = Boundary(r'\Gamma_3', domain)
    
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)
    int_2 = lambda expr: integral(B2, expr)
    # ... bilinear/linear forms
    expr = dot(grad(v), grad(u))
    a1 = BilinearForm((v,u), int_0(expr))

    expr = v*u
    a2 = BilinearForm((v,u), int_0(expr))

    expr = v*u
    a_B1 = BilinearForm((v, u), int_1(expr))

    expr = x*y*v
    l1 = LinearForm(v, int_0(expr))

    expr = cos(x+y)*v
    l2 = LinearForm(v, int_0(expr))

    expr = x*y*v
    l_B2 = LinearForm(v, int_1(expr))
    # ...

    # ...
    equation = find(u, forall=v, lhs=a1(u,v), rhs=l1(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + alpha*a2(v,u))
    equation = find(u, forall=v, lhs=a(u,v), rhs=l1(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    equation = find(u, forall=v, lhs=a(u,v), rhs=l1(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    l = LinearForm(v, l1(v) + l2(v))
    equation = find(u, forall=v, lhs=a(u,v), rhs=l(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    l = LinearForm(v, l1(v) + alpha*l_B2(v))
    equation = find(u, forall=v, lhs=a(u,v), rhs=l(v))
    # ...

    # ... using bc
    equation = find(u, forall=v, lhs=a(u,v), rhs=l(v), bc=EssentialBC(u, 0, B1))
    # ...

##==============================================================================
# TODO shall we have 'find' with Newton?
#def test_find_newton_2d_1():
#
#    # ... abstract model
#    B1 = Boundary(r'\Gamma_1', domain)
#
#    V = ScalarFunctionSpace('V', domain)
#
#    x,y = domain.coordinates
#
#    Un = element_of(V, name='Un')
#
#    v = element_of(V, name='v')
#
#    f  = -4.*exp(-Un)
#    l = LinearForm(v, dot(grad(v), grad(Un)) - f*v )
#
#    eq = NewtonIteration(l, Un, trials='u')
#    # ...
#
#    u = element_of(V, name='u')
#
#    # ...
#    expected =  -4.0*u*v*exp(-Un) + dx(u)*dx(v) + dy(u)*dy(v)
#
#    expr = TerminalExpr(eq.lhs)[0]
#    assert(expr.expr == expected)
#    # ...
#
#    # ...
#    expected = -4.0*v*exp(-Un) - dx(Un)*dx(v) - dy(Un)*dy(v)
#
#    expr = TerminalExpr(eq.rhs)[0]
#    assert(expr.expr == expected)
#    # ...
#
#    # ...
#    bc = EssentialBC(u, 0, B1)
#
#    # can be called with rhs = 0
#    equation = find(u, forall=v, lhs=l, rhs=0, fields=Un, bc=bc)
#    print(equation)
#
#    # or without rhs
#    equation = find(u, forall=v, lhs=l, fields=Un, bc=bc)
#    print(equation)
#    # ...

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
