# coding: utf-8

from sympy import exp

from sympde.calculus.core import grad, dot
from sympde.topology.space import ScalarFunctionSpace
from sympde.topology.space import element_of
from sympde.topology.domain import Domain
from sympde.topology.basic import Boundary
from sympde.expr.expr import LinearForm, integral
from sympde.expr.equation import EssentialBC, NewtonIteration

DIM = 2
domain = Domain('Omega', dim=DIM)

#==============================================================================
def test_newton_2d_1():

    # ... abstract model
    B1 = Boundary(r'\Gamma_1', domain)

    V = ScalarFunctionSpace('V', domain)

    x,y = domain.coordinates

    Un = element_of(V, name='Un')
    v  = element_of(V, name='v')
    
    int_0 = lambda expr: integral(domain , expr)

    f  = -4.*exp(-Un)
    l = LinearForm(v, int_0(dot(grad(v), grad(Un)) - f*v) )

    u  = element_of(V, name='u')
    eq = NewtonIteration(l, Un, trials=u)

    expected = int_0(-4.0*u*v*exp(-Un) + dot(grad(u), grad(v)))

    assert(eq.lhs.expr == expected)

    bc = EssentialBC(u, 0, B1)
    eq = NewtonIteration(l, Un, bc=bc, trials=u)


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
