# coding: utf-8

from sympde.core     import Constant
from sympde.calculus import grad, dot
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import element_of
from sympde.topology import Domain, Boundary, NormalVector
from sympde.expr     import EssentialBC

#==============================================================================
def test_essential_bc_1():
    domain = Domain('Omega', dim=2)

    V = ScalarFunctionSpace('V', domain)
    W = VectorFunctionSpace('W', domain)

    v = element_of(V, name='v')
    w = element_of(W, name='w')

    B1 = Boundary(r'\Gamma_1', domain)
    nn = NormalVector('nn')

    # ... scalar case
    bc = EssentialBC(v, 0, B1)

    assert( bc.variable == v )
    assert( bc.order == 0 )
    assert( bc.normal_component == False )
    assert( bc.index_component == None )
    # ...

    # ... scalar case
    bc = EssentialBC(dot(grad(v), nn), 0, B1)

    assert( bc.variable == v )
    assert( bc.order == 1 )
    assert( bc.normal_component == False )
    assert( bc.index_component == None )
    # ...

    # ... vector case
    bc = EssentialBC(w, 0, B1)
    assert( bc.variable == w )
    assert( bc.order == 0 )
    assert( bc.normal_component == False )
    assert( bc.index_component == [0,1] )
    # ...

    # ... vector case
    bc = EssentialBC(dot(w, nn), 0, B1)

    assert( bc.variable == w )
    assert( bc.order == 0 )
    assert( bc.normal_component == True )
    assert( bc.index_component == None )
    # ...

    # ... vector case
    bc = EssentialBC(w[0], 0, B1)
    assert( bc.variable == w )
    assert( bc.order == 0 )
    assert( bc.normal_component == False )
    assert( bc.index_component == [0] )
    # ...


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
