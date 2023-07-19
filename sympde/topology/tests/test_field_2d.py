# coding: utf-8

from sympy import Matrix

from sympde.calculus import grad, inner
from sympde.topology import Domain
from sympde.topology import dx
from sympde.topology import VectorFunctionSpace
from sympde.topology import element_of

#==============================================================================

DIM = 2
domain = Domain('Omega', dim=DIM)

# ...
def test_field_2d_1():
    print('============ test_field_2d_1 =============')

#    x, y = domain.coordinates

    W = VectorFunctionSpace('W', domain)
    F = element_of(W, 'F')

    assert( dx(F) == Matrix([[dx(F[0]), dx(F[1])]]) )

    # TODO not working yet => check it for VectorFunction also
#    print(dx(x*F))

    expr = inner(grad(F), grad(F))
    print(expr)


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
