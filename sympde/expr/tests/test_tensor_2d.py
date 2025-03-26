# coding: utf-8

# TODO: - add assert to every test

from sympy.core.containers import Tuple

from sympde.core     import Constant
from sympde.calculus import grad, dot, curl, div
#from sympde.calculus import laplace
#from sympde.topology import dx
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import elements_of
from sympde.topology import Boundary
from sympde.topology import Domain
from sympde.topology import Mapping

from sympde.expr.expr import BilinearForm, integral
from sympde.expr.evaluation import TensorExpr


#==============================================================================
def test_tensorize_2d_1():

    domain = Domain('Omega', dim=2)

    mu    = Constant('mu'   , is_real=True)

    V = ScalarFunctionSpace('V', domain)
    u, v = elements_of(V, names='u, v')

    int_0 = lambda expr: integral(domain , expr)

    # ...
#    a = BilinearForm((u,v), u*v)
    a = BilinearForm((u,v), int_0(mu*u*v + dot(grad(u),grad(v))))
#    a = BilinearForm((u,v), dot(grad(u),grad(v)))
#    a = BilinearForm((u,v), dx(u)*v)
#    a = BilinearForm((u,v), laplace(u)*laplace(v))

    expr = TensorExpr(a, domain=domain)
    print(expr)
    # ...

#==============================================================================
def test_tensorize_2d_2():

    domain = Domain('Omega', dim=2)

    V = VectorFunctionSpace('V', domain)
    u, v = elements_of(V, names='u, v')

    int_0 = lambda expr: integral(domain , expr)
    # ...
#    a = BilinearForm((u,v), dot(u,v))
    a = BilinearForm((u,v), int_0(curl(u)*curl(v) + div(u)*div(v)))

    expr = TensorExpr(a, domain=domain)
    print(expr)
    # ...

#==============================================================================
def test_tensorize_2d_1_mapping():

    DIM = 2

    M = Mapping('Map', dim=DIM)

    domain = M(Domain('Omega', dim=DIM))

#    mu    = Constant('mu'   , is_real=True)

    V = ScalarFunctionSpace('V', domain)
    u, v = elements_of(V, names='u, v')

    int_0 = lambda expr: integral(domain , expr)
    # ...
#    a = BilinearForm((u,v), u*v)
#    a = BilinearForm((u,v), mu*u*v + dot(grad(u),grad(v)))
    a = BilinearForm((u,v), int_0(dot(grad(u),grad(v))))
#    a = BilinearForm((u,v), dx(u)*v)
#    a = BilinearForm((u,v), laplace(u)*laplace(v))

    expr = TensorExpr(a, domain=domain)
    print(expr)
    # ...

#==============================================================================
def test_tensorize_2d_2_mapping():

    DIM = 2
    M = Mapping('M', dim=DIM)
    domain = M(Domain('Omega', dim=DIM))

    V = VectorFunctionSpace('V', domain)
    u, v = elements_of(V, names='u, v')

    c = Constant('c')

    int_0 = lambda expr: integral(domain , expr)

    a = BilinearForm((u,v), int_0(c * div(v) * div(u) + curl(v) * curl(u)))
    expr = TensorExpr(a, domain=domain)
    print(expr)

#==============================================================================
def test_tensorize_2d_3():

    domain = Domain('Omega', dim=2)

    V = ScalarFunctionSpace('V', domain)
    u, v = elements_of(V, names='u,v')

    bx = Constant('bx')
    by = Constant('by')
    b  = Tuple(bx, by)

    expr = integral(domain, dot(b, grad(v)) * dot(b, grad(u)))
    a = BilinearForm((u,v), expr)

    print(TensorExpr(a, domain=domain))
    print('')


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
