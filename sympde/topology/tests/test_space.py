# coding: utf-8

from sympde.calculus import grad, div
from sympde.topology import Domain
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
#from sympde.topology import ProductSpace
from sympde.topology import element_of
from sympde.topology import H1Space, HcurlSpace, HdivSpace, L2Space, UndefinedSpace
from sympde.topology import ScalarFunction, VectorFunction
from sympde.topology import Projector

#==============================================================================
def test_space_1d_1():

    DIM = 1
    domain = Domain('Omega', dim=DIM)

    V1 = ScalarFunctionSpace('V1', domain)
    V2 = ScalarFunctionSpace('V2', domain)
    V3 = ScalarFunctionSpace('V3', domain)

    U1 = ScalarFunctionSpace('U1', domain)
    U2 = ScalarFunctionSpace('U2', domain)
#    U3 = ScalarFunctionSpace('U3', domain, shape=2)
#
#    V = ProductSpace(V1, V2, V3)
#    assert(V.ldim == DIM)
#    assert(V.shape == 3)
#    assert(V.name == 'V1V2V3')
#
#    U = ProductSpace(U1, U2, U3)
#    assert(U.ldim == DIM)
#    assert(U.shape == 4)
#    assert(U.name == 'U1U2U3')

#==============================================================================
def test_space_2d_1():

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    V1 = ScalarFunctionSpace('V1', domain)
    V2 = ScalarFunctionSpace('V2', domain)
    V3 = ScalarFunctionSpace('V3', domain)

    U1 = ScalarFunctionSpace('U1', domain)
    U2 = ScalarFunctionSpace('U2', domain)
    U3 = VectorFunctionSpace('U3', domain)

    # ...
    V = V1 * V2 * V3

    assert(V.ldim == DIM)
    assert(V.shape == 3)
    assert(V.name == 'V1V2V3')
    # ...

    # ...
    U = U1 * U2 * U3

    assert(U.ldim == DIM)
    assert(U.shape == 4)
    assert(U.name == 'U1U2U3')
    # ...

    # ...
    v1, v2, v3 = element_of(V, 'v1, v2, v3')

    assert(isinstance(v1, ScalarFunction))
    assert(isinstance(v2, ScalarFunction))
    assert(isinstance(v3, ScalarFunction))
    assert(v1.space is V1)
    assert(v2.space is V2)
    assert(v3.space is V3)
    assert(v1.name == 'v1')
    assert(v2.name == 'v2')
    assert(v3.name == 'v3')
    # ...

    # ...
    u1, u2, u3 = element_of(U, 'u1, u2, u3')

    assert(isinstance(u1, ScalarFunction))
    assert(isinstance(u2, ScalarFunction))
    assert(isinstance(u3, VectorFunction))
    assert(u1.name == 'u1')
    assert(u2.name == 'u2')
    assert(u3.name == 'u3')
    # ...


#==============================================================================
def test_space_3d_1():

    DIM = 3
    domain = Domain('Omega', dim=DIM)

    V1 = ScalarFunctionSpace('V1', domain)
    V2 = ScalarFunctionSpace('V2', domain)
    V3 = ScalarFunctionSpace('V3', domain)

    U1 = ScalarFunctionSpace('U1', domain)
    U2 = ScalarFunctionSpace('U2', domain)
#    U3 = ScalarFunctionSpace('U3', domain, shape=2)
#
#    V = ProductSpace(V1, V2, V3)
#    assert(V.ldim == DIM)
#    assert(V.shape == 3)
#    assert(V.name == 'V1V2V3')
#
#    U = ProductSpace(U1, U2, U3)
#    assert(U.ldim == DIM)
#    assert(U.shape == 4)
#    assert(U.name == 'U1U2U3')

#==============================================================================
def test_space_2d_2():

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    H1    =       ScalarFunctionSpace('V0', domain, kind='H1')
    Hcurl = VectorFunctionSpace('V1', domain, kind='Hcurl')
    Hdiv  = VectorFunctionSpace('V2', domain, kind='Hdiv')
    L2    =       ScalarFunctionSpace('V3', domain, kind='L2')
    V     =       ScalarFunctionSpace('V', domain, kind=None)
    W     = VectorFunctionSpace('W', domain, kind=None)

    assert(H1.kind    == H1Space)
    assert(Hcurl.kind == HcurlSpace)
    assert(Hdiv.kind  == HdivSpace)
    assert(L2.kind    == L2Space)
    assert(V.kind     == UndefinedSpace)
    assert(W.kind     == UndefinedSpace)

    assert(H1.regularity    > L2.regularity)
    assert(H1.regularity    > Hcurl.regularity)
    assert(Hcurl.regularity > L2.regularity)

#==============================================================================
def test_projector_2d_1():

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain, kind=None)
    W = VectorFunctionSpace('W', domain, kind=None)

    v, w = element_of(V*W, ['v', 'w'])

    # ...
    P_V = Projector(V)
    assert(P_V.space == V)

    Pv = P_V(v)
    assert(isinstance(Pv, ScalarFunction))
    assert(Pv == v)
    assert(grad(Pv**2) == 2*v*grad(v))

    Pdiv_w = P_V(div(w))
    assert(isinstance(Pdiv_w, ScalarFunction))
    # ...

    # ...
    P_W = Projector(W)
    assert(P_W.space == W)

    Pw = P_W(w)
    assert(isinstance(Pw, VectorFunction))
    assert(Pw == w)

    Pgrad_v = P_W(grad(v))
    assert(isinstance(Pgrad_v, VectorFunction))
    assert(P_W(Pgrad_v) == Pgrad_v)
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

#test_space_operators_2d_1()
