# coding: utf-8

from sympde.topology import Domain
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import ProductSpace
from sympde.topology import H1Space, HcurlSpace, HdivSpace, L2Space

#==============================================================================
def test_space_1d_1():

    DIM = 1
    domain = Domain('Omega', dim=DIM)

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    V3 = FunctionSpace('V3', domain)

    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
#    U3 = FunctionSpace('U3', domain, shape=2)
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

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    V3 = FunctionSpace('V3', domain)

    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
#    U3 = FunctionSpace('U3', domain, shape=2)
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
def test_space_3d_1():

    DIM = 3
    domain = Domain('Omega', dim=DIM)

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    V3 = FunctionSpace('V3', domain)

    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
#    U3 = FunctionSpace('U3', domain, shape=2)
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

    H1    =       FunctionSpace('V0', domain, kind='H1')
    Hcurl = VectorFunctionSpace('V1', domain, kind='Hcurl')
    Hdiv  = VectorFunctionSpace('V2', domain, kind='Hdiv')
    L2    =       FunctionSpace('V3', domain, kind='L2')

    assert(H1.kind == H1Space)
    assert(Hcurl.kind == HcurlSpace)
    assert(Hdiv.kind == HdivSpace)
    assert(L2.kind == L2Space)


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()

test_space_2d_2()
