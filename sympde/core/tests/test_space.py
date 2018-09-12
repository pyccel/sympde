# coding: utf-8

from sympde.core import Domain
from sympde.core import FunctionSpace
from sympde.core import ProductSpace

# ...
def test_space_1d():
    print('============ test_space_1d ==============')

    DIM = 1
    domain = Domain('Omega', dim=DIM)

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    V3 = FunctionSpace('V3', domain)

    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
    U3 = FunctionSpace('U3', domain, shape=2)

    V = ProductSpace(V1, V2, V3)
    assert(V.ldim == DIM)
    assert(V.shape == 3)
    assert(V.name == 'V1V2V3')

    U = ProductSpace(U1, U2, U3)
    assert(U.ldim == DIM)
    assert(U.shape == 4)
    assert(U.name == 'U1U2U3')
# ...

# ...
def test_space_2d():
    print('============ test_space_2d ==============')

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    V3 = FunctionSpace('V3', domain)

    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
    U3 = FunctionSpace('U3', domain, shape=2)

    V = ProductSpace(V1, V2, V3)
    assert(V.ldim == DIM)
    assert(V.shape == 3)
    assert(V.name == 'V1V2V3')

    U = ProductSpace(U1, U2, U3)
    assert(U.ldim == DIM)
    assert(U.shape == 4)
    assert(U.name == 'U1U2U3')
# ...

# ...
def test_space_3d():
    print('============ test_space_3d ==============')

    DIM = 3
    domain = Domain('Omega', dim=DIM)

    V1 = FunctionSpace('V1', domain)
    V2 = FunctionSpace('V2', domain)
    V3 = FunctionSpace('V3', domain)

    U1 = FunctionSpace('U1', domain)
    U2 = FunctionSpace('U2', domain)
    U3 = FunctionSpace('U3', domain, shape=2)

    V = ProductSpace(V1, V2, V3)
    assert(V.ldim == DIM)
    assert(V.shape == 3)
    assert(V.name == 'V1V2V3')

    U = ProductSpace(U1, U2, U3)
    assert(U.ldim == DIM)
    assert(U.shape == 4)
    assert(U.name == 'U1U2U3')
# ...

# .....................................................
if __name__ == '__main__':

    test_space_1d()
    test_space_2d()
    test_space_3d()
