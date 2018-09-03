# coding: utf-8

from symfe.core import BasicSobolevSpace
from symfe.core import ProductSpace

# ...
def test_space():
    print('============ test_space ==============')

    V1 = BasicSobolevSpace('V1', ldim=2)
    V2 = BasicSobolevSpace('V2', ldim=2)
    V3 = BasicSobolevSpace('V3', ldim=2)

    U1 = BasicSobolevSpace('U1', ldim=2)
    U2 = BasicSobolevSpace('U2', ldim=2)
    U3 = BasicSobolevSpace('U3', ldim=2, shape=2)

    V = ProductSpace(V1, V2, V3)
    assert(V.ldim == 2)
    assert(V.shape == 3)
    assert(V.name == 'V1V2V3')

    U = ProductSpace(U1, U2, U3)
    assert(U.ldim == 2)
    assert(U.shape == 4)
    assert(U.name == 'U1U2U3')
# ...

# .....................................................
if __name__ == '__main__':

    test_space()
