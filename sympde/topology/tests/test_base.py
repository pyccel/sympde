# coding: utf-8

from sympy.tensor import Indexed

from sympde.topology import InteriorDomain, Union
from sympde.topology import Boundary, NormalVector, TangentVector
from sympde.topology import Topology, Edge
from sympde.topology import Domain

# ...
def test_interior_domain():
    D1 = InteriorDomain('D1', dim=2)
    D2 = InteriorDomain('D2', dim=2)

    D = Union(D1, D2)

    assert(D.dim == 2)
    assert(len(D) == 2)
# ...

# ...
def test_topology_1():
    Omega_1 = InteriorDomain('Omega_1', dim=2)
    Omega_2 = InteriorDomain('Omega_2', dim=2)

    topo = Topology()

    B = Edge('B')

    Gamma_11 = Boundary('Gamma_1', Omega_1)
    Gamma_12 = Boundary('Gamma_2', Omega_2)

    topo[B] = (Gamma_11, Gamma_12)
# ...

# ...
def test_topology_2():
    topo = Topology(filename='square_mp_0.h5')
# ...

# ...
def test_domain_1():
    Omega_1 = InteriorDomain('Omega_1', dim=2)
    Omega_2 = InteriorDomain('Omega_2', dim=2)

    Gamma_1 = Boundary('Gamma_1', Omega_1)
    Gamma_2 = Boundary('Gamma_2', Omega_2)
    Gamma_3 = Boundary('Gamma_3', Omega_2)

    Omega = Domain('Omega', [Omega_1, Omega_2], [Gamma_1, Gamma_2, Gamma_3])

    assert( Omega.dim == 2 )
    assert( len(Omega.interior) == 2 )
    assert( len(Omega.boundary) == 3 )
# ...

# ...
def test_domain_2():
    topo = Topology(filename='square_mp_0.h5')
    Omega = Domain('Omega', topology=topo)

    assert( isinstance(Omega.interior, Union) )
    assert( len(Omega.interior) == 2 )
    assert( len(Omega.boundary) == 8 )
# ...

# .....................................................
if __name__ == '__main__':

    test_interior_domain()
    test_topology_1()
    test_topology_2()
    test_domain_1()
    test_domain_2()
