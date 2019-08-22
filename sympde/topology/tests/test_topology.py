# coding: utf-8

from collections import OrderedDict

from sympy.tensor import Indexed

from sympde.topology import InteriorDomain, Union
from sympde.topology import Boundary, NormalVector, TangentVector
from sympde.topology import Connectivity, Edge
from sympde.topology import Domain, ElementDomain
from sympde.topology import Area
from sympde.topology import Interface
from sympde.topology import Line, Square

import os

base_dir = os.path.dirname(os.path.realpath(__file__))
topo_dir = os.path.join(base_dir, 'data')


#==============================================================================
def test_interior_domain():
    D1 = InteriorDomain('D1', dim=2)
    D2 = InteriorDomain('D2', dim=2)

    assert( D1.todict() == OrderedDict([('name', 'D1')]) )
    assert( D2.todict() == OrderedDict([('name', 'D2')]) )

    assert( Union(D2, D1) == Union(D1, D2) )

    D = Union(D1, D2)

    assert(D.dim == 2)
    assert(len(D) == 2)
    assert( D.todict() == [OrderedDict([('name', 'D1')]),
                           OrderedDict([('name', 'D2')])] )

#==============================================================================
def test_topology_1():
    # ... create a domain with 2 subdomains A and B
    A = InteriorDomain('A', dim=2)
    B = InteriorDomain('B', dim=2)

    connectivity = Connectivity()

    bnd_A_1 = Boundary('Gamma_1', A)
    bnd_A_2 = Boundary('Gamma_2', A)
    bnd_A_3 = Boundary('Gamma_3', A)

    bnd_B_1 = Boundary('Gamma_1', B)
    bnd_B_2 = Boundary('Gamma_2', B)
    bnd_B_3 = Boundary('Gamma_3', B)

    connectivity['I'] = (bnd_A_1, bnd_B_2)

    Omega = Domain('Omega',
                   interiors=[A, B],
                   boundaries=[bnd_A_2, bnd_A_3, bnd_B_1, bnd_B_3],
                   connectivity=connectivity)

    interfaces = Omega.interfaces
    assert(isinstance(interfaces, Interface))

    # export
    Omega.export('omega.h5')
    # ...

    # read it again and check that it has the same description as Omega
    D = Domain.from_file('omega.h5')
    assert( D.todict() == Omega.todict() )

#==============================================================================
def test_domain_1():
    Omega_1 = InteriorDomain('Omega_1', dim=2)
    Omega_2 = InteriorDomain('Omega_2', dim=2)

    Gamma_1 = Boundary('Gamma_1', Omega_1)
    Gamma_2 = Boundary('Gamma_2', Omega_2)
    Gamma_3 = Boundary('Gamma_3', Omega_2)

    Omega = Domain('Omega',
                   interiors=[Omega_1, Omega_2],
                   boundaries=[Gamma_1, Gamma_2, Gamma_3])

    assert( Omega.dim == 2 )
    assert( len(Omega.interior) == 2 )
    assert( len(Omega.boundary) == 3 )

#==============================================================================
def test_boundary_1():
    Omega_1 = InteriorDomain('Omega_1', dim=2)

    Gamma_1 = Boundary('Gamma_1', Omega_1)
    Gamma_2 = Boundary('Gamma_2', Omega_1)

    Omega = Domain('Omega',
                   interiors=[Omega_1],
                   boundaries=[Gamma_1, Gamma_2])

    assert(Omega.boundary == Union(Gamma_1, Gamma_2))
    assert(Omega.boundary.complement(Gamma_1) == Gamma_2)
    assert(Omega.boundary - Gamma_1 == Gamma_2)

#==============================================================================
def test_boundary_2():
    Omega_1 = InteriorDomain('Omega_1', dim=2)

    Gamma_1 = Boundary('Gamma_1', Omega_1)
    Gamma_2 = Boundary('Gamma_2', Omega_1)
    Gamma_3 = Boundary('Gamma_3', Omega_1)

    Omega = Domain('Omega',
                   interiors=[Omega_1],
                   boundaries=[Gamma_1, Gamma_2, Gamma_3])

    assert(Omega.boundary == Union(Gamma_1, Gamma_2, Gamma_3))
    assert(Omega.boundary.complement(Gamma_1) == Union(Gamma_2, Gamma_3))
    assert(Omega.boundary - Gamma_1 == Union(Gamma_2, Gamma_3))

#==============================================================================
def test_boundary_3():
    Omega_1 = InteriorDomain('Omega_1', dim=2)

    Gamma_1 = Boundary(r'\Gamma_1', Omega_1, axis=0, ext=-1)
    Gamma_4 = Boundary(r'\Gamma_4', Omega_1, axis=1, ext=1)

    Omega = Domain('Omega',
                   interiors=[Omega_1],
                   boundaries=[Gamma_1, Gamma_4])

    assert(Omega.get_boundary(axis=0, ext=-1) == Gamma_1)
    assert(Omega.get_boundary(axis=1, ext=1) == Gamma_4)

#==============================================================================
def test_element():
    D1 = InteriorDomain('D1', dim=2)
    D2 = InteriorDomain('D2', dim=2)

    D = Union(D1, D2)

    e1 = ElementDomain()

    a = Area(e1)
    print(a)

    a = Area(D1)
    print(a)

    assert(Area(D) ==  Area(D1) + Area(D2))

#==============================================================================
def test_domain_join_line():

    # ... line
    A = Line('A')
    B = Line('B')
    C = Line('C')
    # ...

    # ...
    AB = A.join(B, name = 'AB',
               bnd_minus = A.get_boundary(axis=0, ext=1),
               bnd_plus  = B.get_boundary(axis=0, ext=-1))

    print(AB)
    print(AB.interior)
    print(AB.interfaces)
    print(AB.connectivity)
    print('')
    # ...

    # ...
    AB = A.join(B, name = 'AB',
               bnd_minus = A.get_boundary(axis=0, ext=1),
               bnd_plus  = B.get_boundary(axis=0, ext=-1))
    ABC = AB.join(C, name = 'ABC',
               bnd_minus = A.get_boundary(axis=0, ext=1),
               bnd_plus  = B.get_boundary(axis=0, ext=-1))

    print(ABC)
    print(ABC.interior)
    print(ABC.interfaces)
    print(list(ABC.connectivity.items()))
    print('')
    # ...

#==============================================================================
def test_domain_join_square():

    # ... line
    A = Square('A')
    B = Square('B')
    C = Square('C')
    # ...

    # ...
    AB = A.join(B, name = 'AB',
               bnd_minus = A.get_boundary(axis=0, ext=1),
               bnd_plus  = B.get_boundary(axis=0, ext=-1))

    print(AB)
    print(AB.interior)
    print(AB.interfaces)
    print(AB.connectivity)
    print('')
    # ...

    # ...
    AB = A.join(B, name = 'AB',
               bnd_minus = A.get_boundary(axis=0, ext=1),
               bnd_plus  = B.get_boundary(axis=0, ext=-1))
    ABC = AB.join(C, name = 'ABC',
               bnd_minus = B.get_boundary(axis=0, ext=1),
               bnd_plus  = C.get_boundary(axis=0, ext=-1))

    print(ABC)
    print(ABC.interior)
    print(ABC.interfaces)
    print(list(ABC.connectivity.items()))
    print('')
    # ...



#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
