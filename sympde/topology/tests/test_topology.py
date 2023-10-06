# coding: utf-8

from sympy.tensor import Indexed

from sympde.topology import InteriorDomain, Union
from sympde.topology import Boundary, NormalVector, TangentVector
from sympde.topology import Connectivity, Edge
from sympde.topology import Domain, ElementDomain
from sympde.topology import Area, Mapping
from sympde.topology import Interface
from sympde.topology import Line, Square, Cube
from sympde.topology import IdentityMapping

import os

base_dir = os.path.dirname(os.path.realpath(__file__))
topo_dir = os.path.join(base_dir, 'data')


#==============================================================================
def test_interior_domain():
    D1 = InteriorDomain('D1', dim=2)
    D2 = InteriorDomain('D2', dim=2)

    assert( D1.todict() == {'name': 'D1', 'mapping':'None'} )
    assert( D2.todict() == {'name': 'D2', 'mapping':'None'} )

    assert( Union(D2, D1) == Union(D1, D2) )

    D = Union(D1, D2)

    assert(D.dim == 2)
    assert(len(D) == 2)
    assert( D.todict() == [{'name': 'D1', 'mapping':'None'},
                           {'name': 'D2', 'mapping':'None'}] )


#==============================================================================
def test_topology_1():

    # ... create a domain with 2 subdomains D1 and D2
    A = Square('A')
    B = Square('B')
    # ...

    M1 = Mapping('M1', dim=2)
    M2 = Mapping('M2', dim=2)

    D1 = M1(A)
    D2 = M2(B)

    domains = [D1, D2]
    connectivity = [((0, 0, 1), (1, 0, -1))]
    Omega = Domain.join(domains, connectivity, 'domain')

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

    AB_bnd_minus = A.get_boundary(axis=0, ext=1)
    AB_bnd_plus  = B.get_boundary(axis=0, ext=-1)


    domains = [A, B]
    connectivity = [((0, 0, 1), (1, 0, -1))]
    AB = Domain.join(domains, connectivity, 'AB')

    AB_bnd_minus = A.get_boundary(axis=0, ext=1)
    AB_bnd_plus  = B.get_boundary(axis=0, ext=-1)

    BC_bnd_minus = B.get_boundary(axis=0, ext=1)
    BC_bnd_plus  = C.get_boundary(axis=0, ext=-1)

    assert AB.interior   == Union(A.interior, B.interior)
    assert AB.interfaces == Interface('A|B', AB_bnd_minus, AB_bnd_plus)
    print(AB.connectivity)
    print('')
    # ...

    domains = [A, B, C]
    connectivity = [((0, 0, 1), (1, 0, -1)), ((1, 0, 1), (2, 0, -1))]
    ABC = Domain.join(domains, connectivity, 'ABC')

    assert ABC.interior == Union(A.interior, B.interior, C.interior)
    assert ABC.interfaces == Union(Interface('A|B', AB_bnd_minus, AB_bnd_plus),Interface('B|C', BC_bnd_minus, BC_bnd_plus))
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


    domains = [A, B]
    connectivity = [((0, 0, 1), (1, 0, -1))]
    AB = Domain.join(domains, connectivity, 'AB')

    AB_bnd_minus = A.get_boundary(axis=0, ext=1)
    AB_bnd_plus  = B.get_boundary(axis=0, ext=-1)

    BC_bnd_minus = B.get_boundary(axis=0, ext=1)
    BC_bnd_plus  = C.get_boundary(axis=0, ext=-1)

    print(AB)
    assert AB.interior   == Union(A.interior, B.interior)
    assert AB.interfaces == Interface('A|B', AB_bnd_minus, AB_bnd_plus, ornt=1)
    print(AB.connectivity)
    # ...

    domains = [A, B, C]
    connectivity = [((0, 0, 1),(1, 0, -1)), ((1, 0, 1), (2, 0, -1))]
    ABC = Domain.join(domains, connectivity, 'ABC')


    print(ABC)
    assert ABC.interior == Union(A.interior, B.interior, C.interior)
    assert ABC.interfaces == Union(Interface('A|B', AB_bnd_minus, AB_bnd_plus, ornt=1),
                                   Interface('B|C', BC_bnd_minus, BC_bnd_plus, ornt=1))
    print(list(ABC.connectivity.items()))
    print('')
    # ...

#==============================================================================
def test_get_subdomain():
    A = Square('A')
    B = Square('B')
    C = Square('C')
    # ...

    domains = [A, B]
    connectivity = [((0, 0, 1), (1, 0, -1))]
    AB = Domain.join(domains, connectivity, 'AB')

    # ...

    domains = [A, B, C]
    connectivity = [((0, 0, 1), (1, 0, -1)), ((1, 0, 1), (2, 0, -1))]
    ABC = Domain.join(domains, connectivity, 'ABC')

    A_1 = AB.get_subdomain('A')
    A_2 = ABC.get_subdomain('A')

    assert A_1.boundary == A_2.boundary == A.boundary

    A_pipe_B = ABC.get_subdomain(('A', 'B'))

    assert A_pipe_B.boundary == AB.boundary
    assert A_pipe_B.interfaces == AB.interfaces

    A_pipe_B_pipe_C_1 = ABC.get_subdomain(('A', 'B', 'C'))

    assert A_pipe_B_pipe_C_1 is ABC

    A_pipe_B_pipe_C_2 = ABC.get_subdomain(('A', 'ABC'))

    assert A_pipe_B_pipe_C_2 is ABC

#==============================================================================
def test_2d_domain_without_bnd():

    OmegaLog1 = Square('OmegaLog1', bounds1 = (0,.5), bounds2 = (0,.5))
    mapping_1 = IdentityMapping('M1',2)
    domain_1     = mapping_1(OmegaLog1)
    OmegaLog2 = Square('OmegaLog2', bounds1 = (0,.5), bounds2 = (.5,1.))
    mapping_2 = IdentityMapping('M2',2)
    domain_2     = mapping_2(OmegaLog2)
    OmegaLog3 = Square('OmegaLog3', bounds1 = (.5,1.), bounds2 = (0,.5))
    mapping_3 = IdentityMapping('M3',2)
    domain_3     = mapping_3(OmegaLog3)
    OmegaLog4 = Square('OmegaLog4', bounds1 = (.5,1.), bounds2 = (.5,1.))
    mapping_4 = IdentityMapping('M4',2)
    domain_4     = mapping_4(OmegaLog4)

    domains=[domain_1,domain_2,domain_3,domain_4]

    connectivity = [((0, 0, 1), (2, 0,-1)),
                    ((1, 0, 1), (3, 0,-1)),
                    ((2, 0, 1), (0, 0,-1)),
                    ((3, 0, 1), (1, 0,-1)),
                    ((0, 1, 1), (1, 1,-1)),
                    ((2, 1, 1), (3, 1,-1)),
                    ((1, 1, 1), (0, 1,-1)),
                    ((3, 1, 1), (2, 1,-1))]
    domain = Domain.join(domains, connectivity, 'domain')

    assert len(domain.interior) == 4
    assert len(domain.interfaces) == 8

def test_3d_domain_without_bnd():

    OmegaLog1 = Cube('OmegaLog1', bounds1 = (0,.5), bounds2 = (0,.5), bounds3 = (0,1))
    mapping_1 = IdentityMapping('M1',2)
    domain_1     = mapping_1(OmegaLog1)
    OmegaLog2 = Cube('OmegaLog2', bounds1 = (0,.5), bounds2 = (.5,1.), bounds3 = (0,1))
    mapping_2 = IdentityMapping('M2',2)
    domain_2     = mapping_2(OmegaLog2)
    OmegaLog3 = Cube('OmegaLog3', bounds1 = (.5,1.), bounds2 = (0,.5), bounds3 = (0,1))
    mapping_3 = IdentityMapping('M3',2)
    domain_3     = mapping_3(OmegaLog3)
    OmegaLog4 = Cube('OmegaLog4', bounds1 = (.5,1.), bounds2 = (.5,1.), bounds3 = (0,1))
    mapping_4 = IdentityMapping('M4',2)
    domain_4     = mapping_4(OmegaLog4)

    domains=[domain_1,domain_2,domain_3,domain_4]

    connectivity = [((0, 0, 1),(2, 0, -1),(1,1,1)),
                    ((1,0,1),(3,0,-1),(1,1,1)),
                    ((2,0,1),(0,0,-1),(1,1,1)),
                    ((3,0,1),(1,0,-1),(1,1,1)),
                    ((0,1,1),(1,1,-1),(1,1,1)),
                    ((2,1,1),(3,1,-1),(1,1,1)),
                    ((1,1,1),(0,1,-1),(1,1,1)),
                    ((3,1,1),(2,1,-1),(1,1,1))]
    domain = Domain.join(domains, connectivity, 'domain')

    assert len(domain.interior) == 4
    assert len(domain.interfaces) == 8
#==============================================================================
def test_hash():
    A = Square('A', bounds1=(0, 1), bounds2=(0, 1))
    hash_1 = hash(A)

    A = Square('A', bounds1=(0, 1), bounds2=(1, 2))
    hash_2 = hash(A)

    assert hash_1 != hash_2

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

    # Remove output file generated by test_topology_1()
    fname = 'omega.h5'
    if os.path.exists(fname):
        os.remove(fname)

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
