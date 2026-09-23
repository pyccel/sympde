import os
import copy
import numpy as np
import pytest

from sympde.topology import InteriorDomain, Union
from sympde.topology import Boundary
from sympde.topology import Domain, ElementDomain
from sympde.topology import Area, Mapping, InterfaceMapping, MultiPatchMapping
from sympde.topology import Interface, PatchVertex, SharedVertex
from sympde.topology import Line, Square, Cube
from sympde.topology import IdentityMapping

base_dir = os.path.dirname(os.path.realpath(__file__))
topo_dir = os.path.join(base_dir, 'data')

#==============================================================================
def test_interior_domain():
    D1 = InteriorDomain('D1', dim=2)
    D2 = InteriorDomain('D2', dim=2)

    assert D1.todict() == {'name': 'D1', 'mapping':'None'}
    assert D2.todict() == {'name': 'D2', 'mapping':'None'}

    assert Union(D2, D1) == Union(D1, D2)

    D = Union(D1, D2)

    assert D.dim == 2
    assert len(D) == 2
    assert D.todict() == [{'name': 'D1', 'mapping':'None'},
                          {'name': 'D2', 'mapping':'None'}]

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

    patches = [D1, D2]
    connectivity = [((0, 0, 1), (1, 0, -1), +1)]
    Omega = Domain.join(patches, connectivity, 'domain')

    interfaces = Omega.interfaces
    assert isinstance(interfaces, Interface)

    # export
    Omega.export('omega.h5')
    # ...

    # read it again and check that it has the same description as Omega
    D = Domain.from_file('omega.h5')

    assert D.todict() == Omega.todict()

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

    assert Omega.dim == 2
    assert len(Omega.interior) == 2
    assert len(Omega.boundary) == 3

#==============================================================================
def test_boundary_1():
    Omega_1 = InteriorDomain('Omega_1', dim=2)

    Gamma_1 = Boundary('Gamma_1', Omega_1)
    Gamma_2 = Boundary('Gamma_2', Omega_1)

    Omega = Domain('Omega',
                   interiors=[Omega_1],
                   boundaries=[Gamma_1, Gamma_2])

    assert Omega.boundary == Union(Gamma_1, Gamma_2)
    assert Omega.boundary.complement(Gamma_1) == Gamma_2
    assert Omega.boundary - Gamma_1 == Gamma_2

#==============================================================================
def test_boundary_2():
    Omega_1 = InteriorDomain('Omega_1', dim=2)

    Gamma_1 = Boundary('Gamma_1', Omega_1)
    Gamma_2 = Boundary('Gamma_2', Omega_1)
    Gamma_3 = Boundary('Gamma_3', Omega_1)

    Omega = Domain('Omega',
                   interiors=[Omega_1],
                   boundaries=[Gamma_1, Gamma_2, Gamma_3])

    assert Omega.boundary == Union(Gamma_1, Gamma_2, Gamma_3)
    assert Omega.boundary.complement(Gamma_1) == Union(Gamma_2, Gamma_3)
    assert Omega.boundary - Gamma_1 == Union(Gamma_2, Gamma_3)

#==============================================================================
def test_boundary_3():
    Omega_1 = InteriorDomain('Omega_1', dim=2)

    Gamma_1 = Boundary(r'\Gamma_1', Omega_1, axis=0, ext=-1)
    Gamma_4 = Boundary(r'\Gamma_4', Omega_1, axis=1, ext= 1)

    Omega = Domain('Omega',
                   interiors=[Omega_1],
                   boundaries=[Gamma_1, Gamma_4])

    assert Omega.get_boundary(axis=0, ext=-1) == Gamma_1
    assert Omega.get_boundary(axis=1, ext= 1) == Gamma_4

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

    assert Area(D) ==  Area(D1) + Area(D2)

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
    connectivity = [((0, 0, 1), (1, 0, -1), None)]
    AB = Domain.join(domains, connectivity, 'AB')

    AB_bnd_minus = A.get_boundary(axis=0, ext=1)
    AB_bnd_plus  = B.get_boundary(axis=0, ext=-1)

    BC_bnd_minus = B.get_boundary(axis=0, ext=1)
    BC_bnd_plus  = C.get_boundary(axis=0, ext=-1)

    assert AB.interior   == Union(A.interior, B.interior)
    assert AB.interfaces == Interface(
        'A|B', AB_bnd_minus, AB_bnd_plus, None)
    print(AB.connectivity)
    print('')
    # ...

    domains = [A, B, C]
    connectivity = [((0, 0, 1), (1, 0, -1), None),
                    ((1, 0, 1), (2, 0, -1), None)]
    ABC = Domain.join(domains, connectivity, 'ABC')

    assert ABC.interior == Union(A.interior, B.interior, C.interior)
    assert ABC.interfaces == Union(
        Interface('A|B', AB_bnd_minus, AB_bnd_plus, None),
        Interface('B|C', BC_bnd_minus, BC_bnd_plus, None))
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


    patches = [A, B]
    connectivity = [((0, 0, 1), (1, 0, -1), +1)]
    AB = Domain.join(patches, connectivity, 'AB')

    AB_bnd_minus = A.get_boundary(axis=0, ext= 1)
    AB_bnd_plus  = B.get_boundary(axis=0, ext=-1)

    BC_bnd_minus = B.get_boundary(axis=0, ext= 1)
    BC_bnd_plus  = C.get_boundary(axis=0, ext=-1)

    print(AB)
    assert AB.interior   == Union(A.interior, B.interior)
    assert AB.interfaces == Interface(
        'A|B', AB_bnd_minus, AB_bnd_plus, +1)
    print(AB.connectivity)
    # ...

    patches = [A, B, C]
    connectivity = [((0, 0, 1), (1, 0, -1), +1),
                    ((1, 0, 1), (2, 0, -1), +1)]
    ABC = Domain.join(patches, connectivity, 'ABC')


    print(ABC)
    assert ABC.interior == Union(A.interior, B.interior, C.interior)
    assert ABC.interfaces == Union(
        Interface('A|B', AB_bnd_minus, AB_bnd_plus, +1),
        Interface('B|C', BC_bnd_minus, BC_bnd_plus, +1))
    print(list(ABC.connectivity.items()))
    print('')
    # ...

#==============================================================================
def test_get_subdomain():
    A = Square('A')
    B = Square('B')
    C = Square('C')
    # ...

    patches = [A, B]
    connectivity = [((0, 0, 1), (1, 0, -1), +1)]
    AB = Domain.join(patches, connectivity, 'AB')

    # ...

    patches = [A, B, C]
    connectivity = [((0, 0, 1), (1, 0, -1), +1),
                    ((1, 0, 1), (2, 0, -1), +1)]
    ABC = Domain.join(patches, connectivity, 'ABC')

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
def test_get_subdomain_of_three_from_four_patches():
    patches = [Square(name) for name in 'ABCD']
    connectivity = [
        ((0, 0, +1), (1, 0, -1), +1),
        ((1, 0, +1), (2, 0, -1), +1),
        ((2, 0, +1), (3, 0, -1), +1),
    ]
    domain = Domain.join(patches, connectivity, 'ABCD')
    subdomain = domain.get_subdomain(('A', 'B', 'C'))
    expected = Domain.join(patches[:3], connectivity[:2], 'expected')

    assert subdomain.name == 'A|B|C'
    assert subdomain.patches == expected.patches
    assert subdomain.boundary == expected.boundary
    assert tuple(subdomain.interface_map) == ('A|B', 'B|C')
    assert subdomain.interface_map['A|B'] is domain.interface_map['A|B']
    assert subdomain.interface_map['B|C'] is domain.interface_map['B|C']


#==============================================================================
def test_get_subdomain_preserves_repeated_and_self_interfaces():
    patch_a = Square('A')
    patch_b = Square('B')
    patch_c = Square('C')
    domain = Domain.join(
        [patch_a, patch_b, patch_c],
        [
            ((0, 1, -1), (0, 1, +1), +1),
            ((0, 0, -1), (1, 0, -1), +1),
            ((0, 0, +1), (1, 0, +1), +1),
            ((1, 1, +1), (2, 1, -1), +1),
        ],
        'ABC')

    patch_subdomain = domain.get_subdomain('A')
    pair_subdomain = domain.get_subdomain(('A', 'B'))

    assert tuple(patch_subdomain.interface_map) == ('A|A',)
    assert patch_subdomain.interface_map['A|A'] is \
        domain.interface_map['A|A']
    assert {
        (side.domain.name, side.axis, side.ext)
        for side in patch_subdomain.exterior_sides
    } == {('A', 0, -1), ('A', 0, +1)}

    assert tuple(pair_subdomain.interface_map) == (
        'A|A', 'A|B', 'A|B#2')
    assert all(
        pair_subdomain.interface_map[name] is domain.interface_map[name]
        for name in pair_subdomain.interface_map)
    assert ('B', 1, +1) in {
        (side.domain.name, side.axis, side.ext)
        for side in pair_subdomain.exterior_sides
    }


#==============================================================================
def test_get_subdomain_preserves_mapped_3d_topology():
    logical_patches = [Cube(name) for name in 'ABCD']
    patches = [
        IdentityMapping(f'F{name}', dim=3)(patch)
        for name, patch in zip('ABCD', logical_patches)
    ]
    domain = Domain.join(
        patches,
        [
            ((0, 0, +1), (1, 1, -1), (-1, +1, -1)),
            ((1, 2, +1), (2, 0, -1), (+1, -1, +1)),
            ((2, 1, +1), (3, 2, -1), (-1, -1, +1)),
        ],
        'mapped-ABCD')
    selected_names = tuple(patch.interior.name for patch in patches[:3])
    subdomain = domain.get_subdomain(selected_names)
    retained_interfaces = tuple(domain.interface_map.values())[:2]
    cut_side = tuple(domain.interface_map.values())[2].minus

    assert isinstance(subdomain.mapping, MultiPatchMapping)
    assert subdomain.logical_domain is not None
    assert subdomain.patches == tuple(patch.interior for patch in patches[:3])
    assert tuple(subdomain.interface_map.values()) == retained_interfaces
    assert cut_side in subdomain.exterior_sides
    assert tuple(
        interface.axis_map for interface in subdomain.interface_map.values()
    ) == tuple(interface.axis_map for interface in retained_interfaces)
    assert set(subdomain.mapping.mappings) == {
        patch.interior.logical_domain for patch in patches[:3]}

    for interface in subdomain.interface_map.values():
        logical_interface = interface.logical_domain
        assert logical_interface is not None
        assert subdomain.logical_domain.interface_map[
            str(logical_interface.name)] is logical_interface


#==============================================================================
def test_domain_join_rejects_nested_multipatch_inputs():
    patch_a = Square('A')
    patch_b = Square('B')
    patch_c = Square('C')
    joined = Domain.join(
        [patch_a, patch_b],
        [((0, 0, +1), (1, 0, -1), +1)],
        'AB')

    with pytest.raises(ValueError, match='atomic patch domains'):
        Domain.join([joined, patch_c], [], 'ABC')

#==============================================================================
def test_2d_domain_without_bnd():

    OmegaLog1 = Square('OmegaLog1', bounds1 = (0,.5), bounds2 = (0,.5))
    mapping_1 = IdentityMapping('M1', 2)
    domain_1  = mapping_1(OmegaLog1)
    OmegaLog2 = Square('OmegaLog2', bounds1 = (0,.5), bounds2 = (.5,1.))
    mapping_2 = IdentityMapping('M2', 2)
    domain_2  = mapping_2(OmegaLog2)
    OmegaLog3 = Square('OmegaLog3', bounds1 = (.5,1.), bounds2 = (0,.5))
    mapping_3 = IdentityMapping('M3', 2)
    domain_3  = mapping_3(OmegaLog3)
    OmegaLog4 = Square('OmegaLog4', bounds1 = (.5,1.), bounds2 = (.5,1.))
    mapping_4 = IdentityMapping('M4', 2)
    domain_4  = mapping_4(OmegaLog4)

    patches =  [domain_1, domain_2, domain_3, domain_4]
    connectivity = [((0, 0, 1), (2, 0,-1), +1),
                    ((1, 0, 1), (3, 0,-1), +1),
                    ((2, 0, 1), (0, 0,-1), +1),
                    ((3, 0, 1), (1, 0,-1), +1),
                    ((0, 1, 1), (1, 1,-1), +1),
                    ((2, 1, 1), (3, 1,-1), +1),
                    ((1, 1, 1), (0, 1,-1), +1),
                    ((3, 1, 1), (2, 1,-1), +1)]
    domain = Domain.join(patches, connectivity, 'domain')

    assert len(domain.interior) == 4
    assert len(domain.interfaces) == 8

    assert domain.boundary is None

#==============================================================================
def test_3d_domain_with_doubly_connected_patches():

    OmegaLog1 = Cube('OmegaLog1', bounds1 = (0,.5), bounds2 = (0,.5), bounds3 = (0,1))
    mapping_1 = IdentityMapping('M1', 3)
    domain_1  = mapping_1(OmegaLog1)
    OmegaLog2 = Cube('OmegaLog2', bounds1 = (0,.5), bounds2 = (.5,1.), bounds3 = (0,1))
    mapping_2 = IdentityMapping('M2', 3)
    domain_2  = mapping_2(OmegaLog2)
    OmegaLog3 = Cube('OmegaLog3', bounds1 = (.5,1.), bounds2 = (0,.5), bounds3 = (0,1))
    mapping_3 = IdentityMapping('M3', 3)
    domain_3  = mapping_3(OmegaLog3)
    OmegaLog4 = Cube('OmegaLog4', bounds1 = (.5,1.), bounds2 = (.5,1.), bounds3 = (0,1))
    mapping_4 = IdentityMapping('M4', 3)
    domain_4  = mapping_4(OmegaLog4)

    patches = [domain_1, domain_2, domain_3, domain_4]
    connectivity = [((0, 0, 1), (2, 0,-1), (+1, +1, +1)),
                    ((1, 0, 1), (3, 0,-1), (+1, +1, +1)),
                    ((2, 0, 1), (0, 0,-1), (+1, +1, +1)),
                    ((3, 0, 1), (1, 0,-1), (+1, +1, +1)),
                    ((0, 1, 1), (1, 1,-1), (+1, +1, +1)),
                    ((2, 1, 1), (3, 1,-1), (+1, +1, +1)),
                    ((1, 1, 1), (0, 1,-1), (+1, +1, +1)),
                    ((3, 1, 1), (2, 1,-1), (+1, +1, +1))]
    domain = Domain.join(patches, connectivity, 'domain')

    assert len(domain.interior) == 4
    assert len(domain.interfaces) == 8

    assert isinstance(domain.boundary, Union)
    assert all(isinstance(b, Boundary) for b in domain.boundary)
    assert len(domain.boundary) == 8

#==============================================================================
def test_interface_orientation_axis_map():
    domain = Domain.join(
        [Cube('A'), Cube('B')],
    [((0, 0, +1), (1, 1, -1), (-1, +1, -1))],
    'AB')
    interface = domain.interfaces

    assert interface.orientation == (-1, +1, -1)
    assert interface.axis_map == ((1, 2, +1), (2, 0, -1))


#==============================================================================
@pytest.mark.parametrize('data, permutation, directions', [
    ((+1, +1, +1), (0, 1), (+1, +1)),
    ((+1, +1, -1), (0, 1), (+1, -1)),
    ((+1, -1, +1), (0, 1), (-1, +1)),
    ((+1, -1, -1), (0, 1), (-1, -1)),
    ((-1, +1, +1), (1, 0), (+1, +1)),
    ((-1, +1, -1), (1, 0), (+1, -1)),
    ((-1, -1, +1), (1, 0), (-1, +1)),
    ((-1, -1, -1), (1, 0), (-1, -1)),
])
def test_public_interface_orientation_codec_3d(data, permutation, directions):
    domain = Domain.join(
        [Cube('A'), Cube('B')],
        [((0, 0, +1), (1, 1, -1), data)],
        'AB')
    interface = domain.interfaces
    plus_tangents = (0, 2)

    assert interface.orientation == data
    assert interface.axis_map == tuple(
        (minus_axis, plus_tangents[plus_position], direction)
        for minus_axis, plus_position, direction
        in zip((1, 2), permutation, directions)
    )


#==============================================================================
def test_public_interface_orientation_codec_1d_and_2d():
    line_domain = Domain.join(
        [Line('A'), Line('B')],
        [((0, 0, +1), (1, 0, -1), None)],
        'lines')
    forward = Domain.join(
        [Square('C'), Square('D')],
        [((0, 0, +1), (1, 0, -1), +1)],
        'forward')
    reversed_domain = Domain.join(
        [Square('E'), Square('F')],
        [((0, 0, +1), (1, 0, -1), -1)],
        'reversed')

    assert line_domain.interfaces.orientation is None
    assert line_domain.interfaces.axis_map == ()
    assert forward.interfaces.orientation == +1
    assert reversed_domain.interfaces.orientation == -1

    with pytest.raises(ValueError, match='1D'):
        Domain.join(
            [Line('G'), Line('H')],
            [((0, 0, +1), (1, 0, -1), +1)], 'invalid')
    with pytest.raises(ValueError, match='2D'):
        Domain.join(
            [Square('G'), Square('H')],
            [((0, 0, +1), (1, 0, -1), 0)], 'invalid')
    with pytest.raises(TypeError, match='3D'):
        Domain.join(
            [Cube('G'), Cube('H')],
            [((0, 0, +1), (1, 0, -1), (+1, +1))], 'invalid')
    with pytest.raises(ValueError, match='each 3D'):
        Domain.join(
            [Cube('I'), Cube('J')],
            [((0, 0, +1), (1, 0, -1), (+1, 0, +1))], 'invalid')


#==============================================================================
@pytest.mark.parametrize('orientation', [1.0, True])
def test_2d_interface_orientation_rejects_non_integer_values(orientation):
    with pytest.raises(TypeError, match='2D interface orientation.*integer'):
        Domain.join(
            [Square('A'), Square('B')],
            [((0, 0, +1), (1, 0, -1), orientation)],
            'invalid')


#==============================================================================
@pytest.mark.parametrize('orientation', [
    (1.9, -1.2, 1.0),
    (+1, -1, True),
])
def test_3d_interface_orientation_rejects_non_integer_values(orientation):
    with pytest.raises(TypeError, match='3D.*contain integers'):
        Domain.join(
            [Cube('A'), Cube('B')],
            [((0, 0, +1), (1, 0, -1), orientation)],
            'invalid')


#==============================================================================
@pytest.mark.parametrize('side, expected', [
    ((0, 0.0, +1), 'axis'),
    ((0, False, +1), 'axis'),
    ((0, 0, 1.0), 'extremity'),
    ((0, 0, True), 'extremity'),
])
def test_interface_sides_reject_non_integer_axis_and_extremity(side, expected):
    with pytest.raises(TypeError, match=rf'{expected} must be an integer'):
        Domain.join(
            [Square('A'), Square('B')],
            [(side, (1, 0, -1), +1)],
            'invalid')


#==============================================================================
def test_boundary_lookup_rejects_non_integer_axis_and_extremity():
    patch = Square('A')

    with pytest.raises(TypeError, match='axis must be an integer'):
        patch.get_boundary(axis=0.0, ext=+1)
    with pytest.raises(TypeError, match='extremity must be an integer'):
        patch.get_boundary(axis=0, ext=1.0)
    with pytest.raises(TypeError, match='axis must be an integer'):
        Boundary('invalid', patch.interior, axis=False, ext=+1)


#==============================================================================
def test_domain_join_accepts_omitted_connectivity():
    patch_a = Square('A')
    patch_b = Square('B')

    assert Domain.join([patch_a], name='unused') is patch_a

    domain = Domain.join([patch_a, patch_b], name='AB')
    assert domain.patches == (patch_a.interior, patch_b.interior)
    assert domain.interfaces is None
    assert len(domain.exterior_sides) == 8


#==============================================================================
def test_domain_join_accepts_numpy_integer_patch_indices():
    domain = Domain.join(
        [Square('A'), Square('B')],
        [((np.int64(0), 0, +1), (np.int32(1), 0, -1), +1)],
        'AB')

    interface = domain.interfaces
    assert interface.minus.domain.name == 'A'
    assert interface.plus.domain.name == 'B'


#==============================================================================
def test_domain_join_rejects_duplicate_patch_identity_and_names():
    patch = Square('A')
    with pytest.raises(ValueError, match='same object twice'):
        Domain.join([patch, patch], [], 'duplicate-object')

    with pytest.raises(ValueError, match='interior names must be unique'):
        Domain.join(
            [Square('A', bounds1=(0, 1)),
             Square('A', bounds1=(1, 2))],
            [],
            'duplicate-name')

    logical_patch = Square('logical')
    mapped_patches = [
        IdentityMapping(name, dim=2)(logical_patch)
        for name in ('F', 'G')
    ]
    with pytest.raises(ValueError, match='serialization must be unique'):
        Domain.join(mapped_patches, [], 'duplicate-logical-name')


#==============================================================================
def test_named_multipatch_api_preserves_single_patch_domain_api():
    single_patch = Square('single')

    assert Domain.join([single_patch], [], 'unused') is single_patch
    assert single_patch.patches == (single_patch.interior,)
    assert single_patch.interface_map is single_patch.connectivity
    assert tuple(single_patch.interface_map.values()) == ()
    assert single_patch.exterior_sides == single_patch.boundary.as_tuple()
    assert single_patch.shared_vertices == ()

    patch_a = Square('A')
    patch_b = Square('B')
    domain = Domain.join(
        [patch_a, patch_b],
        interfaces=[(
            (patch_a, 0, +1),
            (patch_b, 1, -1),
            -1,
        )],
        name='AB')

    interface, = domain.interface_map.values()
    assert domain.patches == (patch_a.interior, patch_b.interior)
    assert interface is domain.interfaces
    assert interface.minus_side is interface.minus
    assert interface.plus_side is interface.plus
    assert interface.minus_side.patch is patch_a.interior
    assert interface.minus_side.normal_axis == interface.minus.axis == 0
    assert interface.minus_side.ext == +1
    assert interface.orientation == -1
    assert interface.axis_map == ((1, 0, -1),)
    assert interface.patch_dim == 2
    assert interface.intrinsic_dim == 1
    assert len(domain.exterior_sides) == 6
    assert len(domain.shared_vertices) == 2
    assert all(isinstance(vertex, SharedVertex)
               for vertex in domain.shared_vertices)
    assert all(isinstance(corner, PatchVertex)
               for vertex in domain.shared_vertices
               for corner in vertex.corners)

    with pytest.raises(TypeError, match='either connectivity or interfaces'):
        Domain.join(
            [patch_a, patch_b], connectivity=[], interfaces=[], name='AB')


#==============================================================================
def test_interface_can_be_reconstructed_from_sympy_args():
    A = Cube('A')
    B = Cube('B')
    interface = Interface(
        'A|B', A.get_boundary(0, 1), B.get_boundary(1, -1),
        (-1, +1, -1))

    reconstructed = interface.func(*interface.args)
    copied = copy.copy(interface)

    assert reconstructed == interface
    assert reconstructed.orientation == interface.orientation
    assert copied.orientation == interface.orientation


#==============================================================================
def test_3d_interface_with_different_face_axes_and_shared_corners():
    A = Cube('A')
    B = Cube('B')
    domain = Domain.join(
        [A, B],
        [((0, 0, +1), (1, 1, -1), (-1, +1, -1))],
        'AB')

    interface = domain.interfaces
    assert interface.minus.axis == 0
    assert interface.plus.axis == 1
    assert interface.orientation == (-1, +1, -1)
    assert interface.axis_map == ((1, 2, 1), (2, 0, -1))

    corners = domain.corners.as_tuple()
    assert len(corners) == 4
    coordinate_pairs = {
        tuple(corner.coordinates for corner in shared_corner.corners)
        for shared_corner in corners
    }
    assert ((1, 0, 1), (0, 0, 0)) in coordinate_pairs


#==============================================================================
def test_mapping_multipatch_domain_preserves_interface_orientation():
    logical_domain = Domain.join(
        [Cube('A'), Cube('B')],
        [((0, 0, 1), (1, 1, -1), (-1, +1, -1))],
        'AB')
    mapping = Mapping('F', dim=3)
    mapped_domain = mapping(logical_domain)
    logical_interface = logical_domain.interfaces
    mapped_interface = mapped_domain.interfaces

    assert mapped_interface.orientation == logical_interface.orientation
    assert mapped_interface.axis_map == logical_interface.axis_map
    assert isinstance(mapped_interface.mapping, InterfaceMapping)
    assert mapped_interface.logical_domain is logical_interface
    assert mapped_interface.minus.logical_domain == logical_interface.minus
    assert mapped_interface.plus.logical_domain == logical_interface.plus
    assert mapping(logical_interface) is mapped_interface


#==============================================================================
def test_generated_interface_names_do_not_collide_with_patch_names():
    patch_a = Square('A')
    patch_b = Square('B')
    patch_b2 = Square('B#2')
    domain = Domain.join(
        [patch_a, patch_b, patch_b2],
        [
            ((patch_a, 0, -1), (patch_b, 0, -1), +1),
            ((patch_a, 0, +1), (patch_b, 0, +1), +1),
            ((patch_a, 1, -1), (patch_b2, 1, -1), +1),
        ],
        'collision-safe')

    assert tuple(domain.interface_map) == (
        'A|B', 'A|B#2', 'A|B#2#2')
    assert len(domain.interface_map) == 3
    assert domain.interface_map['A|B#2'].plus.domain == patch_b.interior
    assert domain.interface_map['A|B#2#2'].plus.domain == patch_b2.interior


#==============================================================================
def test_explicit_interface_names_are_reserved_and_unique():
    patch_a = Square('A')
    patch_b = Square('B')
    domain = Domain.join(
        [patch_a, patch_b],
        [
            ((0, 0, -1), (1, 0, -1), +1),
            ((0, 0, +1), (1, 0, +1), +1, 'A|B'),
        ],
        'named')

    assert tuple(domain.interface_map) == ('A|B#2', 'A|B')
    assert str(domain.interface_map['A|B'].name) == 'A|B'

    with pytest.raises(ValueError, match='duplicate explicit interface name'):
        Domain.join(
            [patch_a, patch_b],
            [
                ((0, 0, -1), (1, 0, -1), +1, 'seam'),
                ((0, 0, +1), (1, 0, +1), +1, 'seam'),
            ],
            'duplicate')

    with pytest.raises(TypeError, match='name must be a string'):
        Domain.join(
            [patch_a, patch_b],
            [((0, 0, -1), (1, 0, -1), +1, 1)],
            'invalid-name')


#==============================================================================
def test_interface_requires_orientation():
    with pytest.raises(TypeError, match='minus, plus, orientation'):
        Domain.join(
            [Square('A'), Square('B')],
            [((0, 0, 1), (1, 0, -1))],
            'AB')

    with pytest.raises(TypeError, match='minus, plus, orientation'):
        Domain.join(
            [Square('C'), Square('D')],
            [{'minus': (0, 0, 1),
              'plus': (1, 0, -1),
              'orientation': +1}],
            'CD')

    with pytest.raises(ValueError, match='2D'):
        Domain.join(
            [Square('A'), Square('B')],
            [((0, 0, 1), (1, 0, -1), 0)],
            'AB')


#==============================================================================
def test_single_patch_self_interface():
    patch = Square('A')
    domain = Domain.join(
        [patch],
        [((0, 0, -1), (0, 0, 1), +1)],
        'periodic-A')

    interface = domain.interfaces
    assert interface.minus.domain == patch.interior
    assert interface.plus.domain == patch.interior
    assert interface.minus.axis == interface.plus.axis == 0
    assert interface.orientation == +1
    assert len(domain.corners) == 2


#==============================================================================
def test_mapped_single_patch_self_interface_uses_single_mapping():
    logical_patch = Square('A')
    mapping = IdentityMapping('F', dim=2)
    patch = mapping(logical_patch)
    domain = Domain.join(
        [patch],
        [((0, 0, -1), (0, 0, 1), +1)],
        'periodic-A')

    assert domain.mapping == mapping
    assert domain.logical_domain.interfaces.orientation == +1
    assert domain.interfaces.minus.domain == patch.interior
    assert domain.interfaces.plus.domain == patch.interior


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
