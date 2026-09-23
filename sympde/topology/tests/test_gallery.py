import pytest
import h5py
import numpy as np
import yaml

from sympy.abc import x,y,z
from sympy import Tuple
from sympy import symbols

from sympde.topology import Interval, ProductDomain, InteriorDomain, Domain
from sympde.topology import Line, Square, Cube, NCubeInterior
from sympde.topology import IdentityMapping


x1, x2, x3 = symbols('x1, x2, x3', real=True)

#==============================================================================
def test_interval():
    Ix = Interval('Ix', coordinate=x)
    Iy = Interval('Iy', coordinate=y)
    Iz = Interval('Iz', coordinate=z)

    for I,i in zip([Ix, Iy, Iz], [x, y, z]):
        assert(I.coordinates == i)

    D_xy = ProductDomain(Ix, Iy)
    assert(D_xy.dim == 2)

#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_unit_line():

    # Create 1D domain (Line) from interval [0, 1]
    domain = Line('line')

    assert isinstance(domain, Line)

    # BasicDomain's attributes
    assert domain.dim  == 1
    assert domain.name == 'line'
    assert domain.coordinates == x1

    # Domain's attributes
    assert isinstance(domain.interior, NCubeInterior)
    assert len(domain.boundary) == 2
    assert domain.dtype == {'type': 'Line',
                            'parameters': {'bounds': [0, 1]}}

    # NCube's attributes
    assert domain.min_coords == (0,)
    assert domain.max_coords == (1,)

    # Line's attributes
    assert domain.bounds == (0, 1)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain

#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_generic_line():

    # Create 1D domain (Line) from interval [-3, 4]
    domain = Line('line', bounds=(-3, 4))

    assert isinstance(domain, Line)

    # BasicDomain's attributes
    assert domain.dim  == 1
    assert domain.name == 'line'
    assert domain.coordinates == x1

    # Domain's attributes
    assert isinstance(domain.interior, NCubeInterior)
    assert len(domain.boundary) == 2
    assert domain.dtype == {'type': 'Line',
                            'parameters': {'bounds': [-3, 4]}}

    # NCube's attributes
    assert domain.min_coords == (-3,)
    assert domain.max_coords == ( 4,)

    # Line's attributes
    assert domain.bounds == (-3, 4)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain

#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_unit_square():

    # Create 2D square domain [0, 1]^2
    domain = Square('square')

    assert isinstance(domain, Square)

    # BasicDomain's attributes
    assert domain.dim  == 2
    assert domain.name == 'square'
    assert domain.coordinates == (x1, x2)

    # Domain's attributes
    assert isinstance(domain.interior, InteriorDomain)

    assert len(domain.boundary) == 4
    assert domain.dtype == {'type': 'Square',
                            'parameters': {'bounds1': [0, 1],
                                           'bounds2': [0, 1]}}

    # NCube's attributes
    assert domain.min_coords == (0, 0)
    assert domain.max_coords == (1, 1)

    # Square's attributes
    assert domain.bounds1 == (0, 1)
    assert domain.bounds2 == (0, 1)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain

#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_rectangle():

    # Create 2D rectangular domain [1, 5] X [3, 7]
    domain = Square('rectangle', bounds1=(1, 5), bounds2=(3, 7))

    assert isinstance(domain, Square)

    # BasicDomain's attributes
    assert domain.dim  == 2
    assert domain.name == 'rectangle'
    assert domain.coordinates == (x1, x2)

    # Domain's attributes
    assert isinstance(domain.interior, InteriorDomain)
    assert len(domain.boundary) == 4
    assert domain.dtype == {'type': 'Square',
                            'parameters': {'bounds1': [1, 5],
                                           'bounds2': [3, 7]}}

    # NCube's attributes
    assert domain.min_coords == (1, 3)
    assert domain.max_coords == (5, 7)

    # Square's attributes
    assert domain.bounds1 == (1, 5)
    assert domain.bounds2 == (3, 7)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain

#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_unit_cube():

    # Create 3D cube domain [0, 1]^3
    domain = Cube('cube')

    assert isinstance(domain, Cube)

    # Check object attributes
    assert domain.dim  == 3
    assert domain.name == 'cube'
    assert domain.coordinates == (x1, x2, x3)

    # Domain's attributes
    assert isinstance(domain.interior, InteriorDomain)
    assert len(domain.boundary) == 6
    assert domain.dtype == {'type': 'Cube',
                            'parameters': {'bounds1': [0, 1],
                                           'bounds2': [0, 1],
                                           'bounds3': [0, 1]}}

    # NCube's attributes
    assert domain.min_coords == (0, 0, 0)
    assert domain.max_coords == (1, 1, 1)

    # Cube's attributes
    assert domain.bounds1 == (0, 1)
    assert domain.bounds2 == (0, 1)
    assert domain.bounds3 == (0, 1)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain

#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_orthogonal_hexahedron():

    # Create 3D orthogonal hexahedron [-1, 1] X [0, 10] X [0, 2]
    domain = Cube('hexahedron', bounds1=(-1, 1), bounds2=(0, 10), bounds3=(0, 2))

    assert isinstance(domain, Cube)

    # Check object attributes
    assert domain.dim  == 3
    assert domain.name == 'hexahedron'
    assert domain.coordinates == (x1, x2, x3)

    # Domain's attributes
    assert isinstance(domain.interior, InteriorDomain)
    assert len(domain.boundary) == 6
    assert domain.dtype == {'type': 'Cube',
                            'parameters': {'bounds1': [-1, 1],
                                           'bounds2': [0, 10],
                                           'bounds3': [0,  2]}}

    # NCube's attributes
    assert domain.min_coords == (-1, 0, 0)
    assert domain.max_coords == (1, 10, 2)

    # Cube's attributes
    assert domain.bounds1 == (-1, 1)
    assert domain.bounds2 == (0, 10)
    assert domain.bounds3 == (0,  2)

    # Export to file, read it again and compare with original domain
    domain.export('domain.h5')
    D = Domain.from_file('domain.h5')
    assert D == domain


#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_multipatch_orientation_roundtrip(tmp_path):
    domain = Domain.join(
        [Cube('A'), Cube('B')],
        [((0, 0, 1), (1, 1, -1), (-1, +1, -1))],
        'AB')
    filename = tmp_path / 'multipatch.h5'

    domain.export(str(filename))
    restored = Domain.from_file(str(filename))

    assert restored.interfaces.orientation == (-1, +1, -1)
    assert restored.interfaces.minus.axis == 0
    assert restored.interfaces.plus.axis == 1
    assert restored.interfaces.axis_map == ((1, 2, 1), (2, 0, -1))


#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_multipatch_file_requires_explicit_orientation(tmp_path):
    domain = Domain.join(
        [Square('A'), Square('B')],
        [((0, 0, 1), (1, 0, -1), +1)],
        'AB')
    topology = domain.todict()
    for interface_data in topology['connectivity'].values():
        interface_data.pop()

    filename = tmp_path / 'missing_orientation.h5'
    with h5py.File(filename, mode='w') as h5:
        h5['topology.yml'] = np.array(
            yaml.safe_dump(topology, sort_keys=None), dtype='S')

    with pytest.raises(ValueError, match='explicit orientation'):
        Domain.from_file(str(filename))


#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_multipatch_file_rejects_legacy_orientation_mapping(tmp_path):
    domain = Domain.join(
        [Square('A'), Square('B')],
        [((0, 0, 1), (1, 0, -1), +1)],
        'AB')
    topology = domain.todict()
    for interface_data in topology['connectivity'].values():
        interface_data[2] = {
            'permutation': [0],
            'directions': [1],
        }

    filename = tmp_path / 'legacy_orientation.h5'
    with h5py.File(filename, mode='w') as h5:
        h5['topology.yml'] = np.array(
            yaml.safe_dump(topology, sort_keys=None), dtype='S')

    with pytest.raises(TypeError, match='2D interface orientation'):
        Domain.from_file(str(filename))


#==============================================================================
@pytest.mark.xdist_group('h5py')
@pytest.mark.parametrize('location, message', [
    ('boundary-axis', 'serialized boundary axis must be an integer'),
    ('interface-extremity',
     'serialized interface extremity must be an integer'),
])
def test_multipatch_file_rejects_non_integer_side_metadata(
        tmp_path, location, message):
    domain = Domain.join(
        [Square('A'), Square('B')],
        [((0, 0, +1), (1, 0, -1), +1)],
        'AB')
    topology = domain.todict()

    if location == 'boundary-axis':
        topology['boundary'][0]['axis'] = 0.0
    else:
        interface_data = next(iter(topology['connectivity'].values()))
        interface_data[0]['ext'] = 1.0

    filename = tmp_path / f'{location}.h5'
    with h5py.File(filename, mode='w') as h5:
        h5['topology.yml'] = np.array(
            yaml.safe_dump(topology, sort_keys=None), dtype='S')

    with pytest.raises(TypeError, match=message):
        Domain.from_file(str(filename))


#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_self_interface_roundtrip(tmp_path):
    patch = Square('A')
    domain = Domain.join(
        [patch],
        [((0, 0, -1), (0, 0, 1), +1)],
        'periodic-A')
    filename = tmp_path / 'self_interface.h5'

    domain.export(str(filename))
    restored = Domain.from_file(str(filename))

    assert restored.todict() == domain.todict()
    assert restored.interfaces.minus.domain == restored.interfaces.plus.domain


#==============================================================================
@pytest.mark.xdist_group('h5py')
@pytest.mark.parametrize('mapped', [False, True])
def test_boundary_free_cube_roundtrip(tmp_path, mapped):
    logical_patch = Cube('A')
    patch = (IdentityMapping('F', dim=3)(logical_patch)
             if mapped else logical_patch)
    domain = Domain.join(
        [patch],
        [
            ((0, 0, -1), (0, 0, +1), (+1, +1, +1)),
            ((0, 1, -1), (0, 1, +1), (+1, +1, +1)),
            ((0, 2, -1), (0, 2, +1), (+1, +1, +1)),
        ],
        'periodic-cube')
    filename = tmp_path / f'periodic_cube_{mapped}.h5'

    assert domain.boundary is None
    assert domain.todict()['boundary'] == []

    domain.export(str(filename))
    restored = Domain.from_file(str(filename))

    assert restored.boundary is None
    assert restored.todict() == domain.todict()
    assert tuple(restored.interface_map) == tuple(domain.interface_map)
    assert all((interface.mapping is not None) == mapped
               for interface in restored.interface_map.values())
    assert all((interface.logical_domain is not None) == mapped
               for interface in restored.interface_map.values())

    for legacy_value in ('missing', None):
        topology = domain.todict()
        if legacy_value == 'missing':
            topology.pop('boundary')
        else:
            topology['boundary'] = None
        legacy_filename = tmp_path / f'periodic_cube_{legacy_value}.h5'
        with h5py.File(legacy_filename, mode='w') as h5:
            h5['topology.yml'] = np.array(
                yaml.safe_dump(topology, sort_keys=None), dtype='S')
        legacy = Domain.from_file(str(legacy_filename))
        assert legacy.boundary is None


#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_single_boundary_and_legacy_dictionary_roundtrip(tmp_path):
    patch = Line('A')
    domain = Domain.join(
        [patch],
        [((0, 0, -1), (0, 0, -1), None, 'identified-left')],
        'half-open-line')
    canonical_filename = tmp_path / 'single_boundary.h5'
    legacy_filename = tmp_path / 'legacy_single_boundary.h5'
    topology = domain.todict()

    assert isinstance(topology['boundary'], list)
    assert len(topology['boundary']) == 1

    domain.export(str(canonical_filename))
    canonical = Domain.from_file(str(canonical_filename))
    assert canonical.todict() == domain.todict()

    topology['boundary'] = topology['boundary'][0]
    with h5py.File(legacy_filename, mode='w') as h5:
        h5['topology.yml'] = np.array(
            yaml.safe_dump(topology, sort_keys=None), dtype='S')

    legacy = Domain.from_file(str(legacy_filename))
    assert legacy.todict() == domain.todict()

    for legacy_value in ('missing', None):
        topology = domain.todict()
        if legacy_value == 'missing':
            topology.pop('boundary')
        else:
            topology['boundary'] = None
        incomplete_filename = tmp_path / f'single_boundary_{legacy_value}.h5'
        with h5py.File(incomplete_filename, mode='w') as h5:
            h5['topology.yml'] = np.array(
                yaml.safe_dump(topology, sort_keys=None), dtype='S')
        incomplete = Domain.from_file(str(incomplete_filename))
        assert incomplete.todict() == domain.todict()


#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_explicit_interface_name_roundtrip(tmp_path):
    patch_a = IdentityMapping('FA', dim=2)(Square('A'))
    patch_b = IdentityMapping('FB', dim=2)(Square('B'))
    domain = Domain.join(
        [patch_a, patch_b],
        [
            ((0, 0, -1), (1, 0, -1), +1, 'lower-seam'),
            ((0, 0, +1), (1, 0, +1), -1, 'upper-seam'),
        ],
        'named-seams')
    filename = tmp_path / 'named_interfaces.h5'

    domain.export(str(filename))
    restored = Domain.from_file(str(filename))

    assert set(restored.interface_map) == {'lower-seam', 'upper-seam'}
    assert restored.todict() == domain.todict()
    assert {
        str(interface.logical_domain.name)
        for interface in restored.interface_map.values()
    } == {'A|B', 'A|B#2'}


#==============================================================================
@pytest.mark.xdist_group('h5py')
def test_file_boundary_must_match_connectivity(tmp_path):
    domain = Domain.join(
        [Square('A'), Square('B')],
        [((0, 0, +1), (1, 0, -1), +1)],
        'AB')
    topology = domain.todict()
    topology['boundary'].pop()
    filename = tmp_path / 'inconsistent_boundary.h5'

    with h5py.File(filename, mode='w') as h5:
        h5['topology.yml'] = np.array(
            yaml.safe_dump(topology, sort_keys=None), dtype='S')

    with pytest.raises(ValueError, match='does not match'):
        Domain.from_file(str(filename))

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

    # Remove output file generated by tests
    import os
    fname = 'domain.h5'
    if os.path.exists(fname):
        os.remove(fname)

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
