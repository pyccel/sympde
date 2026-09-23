# coding: utf-8
from __future__ import annotations

import numpy as np
import h5py
import yaml
import os

from collections import OrderedDict
from itertools import product
from typing import Union as TypeUnion, Optional, List, Dict, Iterable, TYPE_CHECKING
# Union clashes with core.basic.Union

from sympy import Integer
from sympy.core.singleton import Singleton
from sympy.core import Basic, symbols
from sympy.core.containers import Tuple
from sympy.tensor import IndexedBase, Indexed
from sympy.core import Add, Mul, Pow
from sympy.core.expr import AtomicExpr

from sympde.old_sympy_utilities import is_sequence, with_metaclass
from sympde.core.basic import CalculusFunction
from .basic            import BasicDomain, InteriorDomain, Boundary, Union, Connectivity
from .basic            import Interval, Interface
from .basic            import CornerBoundary, CornerInterface
from .basic            import ProductDomain
from .basic            import _as_integer

# TODO fix circular dependency between domain and mapping
if TYPE_CHECKING:
    from sympde.topology.mapping import Mapping
# TODO add pdim

iterable_types = (tuple, list, Tuple, Union)

#==============================================================================
class Domain(BasicDomain):
    """
    Represents an undefined domain.
    A domain is defined by at least one interior domain and possible boundaries.
    A domain without a boundary is either infinite or periodic.
    A domain can also be constructed from a connectivity, in which case, only the
    name and connectivity need to be passed.
    """

    def __new__(cls, name : str, *,
            interiors : TypeUnion[Iterable[InteriorDomain], InteriorDomain, None] = None,
            boundaries : TypeUnion[Iterable[Boundary], Boundary, None] = None,
            dim : Optional[int] = None,
            connectivity : Optional[Connectivity] = None,
            mapping : Optional[Mapping] = None,
            logical_domain : Optional[Domain] = None):
        """
        Interiors or connectivity must be given. When the mapping is given 
        then logical_domain must be specified as well.
    
        Parameters
        ----------
        name : str
            Name of the domain
        interiors : Iterable[InteriorDomain], InteriorDomain or None, optional
            Interior domains
        boundaries : Iterable[Boundary], Boundary or None, optional
            The boundaries of the domain
        dim : int, optional
            Dimension of the space of the domain
        connectivity : Connectivity or None, optional
            Connectivity object with the interfaces of the domain
        mapping : Mapping or None, optional
            Maps the logical domain to the physical domain
        logical_domain : Domain or None, optional
            Logical domain that is mapped to the physical domain
        """
        # ...
        if not isinstance(name, str):
            raise TypeError('> name must be a string')
        # ...

        # ...
        if ( ( interiors is None ) and ( connectivity is None ) and ( dim is None) ):
            raise ValueError('> either interiors or connectivity must be given')
        # ...

        # ...
        if not( interiors is None ):
            if not isinstance( interiors, (*iterable_types, InteriorDomain)):
                raise TypeError('> Expecting an iterable or a InteriorDomain')

            if isinstance( interiors, InteriorDomain ):
                interiors = [interiors]

            else:
                new_interiors = []
                for i in interiors:
                    if isinstance(i , iterable_types):
                        new_interiors += list(i)
                    else:
                        new_interiors.append(i)

                interiors = new_interiors

                if not all([isinstance(i, InteriorDomain) for i in interiors]):
                    raise TypeError('> all interiors must be of type InteriorDomain')

            interiors = Tuple(*interiors)
        # ...

        if not( boundaries is None ):
            if not isinstance( boundaries, (*iterable_types, Boundary)):
                raise TypeError('> Expecting an iterable or a Boundary')

            if isinstance( boundaries, Boundary ):
                boundaries = [boundaries]

            else:
                if not all([isinstance(i, Boundary) for i in boundaries]):
                    raise TypeError('> all boundaries must be of type Boundary')

        else:
            boundaries = []

        boundaries = Tuple(*boundaries)

        if not( connectivity is None ):
            if not isinstance( connectivity, Connectivity ):
                raise TypeError('> Expecting a Connectivity')

            # TODO check that patches appearing in connectivity are in interiors
        else:
            connectivity = Connectivity()

        # ...
        if interiors is None and dim:
            interiors = [InteriorDomain(name, dim=dim)]

        if len(interiors) == 0 and dim is None:
            raise TypeError('No interior domain found')

        elif len(interiors) == 1:
            dtype = interiors[0].dtype
            dim   = interiors[0].dim
            interiors = Union(*interiors)
        else:
            dim   = interiors[0].dim
            interiors = Union(*interiors)
            dtype = [i.dtype for i in interiors]

        assert mapping is None and logical_domain is None or \
        mapping is not None and logical_domain  is not None

        # ...
        boundaries = Union(*boundaries)

        obj = Basic.__new__(cls, name, interiors, boundaries, mapping)
        obj._connectivity   = connectivity
        obj._corners        = None
        obj._dtype          = dtype
        obj._dim            = dim
        obj._logical_domain = logical_domain
        return obj

    @property
    def name(self) -> str:
        return self.args[0]

    @property
    def interior(self) -> TypeUnion[Union, InteriorDomain]:
        """Either a Union object containing the interiors or just the interior 
        domain if there is only one"""
        return self.args[1]

    @property
    def boundary(self) -> TypeUnion[Union, Boundary]:
        """Either a Union object containing the boundaries or just a boundary 
        if there is only one"""
        return self.args[2]

    @property
    def mapping(self) -> Optional[Mapping]:
        """The mapping that maps the logical domain to the physical domain"""
        return self.args[3]

    @property
    def subdomains(self) -> tuple:
        """Return top-dimensional interior regions as a tuple."""
        if isinstance( self.interior, iterable_types):
            subs = self.interior
        else:
            subs = [self.interior]
        return tuple(subs)

    @property
    def patches(self) -> tuple:
        """Top-dimensional patch regions as a stable tuple.

        This additive multipatch view leaves the symbolic ``interior`` and
        ``subdomains`` APIs unchanged.  For a single-patch domain it contains
        exactly that domain's interior.
        """
        if isinstance(self.interior, InteriorDomain):
            return (self.interior,)
        return self.interior.as_tuple()

    @property
    def mappings(self) -> OrderedDict:
        return OrderedDict([(P.logical_domain, P.mapping)
                           for P in self.subdomains])

    @property
    def logical_domain(self) -> Domain:
        """The domain is the image of the logical_domain under the mapping"""
        return self._logical_domain

    @property
    def connectivity(self) -> Connectivity:
        """Contains information about the interfaces"""
        return self._connectivity

    @property
    def interface_map(self) -> Connectivity:
        """Interfaces keyed by their unique names.

        Unlike the symbolic :attr:`interfaces` property, this view always has
        mapping semantics, including when the domain has zero or one interface.
        """
        return self.connectivity

    @property
    def exterior_sides(self) -> tuple:
        """Exterior boundary sides as a stable tuple."""
        boundary = self.boundary
        if boundary is None:
            return ()
        if isinstance(boundary, Boundary):
            return (boundary,)
        return boundary.as_tuple()

    @property
    def dim(self) -> int:
        """Dimension of the space"""
        return self._dim

    @property
    def dtype(self) -> dict:
        """Dictionary containing information about domain"""
        return self._dtype

    @property
    def interfaces(self) -> TypeUnion[Union, Interface, None]:
        """
        Union of the interfaces

        The Union constructor is applied to the interfaces. If there is only 
        one interface it returns the interface object and None if there is no 
        interface.
        """
        return self.connectivity.interfaces

    @property
    def corners(self):
        corners = getattr(self,'_corners', None)
        if corners is None:
            corners = self.get_shared_corners()
        self._corners = corners
        return corners

    @property
    def shared_vertices(self) -> tuple:
        """Shared patch-vertex equivalence classes as a stable tuple."""
        corners = self.corners
        if corners is None:
            return ()
        if isinstance(corners, CornerInterface):
            return (corners,)
        return corners.as_tuple()

    def __len__(self):
        if isinstance(self.interior, InteriorDomain):
            return 1

        elif isinstance(self.interior, Union):
            return len(self.interior)

    @property
    def interior_names(self) -> List[str]:
        if isinstance(self.interior, InteriorDomain):
            return [self.interior.name]

        elif isinstance(self.interior, Union):
            return [i.name for i in self.interior.args]

    def set_interfaces(self, *interfaces):
        for i in interfaces:
            self.connectivity[i.name] = i

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

    def get_boundary(self, axis, ext):
        """
        Return the domain boundary at the given extremity of the required axis.

        Parameters
        ----------
        axis : int | None
            Index of the coordinate (0 <= axis < ndim) which has constant value at the boundary.
            In 1D passing `axis=None` is accepted, in which case it is interpreted as 0.
        ext : {-1, +1}
            Extremity identifier:
              * If -1, the boundary is at the minimum value of $x_{axis}$
              * If +1, the boundary is at the maximum value of $x_{axis}$
        
        Returns
        -------
        Boundary (from sympde.topology.basic)
            The domain boundary of interest.
        """
        if axis is None:
            if self.interior.dim != 1:
                raise ValueError('axis may be None only for a 1D domain')
            axis = 0
        else:
            axis = _as_integer(axis, 'boundary axis')
        ext = _as_integer(ext, 'boundary extremity')

        if not 0 <= axis < self.dim:
            raise ValueError(
                f'boundary axis must be between 0 and {self.dim - 1}')
        if ext not in (-1, 1):
            raise ValueError('boundary extremity must be either -1 or 1')

        if isinstance(self.boundary, Union):
            x = [i for i in self.boundary.args if i.ext == ext and i.axis == axis]
            if len(x) == 0:
                raise ValueError(f'> could not find boundary with axis {axis} and ext {ext}')
            return x[0]

        elif isinstance(self.boundary, Boundary):
            if self.boundary.axis == axis and self.boundary.ext == ext:
                return self.boundary

        raise ValueError(f'> could not find boundary with axis {axis} and ext {ext}')

    def get_interface(self, domain1, domain2):
        interfaces = []
        for i in self.interface:
            if i.plus in [domain1, domain2]:
                if i.minus in [domain1, domain2]:
                    interfaces.append(i)
        if interfaces:
            return Union(*interfaces)

        raise ValueError('> could not find the interface of {} and {}'.format(domain1, domain2))

    def get_interior(self, name):
        """return interior by name."""
        if isinstance(self.interior, Union):
            x = [i for i in self.interior.args if i.name == name]
            if len(x) == 0:
                raise ValueError('> could not find interior {}'.format(name))

            return x[0]

        elif isinstance(self.interior, InteriorDomain):
            if self.interior.name == name:
                return self.interior

            else:
                return None

    def todict(self):
        name         = str(self.name)
        dim          = str(self.dim)
        interior     = self.interior.todict()
        boundary     = [side.todict() for side in self.exterior_sides]
        connectivity = self.connectivity.todict()

        dtype = self.dtype
        if dtype is None:
            dtype = 'None'

        d = {'name':         name,
             'dim':          dim,
             'dtype':        dtype,
             'interior':     interior,
             'boundary':     boundary,
             'connectivity': connectivity}

        return dict(sorted(d.items()))

    def export(self, filename):

        yml = self.todict()

        # Dump metadata to string in YAML file format
        geo = yaml.safe_dump(data=yml, sort_keys=None)

        # Create HDF5 file (in parallel mode if MPI communicator size > 1)
        with h5py.File(filename, mode='w') as h5:
            # Write geometry metadata as fixed-length array of ASCII characters
            h5['topology.yml'] = np.array(geo, dtype='S')

    @classmethod
    def from_file(cls, filename):
        """
        Read the "topology.yml" portion of an HDF5 geometry file and create a (mapped)
        multipatch domain using the information therein.

        Parameters
        ----------
        filename : str
            Name of the HDF5 geometry file to be read.

        Returns
        -------
        Domain
            Multipatch domain.
        """
        # ... check extension of the file
        _, ext = os.path.splitext(filename)

        if ext != '.h5':
            raise ValueError('> Only h5 files are supported')
        # ...
        from sympde.topology.mapping import Mapping

        with h5py.File(filename, mode='r') as h5:
            yml = yaml.safe_load(h5['topology.yml'][()])

        domain_name    = yml['name']
        dim            = int(yml['dim'])
        dtype          = yml['dtype']
        d_interior     = yml['interior']
        has_boundary_metadata = (
            'boundary' in yml and yml['boundary'] is not None)
        d_boundary     = yml.get('boundary', [])
        d_connectivity = yml.get('connectivity', {})

        # Boundary metadata was historically cardinality-dependent: a single
        # side was stored as a dictionary and several sides as a list. Accept
        # that representation while making the canonical format always a list.
        if d_boundary is None:
            d_boundary = []
        elif isinstance(d_boundary, dict):
            d_boundary = [d_boundary]
        elif not isinstance(d_boundary, list):
            raise TypeError('boundary topology must be a list or a dictionary')

        if d_connectivity is None:
            d_connectivity = {}
        elif not isinstance(d_connectivity, dict):
            raise TypeError('connectivity topology must be a dictionary')

        if dtype == 'None':
            dtype = None

        assert dtype is not None
        assert all(dtype)
        if isinstance(d_interior, dict):
            d_interior = [d_interior]
            dtype      = [dtype]

        constructors = [globals()[dt['type']] for dt in dtype]
        interiors    = [cs(i['name'], **dt['parameters']) for cs,i,dt in zip(constructors, d_interior, dtype)]
        mappings     = [Mapping(I['mapping'], dim=dim) if I.get('mapping', "None") != "None" else None for I in d_interior]
        domains      = [mapping(i) if mapping else i for i,mapping in zip(interiors, mappings)]
        patch_index  = {I.name:ind for ind,I in enumerate(interiors)}

        def serialized_integer(value, name):
            """Read legacy decimal strings without accepting numeric floats."""
            if isinstance(value, str):
                try:
                    return int(value)
                except ValueError as error:
                    raise TypeError(f'{name} must be an integer') from error
            return _as_integer(value, name)

        stored_boundary_sides = []
        for boundary_data in d_boundary:
            if not isinstance(boundary_data, dict):
                raise TypeError('each serialized boundary must be a dictionary')
            name = boundary_data['patch']
            axis = serialized_integer(
                boundary_data['axis'], 'serialized boundary axis')
            ext  = serialized_integer(
                boundary_data['ext'], 'serialized boundary extremity')
            if name not in patch_index:
                raise ValueError(
                    f'serialized boundary references unknown patch: {name}')
            domains[patch_index[name]].get_boundary(axis=axis, ext=ext)
            stored_boundary_sides.append((str(name), axis, ext))

        if len(stored_boundary_sides) != len(set(stored_boundary_sides)):
            raise ValueError('serialized boundary contains duplicate sides')

        connectivity = []
        for interface_name, interface_data in d_connectivity.items():
            if not isinstance(interface_name, str):
                raise TypeError('serialized interface names must be strings')
            if len(interface_data) != 3:
                raise ValueError(
                    'multipatch topology files must store two interface sides '
                    'and an explicit orientation; regenerate legacy files '
                    'that contain only two interface sides')
            minus, plus, orientation = interface_data

            minus_name = minus['patch']
            minus_axis = serialized_integer(
                minus['axis'], 'serialized interface axis')
            minus_ext  = serialized_integer(
                minus['ext'], 'serialized interface extremity')
            minus_patch_i = patch_index[minus_name]

            plus_name = plus['patch']
            plus_axis = serialized_integer(
                plus['axis'], 'serialized interface axis')
            plus_ext  = serialized_integer(
                plus['ext'], 'serialized interface extremity')
            plus_patch_i = patch_index[plus_name]
            interface = ((minus_patch_i, minus_axis, minus_ext),
                         (plus_patch_i, plus_axis, plus_ext), orientation,
                         interface_name)

            connectivity.append(interface)

        if len(domains) == 1 and not connectivity:
            domain = domains[0]
        else:
            domain = Domain.join(domains, connectivity, domain_name)

        derived_boundary_sides = []
        for boundary in domain.exterior_sides:
            boundary_data = boundary.todict()
            derived_boundary_sides.append((
                str(boundary_data['patch']),
                int(boundary_data['axis']),
                int(boundary_data['ext'])))

        if has_boundary_metadata and \
           set(stored_boundary_sides) != set(derived_boundary_sides):
            raise ValueError(
                'serialized boundary does not match the exterior boundary '
                'derived from connectivity')

        return domain

    @classmethod
    def join(cls, patches, connectivity=None, name=None, *, interfaces=None):
        """
        Create a domain by joining one or more patches in 1D, 2D, or 3D.

        A single patch is returned unchanged when ``connectivity`` is empty.
        With non-empty connectivity, its distinct boundaries may be joined to
        create self-interfaces such as periodic identifications.

        Parameters
        ----------
        patches : sequence of Domain
            Ordered non-empty collection of atomic, unconnected patches. All
            patches must have the same logical dimension. Their positions
            define the integer patch indices accepted by an interface side.
            Patch objects, interior names, and serialized logical names must be
            unique. Joined multipatch domains must be flattened before being
            passed to this method.

        connectivity : sequence, optional
            Interface descriptions. Each description is the tuple
            ``(minus, plus, orientation)`` or
            ``(minus, plus, orientation, name)``. The optional name is
            preserved exactly and must be unique among all interfaces.
            If both this argument and ``interfaces`` are omitted, the
            connectivity is empty.

            Each side is ``(patch, axis, ext)``. ``patch`` is either a patch
            object or its integer position in ``patches``. ``axis`` is the
            zero-based logical axis normal to the face. ``ext=-1`` selects the
            lower-coordinate face and ``ext=+1`` the upper-coordinate face.

            ``minus`` and ``plus`` define the direction in which orientation
            is read. In 1D, orientation is ``None``. In 2D, it is ``+1`` or
            ``-1``. In 3D, it is ``(flag, sign1, sign2)`` as documented by
            :class:`Interface`.

        interfaces : sequence, optional
            Keyword alias for ``connectivity``.  Provide exactly one of
            ``connectivity`` and ``interfaces``.

        name : str
            Name of the domain.

        Returns
        -------
        Domain
            Multipatch domain.

        Notes
        -----
        The compact orientation convention follows the convention used by
        GeoPDEs; see its `multipatch geometry specification
        <https://github.com/rafavzqz/geopdes/blob/master/geopdes/doc/geo_specs_mp_v21.txt#L193-L237>`_
        and T. Dokken, E. Quak, V. Skytt, *Requirements from Isogeometric
        Analysis for Changes in Product Design Ontologies* (2010).

        Connectivity is symbolic. This method does not infer interfaces from
        physical coordinates or verify that mapped faces coincide. It also
        does not impose a mesh, spline space, or trace-coupling policy.

        Examples
        --------
        A four-patch 2D connectivity may reference patch objects directly:

        >>> from sympde.topology import Cube, Domain, Square
        >>> patches = [Square(f'P{i}') for i in range(4)]
        >>> connectivity = [
        ...     ((patches[0], 0, -1), (patches[1], 0, +1), +1),
        ...     ((patches[1], 1, -1), (patches[3], 1, +1), -1),
        ...     ((patches[0], 1, -1), (patches[2], 1, +1), +1),
        ...     ((patches[2], 0, -1), (patches[3], 0, +1), -1),
        ... ]
        >>> omega = Domain.join(patches, connectivity, name='Omega')

        The equivalent connectivity may use positions in ``patches``:

        >>> connectivity = [
        ...     ((0, 0, -1), (1, 0, +1), +1),
        ...     ((1, 1, -1), (3, 1, +1), -1),
        ...     ((0, 1, -1), (2, 1, +1), +1),
        ...     ((2, 0, -1), (3, 0, +1), -1),
        ... ]
        >>> omega = Domain.join(patches, connectivity, name='Omega')

        A tuple also describes a cross-axis 3D interface:

        >>> A = Cube('A')
        >>> B = Cube('B')
        >>> omega = Domain.join(
        ...     patches=[A, B],
        ...     interfaces=[(
        ...         (A, 0, +1),
        ...         (B, 1, -1),
        ...         (-1, +1, -1),
        ...     )],
        ...     name='Omega3D',
        ... )
        >>> interface, = omega.interface_map.values()
        >>> interface.axis_map
        ((1, 2, 1), (2, 0, -1))
        """
        if not isinstance(patches, (tuple, list)):
            raise TypeError('patches must be a list or tuple')
        if not patches or not all(isinstance(patch, Domain) for patch in patches):
            raise TypeError('patches must contain Domain objects')

        if len({id(patch) for patch in patches}) != len(patches):
            raise ValueError('patches must not contain the same object twice')

        if any(not isinstance(patch.interior, InteriorDomain) or
               patch.connectivity for patch in patches):
            raise ValueError(
                'Domain.join expects atomic patch domains without existing '
                'connectivity; flatten joined domains before joining them')

        patch_names = [str(patch.interior.name) for patch in patches]
        if len(set(patch_names)) != len(patch_names):
            raise ValueError('patch interior names must be unique')

        serialized_patch_names = [
            str(patch.interior.logical_domain.name)
            if patch.interior.logical_domain is not None
            else str(patch.interior.name)
            for patch in patches
        ]
        if len(set(serialized_patch_names)) != len(serialized_patch_names):
            raise ValueError(
                'logical patch names used for serialization must be unique')

        if interfaces is not None:
            if connectivity is not None:
                raise TypeError(
                    'provide either connectivity or interfaces, not both')
            connectivity = interfaces
        elif connectivity is None:
            connectivity = ()
        if not isinstance(connectivity, (tuple, list)):
            raise TypeError('connectivity/interfaces must be a list or tuple')
        if not isinstance(name, str):
            raise TypeError('name must be a string')

        if len(patches) == 1 and not connectivity:
            # Preserve the original single-patch domain when there is no
            # topology to add. A non-empty connectivity may contain a
            # self-interface and must follow the normal joining path.
            return patches[0]

        if not all(p.dim == patches[0].dim for p in patches):
            raise ValueError('all patches must have the same logical dimension')
        ldim = int(patches[0].dim)

        normalized_connectivity = []
        explicit_interface_names = set()
        for interface_data in connectivity:
            if not isinstance(interface_data, tuple) or \
               len(interface_data) not in (3, 4):
                raise TypeError(
                    'an interface must be the tuple '
                    '(minus, plus, orientation) or '
                    '(minus, plus, orientation, name)')

            minus_spec, plus_spec, orientation = interface_data[:3]
            interface_name = (
                interface_data[3] if len(interface_data) == 4 else None)
            if interface_name is not None:
                if not isinstance(interface_name, str):
                    raise TypeError('an explicit interface name must be a string')
                if not interface_name:
                    raise ValueError('an explicit interface name cannot be empty')
                if interface_name in explicit_interface_names:
                    raise ValueError(
                        f'duplicate explicit interface name: {interface_name}')
                explicit_interface_names.add(interface_name)

            normalized_connectivity.append(
                (minus_spec, plus_spec, orientation, interface_name))

        from sympde.topology.mapping import MultiPatchMapping
        # ... connectivity
        interfaces = {}
        boundaries = []
        # Reserve every explicit name before allocating generated ones. This
        # makes explicit identity independent of the interface input order.
        physical_interface_names = set(explicit_interface_names)
        logical_interface_names = set()

        def get_unique_interface_name(base, used_names):
            """Return an unused name while preserving the historical pattern."""
            if base not in used_names:
                used_names.add(base)
                return base

            occurrence = 2
            while f'{base}#{occurrence}' in used_names:
                occurrence += 1
            name = f'{base}#{occurrence}'
            used_names.add(name)
            return name

        def get_boundary(boundary_spec):
            if not isinstance(boundary_spec, tuple) or len(boundary_spec) != 3:
                raise TypeError(
                    'an interface side must be the tuple '
                    '(patch, axis, extremity)')
            patch_ref, axis, ext = boundary_spec
            if isinstance(patch_ref, bool):
                raise TypeError(
                    'a patch reference must be a patch or integer index')
            if isinstance(patch_ref, Domain):
                matching_patches = [
                    patch for patch in patches if patch_ref is patch]
                if not matching_patches:
                    raise ValueError(
                        'an interface references a patch not present in patches')
                patch = matching_patches[0]
            else:
                try:
                    patch_index = _as_integer(patch_ref, 'patch index')
                except TypeError as error:
                    raise TypeError(
                        'a patch reference must be a patch or integer index') \
                        from error
                if not 0 <= patch_index < len(patches):
                    raise IndexError(f'patch index {patch_index} is out of range')
                patch = patches[patch_index]
            axis = _as_integer(axis, 'interface side axis')
            ext = _as_integer(ext, 'interface side extremity')
            if not 0 <= axis < ldim:
                raise ValueError(f'axis must be between 0 and {ldim - 1}')
            if ext not in (-1, 1):
                raise ValueError('boundary extremity must be either -1 or 1')
            return patch.get_boundary(axis=axis, ext=ext)

        for minus_spec, plus_spec, orientation, explicit_name \
                in normalized_connectivity:

            bnd_minus = get_boundary(minus_spec)
            bnd_plus  = get_boundary(plus_spec)
            if explicit_name is None:
                physical_base = (
                    f'{bnd_minus.domain.name}|{bnd_plus.domain.name}')
                interface_name = get_unique_interface_name(
                    physical_base, physical_interface_names)
            else:
                interface_name = explicit_name

            logical_name = None
            if bnd_minus.logical_domain and bnd_plus.logical_domain:
                logical_base = (
                    f'{bnd_minus.logical_domain.domain.name}|'
                    f'{bnd_plus.logical_domain.domain.name}')
                logical_name = get_unique_interface_name(
                    logical_base, logical_interface_names)

            interface = bnd_minus.join(
                bnd_plus, orientation, name=interface_name,
                logical_name=logical_name)

            interfaces[str(interface.name)] = interface
            boundaries.append(bnd_minus)
            boundaries.append(bnd_plus)

        connectivity = Connectivity()
        for k,v in interfaces.items():
            connectivity[k] = v

        # ... boundary
        boundaries = Union(*[b for p in patches for b in p.boundary]).complement(Union(*boundaries))
        if boundaries is None:
            boundaries = ()
        elif isinstance(boundaries, Boundary):
            boundaries = (boundaries,)
        else:
            boundaries = boundaries.as_tuple()

        # ... interiors
        interior_patches = [p.interior for p in patches]
        interiors       = Union(*interior_patches)

        if all(e.mapping for e in interior_patches):
            logical_interiors    = Union(*[e.logical_domain for e in interior_patches])
            logical_boundaries   = [e.logical_domain for e in boundaries]
            logical_connectivity = Connectivity()
            for k,v in connectivity.items():
                logical_connectivity[v.logical_domain.name] = v.logical_domain

            patch_mappings = {
                e.logical_domain: e.mapping for e in interior_patches}
            mapping = (
                interior_patches[0].mapping if len(interior_patches) == 1
                else MultiPatchMapping(patch_mappings))
            logical_domain = Domain(name,
                            interiors=logical_interiors,
                            boundaries=logical_boundaries,
                            connectivity=logical_connectivity)
        else:
            mapping        = None
            logical_domain = None

        # ...
        return Domain(name,
                      interiors=interiors,
                      boundaries=boundaries,
                      connectivity=connectivity,
                      mapping=mapping,
                      logical_domain=logical_domain)

    def get_shared_corners(self):
        """Return equivalence classes of patch vertices joined by interfaces.

        Face vertices are related using each interface's signed tangential-axis
        permutation.  Connected components then handle any patch incidence in
        both 2D and 3D, without walking around interfaces in a prescribed order.
        """
        interfaces = self.interfaces
        if interfaces is None:
            return None
        interfaces = (interfaces,) if isinstance(interfaces, Interface) else tuple(interfaces)

        parent = {}

        def find(vertex):
            parent.setdefault(vertex, vertex)
            while parent[vertex] != vertex:
                parent[vertex] = parent[parent[vertex]]
                vertex = parent[vertex]
            return vertex

        def union(vertex_1, vertex_2):
            root_1 = find(vertex_1)
            root_2 = find(vertex_2)
            if root_1 != root_2:
                if root_2.sort_key() < root_1.sort_key():
                    root_1, root_2 = root_2, root_1
                parent[root_2] = root_1

        def make_vertex(face, extents):
            boundaries = [
                face.domain.get_boundary(axis=axis, ext=ext)
                for axis, ext in enumerate(extents)
            ]
            return CornerBoundary(*boundaries)

        for interface in interfaces:
            minus_axis = interface.minus.axis
            plus_axis  = interface.plus.axis

            for tangent_extents in product((-1, 1), repeat=self.dim - 1):
                minus_extents = [None] * self.dim
                plus_extents  = [None] * self.dim
                minus_extents[minus_axis] = interface.minus.ext
                plus_extents[plus_axis]   = interface.plus.ext
                for source_position, (source_axis, target_axis, direction) \
                        in enumerate(interface.axis_map):
                    source_ext = tangent_extents[source_position]
                    minus_extents[source_axis] = source_ext
                    plus_extents[target_axis]  = direction * source_ext

                union(
                    make_vertex(interface.minus, minus_extents),
                    make_vertex(interface.plus, plus_extents))

        groups = {}
        for vertex in parent:
            groups.setdefault(find(vertex), []).append(vertex)

        shared_corners = [
            CornerInterface(*vertices)
            for _, vertices in sorted(groups.items(), key=lambda item: item[0].sort_key())
            if len(vertices) > 1
        ]
        return Union(*shared_corners)

    def get_subdomain(self, names):
        """
        Return the subdomain induced by a selection of patch names.

        Parameters
        ----------
        names : tuple of str or str
            Names of the patches to retain. If a string is given, a one-patch
            subdomain is returned. If a tuple is given, its order determines
            the generated subdomain name.

        Notes
        -----
        Interfaces whose two sides belong to selected patches are retained,
        including repeated and self-interfaces. If exactly one side is
        selected, that side becomes part of the subdomain's exterior boundary.
        """
        if names == ():
            return None

        if isinstance(names, str):
            names = (names,)
        elif not isinstance(names, tuple):
            raise TypeError('names must be a string or a tuple of strings')

        if not all(isinstance(name, str) for name in names):
            raise TypeError('every subdomain name must be a string')
        if len(set(names)) != len(names):
            raise ValueError('subdomain names must be unique')

        interior_names = tuple(self.interior_names)
        unknown_names = [
            name for name in names
            if name not in interior_names and name != self.name
        ]
        if unknown_names:
            unknown = ', '.join(unknown_names)
            raise ValueError(f'unknown subdomain name(s): {unknown}')

        # Selecting the containing domain or every patch is an identity
        # operation, preserving its name and object identity.
        if self.name in names or set(names) == set(interior_names):
            return self

        selected_names = set(names)
        interior_dict = {interior.name: interior for interior in self.patches}
        interiors = [interior_dict[name] for name in names]

        # Existing exterior sides remain exterior when their patch is selected.
        boundaries = [
            boundary for boundary in self.exterior_sides
            if boundary.domain.name in selected_names
        ]
        interfaces = OrderedDict()

        # Connectivity induces both the retained interfaces and the new cut
        # boundary. Iterating the original mapping preserves interface keys and
        # distinguishes repeated interfaces between the same patch pair.
        for key, interface in self.connectivity.items():
            minus_selected = interface.minus.domain.name in selected_names
            plus_selected = interface.plus.domain.name in selected_names

            if minus_selected and plus_selected:
                interfaces[key] = interface
            elif minus_selected:
                boundaries.append(interface.minus)
            elif plus_selected:
                boundaries.append(interface.plus)

        connectivity = Connectivity(interfaces)
        name = '|'.join(names)
        mapping = None
        logical_domain = None

        if all(interior.mapping is not None for interior in interiors):
            logical_interiors = [
                interior.logical_domain for interior in interiors]
            logical_boundaries = [
                boundary.logical_domain for boundary in boundaries]
            logical_interfaces = OrderedDict()

            if any(interior is None for interior in logical_interiors) or \
               any(boundary is None for boundary in logical_boundaries):
                raise ValueError(
                    'mapped patch topology is missing logical-domain metadata')

            for interface in interfaces.values():
                logical_interface = interface.logical_domain
                if logical_interface is None:
                    raise ValueError(
                        'mapped interface is missing logical-domain metadata')
                logical_interfaces[str(logical_interface.name)] = \
                    logical_interface

            logical_domain = Domain(
                name,
                interiors=logical_interiors,
                boundaries=logical_boundaries,
                connectivity=Connectivity(logical_interfaces))

            if len(interiors) == 1:
                mapping = interiors[0].mapping
            else:
                from sympde.topology.mapping import MultiPatchMapping
                mapping = MultiPatchMapping({
                    interior.logical_domain: interior.mapping
                    for interior in interiors
                })

        return Domain(
            name,
            interiors=interiors,
            boundaries=boundaries,
            connectivity=connectivity,
            mapping=mapping,
            logical_domain=logical_domain)


#==============================================================================
class PeriodicDomain(BasicDomain):

    def __init__(self, domain, periods):

        assert isinstance(domain, Domain)
        self._domain = domain
        self._periods = tuple(periods)
        boundary_dict = domain.boundary.todict()

        names = []
        for bd in boundary_dict:
            if periods[int(bd['axis'])] == True:
                names += [bd['name']]

        boundary = [bd for bd in domain.boundary.args if bd.name not in names]

        if len(boundary)>1:
            self._boundary = Union(*boundary)
        else:
            self._boundary = None

    @property
    def domain(self):
        return self._domain

    @property
    def periods(self):
        return self._periods

    @property
    def boundary(self):
        return self._boundary

    @property
    def dim(self):
        return self.domain.dim

    @property
    def coordinates(self):
        return self.domain.coordinates

    def __hash__(self):
        return hash((self._domain, self._periods))


#==============================================================================
class NCubeInterior(InteriorDomain):

    def __new__(cls, name, dim=None, dtype=None, min_coords=None, max_coords=None,
                    mapping=None, logical_domain=None):

        obj = InteriorDomain.__new__(cls, name, dim=dim, dtype=dtype,
                    mapping=mapping, logical_domain=logical_domain)

        obj._min_coords = min_coords
        obj._max_coords = max_coords

        boundaries = []
        i = 1
        for axis in range(dim):
            for ext in [-1, 1]:
                bnd_name = r'\Gamma_{}'.format(i)
                bd_logical_domain = logical_domain
                if bd_logical_domain:
                    bd_logical_domain = bd_logical_domain.get_boundary(axis=axis, ext=ext)
                Gamma = Boundary(bnd_name, obj, axis=axis, ext=ext, mapping=mapping, logical_domain=bd_logical_domain)
                boundaries += [Gamma]
                i += 1
        obj._boundary   = Union(*boundaries)
        return obj

    @property
    def min_coords(self):
        return self._min_coords

    @property
    def max_coords(self):
        return self._max_coords

    @property
    def boundary(self):
        return self._boundary

    def __hash__(self):
        return hash((self.args, self.min_coords, self.max_coords))

    def get_boundary(self, axis=None, ext=None):
        """return boundary by (axis, ext)."""
        # ...
        if ext is None:
            raise ValueError('boundary extremity must be provided')

        if axis is None:
            if self.dim != 1:
                raise ValueError('axis may be None only for a 1D domain')
            axis = 0
        else:
            axis = _as_integer(axis, 'boundary axis')
        ext = _as_integer(ext, 'boundary extremity')

        if not 0 <= axis < self.dim:
            raise ValueError(
                f'boundary axis must be between 0 and {self.dim - 1}')
        if ext not in (-1, 1):
            raise ValueError('boundary extremity must be either -1 or 1')

        if isinstance(self.boundary, Union):
            x = [i for i in self.boundary.args if i.ext == ext and i.axis==axis]
            if x:return x[0]
        raise ValueError('> could not find boundary with axis {} and ext {}'.format(axis, ext))
#==============================================================================
# Ncube's properties (in addition to Domain's properties):
#   . min_coords (default value is tuple of zeros)
#   . max_coords (default value is tuple of ones)
#
class NCube(Domain):

    def __new__(cls, name, dim, min_coords, max_coords):

        assert isinstance(name, str)
        assert isinstance(dim, (int, Integer))
        assert isinstance(min_coords, iterable_types[:-1])
        assert isinstance(max_coords, iterable_types[:-1])

        if not name:
            raise ValueError("Name must be provided")

        if dim < 1:
            raise ValueError("Number of dimensions must be at least 1")

        if not (dim == len(min_coords) == len(max_coords)):
            raise ValueError("Input arguments must have 'dim' components")

        if not all(xmin < xmax for xmin, xmax in zip(min_coords, max_coords)):
            raise ValueError("Min coordinates must be smaller than max")

        coord_names = 'x1:{}'.format(dim + 1)

        coordinates = symbols(coord_names, real=True)

        # Choose which type to use:
        #   a) if dim <= 3, use Line, Square or Cube;
        #   b) if dim <= 4, use a generic 'NCube' type.
        #
        # Moreover, store all initialization parameters in a 'dtype' dictionary.
        # This dictionary will be written to file when exporting the geometry,
        # and it must contain all information necessary for building a new object
        # by calling the appropriate constructor:
        #
        #   cls = globals()[dtype['type']]
        #   domain = cls(name, **dtype['parameters'])
        #
        min_coords = tuple(float(i) for i in min_coords)
        max_coords = tuple(float(i) for i in max_coords)

        if dim == 1:
            cls = Line
            dtype = {'type': 'Line',
                     'parameters': {'bounds': [min_coords[0], max_coords[0]]}}
        elif dim == 2:
            cls = Square
            dtype = {'type': 'Square',
                     'parameters': {'bounds1': [min_coords[0], max_coords[0]],
                                    'bounds2': [min_coords[1], max_coords[1]]}}
        elif dim == 3:
            cls = Cube
            dtype = {'type': 'Cube',
                     'parameters': {'bounds1': [min_coords[0], max_coords[0]],
                                    'bounds2': [min_coords[1], max_coords[1]],
                                    'bounds3': [min_coords[2], max_coords[2]]}}
        else:
            dtype = {'type': 'NCube',
                     'parameters': {'dim'       : dim,
                                    'min_coords': [*min_coords],
                                    'max_coords': [*max_coords]}}

        interior = NCubeInterior(name, dim=dim, dtype=dtype, min_coords=tuple(min_coords), max_coords=tuple(max_coords))

        # Create instance of given type
        obj = super().__new__(cls, name, interiors=[interior], boundaries=interior.boundary)

        # Store attributes in object
        obj._coordinates = tuple(coordinates)

        # Return object
        return obj

    @classmethod
    def from_file(cls, filename):
        msg = "Class method 'from_file' must be called on 'Domain' base class"
        raise TypeError(msg)

    @property
    def min_coords(self):
        return self.interior.min_coords

    @property
    def max_coords(self):
        return self.interior.max_coords
#==============================================================================
class Line(NCube):

    def __new__(cls, name='Line', bounds=(0, 1)):
        dim = 1
        min_coords = (bounds[0],)
        max_coords = (bounds[1],)
        return super().__new__(cls, name, dim, min_coords, max_coords)

    @property
    def bounds(self):
        return (self.min_coords[0], self.max_coords[0])

#==============================================================================
class Square(NCube):

    def __new__(cls, name='Square', bounds1=(0, 1), bounds2=(0, 1)):
        dim = 2
        min_coords = (bounds1[0], bounds2[0])
        max_coords = (bounds1[1], bounds2[1])
        return super().__new__(cls, name, dim, min_coords, max_coords)

    @property
    def bounds1(self):
        return (self.min_coords[0], self.max_coords[0])

    @property
    def bounds2(self):
        return (self.min_coords[1], self.max_coords[1])

#==============================================================================
class Cube(NCube):

    def __new__(cls, name='Cube', bounds1=(0, 1), bounds2=(0, 1), bounds3=(0, 1)):
        dim = 3
        min_coords = (bounds1[0], bounds2[0], bounds3[0])
        max_coords = (bounds1[1], bounds2[1], bounds3[1])
        return super().__new__(cls, name, dim, min_coords, max_coords)

    @property
    def bounds1(self):
        return (self.min_coords[0], self.max_coords[0])

    @property
    def bounds2(self):
        return (self.min_coords[1], self.max_coords[1])

    @property
    def bounds3(self):
        return (self.min_coords[2], self.max_coords[2])

#==============================================================================
class BoundaryVector(IndexedBase):
    is_commutative = False

class NormalVector(BoundaryVector):
    pass

class MinusNormalVector(NormalVector):
    pass

class PlusNormalVector(NormalVector):
    pass

class TangentVector(BoundaryVector):
    pass

#==============================================================================
class ElementDomain(with_metaclass(Singleton, Basic)):
    pass

#==============================================================================
class BasicArea(AtomicExpr):

    def __new__(cls, domain):
        if not isinstance(domain, (BasicDomain, ElementDomain)):
            raise TypeError('expecting a BasicDomain or ElementDomain')

        return Basic.__new__(cls, domain)

    @property
    def domain(self):
        return self.args[0]

class DomainArea(BasicArea):
    pass

class ElementArea(BasicArea):
    pass

#==============================================================================
class BasicGeometryOperator(CalculusFunction):

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

#==============================================================================
class Area(BasicGeometryOperator):

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    @classmethod
    def eval(cls, *args):
        """."""

        if not args:
            return

        if not len(args) == 1:
            raise ValueError('Expecting one argument')

        expr = args[0]
        if isinstance(expr, Union):
            return Add(*[cls.eval(a) for a in expr.args])

        elif isinstance(expr, ElementDomain):
            return ElementArea(expr)

#        elif isinstance(expr, InteriorDomain):
#            return DomainArea(expr)

        return cls(expr, evaluate=False)



#==============================================================================
def split(domain, value):
    if domain.dtype['type'] == 'Line':
        assert(isinstance(value, (int, float)))

        # TODO assert value <- bounds
        bounds = domain.interior.bounds

        # ... left
        bounds = (bounds[0], value)
        I_left = Line(name='{name}_l'.format(name=domain.name),
                      bounds=bounds)
        # ...

        # ... right
        bounds = (value, bounds[1])
        I_right = Line(name='{name}_r'.format(name=domain.name),
                       bounds=bounds)
        # ...

        # ... interiors
        interiors = [I_left.interior, I_right.interior]
        # ...

        # ... external boundaries
        bnd_left = [b for b in I_left.boundary.as_tuple() if b.ext == -1]
        bnd_left = bnd_left[0]

        bnd_right = [b for b in I_right.boundary.as_tuple() if b.ext == 1]
        bnd_right = bnd_right[0]

        boundaries = [bnd_left, bnd_right]
        # ...

        # ... connectivity: internal interfaces
        int_left = [b for b in I_left.boundary.as_tuple() if b.ext == 1]
        int_left = int_left[0]

        int_right = [b for b in I_right.boundary.as_tuple() if b.ext == -1]
        int_right = int_right[0]

        connectivity = Connectivity()
        connectivity['I'] = (int_left, int_right)
        # ...

        return Domain(domain.name,
                      interiors=interiors,
                      boundaries=boundaries,
                      connectivity=connectivity)

    else:
        raise NotImplementedError('TODO')
