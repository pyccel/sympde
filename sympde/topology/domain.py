# coding: utf-8

import numpy as np
import h5py
import yaml
import os

from collections import abc

from sympy import Integer
from sympy.core.singleton import Singleton
from sympy.core.compatibility import with_metaclass, is_sequence
from sympy.core import Basic, symbols
from sympy.core.containers import Tuple
from sympy.tensor import IndexedBase, Indexed
from sympy.core import Add, Mul, Pow
from sympy.core.expr import AtomicExpr

from sympde.core.basic import CalculusFunction
from .basic            import BasicDomain, InteriorDomain, Boundary, Union, Connectivity
from .basic            import Interval, Interface, CornerBoundary, CornerInterface
from .basic            import ProductDomain

# TODO fix circular dependency between domain and mapping

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

    def __new__(cls, name, *, interiors=None, boundaries=None, dim=None,
                connectivity=None, mapping=None, logical_domain=None):
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
    def name(self):
        return self.args[0]

    @property
    def interior(self):
        return self.args[1]

    @property
    def boundary(self):
        return self.args[2]

    @property
    def mapping(self):
        return self.args[3]

    @property
    def logical_domain(self):
        return self._logical_domain

    @property
    def connectivity(self):
        return self._connectivity

    @property
    def dim(self):
        return self._dim

    @property
    def dtype(self):
        return self._dtype

    @property
    def interfaces(self):
        return self.connectivity.interfaces

    @property
    def corners(self):
        corners = getattr(self,'_corners', None)
        if corners is None:
            corners = self.get_shared_corners()
        self._corners = corners
        return corners

    def __len__(self):
        if isinstance(self.interior, InteriorDomain):
            return 1

        elif isinstance(self.interior, Union):
            return len(self.interior)

    @property
    def interior_names(self):
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
        """return boundary by name or (axis, ext)."""
        # ...

        if axis is None:
            assert(self.interior.dim == 1)
            axis = 0
        # ...

        if isinstance(self.boundary, Union):
            x = [i for i in self.boundary.args if i.ext == ext and i.axis == axis]
            if len(x) == 0:
                raise ValueError('> could not find boundary with axis {} and ext {}'.format(axis, ext))

            return x[0]

        elif isinstance(self.boundary, Boundary):
            if self.boundary.axis == axis and self.boundary.ext == ext:
                return self.boundary

        raise ValueError('> could not find boundary with axis {} and ext {}'.format(axis, ext))

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
        boundary     = self.boundary.todict()
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

    def export( self, filename ):

        yml = self.todict()

        # Dump metadata to string in YAML file format
        geo = yaml.dump( data   = yml,
                         sort_keys = None)

        # Create HDF5 file (in parallel mode if MPI communicator size > 1)
        h5 = h5py.File( filename, mode='w' )

        # Write geometry metadata as fixed-length array of ASCII characters
        h5['topology.yml'] = np.array( geo, dtype='S' )

        # Close HDF5 file
        h5.close()

    @classmethod
    def from_file( cls, filename ):

        # ... check extension of the file
        _, ext = os.path.splitext(filename)

        if not(ext == '.h5'):
            raise ValueError('> Only h5 files are supported')
        # ...
        from sympde.topology.mapping import Mapping, MultiPatchMapping

        h5  = h5py.File( filename, mode='r' )
        yml = yaml.load( h5['topology.yml'][()], Loader=yaml.SafeLoader )

        domain_name    = yml['name']
        dim            = int(yml['dim'])
        dtype          = yml['dtype']
        d_interior     = yml['interior']
        d_boundary     = yml['boundary']
        d_connectivity = yml['connectivity']

        if dtype == 'None': dtype = None

        if dtype is not None:
            if isinstance(d_interior, dict):
                d_interior = [d_interior]
                dtype      = [dtype]

        if dtype is not None and all(dtype):
            constructors = [globals()[dt['type']] for dt in dtype]
            interiors    = [cs(i['name'], **dt['parameters']) for cs,i,dt in zip(constructors, d_interior, dtype)]
            mappings     = [Mapping(I['mapping'], dim=dim) if I.get('mapping', "None") != "None" else None for I in d_interior]
            domains      = [mapping(i) if mapping else i for i,mapping in zip(interiors, mappings)]
            patch_index  = {I.name:ind for ind,I in enumerate(interiors)}

            boundaries   = []
            for bd in d_boundary:
                name = bd['patch']
                axis = bd['axis']
                ext  = bd['ext']
                i    = patch_index[name]
                bd   = domains[i].get_boundary(axis=int(axis), ext=int(ext))
                boundaries.append(bd)

            interfaces = []
            for _,(minus, plus) in d_connectivity.items():
                minus_name = minus['patch']
                minus_axis = minus['axis']
                minus_ext  = minus['ext']
                minus_patch_i = patch_index[minus_name]
                minus_bd   = domains[minus_patch_i].get_boundary(axis=int(minus_axis), ext=int(minus_ext))

                plus_name = plus['patch']
                plus_axis = plus['axis']
                plus_ext  = plus['ext']
                plus_patch_i = patch_index[plus_name]
                plus_bd   = domains[plus_patch_i].get_boundary(axis=int(plus_axis), ext=int(plus_ext))

                interfaces.append([minus_bd, plus_bd, 1])

            domain = domains[0]
            for p in domains[1:]:
                domain = domain.join(p, name=domain_name)

            for I in interfaces:
                domain = domain.join(domain, domain.name, bnd_minus=I[0], bnd_plus=I[1], direction=I[2])

            return domain

        # ... create sympde InteriorDomain (s)
        l_interiors = [Domain(i['name'], dim=dim).interior for i in d_interior]
        mappings    = [Mapping(i['mapping'], dim=dim) if i['mapping'] != 'None' else None for i in d_interior]
        l_interiors_index = {I.name:ind for ind,I in enumerate(l_interiors)}

        # create a dict of interiors accessed by name => needed for boundaries
        d_l_interiors  = {}
        d_l_boundaries = {}
        for i in l_interiors:
            d_l_interiors[i.name] = i
            d_l_boundaries[i.name] = []

        # ... create sympde Boundary (s)

        for desc in d_boundary:
            name  = desc['name']
            patch = d_l_interiors[desc['patch']]
            axis  = desc['axis']
            ext   = desc['ext']

            if axis == 'None': axis = None
            if ext  == 'None': ext  = None

            lbnd = Boundary(name, patch, axis=axis, ext=ext)
            d_l_boundaries[patch.name].append(lbnd)
        # ...

        subdomains = []
        for name,M in zip(d_l_interiors, mappings):
            I = d_l_interiors[name]
            B = d_l_boundaries[name]
            D = Domain(I.name, dim=dim, interiors=I, boundaries=B)
            D = M(D) if M else D
            subdomains.append(D)

#         ... create connectivity
        interfaces = []
        for edge,pair in d_connectivity.items():
            bnds = []
            for desc in pair:
                patch = d_l_interiors[desc['patch']]
                patch_index = l_interiors_index[patch.name]
                name  = desc['name']
                bnd   = Boundary(name, patch)
                M     = mappings[patch_index]
                bnd   = M(bnd) if M else bnd
                bnds.append(bnd)

            interfaces.append(bnds)
        # ...

        domain = subdomains[0]
        for p in subdomains[1:]:
            domain = domain.join(p, name=domain_name)

        for I in interfaces:
            domain = domain.join(domain, domain.name, bnd_minus=I[0], bnd_plus=I[1])

        return domain

    def join(self, other, name, bnd_minus=None, bnd_plus=None, direction=None):

        from sympde.topology.mapping import MultiPatchMapping
        # ... connectivity
        connectivity = Connectivity()

        if bnd_minus and bnd_plus:
            interface          = bnd_minus.join(bnd_plus, direction=direction)
            connectivity[interface.name] = interface

        for k,v in self.connectivity.items():
            connectivity[k] = v

        for k,v in other.connectivity.items():
            connectivity[k] = v

        # ... boundary
        boundaries = Union(self.boundary, other.boundary).complement(Union(bnd_minus, bnd_plus))
        boundaries = boundaries.as_tuple()

        # ... interiors
        interiors       = Union(self.interior, other.interior)
        if all(e.mapping for e in interiors):
            logical_interiors    = Union(*[e.logical_domain for e in interiors])
            logical_boundaries   = [e.logical_domain for e in boundaries]
            logical_connectivity = Connectivity()
            for k,v in connectivity.items():
                logical_connectivity[v.logical_domain.name] = v.logical_domain

            mapping        = MultiPatchMapping({e.logical_domain: e.mapping for e in interiors})
            logical_domain = Domain(name,
                            interiors=logical_interiors,
                            boundaries=logical_boundaries,
                            connectivity=logical_connectivity)
        else:
            mapping              = None
            logical_domain       = None

        # ...
        return Domain(name,
                      interiors=interiors,
                      boundaries=boundaries,
                      connectivity=connectivity,
                      mapping=mapping,
                      logical_domain=logical_domain)

    def get_shared_corners(self):
        """ Compute the corners shared by multiple patches in 2D """

        interfaces   = self.interfaces
        interfaces = (interfaces,) if isinstance(interfaces, Interface) else interfaces

        directions   = {i.plus:i.direction for i in interfaces}
        directions.update({i.minus:i.direction for i in interfaces})

        boundaries    = {i.minus:i.plus for i in interfaces}
        boundaries.update({value:key for key, value in boundaries.items()})

        not_treated_corners = set([tuple(set((b, n))) for b in boundaries for n in b.adjacent_boundaries])
        grouped_corners     = []

        while not_treated_corners:
            corner = not_treated_corners.pop()
            grouped_corners.append([corner])
            if not ( corner[0] in boundaries and corner[1] in boundaries):
                while corner[1] in boundaries:
                    bd1     = boundaries[corner[1]]
                    bd2     = bd1.domain.get_boundary(axis=corner[0].axis, ext=corner[0].ext)
                    corner = (bd1, bd2.rotate(directions[bd1]))
                    grouped_corners[-1].append(corner)

                corner = grouped_corners[-1][0]
                while corner[0] in boundaries:
                    bd2     = boundaries[corner[0]]
                    bd1     = bd2.domain.get_boundary(axis=corner[1].axis, ext=corner[1].ext)
                    corner = (bd1.rotate(directions[bd2]), bd2)
                    grouped_corners[-1].insert(0, corner)

            else:
                while corner[1] in boundaries:
                    bd1     = boundaries[corner[1]]
                    bd2     = bd1.domain.get_boundary(axis=corner[0].axis, ext=corner[0].ext)
                    corner = (bd1, bd2.rotate(directions[bd1]))
                    if corner == grouped_corners[-1][0]:
                        break
                    grouped_corners[-1].append(corner)
                else:
                    corner = grouped_corners[-1][0]
                    while corner[0] in boundaries:
                        bd2     = boundaries[corner[0]]
                        bd1     = bd2.domain.get_boundary(axis=corner[1].axis, ext=corner[1].ext)
                        corner = (bd1.rotate(directions[bd2]), bd2)
                        grouped_corners[-1].insert(0, corner)

            grouped_corners[-1] = tuple(tuple(set(c)) for c in grouped_corners[-1])
            not_treated_corners = not_treated_corners.difference(grouped_corners[-1])

        grouped_corners = set(tuple(grouped_corners))
        grouped_corners = Union(*[CornerInterface(*[CornerBoundary(*e) for e in cs]) for cs in grouped_corners])
        return grouped_corners

    def get_subdomain(self, names):
        """
        Returns an individual patch or a Union of patches of a multipatch domain.

        Parameters
        ----------
        names : tuple of str or str
            Names of the patches to join.
            If a string is given, the corresponding patch will be returned.
            If a tuple of strings is given, the Union of the corresponding subdomains will be returned.

        Notes
        -----
        The subdomain is returned as it was before being joined, which means that its boundary includes the
        boundaries that are part of an interface in the multipatch domain.
        """
        if names == ():
            return None

        if isinstance(names, str):
            names = (names,)
            assert names[0] in self.interior_names

        elif isinstance(names, tuple):
            assert all(isinstance(name, str) for name in names)
            assert len(set(names)) == len(names)
            assert all(name in self.interior_names or name == self.name for name in names)

        # Check trivial case of single patch domain
        if isinstance(self.interior, InteriorDomain):
            assert names[0] == self.interior.name
            return self

        # If all patches are joined we get the full domain
        # Same if the full domain is part of the union
        if len(names) == len(self.interior_names) or self.name in names:
            return self

        # Build dictionary of interiors accessed by names
        interior_dict = {i.name: i for i in self.interior.as_tuple()}

        # Build dictionary of boundaries
        if self.boundary is not None:
            if isinstance(self.boundary, Union):
                boundary_dict = {(b.domain.name, b.axis, b.ext): b for b in self.boundary.as_tuple()}
            else:
                b = self.boundary
                boundary_dict = {(b.domain.name, b.axis, b.ext): b}
        else:
            boundary_dict = {}

        # Build dictionary of interfaces
        if self.interfaces is not None:
            if isinstance(self.interfaces, Union):
                interfaces_dict = {(i.minus.domain.name, i.plus.domain.name): i for i in self.interfaces.as_tuple()}
            elif isinstance(self.interfaces, Interface):
                i = self.interfaces
                interfaces_dict = {(i.minus.domain.name, i.plus.domain.name): i}
        else:
            interfaces_dict = {}

        interfaces = []

        for name in names:
            if name == self.name:
                return self
            interior = interior_dict[name]
            boundaries = [boundary_dict.get((name, axis, ext)) for axis in range(self.dim) for ext in [-1, 1]]

            boundaries = [b for b in boundaries if b is not None]

            # Extract boundaries and interfaces from interfaces_dict
            for other_name in self.interior_names:
                if other_name != name:
                    i_minus = interfaces_dict.pop((name, other_name), None)
                    i_plus = interfaces_dict.pop((other_name, name), None)

                    if other_name not in names:
                        if i_minus is not None:
                            boundaries.append(i_minus.minus)
                        if i_plus is not None:
                            boundaries.append(i_plus.plus)
                    else:
                        if i_plus is not None:
                            interfaces.append((i_plus.name, i_plus))
                        if i_minus is not None:
                            interfaces.append((i_minus.name, i_minus))

            # Create domain with name, interior and boundaries
            new_domain = Domain(name=name, interiors=interior, boundaries=boundaries,
                                mapping=interior.mapping, logical_domain=interior.logical_domain)
            try:
                previous_domain = previous_domain.join(new_domain, name=f"{previous_domain.name}|{new_domain.name}")
            except NameError:
                previous_domain = new_domain

        # Add interfaces
        joined_domain = previous_domain
        for k,v in interfaces:
            joined_domain.connectivity[k] = v

        return joined_domain


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
        assert(not( ext  is None ))

        if axis is None:
            assert(self.dim == 1)
            axis = 0

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
        coordinates = symbols(coord_names)

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


