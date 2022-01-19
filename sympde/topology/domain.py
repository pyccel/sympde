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

        else:
            dtype = interiors[0].dtype
            dim   = interiors[0].dim
            interiors = Union(*interiors)

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
        dim            = yml['dim']
        dtype          = yml['dtype']
        d_interior     = yml['interior']
        d_boundary     = yml['boundary']
        d_connectivity = yml['connectivity']
        mapping        = Mapping('{}_mapping'.format(domain_name), dim=int(dim))

        if dtype == 'None': dtype = None

        if dtype is not None:
            constructor = globals()[dtype['type']]
            return mapping(constructor(domain_name, **dtype['parameters']))

        # ... create sympde InteriorDomain (s)
        interior = [InteriorDomain(i['name'], dim=dim) for i in d_interior]

        # create a dict of interiors accessed by name => needed for boundaries
        d_interior = {}
        for i in interior:
            d_interior[i.name] = i

        if len(interior) == 1:
            interior = interior[0]
        # ...

        # ... create sympde Boundary (s)
        boundary = []
        for desc in d_boundary:
            name  = desc['name']
            patch = d_interior[desc['patch']]
            axis  = desc['axis']
            ext   = desc['ext']

            if axis == 'None': axis = None
            if ext  == 'None': ext  = None

            bnd = Boundary(name, patch, axis=axis, ext=ext)
            boundary.append(bnd)

        if len(boundary) == 1:
            boundary = boundary[0]
        # ...

        # ... create connectivity
        connectivity = Connectivity()
        for edge,pair in d_connectivity.items():
            bnds = []
            for desc in pair:
                patch = d_interior[desc['patch']]
                name  = desc['name']
                bnd   = Boundary(name, patch)
                bnds.append(bnd)

            connectivity[edge] = Interface(edge, bnds[0], bnds[1])
        # ...

        logical_domain = Domain(domain_name, dim=dim,
                                interiors=interior,
                                boundaries=boundary,
                                connectivity=connectivity)

        obj = Domain.__new__(cls, domain_name,
                              interiors=interior,
                              boundaries=boundary,
                              connectivity=connectivity,
                              mapping = mapping, logical_domain=logical_domain)
        return obj

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
        raise ValueError('> could not find boundary {}'.format(name))
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
        intervals   = [Interval('{}_{}'.format(name, c.name), coordinate=c, bounds=(xmin, xmax))
                       for c, xmin, xmax in zip(coordinates, min_coords, max_coords)]

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

        if len(intervals) == 1:
            interior = intervals[0]
        else:
            interior = ProductDomain(*intervals, name=name)

        interior = NCubeInterior(interior, dim=dim, dtype=dtype, min_coords=tuple(min_coords), max_coords=tuple(max_coords))

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


