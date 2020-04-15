# coding: utf-8

import numpy as np
import h5py
import yaml
import yamlloader
import os

from collections import OrderedDict
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
from .basic import BasicDomain, InteriorDomain, Boundary, Union, Connectivity
from .basic import Interval
from .basic import ProductDomain

# TODO add pdim
#==============================================================================
class Domain(BasicDomain):
    """
    Represents an undefined domain.
    A domain is defined by at least one interior domain and possible boundaries.
    A domain without a boundary is either infinite or periodic.
    A domain can also be constructed from a connectivity, in which case, only the
    name and connectivity need to be passed.

    """

    def __new__(cls, name, interiors=None, boundaries=None, dim=None,
                connectivity=None):
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
            if ( not isinstance( interiors, (tuple, list, Tuple)) and
                 not isinstance( interiors, InteriorDomain) ):
                raise TypeError('> Expecting an iterable or a InteriorDomain')

            if isinstance( interiors, InteriorDomain ):
                interiors = [interiors]

            else:
                unions    = [i for i in interiors if isinstance(i, Union)]
                interiors = [i for i in interiors if not isinstance(i, Union)]
                for union in unions:
                    interiors += list(union.as_tuple())

                if not all([isinstance(i, InteriorDomain) for i in interiors]):
                    raise TypeError('> all interiors must be of type InteriorDomain')

            interiors = Tuple(*interiors)
        # ...

        # ...
        if not( boundaries is None ):
            if ( not isinstance( boundaries, (tuple, list, Tuple)) and
                 not isinstance( boundaries, Boundary) ):
                raise TypeError('> Expecting an iterable or a Boundary')

            if isinstance( boundaries, Boundary ):
                boundaries = [boundaries]

            else:
                if not all([isinstance(i, Boundary) for i in boundaries]):
                    raise TypeError('> all boundaries must be of type Boundary')

        else:
            boundaries = []

        boundaries = Tuple(*boundaries)
        # ...

        # ...
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
            interiors = interiors[0]
            dtype = interiors.dtype
            dim   = interiors.dim

        elif len(interiors) > 1:
            dtype = interiors[0].dtype
            dim   = interiors[0].dim
            interiors = Union(*interiors)

        # ...
        if len(boundaries) > 1:
            boundaries = Union(*boundaries)
        # ...

        obj = Basic.__new__(cls, name, interiors, boundaries)
        obj._connectivity = connectivity
        obj._dtype        = dtype
        obj._dim          = dim

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


    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

    def get_boundary(self, name=None, axis=None, ext=None):
        """return boundary by name or (axis, ext)."""
        # ...
        if name is None:
            assert(not( ext  is None ))

            if axis is None:
                assert(self.interior.dim == 1)
                axis = 0

            if axis == 0:
                if ext == -1:
                    name = r'\Gamma_1'

                if ext == 1:
                    name = r'\Gamma_2'

            if axis == 1:
                if ext == -1:
                    name = r'\Gamma_3'

                if ext == 1:
                    name = r'\Gamma_4'

            if axis == 2:
                if ext == -1:
                    name = r'\Gamma_5'

                if ext == 1:
                    name = r'\Gamma_6'
        # ...

        if isinstance(self.boundary, Union):
            x = [i for i in self.boundary.args if i.name == name]
            if len(x) == 0:
                raise ValueError('> could not find boundary {}'.format(name))

            return x[0]

        elif isinstance(self.boundary, Boundary):
            if self.boundary.name == name:
                return self.boundary

            else:
                return None

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

        return OrderedDict(sorted(d.items()))

    def export( self, filename ):

        yml = self.todict()

        # Dump metadata to string in YAML file format
        geo = yaml.dump( data   = yml,
                         Dumper = yamlloader.ordereddict.Dumper )

        # Create HDF5 file (in parallel mode if MPI communicator size > 1)
        h5 = h5py.File( filename, mode='w' )

        # Write geometry metadata as fixed-length array of ASCII characters
        h5['topology.yml'] = np.array( geo, dtype='S' )

        # Close HDF5 file
        h5.close()

    @classmethod
    def from_file( cls, filename ):

        # ... check extension of the file
        basename, ext = os.path.splitext(filename)
        if not(ext == '.h5'):
            raise ValueError('> Only h5 files are supported')
        # ...

        h5  = h5py.File( filename, mode='r' )
        yml = yaml.load( h5['topology.yml'][()], Loader=yaml.SafeLoader )

        domain_name    = yml['name']
        dim            = yml['dim']
        dtype          = yml['dtype']
        d_interior     = yml['interior']
        d_boundary     = yml['boundary']
        d_connectivity = yml['connectivity']

        if dtype == 'None': dtype = None

        if dtype is not None:
            constructor = globals()[dtype['type']]
            return constructor(domain_name, **dtype['parameters'])

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

            connectivity[edge] = bnds
        # ...

        return Domain.__new__(cls, domain_name,
                              interiors=interior,
                              boundaries=boundary,
                              connectivity=connectivity)

    def join(self, other, name, bnd_minus, bnd_plus):
        # ... interiors
        interiors = [self.interior, other.interior]
        # ...

        # ... connectivity
        connectivity = Connectivity()
        # TODO be careful with '|' in psydac
        connectivity['{l}|{r}'.format(l=bnd_minus.domain.name, r=bnd_plus.domain.name)] = (bnd_minus, bnd_plus)
        for k,v in self.connectivity.items():
            connectivity[k] = v

        for k,v in other.connectivity.items():
            connectivity[k] = v
        # ...

        # ... boundary
        boundaries_minus = self.boundary.complement(bnd_minus)
        boundaries_plus  = other.boundary.complement(bnd_plus)

        boundaries = Union(boundaries_minus, boundaries_plus)
        boundaries = boundaries.as_tuple()
        # ...
        return Domain(name,
                      interiors=interiors,
                      boundaries=boundaries,
                      connectivity=connectivity)

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
        return self.domain.__hash__() + self._periods.__hash__()


#==============================================================================
class NCubeInterior(InteriorDomain):
    _min_coords = None
    _max_coords = None
    @property
    def min_coords(self):
        return self._min_coords

    @property
    def max_coords(self):
        return self._max_coords

#==============================================================================
# Ncube's properties (in addition to Domain's properties):
#   . min_coords (default value is tuple of zeros)
#   . max_coords (default value is tuple of ones)
#
class NCube(Domain):

    def __new__(cls, name, dim, min_coords, max_coords):

        assert isinstance(name, str)
        assert isinstance(dim, (int, Integer))
        assert isinstance(min_coords, (tuple, list, Tuple))
        assert isinstance(max_coords, (tuple, list, Tuple))

        if not name:
            raise ValueError("Name must be provided")

        if dim < 1:
            raise ValueError("Number of dimensions must be at least 1")

        if not (dim == len(min_coords) == len(max_coords)):
            raise ValueError("Input arguments must have 'dim' components")

        if not all(xmin < xmax for xmin, xmax in zip(min_coords, max_coords)):
            raise ValueError("Min coordinates must be smaller than max")

        # Choose coordinate names. TODO: use unique convention
        if dim <= 3:
            coord_names = ('x', 'y', 'z')[:dim]
        else:
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

        interior             = NCubeInterior(interior, dtype=dtype)
        interior._min_coords = tuple(min_coords)
        interior._max_coords = tuple(max_coords)

        boundaries = []
        i = 1
        for axis in range(dim):
            for ext in [-1, 1]:
                bnd_name = r'\Gamma_{}'.format(i)
                Gamma = Boundary(bnd_name, interior, axis=axis, ext=ext)
                boundaries += [Gamma]
                i += 1

        # Create instance of given type
        obj = super().__new__(cls, name, interiors=[interior],
                boundaries=boundaries)

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
    pass

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

