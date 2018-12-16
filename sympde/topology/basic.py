# coding: utf-8


from collections import OrderedDict
from collections import abc
import h5py
import yaml
import yamlloader
import os
import string
import random

from sympy.core import Basic, Symbol
from sympy.core.containers import Tuple
from sympy.tensor import IndexedBase

#==============================================================================
class BasicDomain(Basic):
    _dim = None
    _coordinates = None

    @property
    def dim(self):
        return self._dim

    @property
    def coordinates(self):
        dim = self.dim
        if self._coordinates is None:
            xyz = ['x', 'y', 'z'][:dim]
            xyz = [Symbol(i) for i in xyz]
            self._coordinates = xyz

        if dim == 1:
            return self._coordinates[0]
        else:
            return self._coordinates

#==============================================================================
class InteriorDomain(BasicDomain):
    """
    Represents an undefined interior domain.

    Examples

    """
    def __new__(cls, name, dim):
        obj = Basic.__new__(cls, name)
        obj._dim = dim
        return obj

    @property
    def name(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

#==============================================================================
class Union(BasicDomain):
    def __new__(cls, *args):
        args = Tuple(*args)
        if not all( [isinstance(i, BasicDomain) for i in args] ):
            raise TypeError('arguments must be of BasicDomain type')

        assert(len(args) > 1)

        dim = args[0].dim
        dims = [a.dim for a in args[1:]]
        if not all( [d == dim for d in dims]):
            raise ValueError('arguments must have the same dimension')

        return Basic.__new__(cls, *args)

    @property
    def dim(self):
        return self._args[0].dim

    def __len__(self):
        return len(self._args)


#==============================================================================
class Boundary(BasicDomain):
    """
    Represents an undefined boundary over a domain.

    Examples

    """
    def __new__(cls, name, domain):
        obj = Basic.__new__(cls, name)
        obj._domain = domain
        return obj

    @property
    def name(self):
        return self._args[0]

    @property
    def domain(self):
        return self._domain

    @property
    def dim(self):
        return self.domain.dim

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

    def __add__(self, other):
        if isinstance(other, ComplementBoundary):
            raise TypeError('> Cannot add complement of boundary')

        return Union(self, other)

    def todict(self):
        domain = self.domain
        name   = self.name
        return OrderedDict( [('patch', domain.name ),
                             ('name' , name        ) ] )


#==============================================================================
class Edge(object):
    def __init__(self, name):
        self._name = name

    @property
    def name(self):
        return self._name


#==============================================================================
class Topology(abc.Mapping):
    _patches  = []
    _boundaries = None

    def __init__(self, data=None, boundaries=None, filename=None):
        # ...
        if data is None:
            data = {}

        else:
            assert( isinstance( data, (dict, OrderedDict)) )

        self._data = data
        # ...

        # ...
        if boundaries is None:
            boundaries = {}

        else:
            assert( isinstance( boundaries, (list, tuple) ) )

        self._boundaries = boundaries
        # ...

        if not( filename is None ):
            self.read(filename)

    @property
    def boundaries(self):
        return self._boundaries

    @property
    def patches(self):
        return self._patches

    def read( self, filename ):
        # ... check extension of the file
        basename, ext = os.path.splitext(filename)
        if not(ext == '.h5'):
            raise ValueError('> Only h5 files are supported')
        # ...

        kwargs = {}

        h5  = h5py.File( filename, mode='r', **kwargs )
        yml = yaml.load( h5['geometry.yml'].value )

        ldim = yml['ldim']
        pdim = yml['pdim']

        n_patches = len( yml['patches'] )

        # ...
        if n_patches == 0:

            h5.close()
            raise ValueError( "Input file contains no patches." )
        # ...

        # ... read patchs
        patches = []
        for i_patch in range( n_patches ):

            item  = yml['patches'][i_patch]
            dtype = item['type']
            name  = item['name']

            domain = InteriorDomain(name, dim=ldim)
            patches.append(domain)
        # ...

        # ... create a dict of patches, needed for the topology
        d_patches = {}
        for patch in patches:
            d_patches[patch.name] = patch

        d_patches = OrderedDict(sorted(d_patches.items()))
        # ...

        # ... read the topology
        # connectivity
        for k,v in yml['topology']['connectivity'].items():
            edge = Edge(k)
            bnds = []
            for desc in v:
                patch = d_patches[desc['patch']]
                name  = desc['name']
                bnd   = Boundary(name, patch)
                bnds.append(bnd)

            self[edge] = bnds

        # boundaries
        bnds = []
        for desc in yml['topology']['boundaries']:
            patch = d_patches[desc['patch']]
            name  = desc['name']
            bnd   = Boundary(name, patch)
            bnds.append(bnd)
        # ...

        # ... close the h5 file
        h5.close()
        # ...

        # ...
        self._ldim       = ldim
        self._pdim       = pdim
        self._patches    = patches
        self._boundaries = bnds
        # ...

    def export(self, filename):
        raise NotImplementedError('TODO')

    def todict(self):
        # ... create the connectivity
        connectivity = {}
        data = OrderedDict(sorted(self._data.items()))
        for edge, pair in data.items():
            connectivity[edge.name] = [bnd.todict() for bnd in pair]

        connectivity = OrderedDict(sorted(connectivity.items()))
        # ...

        # ...
        boundaries = [bnd.todict() for bnd in self.boundaries]
        # ...

        return {'connectivity': connectivity,
                'boundaries': boundaries}


    def __setitem__(self, key, value):
        assert( isinstance( key, Edge ) )
        assert( isinstance( value, (tuple, list)  ) )
        assert( len(value) in [1, 2] )
        assert( all( [isinstance( P, Boundary ) for P in value ] ) )

        self._data[key] = value

    # ==========================================
    #  abstract methods
    # ==========================================
    def __getitem__(self, key):
        return self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)
    # ==========================================
