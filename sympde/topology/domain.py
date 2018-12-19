# coding: utf-8

import numpy as np
import h5py
import yaml
import yamlloader
import os

from collections import OrderedDict
from collections import abc

from sympy.core import Basic, Symbol
from sympy.core.containers import Tuple
from sympy.tensor import IndexedBase

from .basic import BasicDomain, InteriorDomain, Boundary, Union, Connectivity

# TODO add pdim
#==============================================================================
class Domain(BasicDomain):
    """
    Represents an undefined domain.
    A domain is defined by at least one interior domain and possible boundaries.
    A domain without a boundary is either infinite or periodic.
    A domain can also be constructed from a connectivity, in which case, only the
    name and connectivity need to be passed.

    Examples

    """
    def __new__(cls, name, interiors=None, boundaries=None, dim=None,
                connectivity=None, filename=None):
        # ...
        if not isinstance(name, str):
            raise TypeError('> name must be a string')
        # ...

        # ...
        if ( ( interiors is None ) and ( connectivity is None ) and
             (filename is None) and
             ( dim is None) ):
            raise ValueError('> either interiors or connectivity must be given')
        # ...

        # ...
        if not( filename is None ):
            self.read(filename)
        # ...

        # ...
        if not( interiors is None ):
            if ( not isinstance( interiors, (tuple, list, Tuple)) and
                 not isinstance( interiors, InteriorDomain) ):
                raise TypeError('> Expecting an iterable or a InteriorDomain')

            if isinstance( interiors, InteriorDomain ):
                interiors = [interiors]

            else:
                if not all([isinstance(i, InteriorDomain) for i in interiors]):
                    raise TypeError('> all interiors must be of type InteriorDomain')
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
        # ...

        # ...
        if not( connectivity is None ):
            if not isinstance( connectivity, Connectivity ):
                raise TypeError('> Expecting a Connectivity')

            # TODO check that patches appearing in connectivity are in interiors
        # ...

        # ...
        if not dim is None:
            assert(isinstance( dim, int ))

            interiors = [InteriorDomain(name, dim=dim)]
        # ...

        # ...
        if len(interiors) == 0:
            raise TypeError('No interior domain found')

        elif len(interiors) == 1:
            interiors = interiors[0]

        elif len(interiors) > 1:
            interiors = Union(*interiors)
        # ...

        # ...
        if len(boundaries) > 1:
            boundaries = Union(*boundaries)
        # ...

        obj = Basic.__new__(cls, name, interiors, boundaries)
        obj._connectivity = connectivity

        return obj

    @property
    def name(self):
        return self._args[0]

    @property
    def interior(self):
        return self._args[1]

    @property
    def boundary(self):
        return self._args[2]

    @property
    def connectivity(self):
        return self._connectivity

    @property
    def dim(self):
        return self.interior.dim

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
            return [i.name for i in self.interior._args]


    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

    def get_boundary(self, name):
        """return boundary by name."""
        if isinstance(self.boundary, Union):
            x = [i for i in self.boundary._args if i.name == name]
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
            x = [i for i in self.interior._args if i.name == name]
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
        ldim         = str(self.dim)
        pdim         = str(self.dim)
        interior     = self.interior.todict()
        boundary     = self.boundary.todict()
        connectivity = self.connectivity.todict()

        d = {'name':         name,
             'ldim':         ldim,
             'pdim':         pdim,
             'interior':     interior,
             'boundary':     boundary,
             'connectivity': connectivity}

        return OrderedDict(sorted(d.items()))

    def export( self, filename ):
        """
        Parameters
        ----------
        filename : str
          Name of HDF5 output file.

        """
        yml = self.todict()

        # Dump metadata to string in YAML file format
        geo = yaml.dump( data   = yml,
                         Dumper = yamlloader.ordereddict.Dumper )

        # Create HDF5 file (in parallel mode if MPI communicator size > 1)
        h5 = h5py.File( filename, mode='w' )

        # Write geometry metadata as fixed-length array of ASCII characters
        h5['geometry.yml'] = np.array( geo, dtype='S' )

        # Close HDF5 file
        h5.close()

    def read( self, filename ):
        raise NotImplementedError('TODO')

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

        # ... create a dict of patches, needed for the connectivity
        d_patches = {}
        for patch in patches:
            d_patches[patch.name] = patch

        d_patches = OrderedDict(sorted(d_patches.items()))
        # ...

        # ... read the connectivity
        # connectivity
        for k,v in yml['connectivity'].items():
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
        for desc in yml['boundaries']:
            patch = d_patches[desc['patch']]
            name  = desc['name']
            bnd   = Boundary(name, patch)
            bnds.append(bnd)
        # ...

        # ... close the h3 file
        h5.close()
        # ...

        # ...
        self._ldim       = ldim
        self._pdim       = pdim
        self._patches    = patches
        self._boundaries = bnds
        # ...

#==============================================================================
class BoundaryVector(IndexedBase):
    pass

class NormalVector(BoundaryVector):
    pass

class TangentVector(BoundaryVector):
    pass
