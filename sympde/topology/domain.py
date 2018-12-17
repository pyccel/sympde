# coding: utf-8


from collections import OrderedDict
from collections import abc

from sympy.core import Basic, Symbol
from sympy.core.containers import Tuple
from sympy.tensor import IndexedBase

from .basic import BasicDomain, InteriorDomain, Boundary, Union, Topology

#==============================================================================
class Domain(BasicDomain):
    """
    Represents an undefined domain.
    A domain is defined by at least one interior domain and possible boundaries.
    A domain without a boundary is either infinite or periodic.
    A domain can also be constructed from a topology, in which case, only the
    name and topology need to be passed.

    Examples

    """
    def __new__(cls, name, interiors=None, boundaries=None, dim=None,
                topology=None, filename=None):
        # ...
        if not isinstance(name, str):
            raise TypeError('> name must be a string')
        # ...

        # ...
        if ( ( interiors is None ) and ( topology is None ) and
             (filename is None) and
             ( dim is None) ):
            raise ValueError('> either interiors or topology must be given')
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
        if not( filename is None ):
            topology = Topology(filename=filename)
        # ...

        # ...
        if not( topology is None ):
            if not isinstance( topology, Topology ):
                raise TypeError('> Expecting a Topology')

            interiors  = topology.patches
            boundaries = topology.boundaries
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

        return Basic.__new__(cls, name, interiors, boundaries)

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
    def dim(self):
        return self.interior.dim

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))


#==============================================================================
class BoundaryVector(IndexedBase):
    pass

class NormalVector(BoundaryVector):
    pass

class TangentVector(BoundaryVector):
    pass
