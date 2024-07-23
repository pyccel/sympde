from abc import ABC, ABCMeta, abstractmethod
from sympy import IndexedBase, Matrix, ImmutableDenseMatrix
from sympy.core.containers import Tuple
from .basic            import BasicDomain, InteriorDomain, Boundary, Union, Connectivity, Interface
from .domain import Domain, NCubeInterior

__all__ = (
    'MappingMeta',
    'BaseMapping',
)

class MappingMeta(ABCMeta,type(IndexedBase)):
    pass

#==============================================================================
class BaseMapping(IndexedBase):
    """
    Transformation of coordinates, which can be evaluated.

    F: R^l -> R^p
    F(eta) = x

    with l <= p
    """
    
    _name = None
    _ldim = None
    _pdim = None

    
    def __new__(cls, name, dim=None, **kwargs):
        ldim        = kwargs.pop('ldim', cls._ldim)
        pdim        = kwargs.pop('pdim', cls._pdim)

        dims = [dim, ldim, pdim]
        
        for i,d in enumerate(dims):
            if isinstance(d, (tuple, list, Tuple, Matrix, ImmutableDenseMatrix)):
                if not len(d) == 1:
                    raise ValueError('> Expecting a tuple, list, Tuple of length 1')
                dims[i] = d[0]

        dim, ldim, pdim = dims

        if dim is None:
            assert ldim is not None
            assert pdim is not None
            assert pdim >= ldim
        else:
            ldim = dim
            pdim = dim

        obj = IndexedBase.__new__(cls, name, shape=pdim)
        
        obj._name                = name
        obj._ldim                = ldim
        obj._pdim                = pdim
        
        return obj

    
    
    def __call__(self, *args):
        """ Evaluate mapping at either a list of nd-arrays or the full domain."""
        if len(args) == 1 and isinstance(args[0], BasicDomain):
            domain=args[0]
            assert(isinstance(domain, BasicDomain))
            return MappedDomain(self, domain)
        else:
            pass
            
    
    def jacobian_eval( self, *args ):
        """ Compute Jacobian matrix at the list of nd-arrays. """
        pass 

    def jacobian_inv_eval( self, *args ):
        """ Compute inverse Jacobian matrix at the list of nd-arrays.
            An exception should be raised if the matrix is singular.
        """
        pass
    
    def metric_eval( self, *args ):
        """ Compute components of metric tensor at list of nd-arrays. """
        pass
    
    def metric_det_eval( self, *args ):
        """ Compute determinant of metric tensor at the list of nd-arrays. """
        pass
    
    @property
    def ldim( self ):
        """ Number of logical/parametric dimensions in mapping
            (= number of eta components).
        """
        pass
    
    @property
    def pdim( self ):
        """ Number of physical dimensions in mapping
            (= number of x components).""" 
        pass
    
    
#==============================================================================
class MappedDomain(BasicDomain):
    """."""

    @cacheit
    def __new__(cls, mapping, logical_domain):
        assert(isinstance(mapping,BaseMapping))
        assert(isinstance(logical_domain, BasicDomain))
        if isinstance(logical_domain, Domain):
            kwargs = dict(
            dim            = logical_domain._dim,
            mapping        = mapping,
            logical_domain = logical_domain)
            boundaries     = logical_domain.boundary
            interiors      = logical_domain.interior

            if isinstance(interiors, Union):
                kwargs['interiors'] = Union(*[mapping(a) for a in interiors.args])
            else:
                kwargs['interiors'] = mapping(interiors)

            if isinstance(boundaries, Union):
                kwargs['boundaries'] = [mapping(a) for a in boundaries.args]
            elif boundaries:
                kwargs['boundaries'] = mapping(boundaries)

            interfaces =  logical_domain.connectivity.interfaces
            if interfaces:
                print("interfaces")
                if isinstance(interfaces, Union):
                    interfaces = interfaces.args
                else:
                    interfaces = [interfaces]
                connectivity = {}
                for e in interfaces:
                    connectivity[e.name] = Interface(e.name, mapping(e.minus), mapping(e.plus))
                kwargs['connectivity'] = Connectivity(connectivity)

            name = '{}({})'.format(str(mapping.name), str(logical_domain.name))
            print("return Domain(name,**kwargs)")
            return Domain(name, **kwargs)

        elif isinstance(logical_domain, NCubeInterior):
            name  = logical_domain.name
            dim   = logical_domain.dim
            dtype = logical_domain.dtype
            min_coords = logical_domain.min_coords
            max_coords = logical_domain.max_coords
            name = '{}({})'.format(str(mapping.name), str(name))
            return NCubeInterior(name, dim, dtype, min_coords, max_coords, mapping, logical_domain)
        elif isinstance(logical_domain, InteriorDomain):
            name  = logical_domain.name
            dim   = logical_domain.dim
            dtype = logical_domain.dtype
            name = '{}({})'.format(str(mapping.name), str(name))
            return InteriorDomain(name, dim, dtype, mapping, logical_domain)
        elif isinstance(logical_domain, Boundary):
            name   = logical_domain.name
            axis   = logical_domain.axis
            ext    = logical_domain.ext
            domain = mapping(logical_domain.domain)
            return Boundary(name, domain, axis, ext, mapping, logical_domain)
        else:
            raise NotImplementedError('TODO')