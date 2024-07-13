from abc import ABC, ABCMeta, abstractmethod
from sympy import IndexedBase

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
    
    _ldim         = None
    _pdim         = None

    
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

        
        
        self._name                = name
        self._ldim                = ldim
        self._pdim                = pdim


    def __call__(self, *args):
        """ Evaluate mapping at either a list of nd-arrays or the full domain."""
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