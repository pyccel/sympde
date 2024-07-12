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
    @abstractmethod
    def __call__(self, *args):
        """ Evaluate mapping at either a list of nd-arrays or the full domain."""

    @abstractmethod
    def jacobian_eval( self, *args ):
        """ Compute Jacobian matrix at the list of nd-arrays. """

    @abstractmethod
    def jacobian_inv_eval( self, *args ):
        """ Compute inverse Jacobian matrix at the list of nd-arrays.
            An exception should be raised if the matrix is singular.
        """

    @abstractmethod
    def metric_eval( self, *args ):
        """ Compute components of metric tensor at list of nd-arrays. """

    @abstractmethod
    def metric_det_eval( self, *args ):
        """ Compute determinant of metric tensor at the list of nd-arrays. """

    @property
    @abstractmethod
    def ldim( self ):
        """ Number of logical/parametric dimensions in mapping
            (= number of eta components).
        """

    @property
    @abstractmethod
    def pdim( self ):
        """ Number of physical dimensions in mapping
            (= number of x components)."""