# coding: utf-8

import numpy as np
import itertools as it

from sympy                 import lambdify
from sympy.core            import Symbol

from .basic       import BasicDomain
from .base_mapping import BaseMapping, MappedDomain


# TODO fix circular dependency between sympde.topology.domain and sympde.topology.mapping
# TODO fix circular dependency between sympde.expr.evaluation and sympde.topology.mapping

__all__ = (
    'lambdify_sympde',
    'BaseAnalyticMapping',
)


def lambdify_sympde(variables, expr):
    """
    Custom lambify function that covers the
    shortcomings of sympy's lambdify. Most notably,
    this function uses numpy broadcasting rules to
    compute the shape of the output.

    Parameters
    ----------
    variables : sympy.core.symbol.Symbol or list of sympy.core.symbol.Symbol
        variables that appear in the expression
    expr :
        Sympy expression

    Returns
    -------
    lambda_f : callable
        Lambdified function built using numpy.

    Notes
    -----
    Compared to Sympy's lambdify, this function
    is capable of properly handling constant values,
    and array_like structures where not all components
    depend on all variables. See below.

    Examples
    --------
    >>> import numpy as np
    >>> from sympy import symbols,  Matrix
    >>> from sympde.utilities.utils import lambdify_sympde
    >>> x, y = symbols("x,y")
    >>> expr = Matrix([[x, x + y], [0, y]])
    >>> f = lambdify_sympde([x,y], expr)
    >>> f(np.array([[0, 1]]), np.array([[2], [3]]))
    array([[[[0., 1.],
             [0., 1.]],

            [[2., 3.],
             [3., 4.]]],


           [[[0., 0.],
             [0., 0.]],

            [[2., 2.],
             [3., 3.]]]])
    """
    array_expr = np.asarray(expr)
    scalar_shape = array_expr.shape
    if scalar_shape == ():
        f = lambdify(variables, expr, 'numpy')
        def f_vec_sc(*XYZ):
            b = np.broadcast(*XYZ)
            if b.ndim == 0:
                return f(*XYZ)
            temp = np.asarray(f(*XYZ))
            if b.shape == temp.shape:
                return temp

            result = np.zeros(b.shape)
            result[...] = temp
            return result
        return f_vec_sc

    else:
        scalar_functions = {}
        for multi_index in it.product(*tuple(range(s) for s in scalar_shape)):
            scalar_functions[multi_index] = lambdify(variables, array_expr[multi_index], 'numpy')

        def f_vec_v(*XYZ):
            b = np.broadcast(*XYZ)
            result = np.zeros(scalar_shape + b.shape)
            for multi_index in it.product(*tuple(range(s) for s in scalar_shape)):
                result[multi_index] = scalar_functions[multi_index](*XYZ)
            return result
        return f_vec_v


#==============================================================================
class BaseAnalyticMapping(BaseMapping):
    """
    Represents a BaseAnalyticMapping object.

    Examples

    """

    def __new__(cls, name, dim=None, **kwargs):

        obj = super().__new__(cls, name, dim, **kwargs)

        if obj.expressions :
            obj._func_eval = tuple(lambdify_sympde( obj._logical_coordinates, expr) for expr in obj._expressions)
            obj._jac_eval = lambdify_sympde( obj._logical_coordinates, obj._jac)
            obj._inv_jac_eval = lambdify_sympde( obj._logical_coordinates, obj._inv_jac)
            obj._metric_eval = lambdify_sympde( obj._logical_coordinates, obj._metric)
            obj._metric_det_eval = lambdify_sympde( obj._logical_coordinates, obj._metric_det)
        else:
            raise TypeError("BaseAnalyticMapping should have an expression")

        return obj

    #--------------------------------------------------------------------------
    #Abstract Interface :

    def _evaluate_domain( self, domain ):
        assert(isinstance(domain, BasicDomain))
        return MappedDomain(self, domain)

    def _evaluate( self, *Xs ):
        #int, float or numpy arrays
        if self._func_eval is None :
            raise TypeError("not a callable object")
        else :
            assert len(Xs)==self.ldim
            Xshape = np.shape(Xs[0])
            for X in Xs:
                assert np.shape(X) == Xshape
            return tuple( f( *Xs ) for f in self._func_eval)

    def _jacobian_evaluate( self, *Xs ):
        #int, float or numpy arrays
        if self._jac_eval is None:
            raise TypeError("not a callable object")
        else :
            assert len(Xs)==self.ldim
            Xshape = np.shape(Xs[0]) 
            for X in Xs:
                assert np.shape(X) == Xshape
            return self._jac_eval(*Xs)
   
    def _jacobian_inv_evaluate( self, *Xs ):
        #int, float or numpy arrays 
        if self._inv_jac_eval is None:
            raise TypeError("not a callable object")
        else :
            assert len(Xs)==self.ldim
            Xshape = np.shape(Xs[0]) 
            for X in Xs:
                assert np.shape(X) == Xshape
            return  self._inv_jac_eval(*Xs)
        
    def _metric_evaluate( self, *Xs ):
        if self._metric_eval is None:
            raise TypeError("not a callable object")
        else :
            assert len(Xs)==self.ldim
            Xshape = np.shape(Xs[0]) 
            for X in Xs:
                assert np.shape(X) == Xshape
            return  self._metric_eval(*Xs)
        
    def _metric_det_evaluate( self, *Xs ):
        if self._metric_det_eval is None:
            raise TypeError("not a callable object")
        else :
            assert len(Xs)==self.ldim
            Xshape = np.shape(Xs[0]) 
            for X in Xs:
                assert np.shape(X) == Xshape
            return self._metric_det_eval(*Xs)
    
    def __call__( self, *args ):
        if len(args) == 1 and isinstance(args[0], BasicDomain):
            return self._evaluate_domain(args[0])
        elif all(isinstance(arg, (int, float, Symbol, np.ndarray)) for arg in args):
            return self._evaluate(*args)
        else:
            raise TypeError("Invalid arguments for __call__")

    def jacobian_eval( self, *args ):
        if all(isinstance(arg, (int, float, Symbol, np.ndarray)) for arg in args):
            return self._jacobian_evaluate(*args)
        else:
            raise TypeError("Invalid arguments for jacobian_eval")         

    def jacobian_inv_eval( self, *args ):
        if all(isinstance(arg, (int, float, Symbol, np.ndarray)) for arg in args):
            return self._jacobian_inv_evaluate(*args)
        else:
            raise TypeError("Invalid arguments for jacobian_inv_eval")

    def metric_eval( self, *args ):
        if all(isinstance(arg, (int, float, Symbol, np.ndarray)) for arg in args):
            return self._metric_evaluate(*args)
        else:
            raise TypeError("Invalid arguments for metric_eval")

    def metric_det_eval( self, *args ):
        if all(isinstance(arg, (int, float, Symbol, np.ndarray)) for arg in args):
            return self._metric_det_evaluate(*args)
        else:
            raise TypeError("Invalid arguments for metric_det_eval")

