from .mapping import Mapping
from sympy import lambdify, Symbol

class CallableMapping:

    def __init__( self, mapping, **kwargs ):

        # Extract information from class
        assert isinstance(mapping, Mapping)
        
        variables   = mapping.logical_coordinates
        expressions = mapping.expressions
        
        jac         = mapping.jacobian_expr
        inv_jac     = mapping.jacobian_inv_expr

        metric      = mapping.metric_expr
        metric_det  = mapping.metric_det_expr

        free_symbols = mapping.free_symbols
        params       = {a.name:a for a in free_symbols}
        params.update( kwargs )

        for p in params.values():
            assert not isinstance(p, Symbol)

        self._params = params
        self._symbolic_mapping = mapping

        # Callable function: __call__
        self._func_eval = tuple(lambdify( variables, expr, 'numpy' ) for expr in expressions)

        # Callable function: jac_mat
        self._jacobian = lambdify( variables, jac, 'numpy' )

        # Callable function: jac_mat_inv
        self._jacobian_inv = lambdify( variables, inv_jac, 'numpy' )
   
        # Callable function: metric
        self._metric = lambdify( variables, metric, 'numpy' )

        # Callable function: metric_det
        self._metric_det = lambdify( variables, metric_det, 'numpy' )


    #--------------------------------------------------------------------------
    # Abstract interface
    #--------------------------------------------------------------------------
    def __call__( self, eta ):
        return self._func_eval( eta )

    def jacobian( self, eta ):
        return self._jacobian( eta )

    def jacobian_inv( self, eta ):
        return self._jacobian_inv( eta )

    def metric( self, eta ):
        return self._metric( eta )

    def metric_det( self, eta ):
        return self._metric_det( eta )

    @property
    def ldim( self ):
        return type( self ).symbolic.ldim

    @property
    def pdim( self ):
        return type( self ).symbolic.pdim

    #--------------------------------------------------------------------------
    # Symbolic information
    #--------------------------------------------------------------------------
    @property
    def params( self ):
        return self._params
    default_params = dict( x0=0, y0=0, z0=0, k=0.3, D=0.2, c=0.5 )
