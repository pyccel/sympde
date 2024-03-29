from sympy import Symbol

from sympde.utilities.utils import lambdify_sympde
from .mapping import Mapping, BasicCallableMapping

__all__ = ('CallableMapping',)

#==============================================================================
class CallableMapping(BasicCallableMapping):

    def __init__( self, mapping, **kwargs ):

        # Extract information from class
        assert isinstance(mapping, Mapping)

        variables   = mapping.logical_coordinates
        expressions = mapping.expressions

        jac         = mapping.jacobian_expr
        inv_jac     = mapping.jacobian_inv_expr

        metric      = mapping.metric_expr
        metric_det  = mapping.metric_det_expr

        constants    = mapping.constants
        params       = {a.name:a for a in constants}
        params.update( kwargs )
        for p in params.values():
            assert not isinstance(p, Symbol)

        if params:
            subs = {a:params[a.name] for a in constants}
            expressions = expressions.subs(subs)
            jac         = jac.subs(subs)
            inv_jac     = inv_jac.subs(subs)
            metric      = metric.subs(subs)
            metric_det  = metric_det.subs(subs)

        # Callable function: __call__
        self._func_eval = tuple(lambdify_sympde( variables, expr) for expr in expressions)

        # Callable function: jac_mat
        self._jacobian = lambdify_sympde( variables, jac)

        # Callable function: jac_mat_inv
        self._jacobian_inv = lambdify_sympde( variables, inv_jac)

        # Callable function: metric
        self._metric = lambdify_sympde( variables, metric)

        # Callable function: metric_det
        self._metric_det = lambdify_sympde( variables, metric_det)

        # Symbolic information
        self._params           = params
        self._symbolic_mapping = mapping

    #--------------------------------------------------------------------------
    # Abstract interface
    #--------------------------------------------------------------------------
    def __call__(self, *eta):
        return tuple( f( *eta ) for f in self._func_eval)

    def jacobian(self, *eta):
        return self._jacobian( *eta )

    def jacobian_inv(self, *eta):
        """ Compute the inverse Jacobian matrix, if possible."""
        return self._jacobian_inv(*eta)

    def metric(self, *eta):
        return self._metric( *eta )

    def metric_det(self, *eta):
        return self._metric_det( *eta )

    @property
    def ldim(self):
        return self.symbolic_mapping.ldim

    @property
    def pdim(self):
        return self.symbolic_mapping.pdim

    #--------------------------------------------------------------------------
    # Symbolic information
    #--------------------------------------------------------------------------
    @property
    def params( self ):
        return self._params

    @property
    def symbolic_mapping( self ):
        return self._symbolic_mapping

