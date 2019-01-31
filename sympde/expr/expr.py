# coding: utf-8

from itertools import groupby

from sympy.core import Basic
from sympy.core import Symbol
from sympy.core import Function
from sympy.simplify.simplify import simplify
from sympy import collect
from sympy.series.order import Order
from sympy.core import Expr, Add, Mul, Pow
from sympy import S
from sympy.core.containers import Tuple
from sympy import Indexed, IndexedBase, Matrix, ImmutableDenseMatrix
from sympy import expand
from sympy import Integer, Float
from sympy.core.expr import AtomicExpr
from sympy.physics.quantum import TensorProduct
from sympy.series.series import series
from sympy.core.compatibility import is_sequence

from sympde.core.basic import _coeffs_registery
from sympde.core.basic import CalculusFunction
from sympde.core.basic import Constant
from sympde.core.algebra import (Dot_1d,
                 Dot_2d, Inner_2d, Cross_2d,
                 Dot_3d, Inner_3d, Cross_3d)
from sympde.core.utils import random_string

from sympde.calculus import Dot, Inner, Cross
from sympde.calculus import Grad, Rot, Curl, Div
from sympde.calculus import Bracket
from sympde.calculus import Laplace
from sympde.calculus.core import _generic_ops

from sympde.topology import BasicDomain, Domain, MappedDomain, Union, Interval
from sympde.topology import BoundaryVector, NormalVector, TangentVector, Boundary
from sympde.topology.derivatives import _partial_derivatives
from sympde.topology.derivatives import partial_derivative_as_symbol
from sympde.topology.derivatives import sort_partial_derivatives
from sympde.topology.derivatives import get_atom_derivatives
from sympde.topology.derivatives import dx, dy, dz
from sympde.topology.derivatives import (Grad_1d, Div_1d,
                                         Grad_2d, Curl_2d, Rot_2d, Div_2d,
                                         Grad_3d, Curl_3d, Div_3d)
from sympde.topology.derivatives import Bracket_2d
from sympde.topology.derivatives import Laplace_1d, Laplace_2d, Laplace_3d
from sympde.topology.derivatives import Hessian_1d, Hessian_2d, Hessian_3d
from sympde.topology.space import FunctionSpace
from sympde.topology.space import ProductSpace
from sympde.topology.space import TestFunction
from sympde.topology.space import VectorTestFunction
from sympde.topology.space import IndexedTestTrial
from sympde.topology.space import Unknown, VectorUnknown
from sympde.topology.space import Trace
from sympde.topology.space import Field, VectorField, IndexedVectorField
from sympde.topology.measure import CanonicalMeasure
from sympde.topology.measure import CartesianMeasure
from sympde.topology.measure import Measure

from .errors import UnconsistentError
from .errors import UnconsistentLinearExpressionError

# TODO to be moved here
from .form import is_linear_form


#==============================================================================
def _sanitize_arguments(arguments, is_bilinear=False, is_linear=False):

    # ...
    if is_bilinear or is_linear:

        if is_bilinear:
            test_functions = arguments[0]

        elif is_linear:
            test_functions = arguments

        if isinstance(test_functions, (TestFunction, VectorTestFunction)):
            test_functions = [test_functions]

        elif isinstance(test_functions, (tuple, list, Tuple)):
            are_valid = [isinstance(i, (TestFunction, VectorTestFunction)) for i in test_functions]
            if not all(are_valid):
                raise TypeError('> Wrong arguments for test functions')

        else:
            msg = 'Wrong type for test function(s). given {}'.format(type(test_functions))
            raise TypeError(msg)

        test_functions = Tuple(*test_functions)
    # ...

    # ...
    if is_bilinear:

        trial_functions = arguments[1]
        if isinstance(trial_functions, (TestFunction, VectorTestFunction)):
            trial_functions = [trial_functions]

        elif isinstance(trial_functions, (tuple, list, Tuple)):
            are_valid = [isinstance(i, (TestFunction, VectorTestFunction)) for i in trial_functions]
            if not all(are_valid):
                raise TypeError('> Wrong arguments for trial functions')

        else:
            msg = 'Wrong type for trial function(s). given {}'.format(type(trial_functions))
            raise TypeError(msg)

        trial_functions = Tuple(*trial_functions)
    # ...

    if is_bilinear:
        args = [test_functions, trial_functions]
        args = Tuple(*args)

    else:
        args = Tuple(*test_functions)

    return args


#==============================================================================
class BasicExpr(Expr):

    # TODO use .atoms
    @property
    def fields(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, (Field, VectorField))]
        # no redanduncy
        return sorted(list(set(ls)))

    # TODO use .atoms
    @property
    def constants(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Constant)]
        # no redanduncy
        return list(set(ls))

#==============================================================================
class LinearExpr(BasicExpr):
    is_Function = True

    def __new__(cls, arguments, expr, check=False):
        # ...
        if check and not is_linear_form(expr, arguments):
            msg = '> Expression is not linear'
            raise UnconsistentLinearExpressionError(msg)
        # ...

        args = _sanitize_arguments(arguments, is_linear=True)
        return Basic.__new__(cls, args, expr)

    @property
    def variables(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

    def __call__(self, *args):
        args = Tuple(*args)

        return self.expr.xreplace(dict(list(zip(self.variables, args))))

#==============================================================================
class BilinearExpr(BasicExpr):

    def __new__(cls, arguments, expr, check=False):
        # ...
        if not isinstance(arguments, (tuple, list, Tuple)):
            raise TypeError('(test, trial) must be a tuple, list or Tuple')

        if not(len(arguments) == 2):
            raise ValueError('Expecting a couple (test, trial)')
        # ...

        # ...
        if check and not is_bilinear_form(expr, arguments):
            msg = '> Expression is not bilinear'
            raise UnconsistentLinearExpressionError(msg)
        # ...

        args = _sanitize_arguments(arguments, is_bilinear=True)
        return Basic.__new__(cls, args, expr)

    @property
    def variables(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

    def __call__(self, *args):
        assert(len(args) == 2)

        left,right = args
        if not is_sequence(left):
            left = [left]

        if not is_sequence(right):
            right = [right]

        args = Tuple(*left, *right)

        variables = Tuple(*self.variables[0], *self.variables[1])

        return self.expr.xreplace(dict(list(zip(variables, args))))


#==============================================================================
# TODO to be moved
class BasicIntegral(CalculusFunction):

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

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

    @classmethod
    def eval(cls, *_args):
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]

        boundary = list(expr.atoms(Boundary))
        if boundary:
            return BoundaryIntegral(expr)

        else:
            return DomainIntegral(expr)

#==============================================================================
# TODO to be moved
class DomainIntegral(BasicIntegral):

    @classmethod
    def eval(cls, *_args):
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]

        return cls(expr, evaluate=False)

#==============================================================================
# TODO to be moved
class BoundaryIntegral(BasicIntegral):

    @classmethod
    def eval(cls, *_args):
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]

        return cls(expr, evaluate=False)

#==============================================================================
class BasicForm(Expr):

    # TODO use .atoms
    @property
    def fields(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, (Field, VectorField))]
        # no redanduncy
        return sorted(list(set(ls)))

    # TODO use .atoms
    @property
    def constants(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Constant)]
        # no redanduncy
        return list(set(ls))

#==============================================================================
class LinearForm(BasicForm):
    is_Function = True

    def __new__(cls, arguments, expr):
        if not isinstance(expr, LinearExpr):
            expr = LinearExpr(arguments, expr)

        assert(isinstance(expr, LinearExpr))

        obj = Basic.__new__(cls, BasicIntegral(expr))

        obj._expr = expr

        return obj

    @property
    def variables(self):
        return self.expr.variables

    @property
    def expr(self):
        return self._expr

    def __call__(self, *args):
        args = Tuple(*args)

        return self.expr.xreplace(dict(list(zip(self.variables, args))))

#==============================================================================
class BilinearForm(BasicForm):

    def __new__(cls, arguments, expr):
        if not isinstance(expr, BilinearExpr):
            expr = BilinearExpr(arguments, expr)

        assert(isinstance(expr, BilinearExpr))

        variables = expr.variables
        obj = Basic.__new__(cls, BasicIntegral(expr))

        obj._expr = expr

        return obj

    @property
    def variables(self):
        return self.expr.variables

    @property
    def expr(self):
        return self._expr

    def __call__(self, *args):
        assert(len(args) == 2)

        left,right = args
        if not is_sequence(left):
            left = [left]

        if not is_sequence(right):
            right = [right]

        args = Tuple(*left, *right)

        variables = Tuple(*self.variables[0], *self.variables[1])

        return self.expr.xreplace(dict(list(zip(variables, args))))

