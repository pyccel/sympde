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
from sympde.topology.space import BasicFunctionSpace
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

from .errors import UnconsistentLinearExpressionError
from .basic  import BasicExpr, Call
from .basic  import is_linear_form, _sanitize_arguments

#==============================================================================
class LinearExpr(BasicExpr):
    is_linear = True

    def __new__(cls, arguments, expr, check=False):
        # ... this is a hack to avoid having expressions with integrals, a
        #     consequence of using xreplace
        integral = expr.atoms(BasicIntegral)
        if integral:
            expr = expr._args[0].expr
        # ...

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

    def __call__(self, *args, evaluate=False):
        return Call(self, args, evaluate=evaluate)

#==============================================================================
class BilinearExpr(BasicExpr):
    is_bilinear = True

    def __new__(cls, arguments, expr, check=False):
        # ...
        if not isinstance(arguments, (tuple, list, Tuple)):
            raise TypeError('(test, trial) must be a tuple, list or Tuple')

        if not(len(arguments) == 2):
            raise ValueError('Expecting a couple (test, trial)')
        # ...

        # ... this is a hack to avoid having expressions with integrals, a
        #     consequence of using xreplace
        integral = expr.atoms(BasicIntegral)
        if integral:
            expr = expr._args[0].expr
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

    def __call__(self, *args, evaluate=False):
        return Call(self, args, evaluate=evaluate)

#==============================================================================
# TODO: - Union of Domain and Boundary
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

        if not isinstance(expr, BasicExpr):
            raise NotImplementedError('')

        if isinstance(expr, BasicExpr) and expr.is_linear:
            if isinstance(expr.expr, Add):
                args = [cls.eval(LinearExpr(expr.variables, a))
                        for a in expr.expr.args]
                return Add(*args)

        elif isinstance(expr, BasicExpr) and expr.is_bilinear:
            if isinstance(expr.expr, Add):
                args = [cls.eval(BilinearExpr(expr.variables, a))
                        for a in expr.expr.args]
                return Add(*args)

        boundary = list(expr.atoms(Boundary))
        if boundary:
            return BoundaryIntegral(expr)

        else:
            return DomainIntegral(expr)

#==============================================================================
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
class BoundaryIntegral(BasicIntegral):

    @classmethod
    def eval(cls, *_args):
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]

        boundary = list(expr.atoms(Boundary))
        if len(boundary) > 1:
            raise ValueError('Expecting one single boundary')

        return cls(expr, evaluate=False)

#==============================================================================
def _get_domain(a):
    # expr is an integral of BasicExpr or Add of Integral of BasicExpr
    if isinstance(a, BoundaryIntegral):
        boundary = list(a.atoms(Boundary))
        return boundary[0]

    elif isinstance(a, DomainIntegral):
        expr = a._args[0]
        variables = expr.variables

        # ...
        if expr.is_linear:
            spaces = [u.space for u in variables]

        elif expr.is_bilinear:
            spaces  = [u.space for u in variables[0]]
            spaces += [u.space for u in variables[1]]

        else:
            raise TypeError('Expecting linear or bilinear form')
        # ...

        domains = [V.domain for V in spaces]
        return domains[0]

    elif isinstance(a, Add):
        domains = [_get_domain(i) for i in a.args]
        domains = list(set(domains))
        if len(domains) == 1:
            return domains[0]

        else:
            return domains

    else:
        raise TypeError('Expecting a Boundary or Domain integral')

