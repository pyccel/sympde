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

from .basic import BasicForm, Call
from .expr   import LinearExpr, BilinearExpr
from .expr import BasicIntegral, DomainIntegral, BoundaryIntegral
from .expr import _get_domain


#==============================================================================
class LinearForm(BasicForm):
    is_linear = True

    def __new__(cls, arguments, expr):

        # ...
        calls = list(expr.atoms(Call))
        for call in calls:
            args = call._args[1]
            call_expr = call._args[0]
            newexpr = call_expr(*args, evaluate=True)
            expr = expr.subs(call, newexpr)

        # take the BasicExpr
        forms = list(expr.atoms(LinearForm))
        for form in forms:
            expr = expr.subs(form, form.expr._args[0])
        # ...

        if not isinstance(expr, LinearExpr):
            expr = LinearExpr(arguments, expr)

        assert(isinstance(expr, LinearExpr))

        variables = expr.variables
        expr = BasicIntegral(expr)
        obj = Basic.__new__(cls, variables, expr)

        # ...
        domain = _get_domain(expr)
        if is_sequence(domain):
            dim = domain[0].dim

        else:
            dim = domain.dim

        obj._domain = domain
        obj._dim = dim
        # ...

        return obj

    @property
    def variables(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

    def __call__(self, *args, evaluate=False):
        return Call(self, args, evaluate=evaluate)

#==============================================================================
class BilinearForm(BasicForm):
    is_bilinear = True

    def __new__(cls, arguments, expr):

        # ...
        calls = list(expr.atoms(Call))
        for call in calls:
            args = call._args[1]
            call_expr = call._args[0]
            newexpr = call_expr(*args, evaluate=True)
            expr = expr.subs(call, newexpr)

        # take the BasicExpr
        forms = list(expr.atoms(BilinearForm))
        for form in forms:
            expr = expr.subs(form, form.expr._args[0])
        # ...

        if not isinstance(expr, BilinearExpr):
            expr = BilinearExpr(arguments, expr)

        assert(isinstance(expr, BilinearExpr))

        variables = expr.variables
        expr = BasicIntegral(expr)
        obj = Basic.__new__(cls, variables, expr)

        # ...
        domain = _get_domain(expr)
        if is_sequence(domain):
            dim = domain[0].dim

        else:
            dim = domain.dim

        obj._domain = domain
        obj._dim = dim
        # ...

        return obj

    @property
    def variables(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

    def __call__(self, *args, evaluate=False):
        return Call(self, args, evaluate=evaluate)

