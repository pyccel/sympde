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
    is_Function = True
    is_linear   = False
    is_bilinear = False

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
# TODO: - to be moved
#       - Union of Domain and Boundary
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

        if isinstance(expr, LinearExpr):
            if isinstance(expr.expr, Add):
                args = [cls.eval(LinearExpr(expr.variables, a))
                        for a in expr.expr.args]
                return Add(*args)

        elif isinstance(expr, BilinearExpr):
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

        boundary = list(expr.atoms(Boundary))
        if len(boundary) > 1:
            raise ValueError('Expecting one single boundary')

        return cls(expr, evaluate=False)

#==============================================================================
# TODO check unicity of domain in __new__
class BasicForm(Expr):
    is_Function = True
    is_linear   = False
    is_bilinear = False

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

    @property
    def domain(self):
        if isinstance(self.expr, BoundaryIntegral):
            boundary = list(self.expr.atoms(Boundary))
            return boundary[0]

        elif isinstance(self.expr, DomainIntegral):
            # ...
            if self.is_linear:
                spaces = [u.space for u in self.variables]

            elif self.is_bilinear:
                spaces  = [u.space for u in self.variables[0]]
                spaces += [u.space for u in self.variables[1]]

            else:
                raise TypeError('Expecting linear or bilinear form')
            # ...

            domains = [V.domain for V in spaces]
            return domains[0]

        else:
            raise TypeError('Expecting a Boundary or Domain integral')

#==============================================================================
class LinearForm(BasicForm):
    is_linear = True

    def __new__(cls, arguments, expr):
        if not isinstance(expr, LinearExpr):
            expr = LinearExpr(arguments, expr)

        assert(isinstance(expr, LinearExpr))

        variables = expr.variables
        expr = BasicIntegral(expr)
        return Basic.__new__(cls, variables, expr)

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
        if not isinstance(expr, BilinearExpr):
            expr = BilinearExpr(arguments, expr)

        assert(isinstance(expr, BilinearExpr))

        variables = expr.variables
        expr = BasicIntegral(expr)
        return Basic.__new__(cls, variables, expr)

    @property
    def variables(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

    def __call__(self, *args, evaluate=False):
        return Call(self, args, evaluate=evaluate)

#==============================================================================
# default behavior is not to evaluate the call
class Call(CalculusFunction):

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', False):
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

        if not len(_args) == 2:
            raise ValueError('Expecting two arguments')

        expr = _args[0]
        args = _args[1]

        is_linear   = isinstance(expr, (LinearExpr, LinearForm))
        is_bilinear = isinstance(expr, (BilinearExpr, BilinearForm))

        args = _sanitize_arguments(args,
                                   is_linear=is_linear,
                                   is_bilinear=is_bilinear)

        if is_linear:
            args = Tuple(*args)
            variables = expr.variables
            return expr.xreplace(dict(list(zip(variables, args))))

        if is_bilinear:
            left,right = args
            if not is_sequence(left):
                left = [left]

            if not is_sequence(right):
                right = [right]

            args = Tuple(*left, *right)

            variables = expr.variables
            variables = Tuple(*variables[0], *variables[1])

            return expr.xreplace(dict(list(zip(variables, args))))

        return cls(expr, args, evaluate=False)

#==============================================================================
def _get_size_and_starts(ls):
    n = 0
    d_indices = {}
    for x in ls:
        d_indices[x] = n
        if isinstance(x, TestFunction):
            n += 1

        elif isinstance(x, VectorTestFunction):
            for j in range(0, x.shape[0]):
                d_indices[x[j]] = n + j

            n += x.shape[0]

    return n, d_indices

def _init_matrix(expr):
    assert(isinstance(expr, (BasicForm, BasicExpr)))

    if expr.is_bilinear:

        trials = list(expr.variables[0])
        n_rows, trial_indices = _get_size_and_starts(trials)

        tests = list(expr.variables[1])
        n_cols, test_indices = _get_size_and_starts(tests)

        # ...
        lines = []
        for i in range(0, n_rows):
            line = []
            for j in range(0, n_cols):
                line.append(0)
            lines.append(line)

        M = Matrix(lines)
        # ...

    elif expr.is_linear:

        tests = list(expr.variables[0])
        n_rows, test_indices = _get_size_and_starts(tests)

        # ...
        lines = [0 for i in range(0, n_rows)]
        M = Matrix(lines)
        # ...

    else:
        raise TypeError('Expecting BasicExpr or BasicForm')

    return M

#==============================================================================
class TerminalExpr(CalculusFunction):

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
    def eval(cls, *_args, **kwargs):
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]

        n_rows = kwargs.pop('n_rows', None)
        n_cols = kwargs.pop('n_cols', None)
        dim    = kwargs.pop('dim', None)

        if isinstance(expr, Add):
            args = [cls.eval(a) for a in expr.args]
            return Add(*args)

        elif isinstance(expr, Mul):
            args = [cls.eval(a) for a in expr.args]
            return Mul(*args)

        elif isinstance(expr, BasicForm):
            dim = expr.domain.dim
            return cls.eval(expr.expr, dim=dim)

        elif isinstance(expr, BasicIntegral):
            return cls.eval(expr._args[0], dim=dim)

        elif isinstance(expr, BasicExpr):
            M = _init_matrix(expr)
            return cls.eval(expr.expr, dim=dim)

        elif isinstance(expr, _generic_ops):
            # if i = Dot(...) then type(i) is Grad
            op = type(expr)
            new  = eval('{0}_{1}d'.format(op, dim))

            args = [cls.eval(i, dim=dim) for i in expr.args]
            return new(*args)

        elif isinstance(expr, Matrix):
            n,m = expr.shape
            lines = []
            for i in range(0, n):
                line = []
                for j in range(0, m):
                    line.append(cls.eval(expr[i,j], dim=dim))
                lines.append(line)
            return Matrix(lines)

#        elif isinstance(expr, _coeffs_registery):
#            return expr
#
#        elif isinstance(expr, (TestFunction, VectorTestFunction, IndexedTestTrial,
#                               Field, VectorField, IndexedVectorField)):
#
#            return expr
#
#        return cls(expr, evaluate=False)

        return expr

