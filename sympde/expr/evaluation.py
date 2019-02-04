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

from .basic  import BasicExpr, BasicForm
from .expr   import LinearExpr, BilinearExpr
from .expr   import LinearForm, BilinearForm
from .expr import BasicIntegral, DomainIntegral, BoundaryIntegral
from .expr import Functional
from .expr import _get_domain

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

        return  M, test_indices, trial_indices

    elif expr.is_linear:

        tests = list(expr.variables)
        n_rows, test_indices = _get_size_and_starts(tests)

        # ...
        lines = [0 for i in range(0, n_rows)]
        M = Matrix(lines)
        # ...

        return  M, test_indices, None

    elif expr.is_functional:

        # ...
        M = Matrix([0.])
        # ...

        return  M, None, None

    else:
        raise TypeError('Expecting BasicExpr or BasicForm')

    return M

def _to_matrix_bilinear_form(expr, M, test_indices, trial_indices):

    # ...
    def treat_form(arg, M):
        atoms  = list(arg.atoms(TestFunction))
        atoms += list(arg.atoms(VectorTestFunction))
        atoms += list(arg.atoms(IndexedTestTrial))

        for atom in atoms:
            if atom in test_indices:
                i_row = test_indices[atom]

            elif atom in trial_indices:
                i_col = trial_indices[atom]

            else:
                raise ValueError('> Could not find {}'.format(atom))

        M[i_row, i_col] += arg
        return M
    # ...

    # ...
    if isinstance(expr, Add):
        args = expr.args
        for arg in args:
            if isinstance(arg, Mul):
                M = treat_form(arg, M)

    elif isinstance(expr, Mul):
        M = treat_form(expr, M)

    else:
        raise TypeError('> wrong type, given {}'.format(type(expr)))
    # ...

    return M

def _to_matrix_linear_form(expr, M, test_indices):
    # ...
    def treat_form(arg, M):
        atoms  = list(arg.atoms(TestFunction))
        atoms += list(arg.atoms(VectorTestFunction))
        atoms += list(arg.atoms(IndexedTestTrial))

        for atom in atoms:
            if atom in test_indices:
                i_row = test_indices[atom]

            else:
                raise ValueError('> Could not find {}'.format(atom))

        M[i_row] += arg
        return M
    # ...

    # ...
    if isinstance(expr, Add):
        args = expr.args
        for arg in args:
            M = treat_form(arg, M)

    elif isinstance(expr, Mul):
        M = treat_form(expr, M)

    else:
        raise TypeError('> wrong type, given {}'.format(type(expr)))
    # ...

    return M

def _to_matrix_functional_form(expr, M):
    M[0] += expr

    return M

def _to_matrix_form(expr, M, test_indices, trial_indices):
    if not(test_indices is None) and not(trial_indices is None):
        return _to_matrix_bilinear_form(expr, M, test_indices, trial_indices)

    if not(test_indices is None) and trial_indices is None:
        return _to_matrix_linear_form(expr, M, test_indices)

    if test_indices is None and trial_indices is None:
        return _to_matrix_functional_form(expr, M)

#==============================================================================
class KernelExpression(Basic):
    def __new__(cls, target, expr):
        return Basic.__new__(cls, target, expr)

    @property
    def target(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

#==============================================================================
class DomainExpression(KernelExpression):
    pass

#==============================================================================
class BoundaryExpression(KernelExpression):
    pass


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
            args = [cls.eval(a, dim=dim) for a in expr.args]
            return Add(*args)

        elif isinstance(expr, Mul):
            args = [cls.eval(a, dim=dim) for a in expr.args]
            return Mul(*args)

        elif isinstance(expr, BasicForm):
            # ...
            dim = expr.dim
            is_bilinear = expr.is_bilinear

            domain = expr.domain
            if isinstance(domain, Union):
                domain = list(domain._args)

            elif not is_sequence(domain):
                domain = [domain]
            # ...

            # ...
            d_expr = {}
            for d in domain:
                d_expr[d] = S.Zero
            # ...

            # ...
            if isinstance(expr.expr, Add):
                for a in expr.expr.args:
                    newexpr = cls.eval(a, dim=dim)
                    newexpr = expand(newexpr)

                    # ...
                    domain = _get_domain(a)
                    if isinstance(domain, Union):
                        domain = list(domain._args)

                    elif not is_sequence(domain):
                        domain = [domain]
                    # ...

                    # ...
                    for d in domain:
                        d_expr[d] += newexpr
                    # ...

            else:
                newexpr = cls.eval(expr.expr, dim=dim)
                newexpr = expand(newexpr)

                # ...
                if isinstance(expr, Functional):
                    domain = expr.domain

                else:
                    domain = _get_domain(expr.expr)

                if isinstance(domain, Union):
                    domain = list(domain._args)

                elif not is_sequence(domain):
                    domain = [domain]
                # ...

                # ...
                for d in domain:
                    d_expr[d] += newexpr
                # ...
            # ...

            # ...
            d_new = {}
            for domain, newexpr in d_expr.items():

                M, test_indices, trial_indices = _init_matrix(expr)
                M = _to_matrix_form(newexpr, M, test_indices, trial_indices)

                n,m = M.shape
                if n*m == 1: M = M[0,0]

                d_new[domain] = M
            # ...

            # ...
            ls = []
            for domain, newexpr in d_new.items():
                if isinstance(domain, Boundary):
                    ls += [BoundaryExpression(domain, newexpr)]

                elif isinstance(domain, BasicDomain):
                    ls += [DomainExpression(domain, newexpr)]
                else:
                    raise TypeError('')
            # ...

            return ls

        elif isinstance(expr, BasicIntegral):
            if dim is None:
                domain = _get_domain(expr)
                dim = domain.dim

            return cls.eval(expr._args[0], dim=dim)

        elif isinstance(expr, BasicExpr):
            return cls.eval(expr.expr, dim=dim)

        elif isinstance(expr, _generic_ops):
            # if i = Dot(...) then type(i) is Grad
            op = type(expr)
            new  = eval('{0}_{1}d'.format(op, dim))

            args = [cls.eval(i, dim=dim) for i in expr.args]
            return new(*args)

        elif isinstance(expr, Trace):
            # TODO treate different spaces
            if expr.order == 0:
                return cls.eval(expr.expr, dim=dim)

            elif expr.order == 1:
                # TODO give a name to normal vector
                normal_vector_name = 'n'
                n = NormalVector(normal_vector_name)
                M = cls.eval(expr.expr, dim=dim)
                if dim == 1:
                    return M
                else:
                    if isinstance(M, (Add, Mul)):
                        ls = M.atoms(Tuple)
                        for i in ls:
                            M = M.subs(i, Matrix(i))
                        M = simplify(M)

                    e = 0
                    for i in range(0, dim):
                        e += M[i] * n[i]
                    return e

            else:
                raise ValueError('> Only traces of order 0 and 1 are available')


        elif isinstance(expr, Matrix):
            n,m = expr.shape
            lines = []
            for i in range(0, n):
                line = []
                for j in range(0, m):
                    line.append(cls.eval(expr[i,j], dim=dim))
                lines.append(line)
            return Matrix(lines)

        return expr
