# coding: utf-8

from itertools import product
from collections import OrderedDict
from sympy import Abs, S, cacheit
from sympy import Indexed, Matrix, ImmutableDenseMatrix
from sympy import expand
from sympy.core import Basic
from sympy.core import Add, Mul, Pow
from sympy.core.expr import AtomicExpr
from sympy.core.containers import Tuple
from sympy.simplify.simplify import simplify

from sympde.core.basic import _coeffs_registery
from sympde.core.basic import CalculusFunction
from sympde.core.algebra import (Dot_1d,
                                 Dot_2d, Inner_2d, Cross_2d,
                                 Dot_3d, Inner_3d, Cross_3d)
from sympde.core.utils import random_string

from sympde.calculus import jump, avg, minus, plus
from sympde.calculus import Jump
from sympde.calculus.core import _generic_ops, _diff_ops

from sympde.calculus.matrices import SymbolicDeterminant, Inverse, Transpose
from sympde.calculus.matrices import MatSymbolicPow, MatrixElement, SymbolicTrace

from sympde.topology.mapping import JacobianSymbol, InterfaceMapping, MultiPatchMapping, JacobianInverseSymbol

from sympde.topology.basic   import BasicDomain, Union, Interval
from sympde.topology.domain  import NormalVector, TangentVector
from sympde.topology.basic   import Boundary, Interface
from sympde.topology.basic   import InteriorDomain
from sympde.topology.mapping import LogicalExpr, PullBack

# TODO fix circular dependency between sympde.expr.evaluation and sympde.topology.mapping

from sympde.topology.space import ScalarTestFunction
from sympde.topology.space import VectorTestFunction
from sympde.topology.space import IndexedTestTrial
from sympde.topology.space import Trace
from sympde.topology.space import element_of
from sympde.topology.space import VectorField

from sympde.topology.derivatives import _partial_derivatives
from sympde.topology.derivatives import _logical_partial_derivatives
from sympde.topology.derivatives import get_atom_derivatives
from sympde.topology.derivatives import get_atom_logical_derivatives
from sympde.topology.derivatives import dx, dy, dz
from sympde.topology.derivatives import dx1, dx2, dx3
from sympde.topology.derivatives import (Grad_1d, Div_1d,
                                         Grad_2d, Curl_2d, Rot_2d, Div_2d,
                                         Grad_3d, Curl_3d, Div_3d)
from sympde.topology.derivatives import Bracket_2d
from sympde.topology.derivatives import Laplace_1d, Laplace_2d, Laplace_3d
from sympde.topology.derivatives import Hessian_1d, Hessian_2d, Hessian_3d

from sympde.topology.derivatives import (LogicalGrad_1d, LogicalDiv_1d,
                                         LogicalGrad_2d, LogicalCurl_2d, LogicalRot_2d, LogicalDiv_2d,
                                         LogicalGrad_3d, LogicalCurl_3d, LogicalDiv_3d)

from sympde.topology.derivatives import LogicalBracket_2d
from sympde.topology.derivatives import LogicalLaplace_1d, LogicalLaplace_2d, LogicalLaplace_3d
from sympde.topology.derivatives import LogicalHessian_1d, LogicalHessian_2d, LogicalHessian_3d

from sympde.topology.space       import ScalarFunctionSpace

from .basic import BasicExpr, BasicForm
from .expr  import BilinearForm
from .expr  import Integral
from .expr  import Functional
from .expr  import _get_domain

import numpy as np
#==============================================================================
def is_sequence(a):
    return isinstance(a, (list,tuple,Tuple))

#==============================================================================
def _unpack_functions(ls):

    funcs = []

    for x in ls:
        if isinstance(x, ScalarTestFunction):
            funcs.append(x)
        elif isinstance(x, VectorTestFunction):
            funcs.extend(x[j] for j in range(x.shape[0]))
        else:
            raise TypeError('Can only accept ScalarTestFunction and VectorTestFunction')

    return tuple(funcs)

#==============================================================================
def _get_trials_tests_flattened(expr):
    """
    Get all scalar trial/test functions in the given expression, with vector
    functions decomposed into their scalar components.

    Parameters
    ----------
    expr : BasicForm | BasicExpr
        Symbolic expression.

    Returns
    -------
    trials : tuple
        All scalar trial functions in the given expression, with vector
        functions decomposed into their scalar components.

    tests : tuple
        All scalar test functions in the given expression, with vector
        functions decomposed into their scalar components.
    """

    if not isinstance(expr, (BasicForm, BasicExpr)):
        raise TypeError("Expression must be of type BasicForm or BasicExpr, got '{}' instead".format(type(expr)))

    if expr.is_bilinear:
        trials = _unpack_functions(expr.variables[0])
        tests  = _unpack_functions(expr.variables[1])

    elif expr.is_linear:
        trials = None
        tests  = _unpack_functions(expr.variables)

    elif expr.is_functional:
        trials = None
        tests  = None

    else:
        ValueError('Could not interpret expression as bilinear form, linear form, or functional')

    return trials, tests

#==============================================================================
def _to_matrix_form(expr, *, trials=None, tests=None):
    """
    Create a matrix representation of input expression, based on trial and
    test functions. We have three options:

    1. if both the trial and test functions are given, we treat the expression
       as a bilinear form and convert it to a (n_rows x n_cols) rectangular
       matrix with n_rows = len(tests) and n_cols = len(trials);

    2. if only the test functions are given, we treat the expression as a
       linear form and convert it to a (n_rows x 1) matrix (column vector) with
       n_rows = len(tests);

    3. if neither the trial nor the test functions are given, we treat the
       expression as a scalar functional and convert it to a 1x1 matrix.

    Parameters
    ----------
    expr : sympy.Expr
        Expression corresponding to a bilinear/linear form or functional.

    trials : iterable
        List of all scalar trial functions (after unpacking vector functions).

    tests : iterable
        List of all scalar test functions (after unpacking vector functions).

    Returns
    -------
    M : sympy.matrices.immutable.ImmutableDenseMatrix
        Matrix representation of input expression.

    """
    # Bilinear form
    if trials and tests:
        M = Matrix.zeros(len(tests), len(trials))
        for i, test in enumerate(tests):
            subs_i = {v:0 for v in tests if v != test}
            expr_i = expr.subs(subs_i)
            for j, trial in enumerate(trials):
                subs_j  = {u:0 for u in trials if u != trial}
                M[i, j] = expr_i.subs(subs_j)

    # Linear form
    elif tests:
        M = Matrix.zeros(len(tests), 1)
        for i, test in enumerate(tests):
            subs_i = {v:0 for v in tests if v != test}
            M[i, 0] = expr.subs(subs_i)

    # Functional
    else:
        M = [[expr]]

    return ImmutableDenseMatrix(M)

#==============================================================================
def _split_expr_over_interface(expr, interface, tests=None, trials=None):
    """
    Splits an expression defined on an interface, into
    expressions where the test and trial functions are defined on each side of
    the interface.

    Parameters:
        expr: sympde expression

        interface: interface of a connectivity

        tests: tests functions as given from linear or bilinear forms

        trials: trials functions as given from linear or bilinear forms

    Returns: sympde expression
    """
    # ...
    is_bilinear = not( trials is None ) and not( tests is None )
    is_linear   =    ( trials is None ) and not( tests is None )

    if trials is None: trials = []
    if tests  is None: tests  = []
    # ...

    int_expressions = []
    bnd_expressions = OrderedDict()

    # ...
    # we replace all jumps
    jumps = expr.atoms(Jump)
    args = [j._args[0] for j in jumps]

    for a in args:
        expr = expr.subs({jump(a): minus(a) - plus(a)})
    # ...

    # ...
    d_trials = OrderedDict()
    for u in trials:
        u_minus = minus(u)
        u_plus  = plus(u)
        d_trials[u] = {'-': u_minus, '+': u_plus}

#        # TODO add sub for avg
#        expr = expr.subs({jump(u): u_minus - u_plus})

    d_tests  = OrderedDict()
    for v in tests:
        v_minus = minus(v)
        v_plus  = plus(v)
        d_tests[v] = {'-': v_minus, '+': v_plus}

#        # TODO add sub for avg
#        expr = expr.subs({jump(v): v_minus - v_plus})
    # ...

    # ...
    trials = []
    for u in d_trials.keys():
        u_minus = d_trials[u]['-']
        u_plus  = d_trials[u]['+']
        trials += [u_minus, u_plus]

    tests = []
    for u in d_tests.keys():
        u_minus = d_tests[u]['-']
        u_plus  = d_tests[u]['+']
        tests += [u_minus, u_plus]
    # ...

    # ...
    def _nullify(expr, u, us):
        """nullifies all symbols in us except u."""
        others = list(set(us) - set([u]))
        for other in others:
            expr = expr.subs({other: 0})
        return expr
    # ...
    if is_bilinear:
        for u in d_trials.keys():
            u_minus = d_trials[u]['-']
            u_plus  = d_trials[u]['+']
            for v in d_tests.keys():
                v_minus = d_tests[v]['-']
                v_plus  = d_tests[v]['+']

                # ...
                newexpr = _nullify(expr, u_minus, trials)
                newexpr = _nullify(newexpr, v_minus, tests)
                newexpr = newexpr.subs({u_minus: u, v_minus: v})
                mapping = newexpr.atoms(InterfaceMapping)
                if mapping and not newexpr.is_zero:
                    mapping = list(mapping)[0]
                    newexpr = newexpr.subs(mapping, mapping.minus)

                if not newexpr.is_zero:
                    bnd_expressions[interface.minus] = newexpr
                # ...
                # TODO must call InterfaceExpression afterward
                newexpr = _nullify(expr, u_minus, trials)
                newexpr = _nullify(newexpr, v_plus, tests)
                mapping = newexpr.atoms(InterfaceMapping)
                if mapping:
                    mapping = list(mapping)[0]
                    newexpr = newexpr.subs(mapping, mapping.plus)
                if not newexpr.is_zero:
                    int_expressions += [InterfaceExpression(interface, u_minus, v_plus, newexpr)]
                # ...
                # TODO must call InterfaceExpression afterward
                newexpr = _nullify(expr, u_plus, trials)
                newexpr = _nullify(newexpr, v_minus, tests)
                mapping = newexpr.atoms(InterfaceMapping)
                if mapping:
                    mapping = list(mapping)[0]
                    newexpr = newexpr.subs(mapping, mapping.plus)
                if not newexpr.is_zero:
                    int_expressions += [InterfaceExpression(interface, u_plus, v_minus, newexpr)]
                # ...
                newexpr = _nullify(expr, u_plus, trials)
                newexpr = _nullify(newexpr, v_plus, tests)
                newexpr = newexpr.subs({u_plus: u, v_plus: v})
                mapping = newexpr.atoms(InterfaceMapping)
                if mapping:
                    mapping = list(mapping)[0]
                    newexpr = newexpr.subs(mapping, mapping.plus)
                if not newexpr.is_zero:
                    bnd_expressions[interface.plus] = newexpr
                # ...

    elif is_linear:
        for v in d_tests.keys():
            v_minus = d_tests[v]['-']
            v_plus  = d_tests[v]['+']

            # ...
            newexpr = _nullify(expr, v_minus, tests)
            newexpr = newexpr.subs({v_minus: v})
            mapping = newexpr.atoms(InterfaceMapping)
            if mapping:
                mapping = list(mapping)[0]
                newexpr = newexpr.subs(mapping, mapping.minus)
            if not newexpr.is_zero:
                bnd_expressions[interface.minus] = newexpr
            # ...

            # ...
            newexpr = _nullify(expr, v_plus, tests)
            newexpr = newexpr.subs({v_plus: v})
            mapping = newexpr.atoms(InterfaceMapping)
            if mapping:
                mapping = list(mapping)[0]
                newexpr = newexpr.subs(mapping, mapping.plus)
            if not newexpr.is_zero:
                bnd_expressions[interface.plus] = newexpr
            # ...

    return int_expressions, bnd_expressions


#==============================================================================
class KernelExpression(Basic):
    def __new__(cls, target, expr):
        if isinstance(expr, (Matrix, ImmutableDenseMatrix)):
            n,m = expr.shape
            if n*m == 1: expr = expr[0,0]

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
class InterfaceExpression(KernelExpression):
    def __new__(cls, target, u, v, expr):
        obj = KernelExpression.__new__(cls, target, expr)
        obj._trial = u
        obj._test  = v
        return obj

    @property
    def test(self):
        return self._test

    @property
    def trial(self):
        return self._trial

#==============================================================================
class TerminalExpr(CalculusFunction):

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args, **options)
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
    @cacheit
    def eval(cls, *_args, **kwargs):
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr    = _args[0]
        n_rows  = kwargs.pop('n_rows', None)
        n_cols  = kwargs.pop('n_cols', None)
        dim     = kwargs.pop('dim', None)
        logical = kwargs.pop('logical', None)

        if isinstance(expr, Add):
            args = [cls.eval(a, dim=dim, logical=logical) for a in expr.args]
            o = args[0]
            for arg in args[1:]:
                o = o + arg
            return o

        elif isinstance(expr, Mul):
            args = [cls.eval(a, dim=dim, logical=logical) for a in expr.args]
            o = args[0]
            for arg in args[1:]:
                o = o * arg
            return o

        elif isinstance(expr, Abs):
            return Abs(cls.eval(expr.args[0], dim=dim, logical=logical))

        elif isinstance(expr, Pow):
            base = cls.eval(expr.base, dim=dim, logical=logical)
            exp  = cls.eval(expr.exp, dim=dim, logical=logical)
            return base**exp

        elif isinstance(expr, JacobianSymbol):
            axis = expr.axis
            J    = expr.mapping.jacobian_expr
            if axis is None:
                return J
            elif expr.mapping.ldim > 1:
                return J.col_del(axis)
            elif expr.mapping.ldim == 1:
                return J.eye(1)

        elif isinstance(expr, JacobianInverseSymbol):
            axis = expr.axis
            J    = expr.mapping.jacobian_inv_expr
            if axis is None:
                return J
            else:
                return J.col_del(axis)

        elif isinstance(expr, SymbolicDeterminant):
            return cls.eval(expr.arg, dim=dim, logical=logical).det().factor()

        elif isinstance(expr, SymbolicTrace):
            return cls.eval(expr.arg, dim=dim, logical=logical).trace()

        elif isinstance(expr, Transpose):
            return cls.eval(expr.arg, dim=dim, logical=logical).T

        elif isinstance(expr, Inverse):
            return cls.eval(expr.arg, dim=dim, logical=logical).inv()

        elif isinstance(expr, ScalarTestFunction):
            return expr

        elif isinstance(expr, VectorTestFunction):
            return ImmutableDenseMatrix([[expr[i]] for i in range(dim)])

        elif isinstance(expr, PullBack):
            return cls.eval(expr.expr, dim=dim, logical=True)

        elif isinstance(expr, MatrixElement):
            base = cls.eval(expr.base, dim=dim, logical=logical)
            return base[expr.indices]

        elif isinstance(expr, BasicForm):
            # ...
            dim     = expr.ldim
            domain  = expr.domain

            if not isinstance(domain, Union):
                logical = domain.mapping is None
                domain  = (domain,)
            # ...
            d_expr = OrderedDict()
            for d in domain:
                d_expr[d] = S.Zero
            # ...
            if isinstance(expr.expr, Add):
                for a in expr.expr.args:
                    newexpr = cls.eval(a, dim=dim, logical=logical)

                    # ...
                    try:
                        domain = _get_domain(a)
                        if not isinstance(domain, Union):
                            domain = (domain, )
                    except:
                        pass
                    # ...
                    for d in domain:
                        d_expr[d] += newexpr
                    # ...
            else:
                newexpr = cls.eval(expr.expr, dim=dim, logical=logical)
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
                    d_expr[d] = newexpr
                # ...

            trials, tests = _get_trials_tests_flattened(expr)

            d_new = OrderedDict()
            for domain, newexpr in d_expr.items():

                if newexpr != 0:

                    # TODO ARA make sure thre is no problem with psydac
                    #      we should always take the interior of a domain
                    if not isinstance(domain, (Boundary, Interface, InteriorDomain)):
                        domain = domain.interior

                    d_new[domain] = _to_matrix_form(newexpr, trials=trials, tests=tests)

            # ...
            ls = []
            d_all = OrderedDict()

            # ... treating interfaces
            keys = [k for k in d_new.keys() if isinstance(k, Interface)]
            for interface in keys:
                # ...
                trials = None
                tests  = None
                if expr.is_bilinear:
                    trials = list(expr.variables[0])
                    tests  = list(expr.variables[1])

                elif expr.is_linear:
                    tests  = list(expr.variables)
                # ...

                # ...
                newexpr = d_new[interface]
                ls_int, d_bnd = _split_expr_over_interface(newexpr, interface,
                                                       tests=tests,
                                                       trials=trials)
                # ...

                ls += ls_int
                # ...
                for k, v in d_bnd.items():
                    if k in d_all.keys():
                        d_all[k] += v

                    else:
                        d_all[k] = v
            # ...

            # ... treating subdomains
            keys = [k for k in d_new.keys() if isinstance(k, Union)]
            for domain in keys:

                newexpr = d_new[domain]
                d       = OrderedDict((interior, newexpr) for interior in domain.as_tuple())
                            # ...
                for k, v in d.items():
                    if k in d_all.keys():
                        d_all[k] += v

                    else:
                        d_all[k] = v

            d = OrderedDict()

            for k, v in d_new.items():
                if not isinstance( k, (Interface, Union)):
                    d[k] = d_new[k]

            for k, v in d_all.items():
                if k in d.keys():
                    d[k] += v
                else:
                    d[k] = v

            d_new = d
            # ...
            for domain, newexpr in d_new.items():
                if isinstance(domain, Boundary):
                    ls += [BoundaryExpression(domain, newexpr)]

                elif isinstance(domain, BasicDomain):
                    ls += [DomainExpression(domain, newexpr)]

                else:
                    raise TypeError('not implemented for {}'.format(type(domain)))
            # ...
            return tuple(ls)

        elif isinstance(expr, Integral):
            dim     = expr.domain.dim if dim is None else dim
            logical = expr.domain.mapping is None
            expr    = cls.eval(expr._args[0], dim=dim, logical=logical)
            return expr

        elif isinstance(expr, NormalVector):
            lines = [[expr[i] for i in range(dim)]]
            return ImmutableDenseMatrix(lines)

        elif isinstance(expr, TangentVector):
            lines = [[expr[i] for i in range(dim)]]
            return ImmutableDenseMatrix(lines)

        elif isinstance(expr, BasicExpr):
            return cls.eval(expr.expr, dim=dim, logical=logical)

        elif isinstance(expr, _diff_ops):
            op   = type(expr)
            if logical:
                new  = eval('Logical{0}_{1}d'.format(op, dim))
            else:
                new  = eval('{0}_{1}d'.format(op, dim))

            args = [cls.eval(i, dim=dim, logical=logical) for i in expr.args]
            return new(*args)
        elif isinstance(expr, _generic_ops):
            # if i = Dot(...) then type(i) is Grad
            op = type(expr)
            new  = eval('{0}_{1}d'.format(op, dim))
            args = [cls.eval(i, dim=dim, logical=logical) for i in expr.args]
            return new(*args)

        elif isinstance(expr, Trace):
            # TODO treate different spaces
            if expr.order == 0:
                return cls.eval(expr.expr, dim=dim, logical=logical)

            elif expr.order == 1:
                # TODO give a name to normal vector
                normal_vector_name = 'n'
                n = NormalVector(normal_vector_name)
                M = cls.eval(expr.expr, dim=dim, logical=logical)

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

        elif isinstance(expr, (Matrix, ImmutableDenseMatrix)):
            n,m = expr.shape
            lines = []
            for i in range(0, n):
                line = []
                for j in range(0, m):
                    line.append(cls.eval(expr[i,j], dim=dim, logical=logical))
                lines.append(line)
            return ImmutableDenseMatrix(lines)

        elif isinstance(expr, LogicalExpr):
            M         = expr.mapping
            expr      = cls(expr.expr, dim=expr.dim)
            return LogicalExpr(expr, mapping=M, dim=M.ldim)
        return expr


#==============================================================================
# TODO use random_string for space name and use it for 1d test function
def _split_test_function(expr):

    if isinstance(expr, ScalarTestFunction):

        dim = expr.space.ldim
        name = expr.name

        ls = []
        for i in range(dim):
            Di = Interval()
            Vi = ScalarFunctionSpace('tmp_V_{}'.format(i), domain=Di)

            ai = ScalarTestFunction(Vi, '{name}_{i}'.format(name=name, i=i+1))
            ls += [ai]

        return {expr:tuple(ls)}
    elif isinstance(expr, VectorTestFunction):

        dim = expr.space.ldim
        name = expr.name

        ls = OrderedDict()
        Di = Interval()
        for i in range(dim):
            Vi = ScalarFunctionSpace('tmp_V_{}'.format(i), domain=Di)
            ls[expr[i]] = tuple(ScalarTestFunction(Vi, '{name}_{i}_{j}'.format(name=name, i=i,j=j+1)) for j in range(dim))

        return ls

    elif isinstance(expr, IndexedTestTrial):

        dim = expr.base.space.ldim
        index = expr.indices[0]
        name = expr.base.name

        Di = Interval()

        ls = []
        for i in range(dim):
            Vi = ScalarFunctionSpace('tmp_V_{}'.format(i), domain=Di)
            ai = ScalarTestFunction(Vi, '{name}_{j}_{i}'.format(name=name, j=index,i=i+1))
            ls += [ai]
        return {expr:tuple(ls)}

    else:
        msg = 'Expecting ScalarTestFunction or IndexedTestTrial, given {}'.format(type(expr))
        raise TypeError(msg)

#==============================================================================
def _tensorize_atomic_expr(expr, d_atoms):

    if (isinstance(expr, _partial_derivatives) or
        isinstance(expr,_logical_partial_derivatives)):

        if isinstance(expr, _partial_derivatives):
            atom = get_atom_derivatives(expr)

        else:
            atom = get_atom_logical_derivatives(expr)

        args = d_atoms[atom]

        op = eval(expr.__class__.__name__)
        expr = _tensorize_atomic_expr(expr._args[0], d_atoms)

        u = args[op.grad_index]
        expr = expr.subs(u, op(u))
        return expr

    elif isinstance(expr, (ScalarTestFunction, IndexedTestTrial)):
        for k,e in d_atoms.items():
            if k == expr:
                atoms = e
                break
        return Mul(*atoms)

    else:
        raise TypeError('{}'.format(type(expr)))

#==============================================================================
class Basic1dForm(AtomicExpr):

    def __new__(cls, name, axis, weight=S.One):
        return Basic.__new__(cls, name, axis, weight)

    @property
    def name(self):
        return self._args[0]

    @property
    def axis(self):
        return self._args[1]

    @property
    def weight(self):
        return self._args[2]

    def _sympystr(self, printer):
        sstr = printer.doprint
        name = sstr(self.name)
        axis = sstr(self.axis)
        name = '{name}_{axis}'.format(name=name, axis=axis)

        if self.weight == S.One:
            return name

        else:
            weight = sstr(self.weight)
            return '{name}({weight})'.format(name=name, weight=weight)


class Mass(Basic1dForm):

    def __new__(cls, axis, weight=S.One):
#        name = 'Mass'
        name = 'M'
        return Basic.__new__(cls, name, axis, weight)

class Stiffness(Basic1dForm):

    def __new__(cls, axis, weight=S.One):
#        name = 'Stiffness'
        name = 'S'
        return Basic.__new__(cls, name, axis, weight)

class Advection(Basic1dForm):

    def __new__(cls, axis, weight=S.One):
#        name = 'Advection'
        name = 'A'
        return Basic.__new__(cls, name, axis, weight)

class AdvectionT(Basic1dForm):

    def __new__(cls, axis, weight=S.One):
#        name = 'AdvectionT'
        name = 'AT'
        return Basic.__new__(cls, name, axis, weight)

class Bilaplacian(Basic1dForm):

    def __new__(cls, axis, weight=S.One):
#        name = 'Bilaplacian'
        name = 'B'
        return Basic.__new__(cls, name, axis, weight)


Mass_0 = Mass(0)
Mass_1 = Mass(1)
Mass_2 = Mass(2)

Stiffness_0 = Stiffness(0)
Stiffness_1 = Stiffness(1)
Stiffness_2 = Stiffness(2)

Advection_0 = Advection(0)
Advection_1 = Advection(1)
Advection_2 = Advection(2)

AdvectionT_0 = AdvectionT(0)
AdvectionT_1 = AdvectionT(1)
AdvectionT_2 = AdvectionT(2)

Bilaplacian_0 = Bilaplacian(0)
Bilaplacian_1 = Bilaplacian(1)
Bilaplacian_2 = Bilaplacian(2)

#==============================================================================
def _replace_atomic_expr(expr, trials, tests, d_atoms, logical=False, expand=False):

    expr   = expr.expand(deep=expand)
    masses = [Mass_0, Mass_1, Mass_2]
    stiffs = [Stiffness_0, Stiffness_1, Stiffness_2]
    advs   = [Advection_0, Advection_1, Advection_2]
    advts  = [AdvectionT_0, AdvectionT_1, AdvectionT_2]
    bils   = [Bilaplacian_0, Bilaplacian_1, Bilaplacian_2]

    if not logical:
        ops = [dx, dy, dz]

    else:
        ops = [dx1, dx2, dx3]

    for u in trials:
        u_atoms = d_atoms[u]
        for v in tests:
            v_atoms = d_atoms[v]

            for i,(ui,vi) in enumerate(zip(u_atoms, v_atoms)):

                mass  = masses[i]
                stiff = stiffs[i]
                adv   = advs[i]
                advt  = advts[i]
                bil   = bils[i]

                d = ops[i]
                # Mass
                expr = expr.subs(vi*ui, mass)

                # Stiffness
                expr = expr.subs(d(vi)*d(ui), stiff)

                # Advection
                expr = expr.subs(vi*d(ui), adv)

                # Transpose of Advection
                expr = expr.subs(d(vi)*ui, advt)

                # Bilaplacian
                expr = expr.subs(d(d(vi))*d(d(ui)), bil)
    return expr

#==============================================================================
class TensorExpr(CalculusFunction):

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args, **options)
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
        d_atoms = kwargs.pop('d_atoms', OrderedDict())
        mapping = kwargs.pop('mapping', None)
        expand    = kwargs.pop('expand', False)

        if isinstance(expr, Add):
            args = [cls.eval(a, d_atoms=d_atoms, mapping=mapping) for a in expr.args]
            return Add(*args)

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if isinstance(a, _coeffs_registery)]
            stests  = [a for a in expr.args if not(a in coeffs) and a.atoms(ScalarTestFunction)]
            vtests  = [a for a in expr.args if not(a in coeffs) and a.atoms(VectorTestFunction)]
            vectors = [a for a in expr.args if (not(a in coeffs) and
                                                not(a in stests) and
                                                not(a in vtests))]

            a = S.One
            if coeffs:
                a = Mul(*coeffs)

            b = S.One
            if vectors:

                args = [cls(i, evaluate=False) for i in vectors]
                b = Mul(*args)

            sb = S.One
            if stests:
                args = [cls.eval(i, d_atoms=d_atoms, mapping=mapping) for i in stests]
                sb = Mul(*args)

            vb = S.One
            if vtests:
                args = [cls.eval(i, d_atoms=d_atoms, mapping=mapping) for i in vtests]
                vb = Mul(*args)

            return Mul(a, b, sb, vb)

        elif isinstance(expr, _coeffs_registery):
            return expr

        elif isinstance(expr, Pow):
            b = expr.base
            e = expr.exp
            return Pow(cls.eval(b), e)

        elif isinstance(expr, (Matrix, ImmutableDenseMatrix)):

            n_rows, n_cols = expr.shape

            lines = []
            for i_row in range(0, n_rows):
                line = []
                for i_col in range(0, n_cols):
                    line.append(cls.eval(expr[i_row,i_col], d_atoms=d_atoms, mapping=mapping))

                lines.append(line)

            return ImmutableDenseMatrix(lines)

        elif isinstance(expr, DomainExpression):
            # TODO to be removed
            target = expr.target
            expr   = expr.expr
            return cls.eval(expr, d_atoms=d_atoms, mapping=mapping)

        elif isinstance(expr, BilinearForm):
            trials = list(expr.variables[0])
            tests  = list(expr.variables[1])

            # ... # TODO improve
            terminal_expr = TerminalExpr(expr)[0]
            # ...

            # ...
            variables  = list(terminal_expr.atoms(ScalarTestFunction))
            variables += list(terminal_expr.atoms(IndexedTestTrial))
            # ...

            # ...
            if not(mapping is None):
                dim           = expr.ldim
                terminal_expr = LogicalExpr(terminal_expr.expr, mapping=mapping, dim=dim)
                variables     = [LogicalExpr(e, mapping=mapping, dim=dim) for e in variables ]
                trials        = [LogicalExpr(e, mapping=mapping, dim=dim) for e in trials ]
                tests         = [LogicalExpr(e, mapping=mapping, dim=dim) for e in tests ]
            # ...

            d_atoms = OrderedDict()
            for a in variables:
                new = _split_test_function(a)
                d_atoms[a] = new[a]

            # ...
            expr = cls.eval(terminal_expr, d_atoms=d_atoms, mapping=mapping)
            # ...

            # ...
            trials = [a for a in variables if ((isinstance(a, ScalarTestFunction) and a in trials) or
                                               (isinstance(a, IndexedTestTrial) and a.base in trials))]

            tests = [a for a in variables if ((isinstance(a, ScalarTestFunction) and a in tests) or
                                              (isinstance(a, IndexedTestTrial) and a.base in tests))]
            # ...

            expr = _replace_atomic_expr(expr, trials, tests, d_atoms,
                                        logical=True, expand=expand)

            return expr

        if expr.atoms(ScalarTestFunction) or expr.atoms(IndexedTestTrial):
            return _tensorize_atomic_expr(expr, d_atoms)

        return cls(expr, evaluate=False)
