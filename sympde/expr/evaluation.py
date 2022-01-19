# coding: utf-8

from itertools import product


import numpy as np

from sympy import Abs, S, cacheit
from sympy import Indexed, Matrix, ImmutableDenseMatrix
from sympy import expand
from sympy.core import Basic, Symbol
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
from sympde.calculus import Jump, is_zero
from sympde.calculus.core import _generic_ops, _diff_ops
from sympde.calculus.matrices import SymbolicDeterminant, Inverse, Transpose
from sympde.calculus.matrices import MatSymbolicPow, MatrixElement, SymbolicTrace

from sympde.topology.basic   import BasicDomain, Union, Interval
from sympde.topology.basic   import Boundary, Interface
from sympde.topology.basic   import InteriorDomain
from sympde.topology.domain  import NormalVector, TangentVector, NCube, NCubeInterior
from sympde.topology.mapping import JacobianSymbol, InterfaceMapping, MultiPatchMapping, JacobianInverseSymbol
from sympde.topology.mapping import LogicalExpr, PullBack

# TODO fix circular dependency between sympde.expr.evaluation and sympde.topology.mapping

from sympde.topology.space import ScalarFunction
from sympde.topology.space import VectorFunction
from sympde.topology.space import IndexedVectorFunction
from sympde.topology.space import Trace
from sympde.topology.space import element_of
from sympde.topology.space import ScalarFunctionSpace

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

from .basic import BasicExpr, BasicForm
from .expr  import BilinearForm
from .expr  import Integral
from .expr  import Functional
from .expr  import _get_domain

__all__ = (
    'Advection',
    'AdvectionT',
    'Basic1dForm',
    'Bilaplacian',
    'BoundaryExpression',
    'DomainExpression',
    'InterfaceExpression',
    'KernelExpression',
    'Mass',
    'Stiffness',
    'TensorExpr',
    'TerminalExpr',
    '_get_trials_tests',
    '_replace_atomic_expr',
    '_split_expr_over_interface',
    '_split_test_function',
    '_tensorize_atomic_expr',
    '_to_matrix_form',
    '_unpack_functions',
    'is_sequence',
)

#==============================================================================
def is_sequence(a):
    return isinstance(a, (list,tuple,Tuple))

#==============================================================================
def _unpack_functions(ls):

    funcs = []

    for x in ls:
        if isinstance(x, ScalarFunction):
            funcs.append(x)
        elif isinstance(x, VectorFunction):
            funcs.extend(x[j] for j in range(x.shape[0]))
        else:
            raise TypeError('Can only accept ScalarFunction and VectorFunction')

    return tuple(funcs)

#==============================================================================
def _get_trials_tests(expr, *, flatten=False):
    """
    Get all scalar trial/test functions in the given expression.

    Parameters
    ----------
    expr : BasicForm | BasicExpr
        Symbolic expression.

    flatten: Boolean, optional
        Decompose the vector trial/test functions into their scalar components.

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
        trials = _unpack_functions(expr.variables[0]) if flatten else expr.variables[0]
        tests  = _unpack_functions(expr.variables[1]) if flatten else expr.variables[1]

    elif expr.is_linear:
        trials = None
        tests  = _unpack_functions(expr.variables) if flatten else expr.variables

    elif expr.is_functional:
        trials = None
        tests  = None

    else:
        ValueError('Could not interpret expression as bilinear form, linear form, or functional')

    return trials, tests

#==============================================================================
def _to_matrix_form(expr, *, trials=None, tests=None, domain=None):
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

    The domain is needed when the expression is defined over an interface,
    to delete the plus/minus operators and the InterfaceMapping object.

    Parameters
    ----------
    expr : sympy.Expr
        Expression corresponding to a bilinear/linear form or functional.

    trials : iterable
        List of all scalar trial functions (after unpacking vector functions).

    tests : iterable
        List of all scalar test functions (after unpacking vector functions).

    domain : Domain
        domain of expr

    Returns
    -------
    M : sympy.matrices.immutable.ImmutableDenseMatrix
        Matrix representation of input expression.

    """

    if not isinstance(domain, Interface):
        atoms     = expr.atoms(minus, plus)
        new_atoms = [e.args[0] for e in atoms]
        subs      = tuple(zip(atoms, new_atoms))
        expr      = expr.subs(subs)
        mapping   = expr.atoms(InterfaceMapping)
        if mapping:
            mapping = list(mapping)[0]
            expr    = expr.subs(mapping, mapping.minus)

    # Bilinear form
    if trials and tests:
        M = [[None for j in trials] for i in tests]
        for i, test in enumerate(tests):
            subs_i = {v:0 for v in tests if v != test}
            expr_i = expr.subs(subs_i)
            for j, trial in enumerate(trials):
                subs_j  = {u:0 for u in trials if u != trial}
                M[i][j] = expr_i.subs(subs_j)
        M = Matrix(M)

    # Linear form
    elif tests:
        M = [[None] for i in tests]
        for i, test in enumerate(tests):
            subs_i = {v:0 for v in tests if v != test}
            M[i][0] = expr.subs(subs_i)
        M = Matrix(M)
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

    int_expressions = {}
    bnd_expressions = {}

    # ...
    # we replace all jumps
    jumps = expr.atoms(Jump)
    args = [j._args[0] for j in jumps]

    for a in args:
        expr = expr.subs({jump(a): minus(a) - plus(a)})
    # ...

    # ...
    d_trials = {}
    for u in trials:
        u_minus = minus(u)
        u_plus  = plus(u)
        d_trials[u] = {'-': u_minus, '+': u_plus}

#        # TODO add sub for avg
#        expr = expr.subs({jump(u): u_minus - u_plus})

    d_tests  = {}
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
            for v in d_tests.keys():
                u_minus = d_trials[u]['-']
                u_plus  = d_trials[u]['+']
                v_minus = d_tests[v]['-']
                v_plus  = d_tests[v]['+']

                # ...
                newexpr = _nullify(expr, u_minus, trials)
                newexpr = _nullify(newexpr, v_minus, tests)
                newexpr = newexpr.subs({u_minus: u, v_minus: v})
                mapping = newexpr.atoms(InterfaceMapping)
                if mapping and not is_zero(newexpr):
                    mapping = list(mapping)[0]
                    newexpr = newexpr.subs(mapping, mapping.minus)

                if not is_zero(newexpr):
                    if interface.minus in bnd_expressions:
                        newexpr += bnd_expressions[interface.minus]
                    bnd_expressions[interface.minus] = newexpr.subs(interface, interface.minus)

                # ...
                newexpr = _nullify(expr, u_plus, trials)
                newexpr = _nullify(newexpr, v_plus, tests)
                newexpr = newexpr.subs({u_plus: u, v_plus: v})
                mapping = newexpr.atoms(InterfaceMapping)

                if mapping and not is_zero(newexpr):
                    mapping = list(mapping)[0]
                    newexpr = newexpr.subs(mapping, mapping.plus)

                for nn in newexpr.atoms(NormalVector):
                    newexpr = newexpr.subs(nn, -nn)

                if not is_zero(newexpr):
                    if interface.plus in bnd_expressions:
                        newexpr += bnd_expressions[interface.plus]
                    bnd_expressions[interface.plus] = newexpr.subs(interface, interface.plus)

                # ...
                # TODO must call InterfaceExpression afterward
                newexpr = _nullify(expr, u_minus, trials)
                newexpr = _nullify(newexpr, v_plus, tests)
                mapping = newexpr.atoms(InterfaceMapping)
                if mapping:
                    mapping = list(mapping)[0]
                    for det in newexpr.atoms(SymbolicDeterminant):
                        if det.atoms(InterfaceMapping):
                            newdet = det.subs(mapping, mapping.minus)
                            newexpr = newexpr.subs(det, newdet)

                if not is_zero(newexpr):
                    if isinstance(u, IndexedVectorFunction):
                        u_minus = minus(u.base)
                    if isinstance(v, IndexedVectorFunction):
                        v_plus = plus(v.base)
                    if (u_minus, v_plus) in int_expressions:
                        newexpr += int_expressions[u_minus, v_plus].expr

                    int_expressions[u_minus, v_plus] = InterfaceExpression(interface, u_minus, v_plus, newexpr)

                # ...
                # TODO must call InterfaceExpression afterward
                newexpr = _nullify(expr, u_plus, trials)
                newexpr = _nullify(newexpr, v_minus, tests)
                mapping = newexpr.atoms(InterfaceMapping)
                if mapping:
                    mapping = list(mapping)[0]
                    for det in newexpr.atoms(SymbolicDeterminant):
                        if det.atoms(InterfaceMapping):
                            newdet = det.subs(mapping, mapping.minus)
                            newexpr = newexpr.subs(det, newdet)
                if not is_zero(newexpr):
                    if isinstance(u, IndexedVectorFunction):
                        u_plus = plus(u.base)
                    if isinstance(v, IndexedVectorFunction):
                        v_minus = minus(v.base)
                    if (u_plus, v_minus) in int_expressions:
                        newexpr += int_expressions[u_plus, v_minus].expr

                    int_expressions[u_plus, v_minus] = InterfaceExpression(interface, u_plus, v_minus, newexpr)

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

            if not is_zero(newexpr):
                if interface.minus in bnd_expressions:
                    newexpr += bnd_expressions[interface.minus]

                bnd_expressions[interface.minus] = newexpr
            # ...
            # ...
            newexpr = _nullify(expr, v_plus, tests)
            newexpr = newexpr.subs({v_plus: v})
            mapping = newexpr.atoms(InterfaceMapping)
            if mapping:
                mapping = list(mapping)[0]
                newexpr = newexpr.subs(mapping, mapping.plus)

            if not is_zero(newexpr):
                if interface.plus in bnd_expressions:
                     newexpr += bnd_expressions[interface.plus]

                bnd_expressions[interface.plus] = newexpr
            # ...

    int_expressions = list(int_expressions.values())
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
    """
     This class takes a SymPDE expression and transforms its vector operators into atomic ones
     for a specified domain.

     Parameters
     ----------
     expr: sympy.Expr
      The mathematical expression.

     domain: Domain
      The domain in which the expression is defined.

     Returns
     -------
      r: sympy.Expr or tuple of sympy.Expr
       The atomized expression.
    """
    def __new__(cls, expr, domain, **options):

        if options.pop('evaluate', True):
            r = cls.eval(expr, domain)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, expr, domain)
        else:
            return r

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

    @property
    def expr(self):
        return self._args[0]

    @property
    def domain(self):
        return self._args[1]

    @classmethod
    @cacheit
    def eval(cls, expr, domain):
        """."""

        dim = domain.dim
        if isinstance(expr, Add):
            args = [cls.eval(a, domain=domain) for a in expr.args]
            o = args[0]
            for arg in args[1:]:
                o = o + arg
            return o

        elif isinstance(expr, Mul):
            args = [cls.eval(a, domain=domain) for a in expr.args]
            o = args[0]
            for arg in args[1:]:
                o = o * arg
            return o

        elif isinstance(expr, Abs):
            return Abs(cls.eval(expr.args[0], domain=domain))

        elif isinstance(expr, Pow):
            base = cls.eval(expr.base, domain=domain)
            exp  = cls.eval(expr.exp, domain=domain)
            return base**exp

        elif isinstance(expr, JacobianSymbol):
            axis = expr.axis
            J    = expr.mapping.jacobian_expr
            if not axis is None:
                if expr.mapping.ldim > 1:
                    J = J.col_del(axis)
                if expr.mapping.ldim == 1:
                    J = J.eye(1)
            if isinstance(domain, Interface):
                mapping = expr.mapping
                assert mapping.is_plus is not mapping.is_minus
                if mapping.is_plus and mapping.is_analytical:
                    domain      = domain.plus
                    axis        = domain.axis
                    ext         = domain.ext
                    domain      = domain.domain
                    coordinates = domain.coordinates
                    if isinstance(domain, (NCube, NCubeInterior)):
                        bounds      = domain.min_coords if ext == -1 else domain.max_coords
                        J           = J.subs(coordinates[axis], bounds[axis])
                    newcoords = [Symbol(c.name+"_plus") for c in coordinates]
                    subs      = list(zip(coordinates, newcoords))
                    J = J.subs(subs)

            return J

        elif isinstance(expr, JacobianInverseSymbol):
            axis = expr.axis
            J    = expr.mapping.jacobian_inv_expr
            if not axis is None:
                J = J.col_del(axis)

            if isinstance(domain, Interface):
                mapping = expr.mapping
                assert mapping.is_plus is not mapping.is_minus
                if mapping.is_plus and mapping.is_analytical:
                    domain      = domain.plus
                    axis        = domain.axis
                    ext         = domain.ext
                    domain      = domain.domain
                    coordinates = domain.coordinates
                    if isinstance(domain, (NCube, NCubeInterior)):
                        bounds      = domain.min_coords if ext == -1 else domain.max_coords
                        J = J.subs(coordinates[axis], bounds[axis])
                    newcoords = [Symbol(c.name+"_plus") for c in coordinates]
                    subs      = list(zip(coordinates, newcoords))
                    J = J.subs(subs)

            return J

        elif isinstance(expr, SymbolicDeterminant):
            return cls.eval(expr.arg, domain=domain).det().factor()

        elif isinstance(expr, SymbolicTrace):
            return cls.eval(expr.arg, domain=domain).trace()

        elif isinstance(expr, Transpose):
            return cls.eval(expr.arg, domain=domain).T

        elif isinstance(expr, Inverse):
            return cls.eval(expr.arg, domain=domain).inv()

        elif isinstance(expr, ScalarFunction):
            return expr

        elif isinstance(expr, VectorFunction):
            return ImmutableDenseMatrix([[expr[i]] for i in range(dim)])

        elif isinstance(expr, (minus, plus)):
            newexpr = cls.eval(expr.args[0], domain=domain)
            return type(expr)(newexpr)

        elif isinstance(expr, PullBack):
            return cls.eval(expr.expr, domain=domain)

        elif isinstance(expr, MatrixElement):
            base = cls.eval(expr.base, domain=domain)
            return base[expr.indices]

        elif isinstance(expr, BasicForm):
            # ...
            dim     = expr.ldim
            domains  = expr.domain

            if not isinstance(domains, Union):
                domains  = (domains,)
            # ...
            d_expr = {}
            d_int  = {}
            for d in domains:
                d_expr[d] = S.Zero
            # ...
            if isinstance(expr.expr, Add):
                for a in expr.expr.args:
                    domains = _get_domain(a)
                    if not isinstance(domains, Union):
                        domains = (domains, )

                    # ...
                    for d in domains:
                        d_expr[d] += a
                    # ...
            else:
                # ...
                if isinstance(expr, Functional):
                    domains = expr.domain

                else:
                    domains = _get_domain(expr.expr)

                if isinstance(domains, Union):
                    domains = list(domains._args)

                elif not is_sequence(domains):
                    domains = (domains,)
                # ...

                # ...
                for d in domains:
                    d_expr[d] = expr.expr
                # ...

            trials, tests = _get_trials_tests(expr, flatten=False)
            # ... treating interfaces
            keys = [k for k in d_expr.keys() if isinstance(k, Interface)]
            for interface in keys:
                newexpr = d_expr.pop(interface)
                ls_int, d_bnd = _split_expr_over_interface(newexpr, interface,
                                                       tests=tests,
                                                       trials=trials)
                # ...
                for a in ls_int:
                    if (a.target, a.trial, a.test) in d_int:
                        d_int[a.target, a.trial, a.test] += a.expr
                    else:
                        d_int[a.target, a.trial, a.test] = a.expr

                for k, v in d_bnd.items():
                    if k in d_expr.keys():
                        d_expr[k] += v

                    else:
                        d_expr[k] = v
            # ...
            trials, tests = _get_trials_tests(expr, flatten=True)
            d_new = {}
            for domain, a in d_expr.items():
                newexpr = cls.eval(a, domain=domain)
                if newexpr != 0:
                    # TODO ARA make sure thre is no problem with psydac
                    #      we should always take the interior of a domain
                    if not isinstance(domain, (Boundary, Interface, InteriorDomain)):
                        domain = domain.interior
                    d_new[domain] = _to_matrix_form(newexpr, trials=trials, tests=tests, domain=domain)

            if len(d_new) == 0:
                # ... corner case where the expression is zero
                for domain, a in d_expr.items():
                    if not isinstance(domain, (Boundary, Interface, InteriorDomain)):
                        domain = domain.interior
                    d_new[domain] = _to_matrix_form(S.Zero, trials=trials, tests=tests, domain=domain)

                    if isinstance(domain, Boundary):
                        return (BoundaryExpression(domain, d_new[domain]),)
                    elif isinstance(domain, BasicDomain):
                        return (DomainExpression(domain, d_new[domain]),)

            # ... treating subdomains
            keys = [k for k in d_new.keys() if isinstance(k, Union)]
            for domain in keys:
                newexpr = d_new.pop(domain)
                d       = dict((interior, newexpr) for interior in domain.as_tuple())
                            # ...
                for k, v in d.items():
                    if k in d_new.keys():
                        d_new[k] += v

                    else:
                        d_new[k] = v

            ls = []
            for key in d_int:
                domain, u,v = key
                expr    = d_int[domain, u, v]

                newexpr = cls.eval(expr, domain=domain)
                if newexpr != 0:
                    newexpr = _to_matrix_form(newexpr, trials=trials, tests=tests, domain=domain)
                    ls += [InterfaceExpression(domain, u, v, newexpr)]

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
            domain = expr.domain
            expr   = cls.eval(expr._args[0], domain=domain)
            return expr

        elif isinstance(expr, NormalVector):
            lines = [[expr[i] for i in range(dim)]]
            return ImmutableDenseMatrix(lines)

        elif isinstance(expr, TangentVector):
            lines = [[expr[i] for i in range(dim)]]
            return ImmutableDenseMatrix(lines)

        elif isinstance(expr, BasicExpr):
            return cls.eval(expr.expr, domain=domain)

        elif isinstance(expr, _diff_ops):
            op   = type(expr)
            if domain.mapping is None:
                new  = eval('Logical{0}_{1}d'.format(op, dim))
            else:
                new  = eval('{0}_{1}d'.format(op, dim))

            args = [cls.eval(i, domain=domain) for i in expr.args]
            return new(*args)

        elif isinstance(expr, _generic_ops):
            # if i = Dot(...) then type(i) is Grad
            op = type(expr)
            new  = eval('{0}_{1}d'.format(op, dim))
            args = [cls.eval(i, domain=domain) for i in expr.args]
            return new(*args)

        elif isinstance(expr, Trace):
            # TODO treate different spaces
            if expr.order == 0:
                return cls.eval(expr.expr, domain=domain)

            elif expr.order == 1:
                # TODO give a name to normal vector
                normal_vector_name = 'n'
                n = NormalVector(normal_vector_name)
                M = cls.eval(expr.expr, domain=domain)

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
                    line.append(cls.eval(expr[i,j], domain=domain))
                lines.append(line)
            return ImmutableDenseMatrix(lines)

        elif isinstance(expr, LogicalExpr):
            domain    = expr.domain
            expr      = cls(expr.expr, domain=domain)
            return LogicalExpr(expr, domain=domain)

        return expr


#==============================================================================
# TODO use random_string for space name and use it for 1d test function
def _split_test_function(expr):

    if isinstance(expr, ScalarFunction):

        dim = expr.space.ldim
        name = expr.name

        ls = []
        for i in range(dim):
            Di = Interval()
            Vi = ScalarFunctionSpace('tmp_V_{}'.format(i), domain=Di)

            ai = ScalarFunction(Vi, '{name}_{i}'.format(name=name, i=i+1))
            ls += [ai]

        return {expr:tuple(ls)}
    elif isinstance(expr, VectorFunction):

        dim = expr.space.ldim
        name = expr.name

        ls = {}
        Di = Interval()
        for i in range(dim):
            Vi = ScalarFunctionSpace('tmp_V_{}'.format(i), domain=Di)
            ls[expr[i]] = tuple(ScalarFunction(Vi, '{name}_{i}_{j}'.format(name=name, i=i,j=j+1)) for j in range(dim))

        return ls

    elif isinstance(expr, IndexedVectorFunction):

        dim = expr.base.space.ldim
        index = expr.indices[0]
        name = expr.base.name

        Di = Interval()

        ls = []
        for i in range(dim):
            Vi = ScalarFunctionSpace('tmp_V_{}'.format(i), domain=Di)
            ai = ScalarFunction(Vi, '{name}_{j}_{i}'.format(name=name, j=index,i=i+1))
            ls += [ai]
        return {expr:tuple(ls)}

    else:
        msg = 'Expecting ScalarFunction or IndexedVectorFunction, given {}'.format(type(expr))
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

    elif isinstance(expr, (ScalarFunction, IndexedVectorFunction)):
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
        d_atoms = kwargs.pop('d_atoms', {})
        domain  = kwargs.pop('domain', None)
        expand  = kwargs.pop('expand', False)

        if isinstance(expr, Add):
            args = [cls.eval(a, d_atoms=d_atoms, domain=domain) for a in expr.args]
            return Add(*args)

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if isinstance(a, _coeffs_registery)]
            stests  = [a for a in expr.args if not(a in coeffs) and a.atoms(ScalarFunction)]
            vtests  = [a for a in expr.args if not(a in coeffs) and a.atoms(VectorFunction)]
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
                args = [cls.eval(i, d_atoms=d_atoms, domain=domain) for i in stests]
                sb = Mul(*args)

            vb = S.One
            if vtests:
                args = [cls.eval(i, d_atoms=d_atoms, domain=domain) for i in vtests]
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
                    line.append(cls.eval(expr[i_row,i_col], d_atoms=d_atoms, domain=domain))

                lines.append(line)

            return ImmutableDenseMatrix(lines)

        elif isinstance(expr, DomainExpression):
            # TODO to be removed
            return cls.eval(expr.expr, d_atoms=d_atoms, domain=domain)

        elif isinstance(expr, BilinearForm):
            trials = expr.variables[0]
            tests  = expr.variables[1]
            fields = expr.fields

            # ... # TODO improve
            terminal_expr = TerminalExpr(expr, domain)[0]
            # ...

            # ...
            # Collect all variables in the expression (fields must be excluded)
            variables  = [e for e in terminal_expr.atoms(ScalarFunction) if e not in fields]
            variables += [e for e in terminal_expr.atoms(IndexedVectorFunction) if e.base not in fields]
            # ...

            # ...
            if domain is not None and domain.mapping is not None:
                terminal_expr = LogicalExpr(terminal_expr.expr, domain)
                variables     = [LogicalExpr(e, domain) for e in variables ]
                trials        = [LogicalExpr(e, domain) for e in trials ]
                tests         = [LogicalExpr(e, domain) for e in tests ]
            # ...

            # Prepare dictionary for '_tensorize_atomic_expr', which should
            # process all variables but leave fields unchanged:
            d_atoms = dict(
                [(v,_split_test_function(v)[v]) for v in variables] + [(f, (f,)) for f in fields])

            # ...
            expr = cls.eval(terminal_expr, d_atoms=d_atoms, domain=domain)
            # ...

            # ...
            trials = [a for a in variables if ((isinstance(a, ScalarFunction) and a in trials) or
                                               (isinstance(a, IndexedVectorFunction) and a.base in trials))]

            tests = [a for a in variables if ((isinstance(a, ScalarFunction) and a in tests) or
                                              (isinstance(a, IndexedVectorFunction) and a.base in tests))]
            # ...

            expr = _replace_atomic_expr(expr, trials, tests, d_atoms,
                                        logical=True, expand=expand)

            return expr

        if expr.atoms(ScalarFunction) or expr.atoms(IndexedVectorFunction):
            return _tensorize_atomic_expr(expr, d_atoms)

        return cls(expr, evaluate=False)
