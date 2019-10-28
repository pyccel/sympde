# coding: utf-8

from itertools import product

from sympy import S
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
from sympde.calculus.core import _generic_ops

from sympde.topology import BasicDomain, Union, Interval
from sympde.topology import NormalVector, TangentVector
from sympde.topology import Boundary, Interface
from sympde.topology import InteriorDomain
from sympde.topology import DetJacobian
from sympde.topology import SymbolicDeterminant
from sympde.topology import LogicalExpr
from sympde.topology.space import ScalarFunctionSpace
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

from .basic import BasicExpr, BasicForm
from .expr  import BilinearForm
from .expr  import DomainIntegral, BoundaryIntegral, InterfaceIntegral
from .expr  import Functional
from .expr  import _get_domain

#==============================================================================
def is_sequence(a):
    return isinstance(a, (list,tuple,Tuple))

#==============================================================================
# TODO use this function everywhere it is needed
def zero_matrix(n_rows, n_cols):
    lines = []
    for i in range(0, n_rows):
        line = []
        for j in range(0, n_cols):
            line.append(0)
        lines.append(line)

    return Matrix(lines)

#==============================================================================
def _get_size_and_starts(ls):
    n = 0
    d_indices = {}
    for x in ls:
        d_indices[x] = n
        if isinstance(x, ScalarTestFunction):
            n += 1

        elif isinstance(x, VectorTestFunction):
            for j in range(0, x.shape[0]):
                d_indices[x[j]] = n + j

            n += x.shape[0]

    return n, d_indices

#==============================================================================
def _init_matrix(expr):
    assert(isinstance(expr, (BasicForm, BasicExpr)))

    if expr.is_bilinear:

        trials = list(expr.variables[0])
        n_cols, trial_indices = _get_size_and_starts(trials)

        tests = list(expr.variables[1])
        n_rows, test_indices = _get_size_and_starts(tests)

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

#==============================================================================
def _to_matrix_bilinear_form(expr, M, test_indices, trial_indices):

    # ...
    def treat_form(arg, M):
        atoms  = list(arg.atoms(ScalarTestFunction))
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

    elif isinstance(expr, (ScalarTestFunction, VectorTestFunction, IndexedTestTrial)):
        M = treat_form(expr, M)
    else:
        raise TypeError('> wrong type, given {}'.format(type(expr)))
    # ...

    return M

#==============================================================================
def _to_matrix_linear_form(expr, M, test_indices):
    # ...
    def treat_form(arg, M):
        atoms  = list(arg.atoms(ScalarTestFunction))
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

    elif isinstance(expr, (ScalarTestFunction, VectorTestFunction, IndexedTestTrial)):
        M = treat_form(expr, M)

    else:
        raise TypeError('> wrong type, given {}'.format(type(expr)))
    # ...

    return M

#==============================================================================
def _to_matrix_functional_form(expr, M):
    M[0] += expr

    return M

#==============================================================================
def _to_matrix_form(expr, M, test_indices, trial_indices):
    if not(test_indices is None) and not(trial_indices is None):
        return _to_matrix_bilinear_form(expr, M, test_indices, trial_indices)

    if not(test_indices is None) and trial_indices is None:
        return _to_matrix_linear_form(expr, M, test_indices)

    if test_indices is None and trial_indices is None:
        return _to_matrix_functional_form(expr, M)

#==============================================================================
def _split_expr_over_subdomains(expr, interiors, tests=None, trials=None):
    """
    Splits an expression defined on a domain, having multiple interiors, into
    expressions where the test and trial functions are defined on each side of
    the subdomain.

    Parameters:
        expr: sympde expression

        interiors: an interior or union of interiors

        tests: tests functions as given from linear or bilinear forms

        trials: trials functions as given from linear or bilinear forms

    Returns: sympde expression
    """
    # ...
    def _new_atom(v, interior):
        new = '{v}_{domain}'.format( v      = v.name,
                                     domain = interior.name )
        return element_of(v.space, name=new)
    # ...

    # ...
    is_bilinear = not( trials is None ) and not( tests is None )
    is_linear   =    ( trials is None ) and not( tests is None )

    if trials is None: trials = []
    if tests  is None: tests  = []
    # ...

    # ...
    d_expr = {}
    for interior in interiors:
        d_expr[interior] = 0
    # ...

    # ... create trial/test functions on each subdomain
    d_trials = {}
    d_tests  = {}
    for interior in interiors:
        # ...
        news = []
        for v in trials:
            new = _new_atom(v, interior)
            news.append(new)
        d_trials[interior] = news
        # ...

        # ...
        news = []
        for v in tests:
            new = _new_atom(v, interior)
            news.append(new)
        d_tests[interior] = news
        # ...
    # ...

    # ...
    for interior in interiors:
        d_expr[interior] = expr

        # use trial functions on each subdomain
        olds = trials
        news = d_trials[interior]
        for old,new in zip(olds, news):
            d_expr[interior] = d_expr[interior].subs({old: new})

        # use test functions on each subdomain
        olds = tests
        news = d_tests[interior]
        for old,new in zip(olds, news):
            d_expr[interior] = d_expr[interior].subs({old: new})
    # ...

    return d_expr

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
    def _new_atom(v, label):
        new = '{v}_{label}'.format( v     = v.name,
                                    label = label )
        return element_of(v.space, name=new)
    # ...

    # ...
    is_bilinear = not( trials is None ) and not( tests is None )
    is_linear   =    ( trials is None ) and not( tests is None )

    if trials is None: trials = []
    if tests  is None: tests  = []
    # ...

    # ...
    B_minus = interface.minus
    B_plus  = interface.plus
    boundaries = (B_minus, B_plus)
    labels     = ('minus', 'plus')
    # ...

    int_expressions = []
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

    expr = expand(expr)
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

    # ...
    def _not_zero_matrix(M):
        n,m = expr.shape
        return any([M[i,j] != 0 for i,j in product(range(n), range(m))])
    # ...

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
                if _not_zero_matrix(newexpr):
                    bnd_expressions[interface.minus] = newexpr
                # ...

                # ...
                # TODO must call InterfaceExpression afterward
                newexpr = _nullify(expr, u_minus, trials)
                newexpr = _nullify(newexpr, v_plus, tests)
                if _not_zero_matrix(newexpr):
                    int_expressions += [InterfaceExpression(interface, newexpr)]
                # ...

                # ...
                # TODO must call InterfaceExpression afterward
                newexpr = _nullify(expr, u_plus, trials)
                newexpr = _nullify(newexpr, v_minus, tests)
                if _not_zero_matrix(newexpr):
                    int_expressions += [InterfaceExpression(interface, newexpr)]
                # ...

                # ...
                newexpr = _nullify(expr, u_plus, trials)
                newexpr = _nullify(newexpr, v_plus, tests)
                newexpr = newexpr.subs({u_plus: u, v_plus: v})
                if _not_zero_matrix(newexpr):
                    bnd_expressions[interface.plus] = newexpr
                # ...

    elif is_linear:
        for v in d_tests.keys():
            v_minus = d_tests[v]['-']
            v_plus  = d_tests[v]['+']

            # ...
            newexpr = _nullify(expr, v_minus, tests)
            newexpr = newexpr.subs({v_minus: v})
            if _not_zero_matrix(newexpr):
                bnd_expressions[interface.minus] = newexpr
            # ...

            # ...
            newexpr = _nullify(expr, v_plus, tests)
            newexpr = newexpr.subs({v_plus: v})
            if _not_zero_matrix(newexpr):
                bnd_expressions[interface.plus] = newexpr
            # ...
    # ...

    return int_expressions, bnd_expressions


#==============================================================================
class KernelExpression(Basic):
    def __new__(cls, target, expr):
        assert(isinstance(expr, (Matrix, ImmutableDenseMatrix)))
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
    pass


#==============================================================================
class TerminalExpr(CalculusFunction):

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):

            args = cls._annotate(*args)
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

    # TODO should we keep it?
    def _annotate(*args):
        args = list(args)
        expr = args[0]
        if isinstance(expr, BasicForm):
            if not expr.is_annotated:
                expr = expr.annotate()
        else:
            if isinstance(expr, (Add, Mul)):
                indexed_fields = list(expr.atoms(IndexedTestTrial))
                new_indexed_fields = [VectorField(F.base.space,F.base.name) for F in indexed_fields]
                new_indexed_fields = [new_F[F.indices[0]] for new_F,F in zip(new_indexed_fields, indexed_fields)]
                expr = expr.subs(zip(indexed_fields, new_indexed_fields))
                fields = list(expr.atoms(ScalarTestFunction,VectorTestFunction).difference(indexed_fields))
                new_fields = [f.space.field(f.name) for f in fields]
                expr = expr.subs(zip(fields, new_fields))


        args[0] = expr
        return args

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
            o = args[0]
            for arg in args[1:]:
                o = o + arg
            return o

        elif isinstance(expr, Mul):
            args = [cls.eval(a, dim=dim) for a in expr.args]
            o = args[0]
            for arg in args[1:]:
                o = o * arg
            return o

        elif isinstance(expr, (ScalarTestFunction, VectorTestFunction)):
            return expr

        elif isinstance(expr, BasicForm):
            # ...
            dim = expr.ldim
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

            if isinstance(expr.expr, Add):
                for a in expr.expr.args:
                    newexpr = cls.eval(a, dim=dim)
                    newexpr = expand(newexpr)

                    # ...
                    try:
                        domain = _get_domain(a)
                        if isinstance(domain, Union):
                            domain = list(domain._args)

                        elif not is_sequence(domain):
                            domain = [domain]
                    except:
                        pass
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

                if newexpr != 0:
                    M, test_indices, trial_indices = _init_matrix(expr)
                    M = _to_matrix_form(newexpr, M, test_indices, trial_indices)

                    # TODO ARA make sure thre is no problem with psydac
                    #      we should always take the interior of a domain
                    if not isinstance(domain, (Boundary, Interface, InteriorDomain)):
                        domain = domain.interior

                    d_new[domain] = M
            # ...

            # ...
            ls = []
            d_all = {}
            # ...
#            print(d_new)
#            print([type(i) for i in d_new.keys()])

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
                ls_int, d = _split_expr_over_interface(newexpr, interface,
                                                       tests=tests,
                                                       trials=trials)
                # ...

                # ...
                ls += ls_int
                # ...

                # ...
                for k, v in d.items():
                    if k in d_all.keys():
                        d_all[k] += v

                    else:
                        d_all[k] = v
                # ...
            # ...

            # ... treating subdomains
            keys = [k for k in d_new.keys() if isinstance(k, Union)]
            for domain in keys:
                # ...
                trials = None
                tests  = None
                if expr.is_bilinear:
                    trials = list(expr.variables[0])
                    tests  = list(expr.variables[1])

                elif expr.is_linear:
                    tests  = list(expr.variables)

                else:
                    raise NotImplementedError('Only Bilinear and Linear forms are available')
                # ...

                # ...
                newexpr = d_new[domain]
                d = _split_expr_over_subdomains(newexpr, domain.as_tuple(),
                                                tests=tests, trials=trials)
                # ...

                # ...
                for k, v in d.items():
                    if k in d_all.keys():
                        d_all[k] += v

                    else:
                        d_all[k] = v
                # ...
            # ...

            # ...
            d = {}

            for k, v in d_new.items():
                if not isinstance( k, (Interface, Union) ):
                    d[k] = d_new[k]

            for k, v in d_all.items():
                if k in d.keys():
                    d[k] += v

                else:
                    d[k] = v

            d_new = d
            # ...

            # ...
            for domain, newexpr in d_new.items():
                if isinstance(domain, Boundary):
                    ls += [BoundaryExpression(domain, newexpr)]

                elif isinstance(domain, Interface):
                    ls += [InterfaceExpression(domain, newexpr)]

                elif isinstance(domain, BasicDomain):
                    ls += [DomainExpression(domain, newexpr)]

                else:
                    raise TypeError('not implemented for {}'.format(type(domain)))
            # ...
            return ls

        elif isinstance(expr, (DomainIntegral, BoundaryIntegral, InterfaceIntegral)):
            if dim is None:
                domain = expr.domain
                dim = domain.dim

            return cls.eval(expr._args[0], dim=dim)

        elif isinstance(expr, NormalVector):
            lines = [[expr[i] for i in range(dim)]]
            return Matrix(lines)

        elif isinstance(expr, TangentVector):
            lines = [[expr[i] for i in range(dim)]]
            return Matrix(lines)

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


#==============================================================================
# TODO use random_string for space name and use it for 1d test function
def _split_test_function(expr):

    if isinstance(expr, ScalarTestFunction):

        dim = expr.space.ldim
        name = expr.name

        ls = []
        for i in range(0, dim):
            Di = Interval()
            Vi = ScalarFunctionSpace('tmp_V_{}'.format(i), domain=Di)

            ai = ScalarTestFunction(Vi, '{name}{i}'.format(name=name, i=i+1))
            ls += [ai]

        return ls

    elif isinstance(expr, IndexedTestTrial):

        i = expr.indices
        assert(len(i) == 1)
        i = i[0]

        V = expr.base.space
        Vi = ScalarFunctionSpace('tmpV_{}'.format(i), V.domain)
        vi = ScalarTestFunction(Vi, '{test}{i}'.format(test=expr.base.name, i=i))

        return _split_test_function(vi)

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
        atoms = d_atoms[expr]
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
def _replace_atomic_expr(expr, trials, tests, d_atoms, logical=False):

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
        mapping = kwargs.pop('mapping', None)

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

            return Matrix(lines)

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

            d_atoms = {}
            for a in variables:
                new = _split_test_function(a)
                d_atoms[a] = new
            # ...

            # ...
            logical = False
            if not(mapping is None):
                logical = True
                terminal_expr = LogicalExpr(mapping, terminal_expr.expr)

                det_M = DetJacobian(mapping)
                det   = SymbolicDeterminant(mapping)
                terminal_expr = terminal_expr.subs(det_M, det)
                terminal_expr = expand(terminal_expr)
            # ...

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
                                        logical=logical)

            return expr

        if expr.atoms(ScalarTestFunction) or expr.atoms(IndexedTestTrial):
            return _tensorize_atomic_expr(expr, d_atoms)

        return cls(expr, evaluate=False)
