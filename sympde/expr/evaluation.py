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
from sympde.topology.derivatives import _logical_partial_derivatives
from sympde.topology.derivatives import partial_derivative_as_symbol
from sympde.topology.derivatives import sort_partial_derivatives
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
from sympde.topology.space import BasicFunctionSpace
from sympde.topology.space import FunctionSpace
from sympde.topology.space import ProductSpace
from sympde.topology.space import ScalarTestFunction
from sympde.topology.space import VectorTestFunction
from sympde.topology.space import Element,IndexedElement
from sympde.topology.space import IndexedTestTrial
from sympde.topology.space import Unknown, VectorUnknown
from sympde.topology.space import Trace
from sympde.topology.space import ScalarField, VectorField, IndexedVectorField
from sympde.topology.measure import CanonicalMeasure
from sympde.topology.measure import CartesianMeasure
from sympde.topology.measure import Measure
from sympde.topology import Mapping, DetJacobian
from sympde.topology import SymbolicDeterminant, SymbolicCovariant, SymbolicContravariant

from sympde.topology import LogicalExpr

from .basic  import BasicExpr, BasicForm
from .expr   import LinearExpr, BilinearExpr
from .expr   import LinearForm, BilinearForm, Norm
from .equation import Equation
from .expr import BasicIntegral, DomainIntegral, BoundaryIntegral
from .expr import Functional
from .expr import _get_domain

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

    else:
        raise TypeError('> wrong type, given {}'.format(type(expr)))
    # ...

    return M

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
            # ...z
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

            # ...
            if isinstance(expr, (Equation, BilinearForm)):
            
                test_functions  = expr.test_functions
                trial_functions = expr.trial_functions
                new_test_functions  = [f.space.element(f.name) for f in test_functions]
                new_trial_functions = [f.space.element(f.name) for f in trial_functions]
                expr = expr.subs(zip(test_functions, new_test_functions))
                expr = expr.subs(zip(trial_functions, new_trial_functions))
            elif isinstance(expr, LinearForm):
            
                test_functions  = expr.test_functions
                new_test_functions  = [f.space.element(f.name) for f in test_functions]
                expr = expr.subs(zip(test_functions, new_test_functions))
            elif isinstance(expr, Functional):
            
                indexed_fields = list(expr.expr.atoms(IndexedElement))
                new_indexed_fields = [VectorField(F.base.space,F.base.name) for F in indexed_fields]
                new_indexed_fields = [new_F[F.indices[0]] for new_F,F in zip(new_indexed_fields, indexed_fields)]
                expr = expr.subs(zip(indexed_fields, new_indexed_fields))
                fields = list(expr.expr.atoms(Element).difference(indexed_fields))
                new_fields = [f.space.field(f.name) for f in fields]
                expr = expr.subs(zip(fields, new_fields))
            else:
                raise NotImplementedError('TODO')

            
            # ...
            
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


#==============================================================================
# TODO use random_string for space name and use it for 1d test function
def _split_test_function(expr):

    if isinstance(expr, ScalarTestFunction):

        dim = expr.space.ldim
        name = expr.name

        ls = []
        for i in range(0, dim):
            Di = Interval()
            Vi = FunctionSpace('tmp_V_{}'.format(i), domain=Di)

            ai = ScalarTestFunction(Vi, '{name}{i}'.format(name=name, i=i))
            ls += [ai]

        return ls

    elif isinstance(expr, IndexedTestTrial):

        i = expr.indices
        assert(len(i) == 1)
        i = i[0]

        V = expr.base.space
        Vi = FunctionSpace('tmpV_{}'.format(i), V.domain)
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
        print(expr)
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
