# coding: utf-8

# TODO - transpose of BilinearForm
#      - add unknown status if only one space is given to the BilinearForm
#        => can not be evaluated except if it is called on a couple test/trial
#      - check that a BilinearForm is bilinear (using Constant)
#      - check that a LinearForm is linear
#      - add is_symmetric property for BilinearForm

from numpy import zeros
from itertools import groupby

from sympy.core import Basic
from sympy.core import Symbol
from sympy.core import Function
from sympy.core import Expr, Add, Mul, Pow
from sympy import S
from sympy.core.containers import Tuple
from sympy import preorder_traversal
from sympy import Indexed, IndexedBase, Matrix, ImmutableDenseMatrix
from sympy.matrices.dense import MutableDenseMatrix
from sympy.physics.quantum import TensorProduct
from sympy import expand
from sympy import Integer, Float
from sympy.core.expr import AtomicExpr
from sympy.physics.quantum import TensorProduct

from .measure import CanonicalMeasure
from .measure import CartesianMeasure
from .measure import Measure
from .geometry import BasicDomain, Domain, Line
from .geometry import BoundaryVector, NormalVector, TangentVector, Boundary
from .derivatives import _partial_derivatives
from .derivatives import partial_derivative_as_symbol
from .derivatives import sort_partial_derivatives
from .derivatives import get_atom_derivatives
from .derivatives import dx, dy, dz
from .derivatives import (Grad_1d, Div_1d,
                          Grad_2d, Curl_2d, Rot_2d, Div_2d,
                          Grad_3d, Curl_3d, Div_3d)

from .basic import _coeffs_registery
from .basic import CalculusFunction
from .basic import Field, Constant

from .generic import Dot, Inner, Cross
from .generic import Grad, Rot, Curl, Div
from .generic import _generic_ops

from .algebra import (Dot_1d,
                      Dot_2d, Inner_2d, Cross_2d,
                      Dot_3d, Inner_3d, Cross_3d)

from .space import FunctionSpace
from .space import ProductSpace
from .space import TestFunction
from .space import VectorTestFunction
from .space import IndexedTestTrial
from .space import Unknown, VectorUnknown
from .space import Trace

from .errors import UnconsistentError


def _initialize_measure(measure, coordinates):
    if not( measure is None ):
        return measure

    if not( coordinates is None ):
        return Measure(coordinates)

    else:
        raise NotImplementedError('')

def _initialize_boundary(expr):

    traces = expr.atoms(Trace)
    boundaries = [trace.boundary for trace in traces]
    boundaries = list(set(boundaries)) # remove redanduncy

    boundary = None
    if len(boundaries) == 0:
        boundary = None

    elif len(boundaries) == 1:
        boundary = boundaries[0]

    elif (len(boundaries) > 1):
        if not is_sum_of_form_calls(expr):
            msg = '> BilinearForm can not be defined on different boundaries'
            raise UnconsistentError(msg)

        boundary = True

    # ...
    if isinstance(boundary, Boundary):
        # ...
        if isinstance(expr, Add) and not is_sum_of_form_calls(expr):
            args = expr.args
            args = [a for a in args if not a.atoms(Trace)]
            if args:
                msg = '> Only boundary terms, using traces, are allowed'
                raise UnconsistentError(msg)
        # ...
    # ...

    return boundary


class BasicForm(Expr):
    _name = None
    _boundary = None

    @property
    def fields(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Field)]
        # no redanduncy
        return sorted(list(set(ls)))

    @property
    def constants(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Constant)]
        # no redanduncy
        return list(set(ls))

    @property
    def domain(self):
        return self._domain

    @property
    def boundary(self):
        return self._boundary

    @property
    def measure(self):
        return self._measure

    @property
    def mapping(self):
        return self._mapping

    @property
    def name(self):
        return self._name


# TODO we should check that the only free symbols are fields, constants or coordinates
class Integral(BasicForm):
    """

    Examples

    """
    _ldim = None
    _coordinates = None
    # TODO we must remove space from here
    def __new__(cls, expr, domain, coordinates=None, measure=None,
                mapping=None, name=None):

        if not isinstance(domain, BasicDomain):
            raise TypeError('> Expecting a BasicDomain object for domain')

        # ... check that there are no test functions in the expression
        ls = [a for a in expr.free_symbols if isinstance(a, (TestFunction, VectorTestFunction))]
        if not(len(ls) == 0):
            raise TypeError('Cannot use test functions in Integral')
        # ...

        # compute dim from fields if available
        ls = [a for a in expr.free_symbols if isinstance(a, Field)]
        if ls:
            F = ls[0]
            space = F.space
            ldim = F.space.ldim

            if coordinates is None:
                coordinates = F.space.coordinates

        else:
            if coordinates is None:
                raise ValueError('> Coordinates must be provided if the expression has no fields')

            coordinates = [str(i) for i in coordinates]
            ldim = len(coordinates)

            ID = abs(hash(expr))
            space = FunctionSpace('space_{}'.format(ID), domain,
                                  coordinates=coordinates)

            if ldim == 1:
                coordinates = coordinates[0]

        measure = _initialize_measure(measure, coordinates)

        # get boundary terms
        boundary = _initialize_boundary(expr)

        # so that we may pass mapping as False
        if not mapping: mapping=None

        obj = Basic.__new__(cls, expr)
        obj._ldim = ldim
        obj._coordinates = coordinates
        obj._boundary = boundary
        obj._domain = domain
        obj._measure = measure
        obj._mapping = mapping
        obj._space = space
        obj._name = name

        return obj

    @property
    def expr(self):
        return self._args[0]

    @property
    def ldim(self):
        return self._ldim

    @property
    def space(self):
        return self._space

    @property
    def coordinates(self):
        return self._coordinates

    def _sympystr(self, printer):
        sstr = printer.doprint
        expr = self.expr
        return sstr(expr)

    # TODO how to implement this?
    def __call__(self, *args):
        raise NotImplementedError('')


class LinearForm(BasicForm):
    """

    Examples

    """
    def __new__(cls, arguments, expr, measure=None, mapping=None, name=None):

        args = _sanitize_form_arguments(arguments, expr, is_linear=True)
        obj = Basic.__new__(cls, args, expr)

        # TODO must check that all domains are the same
        domain = obj.test_spaces[0].domain
        measure = _initialize_measure(measure, obj.coordinates)

        # get boundary terms
        boundary = _initialize_boundary(expr)

        # so that we may pass mapping as False
        if not mapping: mapping=None

        obj._domain = domain
        obj._boundary = boundary
        obj._measure = measure
        obj._mapping = mapping
        obj._name = name

        return obj

    @property
    def variables(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

    @property
    def test_functions(self):
        return self.variables

    @property
    def ldim(self):
        return self.test_spaces[0].ldim

    @property
    def coordinates(self):
        return self.test_spaces[0].coordinates

    @property
    def test_spaces(self):
        return [u.space for u in self.test_functions]

    def _sympystr(self, printer):
        sstr = printer.doprint
        expr = self.expr
        return sstr(expr)

    def __call__(self, *args):
        args = Tuple(*args)
        return FormCall(self, args)


class BilinearForm(BasicForm):
    """

    Examples

    """
    def __new__(cls, arguments, expr, measure=None, mapping=None, name=None):

        # ...
        if not isinstance(arguments, (tuple, list, Tuple)):
            raise TypeError('(test, trial) must be a tuple, list or Tuple')

        if not(len(arguments) == 2):
            raise ValueError('Expecting a couple (test, trial)')

        args = _sanitize_form_arguments(arguments, expr, is_bilinear=True)
        obj = Basic.__new__(cls, args, expr)

        # TODO must check that all domains are the same
        domain = obj.test_spaces[0].domain
        measure = _initialize_measure(measure, obj.coordinates)

        # get boundary terms
        boundary = _initialize_boundary(expr)

        # so that we may pass mapping as False
        if not mapping: mapping=None

        obj._domain = domain
        obj._boundary = boundary
        obj._measure = measure
        obj._mapping = mapping
        obj._name = name

        return obj

    @property
    def variables(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

    @property
    def test_functions(self):
        return self.variables[0]

    @property
    def trial_functions(self):
        return self.variables[1]

    @property
    def ldim(self):
        return self.test_spaces[0].ldim

    @property
    def coordinates(self):
        return self.test_spaces[0].coordinates

    @property
    def test_spaces(self):
        return [u.space for u in self.test_functions]

    @property
    def trial_spaces(self):
        return [u.space for u in self.trial_functions]

    def _sympystr(self, printer):
        sstr = printer.doprint
        expr = self.expr
        return sstr(expr)

    def __call__(self, *args):
        if not(len(args) == 2):
            raise ValueError('Expecting a couple (test, trial)')

        return FormCall(self, args)


class BilinearAtomicForm(BilinearForm, AtomicExpr):
    """

    Examples

    """
    _name = None

    @property
    def name(self):
        return self._name

    def _sympystr(self, printer):
        sstr = printer.doprint
        name = sstr(self.name)

        test = [sstr(i) for i in self.test_functions]
        test = ','.join(i for i in test)

        trial = [sstr(i) for i in self.trial_functions]
        trial = ','.join(i for i in trial)

        return '{name}({test},{trial})'.format(name=name, trial=trial, test=test)

class Mass(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'Mass'
    def __new__(cls, test, trial):

        test_trial = [test, trial]
        expr = test * trial

        return BilinearForm.__new__(cls, test_trial, expr)

class Stiffness(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'Stiffness'
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = d(test) * d(trial)

        return BilinearForm.__new__(cls, test_trial, expr)

class Advection(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'Advection'
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = test * d(trial)

        return BilinearForm.__new__(cls, test_trial, expr)

class AdvectionT(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'AdvectionT'
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = d(test) * trial

        return BilinearForm.__new__(cls, test_trial, expr)


class Kron(BilinearAtomicForm):

    """."""

    def __new__(cls, ls):

        ln = len(ls)
        dim = len(ls[0])
        obj = Basic.__new__(cls, ls)
        if dim >= 3:
            raise NotImplementedError('TODO')

        args1 = ['A%s'%i for i in range(ln)]
        args2 = ['B%s'%i for i in range(ln)]
        try:
            import importlib
            package = importlib.import_module("sympde.codegen.templates.kron")
        except:
            raise ImportError('could not import kron_dot')
        name = 'kron_dot_2d'
        template = getattr(package,name)
        body = [template['body'].format(MAT1=args1[i],MAT2=args2[i]) for i in range(ln)]
        body = '\n'.join(i for i in body)
        args = ','.join(arg for arg in args1+args2)
        function = template['function'].format(___MAT_ARGS___=args,__BODY__=body)
        args_types = ','.join('double[:,:]' for i in range(2*ln))
        header = template['header'].format(__ARGS_TYPES__=args_types)
        dot = compile(function,'','single')
        dic = {}
        eval(dot,dic)
        _dot = dic[name]
        setattr(obj, '_dot',_dot)
        return obj

    @property
    def args(self):
        return self._args[0]

    def dot(self, x):
        space = x.space
        args = list(zip(*self.args))
        args1 = args[0]
        args2 = args[1]
        args1 = [arg._data for arg in args1]
        args2 = [arg._data for arg in args2]
        #args1 = [arg._data.T for arg in args1]
        #args2 = [arg._data.T for arg in args2]
        starts = space.starts
        ends   = space.ends
        pads   = space.pads

        from spl.linalg.stencil import StencilVector
        Y      = StencilVector(space)
        X_tmp  = StencilVector(space)
        #self._dot(starts,ends,pads,x._data.T,Y._data.T,X_tmp._data.T,*args1,*args2)
        args = list(args1) + list(args2)
        self._dot(starts,ends,pads,x._data,Y._data,X_tmp._data,*args)
        return Y


    def __str__(self):
        return 'Kron'

    def _sympystr(self, printer):
        return 'Kron'


class FormCall(AtomicExpr):

    is_commutative = False

    def __new__(cls, expr, args, name=None):

        if not isinstance(expr, (BilinearForm, LinearForm)):
            raise TypeError('> Expecting BilinearForm, LinearForm')

        if not isinstance(args, (list, tuple, Tuple)):
            args = [args]

        if isinstance(expr, BilinearForm):
            expr = subs_bilinear_form(expr, args)

        if isinstance(expr, LinearForm):
            expr = subs_linear_form(expr, args)

        # ...
        if not name:
            if expr.name:
                name = expr.name

            else:
                name = 'FormCall'
#                raise ValueError('Callable Bilinear/Linear form must have a name')

        args = Tuple(*args)
        obj = Basic.__new__(cls, expr, args, name)
        # ...

        return obj

    @property
    def expr(self):
        return self._args[0]

    @property
    def arguments(self):
        return self._args[1]

    @property
    def name(self):
        return self._args[2]

    @property
    def form_name(self):
        return self.name

    def _sympystr(self, printer):
        sstr = printer.doprint

        expr = self.expr

        name = sstr(self.name)

        # ...
        test = [sstr(i) for i in expr.test_functions]
        test_str = ','.join(i for i in test)

        if len(test) == 1:
            test = test_str

        else:
            test = '({})'.format(test_str)
        # ...

        if isinstance(expr, BilinearForm):
            # ...
            trial = [sstr(i) for i in expr.trial_functions]
            trial_str = ','.join(i for i in trial)

            if len(trial) == 1:
                trial = trial_str

            else:
                trial = '({})'.format(trial_str)
            # ...

            return '{name}({test},{trial})'.format(name=name,
                                                    trial=trial,
                                                    test=test)

        if isinstance(expr, LinearForm):
            return '{name}({test})'.format(name=name, test=test)

def is_mul_of_form_call(expr):
    if not isinstance(expr, Mul):
        return False

    are_calls = [isinstance(i, FormCall) for i in expr.args]
    any_are_calls = any(are_calls)
    if not any_are_calls:
        return False

    if (isinstance(any_are_calls, (list, tuple, Tuple)) and
        (len(any_are_calls) > 1)):
        raise TypeError('> Cannot multiply calls of Bilinear/Linear forms')

    return True

def is_sum_of_form_calls(expr):
    if not isinstance(expr, Add):
        return False

    are_valid = [isinstance(i, FormCall) or is_mul_of_form_call(i) for i in expr.args]

    all_are_valid = all(are_valid)
    if any(are_valid) and not(all_are_valid):
        raise TypeError('> Invalid expression')

    return all_are_valid


def _sanitize_form_arguments(arguments, expr, is_bilinear=False, is_linear=False):

    is_linear = is_linear or (len(expr.atoms(LinearForm)) > 0)
    is_bilinear = is_bilinear or (len(expr.atoms(BilinearForm)) > 0)

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


# ...
def atomize(expr, dim=None):
    """
    """
    if not isinstance(expr, (Expr,
                             _partial_derivatives, _generic_ops,
                             TestFunction, VectorTestFunction, Indexed,
                             Field, Constant, Symbol, Function,
                             BoundaryVector, Trace,
                             Integer, Float, Matrix, ImmutableDenseMatrix,
                             list, tuple, Tuple)):
        msg = ('> Wrong input type.')

        raise TypeError(msg, ', given ', expr, type(expr))

    # ... replace a FormCall by its expression
    calls = expr.atoms(FormCall)
    for call in calls:
        expr = expr.subs(call, call.expr)
    # ...

#    print('> expr [atomize] = ', expr, type(expr))

    # ... compute dim if None
    if dim is None:
        ls = [i for i in expr.free_symbols if isinstance(i, (TestFunction,
                                                             VectorTestFunction,
                                                             Field))]

#        ls = expr.atoms((TestFunction, VectorTestFunction, Field))
#        ls = list(ls)
        if ls:
            atom = ls[0]
            if atom.space is None:
                raise ValueError('Expecting atom to be associated to a space')

            dim = atom.space.ldim
    # ...

    if isinstance(expr, (list, tuple, Tuple)):
        args = [atomize(i, dim=dim) for i in expr]
        return Tuple(*args)

    elif isinstance(expr, Add):
        args = [atomize(i, dim=dim) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [atomize(i, dim=dim) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)

    elif isinstance(expr, Pow):

        b = atomize(expr.base, dim=dim)
        e = expr.exp

        return Pow(b, e)

    elif isinstance(expr, BasicForm):

        return atomize(expr.expr, dim=dim)

    elif isinstance(expr, Trace):
        # TODO treate different spaces
        if expr.order == 0:
            return atomize(expr.expr, dim=dim)

        elif expr.order == 1:
            # TODO must be passed as key word to atomize
            normal_vector_name = 'n'
            n = NormalVector(normal_vector_name)
            M = atomize(expr.expr, dim=dim)
            if dim == 1:
                return M
            else:
                e = 0
                for i in range(0, dim):
                    e += M[i] * n[i]
                return e

        else:
            raise ValueError('> Only traces of order 0 and 1 are available')

    elif isinstance(expr, (Dot, Inner, Cross, Grad, Rot, Curl, Div)):
        # if i = Dot(...) then type(i) is Grad
        op = type(expr)
        new  = eval('{0}_{1}d'.format(op, dim))

        args = [atomize(i, dim=dim) for i in expr.args]
        return new(*args)

    elif isinstance(expr, Matrix):
        n,m = expr.shape
        lines = []
        for i in range(0, n):
            line = []
            for j in range(0, m):
                line.append(atomize(expr[i,j], dim=dim))
            lines.append(line)
        return Matrix(lines)

    return expr
# ...

# ...
def _evaluate_core(a, verbose=False, variables=None, M=None):

    # ...
    if not isinstance(a, (BasicForm, Add, Mul)):
        msg = 'Expecting a BasicForm, Add or Mul. Given {}'.format(type(a))
        raise TypeError(msg)
    # ...

    # ...
    variables = []
    if isinstance(a, (BilinearForm, LinearForm)):
        variables = a.variables

    elif isinstance(a, FormCall):
        variables = a.variables
    # ...

    # ...
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
    # ...

    # ...
    tests = []
    trials = []
    # ...

    # ...
    if isinstance(a, BilinearForm):
        tests = list(a.test_functions)
        n_rows, test_indices = _get_size_and_starts(a.test_functions)

        trials = list(a.trial_functions)
        n_cols, trial_indices = _get_size_and_starts(a.trial_functions)

        lines = []
        for i in range(0, n_rows):
            line = []
            for j in range(0, n_cols):
                line.append(0)
            lines.append(line)

        M = Matrix(lines)

    elif isinstance(a, LinearForm):
        tests = list(a.test_functions)
        n_rows, test_indices = _get_size_and_starts(a.test_functions)

        lines = [0 for i in range(0, n_rows)]
        M = Matrix(lines)
    # ...

    # ...
    if isinstance(a, Add):
        args = [_evaluate_core(i, verbose=verbose, variables=variables, M=M)
                for i in a.args]

        return Add(*args)

    elif isinstance(a, Mul):
        # a coeff can be a symbol, otherwise the expression c1 * a
        # raises an error
        coeffs  = [i for i in a.args if isinstance(i, _coeffs_registery) or isinstance(i, Symbol)]
        vectors = [i for i in a.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [_evaluate_core(i, verbose=verbose, variables=variables, M=M)
                    for i in vectors]
            j = Mul(*args)

        return Mul(i, j)
    # ...

    dim = a.ldim
    expr = a.expr

    # convert generic operators to atomic ones
    expr = atomize(expr, dim=dim)

    # we need to expand the expression so that we have a sum of product
    expr = expand(expr)

    if verbose:
        print('> atomized   >>> {0}'.format(expr))

    # ...
    def __evaluate_core_LinearForm(expr, M):
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
    # ...

    # ...
    def __evaluate_core_BilinearForm(expr, M):

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
    # ...

    # ...
    if isinstance(a, BilinearForm):
        M = __evaluate_core_BilinearForm(expr, M)

        # returning scalars when possibl
        if (n_rows == 1) and (n_cols == 1):
            return M[0, 0]

    elif isinstance(a, LinearForm):
        M = __evaluate_core_LinearForm(expr, M)

        # returning scalars when possibl
        if (n_rows == 1):
            return M[0]

    elif isinstance(a, Integral):
        return expr
    # ...

    return M
# ...

# ...
def _evaluate_bnd(a, bnd_calls, verbose=False):
    if verbose:
        print('> bnd calls = ', bnd_calls)

    a_expr = a
    if isinstance(a, BasicForm) and is_sum_of_form_calls(a.expr):
        a_expr = a.expr

    # ...
    boundaries = []
    groups = []
    keyfunc = lambda call: call.expr.boundary.name
    bnd_calls = sorted(bnd_calls, key=keyfunc)
    for k, g in groupby(bnd_calls, keyfunc):
        list_g = list(g)
        a_bnd = _extract_linear_combination(a_expr, list_g)
        groups.append(a_bnd)
        bnd = list_g[0].expr.boundary
        boundaries.append(bnd)

    if verbose:
        print('> groups     = ', groups)
        print('> boundaries = ', boundaries)
    # ...

    # ...
    groups_M = []
    for bnd, ai in zip(boundaries, groups):
        a_bnd = ai
        for call in ai.atoms(FormCall):
            a_bnd = a_bnd.subs(call, call.expr)

        M_bnd = _evaluate_core(a_bnd, verbose=verbose)

        groups_M.append(BoundaryExpression(bnd, M_bnd))

    if verbose:
        print('> groups_M = ', groups_M)
    # ...

    return groups_M
# ...

def _extract_linear_combination(expr, ls):
    """returns a new expression for terms that are in ls only."""
    # something like a1 + a2 or a1 + alpha * a2
    if isinstance(expr, Add):
        args = []
        for arg in expr.args:
            # somthing like alpha*a4
            if isinstance(arg, Mul):
                m_args = [i for i in arg.args if i in ls]
                if m_args:
                    args += [arg]

            elif arg in ls:
                args += [arg]

        expr = Add(*args)
    return expr

# TODO check that a is a Form, FormCall or linear combination of them
def evaluate(a, verbose=False):
    calls = a.atoms(FormCall)

    bnd_calls = []
    if calls:
        bnd_calls = [a for a in calls if a.expr.boundary]
        calls = [a for a in calls if not(a in bnd_calls)]

        if verbose:
            print('> calls = ', calls)

    expr_bnd = []
    if bnd_calls:
        expr_bnd = _evaluate_bnd(a, bnd_calls, verbose=verbose)

    expr_domain = []
    if calls:
        # TODO - must check that calls have the same domein
        #      - shall we need to add a groupby here too?
        domain = calls[0].expr.domain

        a_expr = a
        if (isinstance(a, BasicForm) and is_sum_of_form_calls(a.expr) and
            bnd_calls):
            a_expr = a.expr

        a = _extract_linear_combination(a_expr, calls)
        if verbose:
            print('> a = ', a)

        # ... replace a FormCall by its expression
        for call in calls:
            a = a.subs(call, call.expr)
        # ...

        expr = _evaluate_core(a, verbose=verbose)
        boundary = list(a.atoms(Boundary))
        if not boundary:
            expr_domain = [DomainExpression(domain, expr)]

        else:
            boundary = boundary[0]
            expr_domain = [BoundaryExpression(boundary, expr)]

    elif isinstance(a, BasicForm):
        domain = a.domain
        expr = _evaluate_core(a, verbose=verbose)
        boundary = list(a.atoms(Boundary))

        expr = _evaluate_core(a, verbose=verbose)
        boundary = list(a.atoms(Boundary))
        if not boundary:
            expr_domain = [DomainExpression(domain, expr)]
        else:
            boundary = boundary[0]
            expr_domain = [BoundaryExpression(boundary, expr)]

    return expr_bnd + expr_domain


class KernelExpression(Basic):
    def __new__(cls, target, expr):
        return Basic.__new__(cls, target, expr)

    @property
    def target(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

class DomainExpression(KernelExpression):
    pass

class BoundaryExpression(KernelExpression):
    pass


# TODO - get dim from atoms
#      - check coefficinets/functions
def _tensorize_core(expr, dim, tests, trials):

    if isinstance(expr, Add):
        args = [_tensorize_core(i, dim, tests, trials) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        args = [i for i in expr.args if not(i in coeffs)]

        d_atoms = {}
        _coordinates = ['x', 'y', 'z']
        test_trial = list(tests) + list(trials)
        for a in test_trial:
            d_atoms[a] = []

            new = S.One
            for i in range(0, dim):
                coord = _coordinates[i]
                Vi = FunctionSpace('V_{}'.format(i),
                                       domain=Line(),
                                       coordinates=[coord])

                ai = TestFunction(Vi, '{test}{i}'.format(test=a.name, i=i))
                d_atoms[a].append(ai)

                new *= ai
            expr = expr.subs({a: new})

        # make sure we have sum of products
        expr = expand(expr)

        # ...
        # TODO - improve this later
        #      - must distinguish between test/trial
        assert(len(tests) == 1)
        assert(len(trials) == 1)

        v = tests[0]
        u = trials[0]

        ops = {'x': dx, 'y': dy, 'z': dz}

        for ui,vi in zip(d_atoms[u], d_atoms[v]):
            coord = ui.space.coordinates.name
            d = ops[coord]

            # ... Mass
            old = vi*ui
            new = Mass(vi,ui)

            expr = expr.subs({old: new})
            # ...

            # ... Stiffness
            old = d(vi)*d(ui)
            new = Stiffness(vi,ui)

            expr = expr.subs({old: new})
            # ...

            # ... Advection
            old = vi*d(ui)
            new = Advection(vi,ui)

            expr = expr.subs({old: new})
            # ...

            # ... Transpose of Advection
            old = d(vi)*ui
            new = AdvectionT(vi,ui)

            expr = expr.subs({old: new})
            # ...


        expr = subs_mul(expr)
        # ...


    return expr


def tensorize(a):

    if not isinstance(a, BilinearForm):
        raise TypeError('Expecting a BilinearForm')

    dim = a.ldim
    expr = a.expr
    tests = a.test_functions
    trials = a.trial_functions

    # ... TODO this is not the best thing to do
    #          we should modify matricize and call it here
    e = evaluate(a)
    if isinstance(e, (Matrix, ImmutableDenseMatrix)):

        assert(len(a.test_spaces) == 1)
        assert(len(a.trial_spaces) == 1)

        assert(len(tests) == 1)
        assert(len(trials) == 1)

        V = a.test_spaces[0]
        U = a.trial_spaces[0]

        v = tests[0]
        u = trials[0]

        from sympde.core import FunctionSpace
        V = FunctionSpace(V.name, V.domain, coordinates=[i.name for i in V.coordinates])
        U = FunctionSpace(U.name, U.domain, coordinates=[i.name for i in U.coordinates])

        vv = TestFunction(V, name=v.name*2)
        uu = TestFunction(U, name=u.name*2)

        tests = [TestFunction(V, name=i.name*2) for i in tests]
        trials = [TestFunction(U, name=i.name*2) for i in trials]

        print(e)
        raise NotImplementedError('inv_normalize has been removed')
#        expr = inv_normalize(e, {v: vv, u: uu})
        print(expr)
        n_rows, n_cols = expr.shape
        lines = []
        for i in range(0, n_rows):
            line = []
            for j in range(0, n_cols):
                eij = _tensorize_core(expr[i,j], dim, tests, trials)
                line.append(eij)
            lines.append(line)
        expr = Matrix(lines)

    else:

        expr = atomize(expr)
        expr = _tensorize_core(expr, dim, tests, trials)
    # ...

    return expr

def subs_mul(expr):
    """substitute Mul with TensorProduct"""

    if isinstance(expr,(Add, Mul)):
        args = expr.args
        args = [subs_mul(arg) for arg in args]

        if isinstance(expr, Mul):

            return TensorProduct(*args)
        elif isinstance(expr, Add):

            return Add(*args)
    else:

        return expr


def subs_bilinear_form(form, newargs):
    # ...
    test_trial = _sanitize_form_arguments(newargs, form, is_bilinear=True)

    if not isinstance(test_trial, (tuple, list, Tuple)):
        raise TypeError('(test, trial) must be a tuple, list or Tuple')

    if not(len(test_trial) == 2):
        raise ValueError('Expecting a couple (test, trial)')
    # ...

    # ...
    test_functions = test_trial[0]
    if isinstance(test_functions, (TestFunction, VectorTestFunction)):
        test_functions = [test_functions]

    elif isinstance(test_functions, (tuple, list, Tuple)):
        test_functions = Tuple(*test_functions)
    # ...

    # ...
    trial_functions = test_trial[1]
    if isinstance(trial_functions, (TestFunction, VectorTestFunction)):
        trial_functions = [trial_functions]

    elif isinstance(trial_functions, (tuple, list, Tuple)):
        trial_functions = Tuple(*trial_functions)
    # ...

    # in order to avoid problems when swapping indices, we need to create
    # temp symbols

    domain = form.domain

    # ...
    d_tmp = {}
    for x in trial_functions:
        name = 'trial_{}'.format(abs(hash(x)))
        if isinstance(x, VectorTestFunction):
            X = VectorUnknown(name, domain, shape=x.shape)

        else:
            X = Unknown(name, domain)

        d_tmp[X] = x
    # ...

    expr = form.expr

    # ... replacing trial functions by tmp symbols
    d = {}
    for k,v in zip(form.trial_functions, d_tmp):
        d[k] = v
    expr = expr.subs(d)
    # ...

    # ... replacing test functions
    d = {}
    for k,v in zip(form.test_functions, test_functions):
        d[k] = v
    expr = expr.subs(d)
    # ...

    # ... replacing trial functions from tmp symbols
    expr = expr.subs(d_tmp)
    # ...

    # ...
    if len(test_functions) == 1: test_functions = test_functions[0]
    if len(trial_functions) == 1: trial_functions = trial_functions[0]

    test_trial = (test_functions, trial_functions)
    # ...

    return BilinearForm(test_trial, expr, measure=form.measure,
                        mapping=form.mapping, name=form.name)


def subs_linear_form(form, newargs):
    # ...
    test_functions = _sanitize_form_arguments(newargs, form, is_linear=True)
    test_functions = test_functions[0]

    if isinstance(test_functions, (TestFunction, VectorTestFunction)):
        test_functions = [test_functions]

    elif isinstance(test_functions, (tuple, list, Tuple)):
        test_functions = list(*test_functions)
    # ...

    expr = form.expr

    # ... replacing test functions
    d = {}
    for k,v in zip(form.test_functions, test_functions):
        d[k] = v
    expr = expr.subs(d)
    # ...

    if len(test_functions) == 1: test_functions = test_functions[0]

    return LinearForm(test_functions, expr, measure=form.measure,
                      mapping=form.mapping, name=form.name)
