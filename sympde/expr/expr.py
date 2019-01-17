# coding: utf-8

# TODO - transpose of BilinearForm
#      - add unknown status if only one space is given to the BilinearForm
#        => can not be evaluated except if it is called on a couple test/trial
#      - check that a BilinearForm is bilinear (using Constant)
#      - check that a LinearForm is linear
#      - add is_symmetric property for BilinearForm

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
from .errors import UnconsistentLinearFormError
from .errors import UnconsistentBilinearFormError


#==============================================================================
def _initialize_measure(measure, coordinates):
    if not( measure is None ):
        return measure

    if not( coordinates is None ):
        return Measure(coordinates)

    else:
        raise NotImplementedError('')

#==============================================================================
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

        boundary = Union(*boundaries)

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


#==============================================================================
class BasicForm(Expr):
    _name = None
    _boundary = None

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


#==============================================================================
# TODO we should check that the only free symbols are fields, constants or coordinates
class Integral(BasicForm):
    """

    Examples

    """
    _ldim = None
    _coordinates = None
    def __new__(cls, expr, domain, measure=None, name=None):
        # ... treat union of domains
        #     TODO improve
        unions = expr.atoms(Union)
        if unions:
            if not( len(unions) == 1 ):
                raise NotImplementedError('only one union is available for the moment')

            domains = []
            for i in list(unions):
                domains += list(i._args)
            domains = set(domains)
            if len(domains) > 1:
                forms = []
                for domain in domains:
                    i = list(unions)[0]
                    _expr = expr.replace(i, domain)
                    form  = Integral(_expr, domain, measure=measure,
                                      name=None)

                    forms.append(form)

                expr = Add(*forms)
        # ...


        if not isinstance(domain, BasicDomain):
            raise TypeError('> Expecting a BasicDomain object for domain')

        # ... check that there are no test functions in the expression
        ls = [a for a in expr.free_symbols if isinstance(a, (TestFunction, VectorTestFunction))]
        if not(len(ls) == 0):
            raise TypeError('Cannot use test functions in Integral')
        # ...

        # ...
        coordinates = domain.coordinates
        ldim = domain.dim
        # ...

        # compute dim from fields if available
        ls = list(expr.atoms((Field, VectorField)))
        if ls:
            F = ls[0]
            space = F.space

        else:
            tag = random_string( 3 )
            space_name = 'space_{}'.format(tag)
            space = FunctionSpace(space_name, domain)
            # TODO vector case

        # check if we are using a mapping
        mapping = None
        if isinstance( domain, MappedDomain ):
            mapping = domain.mapping

        measure = _initialize_measure(measure, coordinates)

        # get boundary terms
        boundary = _initialize_boundary(expr)

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


#==============================================================================
class LinearForm(BasicForm):
    """

    Examples

    """
    def __new__(cls, arguments, expr, measure=None, name=None, check=False):
        # ... treat union of domains
        #     TODO improve
        unions = expr.atoms(Union)
        if unions:
            if not( len(unions) == 1 ):
                raise NotImplementedError('only one union is available for the moment')

            domains = []
            for i in list(unions):
                domains += list(i._args)
            domains = set(domains)
            if len(domains) > 1:
                forms = []
                for domain in domains:
                    i = list(unions)[0]
                    _expr = expr.replace(i, domain)
                    _name = None
                    if not( name is None ):
                        _name = '{name}_{domain}'.format(name=name, domain=domain.name)

                    form  = LinearForm(arguments, _expr, measure=measure,
                                       name=_name)

                    forms.append(form(arguments))

                expr = Add(*forms)
        # ...

        # ...
        calls = list(expr.atoms(FormCall))
        if check and not calls:
            if not is_linear_form(expr, arguments):
                msg = '> Expression is not bilinear'
                raise UnconsistentLinearFormError(msg)
        # ...

        args = _sanitize_form_arguments(arguments, expr, is_linear=True)
        obj = Basic.__new__(cls, args, expr)

        # TODO must check that all domains are the same
        domain = obj.test_spaces[0].domain

        # check if we are using a mapping
        mapping = None
        if isinstance( domain, MappedDomain ):
            mapping = domain.mapping

        measure = _initialize_measure(measure, obj.coordinates)

        # get boundary terms
        boundary = _initialize_boundary(expr)

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

        if isinstance(self.expr, FormCall):
            call = self.expr
            return FormCall(call.expr, args)

        return FormCall(self, args)



#==============================================================================
class BilinearForm(BasicForm):
    """

    Examples

    """

    def __new__(cls, arguments, expr, measure=None, name=None, check=False):
        # ... treat union of domains
        #     TODO improve
        unions = expr.atoms(Union)
        if unions:
            if not( len(unions) == 1 ):
                raise NotImplementedError('only one union is available for the moment')

            domains = []
            for i in list(unions):
                domains += list(i._args)
            domains = set(domains)
            if len(domains) > 1:
                forms = []
                for domain in domains:
                    i = list(unions)[0]
                    _expr = expr.replace(i, domain)
                    form  = BilinearForm(arguments, _expr, measure=measure,
                                         name=None)

                    forms.append(form(*arguments))

                expr = Add(*forms)
        # ...

        # ...
        if not isinstance(arguments, (tuple, list, Tuple)):
            raise TypeError('(test, trial) must be a tuple, list or Tuple')

        if not(len(arguments) == 2):
            raise ValueError('Expecting a couple (test, trial)')
        # ...

        # ...
        calls = list(expr.atoms(FormCall))
        if check and not calls:
            if not is_bilinear_form(expr, arguments):
                msg = '> Expression is not bilinear'
                raise UnconsistentBilinearFormError(msg)
        # ...

        args = _sanitize_form_arguments(arguments, expr, is_bilinear=True)
        obj = Basic.__new__(cls, args, expr)

        # TODO must check that all domains are the same
        domain = obj.test_spaces[0].domain

        # check if we are using a mapping
        mapping = None
        if isinstance( domain, MappedDomain ):
            mapping = domain.mapping

        measure = _initialize_measure(measure, obj.coordinates)

        # get boundary terms
        boundary = _initialize_boundary(expr)

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

        # ...
        expr = self
        if isinstance(self.expr, FormCall):
            expr = self.expr.expr
        # ...

        are_dummy = lambda xs: [isinstance(i, (TestFunction, VectorTestFunction))
                                for i in xs]

        # ...
        test_functions  = args[0]
        if not isinstance(test_functions, (tuple, list, Tuple)):
            test_functions = [test_functions]
        test_functions = Tuple(*test_functions)
        # ...

        # ...
        trial_functions = args[1]
        if not isinstance(trial_functions, (tuple, list, Tuple)):
            trial_functions = [trial_functions]
        trial_functions = Tuple(*trial_functions)
        # ...

        # ...
        if not all(are_dummy(test_functions)) and not all(are_dummy(trial_functions)):
            # This should return an Integral
            raise NotImplementedError('')

        else:
            free_variables = None
            target = None
            if not any(are_dummy(test_functions)):
                assert(all(are_dummy(trial_functions)))

                args           = trial_functions
                free_variables = test_functions
                target         = 'test'

            elif not any(are_dummy(trial_functions)):
                assert(all(are_dummy(test_functions)))

                args           = test_functions
                free_variables = trial_functions
                target         = 'trial'

            if not( free_variables is None ):
                if target == 'test':
                    target = expr.variables[0]

                elif target == 'trial':
                    target = expr.variables[1]

                expr = expr.expr

                for old, new in zip(target, free_variables):
                    expr = expr.subs(old, new)

                expr = LinearForm(args, expr, check=False)
        # ...

        return FormCall(expr, args)


#==============================================================================
class Norm(Integral):
    def __new__(cls, expr, domain, kind='l2', measure=None,
                name=None):
        # ...
        tests = expr.atoms((TestFunction, VectorTestFunction))
        if tests:
            msg = '> Expecting an Expression without test functions'
            raise UnconsistentArgumentsError(msg)

        if not isinstance(expr, (Expr, Matrix, ImmutableDenseMatrix)):
            msg = '> Expecting Expr, Matrix, ImmutableDenseMatrix'
            raise UnconsistentArgumentsError(msg)

        # ...

        # ...
        if not(kind in ['l2', 'h1']):
            raise ValueError('> Only L2, H1 norms are available')
        # ...

        # ...
        if name is None:
            name = random_string( 3 )

        name = '{kind}norm_{name}'.format(kind=kind, name=name)
        # ...

        # ...
        is_vector = isinstance(expr, (Matrix, Tuple, list, tuple))
        if is_vector:
            expr = Matrix(expr)
        # ...

        # ...
        exponent = None
        if kind == 'l2':
            exponent = 2

            if not is_vector:
                expr = expr*expr

            else:
                if not( expr.shape[1] == 1 ):
                    raise ValueError('Wrong expression for Matrix. must be a row')

                v = Tuple(*expr[:,0])
                expr = Dot(v, v)

        elif kind == 'h1':
            exponent = 2

            if not is_vector:
                expr = Dot(Grad(expr), Grad(expr))

            else:
                if not( expr.shape[1] == 1 ):
                    raise ValueError('Wrong expression for Matrix. must be a row')

                v = Tuple(*expr[:,0])
                expr = Inner(Grad(v), Grad(v))
        # ...

        obj = Integral.__new__(cls, expr, domain, measure=name, name=name)

        obj._exponent = exponent

        return obj

    @property
    def exponent(self):
        return self._exponent


#==============================================================================
class BilinearAtomicForm(BilinearForm, AtomicExpr):
    """

    Examples

    """

    def _sympystr(self, printer):
        sstr = printer.doprint
        name = sstr(self.name)

        test = [sstr(i) for i in self.test_functions]
        test = ','.join(i for i in test)

        trial = [sstr(i) for i in self.trial_functions]
        trial = ','.join(i for i in trial)

        return '{name}({test},{trial})'.format(name=name, trial=trial, test=test)

#==============================================================================
class Mass(BilinearAtomicForm):
    """

    Examples

    """
    def __new__(cls, test, trial):

        test_trial = [test, trial]
        expr = test * trial

        return BilinearForm.__new__(cls, test_trial, expr, name='Mass')

#==============================================================================
class Stiffness(BilinearAtomicForm):
    """

    Examples

    """
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = d(test) * d(trial)

        return BilinearForm.__new__(cls, test_trial, expr, name='Stiffness')

#==============================================================================
class Advection(BilinearAtomicForm):
    """

    Examples

    """
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = test * d(trial)

        return BilinearForm.__new__(cls, test_trial, expr, name='Advection')

#==============================================================================
class AdvectionT(BilinearAtomicForm):
    """

    Examples

    """
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = d(test) * trial

        return BilinearForm.__new__(cls, test_trial, expr, name='AdvectionT')

#==============================================================================
class Bilaplacian(BilinearAtomicForm):
    """

    Examples

    """
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = d(d(test)) * d(d(trial))

        return BilinearForm.__new__(cls, test_trial, expr, name='Bilaplacian')


#==============================================================================
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


#==============================================================================
class FormCall(AtomicExpr):

    is_commutative = False

    def __new__(cls, expr, args):

        # ...
        if not isinstance(args, (list, tuple, Tuple)):
            args = [args]
        # ...

        # ...
        name = expr.name
        if isinstance( expr, BilinearForm ):
            expr = subs_form(expr, args)

            if not isinstance(expr, BilinearForm):
                expr = BilinearForm(args, expr, name=name, check=False)

        elif isinstance( expr, LinearForm ):
            expr = subs_form(expr, args)

            if not isinstance(expr, LinearForm):
                expr = LinearForm(args, expr, name=name, check=False)

        else:
            raise TypeError('> Expecting BilinearForm, LinearForm')
        # ...

        # ...
        args = Tuple(*args)
        obj = Basic.__new__(cls, expr, args)
        # ...

        return obj

    @property
    def expr(self):
        return self._args[0]

    @property
    def arguments(self):
        return self._args[1]


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
    if isinstance(expr, FormCall):  return True
    if not isinstance(expr, Add):   return False

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
                             VectorField,
                             BoundaryVector, Trace,
                             Integer, Float, Matrix, ImmutableDenseMatrix,
                             list, tuple, Tuple)):
        msg = ('> Wrong input type.')

        raise TypeError(msg, ', given ', expr, type(expr))

    # ... replace a FormCall by its expression
    calls = expr.atoms(FormCall)
    for call in calls:
        expr = expr.subs(call, call.expr)
#        expr = expr.subs(call, call.expr, simultaneous=True)
#        expr = expr.xreplace({call: call.expr})
#        expr = expr.replace(call, call.expr)
    # ...

#    print('> expr [atomize] = ', expr, type(expr))

    # ... compute dim if None
    if dim is None:
        ls = [i for i in expr.free_symbols if isinstance(i, (TestFunction,
                                                             VectorTestFunction,
                                                             Field,
                                                             VectorField))]

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

    elif isinstance(expr, _generic_ops):
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
        # TODO treat Mul node
        if isinstance(a.boundary, Union):
            if isinstance(a.expr, Add):
                newargs = []
                for i in a.expr.args:

                    expr = i.expr
                    if isinstance(expr, BasicForm):
                        expr = expr.expr

                    if is_sum_of_form_calls(expr):
                        newargs += expr.args

                    else:
                        newargs.append(i)

                a_expr = Add(*newargs)

            elif isinstance(a.expr, Mul):
                raise NotImplementedError('')

            else:
                a_expr = a.expr

        else:
            a_expr = a.expr

    if isinstance(a_expr, FormCall):
        a_expr = a_expr.expr.expr

    # ...
    boundaries = []
    groups = []
    keyfunc = lambda call: call.expr.boundary
    for bnd, g in groupby(bnd_calls, keyfunc):
        ls = list(g)
        a_bnd = _extract_linear_combination(a_expr, ls)
        groups.append(a_bnd)
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
    # ...
    _calls = a.atoms(FormCall)
    calls = []
    for call in _calls:
        mycalls = call.expr.atoms(FormCall)
        if mycalls:
            calls += list(mycalls)
        else:
            calls += [call]

    # remove redundancy
    calls = list(set(calls))
    # ...

    bnd_calls = []
    if calls:
        bnd_calls = [a for a in calls if a.expr.boundary]
        calls = [a for a in calls if not(a in bnd_calls)]

        if verbose:
            print('> calls = ', calls)

    expr_bnd = []
    if bnd_calls:
        expr_bnd = _evaluate_bnd(a, bnd_calls, verbose=verbose)

    bnd_done = [i.target for i in expr_bnd]

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
        boundary = list(a.atoms(Boundary))

        if not boundary:
            domain = a.domain
            expr = _evaluate_core(a, verbose=verbose)
            expr_domain = [DomainExpression(domain, expr)]

        # TODO nitsch case
#        else:
#            expr_domain = []
#            boundary = [i for i in boundary if not(i in bnd_done)]
#            for bnd in boundary:
#                expr_domain += [BoundaryExpression(bnd, expr)]

    return expr_bnd + expr_domain


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
        _coordinates = [Symbol(i) for i in _coordinates]
        test_trial = list(tests) + list(trials)
        for a in test_trial:
            d_atoms[a] = []

            new = S.One
            for i in range(0, dim):
                coord = _coordinates[i]
                Di = Interval(coordinate=coord)
                Vi = FunctionSpace('V_{}'.format(i), domain=Di)

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

            # ... Bilaplacian
            old = d(d(vi))*d(d(ui))
            new = Bilaplacian(vi,ui)

            expr = expr.subs({old: new})
            # ...

        expr = subs_mul(expr)
        # ...

    return expr

#==============================================================================
def _tensorize_weights(expr):

    if isinstance(expr, Add):
        args = []
        for term in expr.args:
            #print('> ', term, type(term))
            arg = _tensorize_weights(term)
            args.append(arg)
        expr = Add(*args)

    elif isinstance(expr, Mul):
        args = []
        for term in expr.args:
#            print('>> ', term, type(term))
            arg = _tensorize_weights(term)
            args.append(arg)

        tensor = [a for a in args if isinstance(a, TensorProduct)]
        weights = [a for a in args if not( a in tensor )]

        if tensor:

            tensor = tensor[0]
            forms = tensor.args

            coords = [a.coordinates for a in forms]

    #        print(forms)
    #        print(coords)
    #        print(weights)

    #        # ...
    #        d_args = {}
    #        for x in coords:
    #            d_args[x] = []
    #
    #        for x in coords:
    #            for a in weights:
    #                # TODO improve for functions => separability
    #                ls = a.atoms(Symbol)
    #                if x in ls:
    #                    print('found ', x, ' in ', a)
    #        # ...

        expr = Mul(*args)

    elif isinstance(expr, TensorProduct):
        args = []
        for term in expr.args:
#            print('>>> ', term, type(term))
            arg = _tensorize_weights(term)
#            if not( arg is S.One ):
#                args.append(arg)
            if isinstance(term, BilinearAtomicForm):
                coords = term.domain.coordinates
                #print(coords)
        expr = TensorProduct(*args)

    return expr

#==============================================================================
def tensorize(a):

    if not isinstance(a, BilinearForm):
        raise TypeError('Expecting a BilinearForm')

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

    if is_sum_of_form_calls(a.expr):
        # ...
        n_rows, test_indices = _get_size_and_starts(a.test_functions)
        n_cols, trial_indices = _get_size_and_starts(a.trial_functions)

        lines = []
        for i in range(0, n_rows):
            line = []
            for j in range(0, n_cols):
                line.append(0)
            lines.append(line)

        M = Matrix(lines)
        # ...

        calls = a.atoms(FormCall)
        for call in calls:
            t = tensorize(call.expr)

            atoms = call.arguments
            i_row = None
            i_col = None
            l_row = 1
            l_col = 1
            for atom in atoms:
                if atom in test_indices:
                    i_row = test_indices[atom]

                    if isinstance(atom, VectorTestFunction):
                        l_row = atom.shape[0]

                elif atom in trial_indices:
                    i_col = trial_indices[atom]

                    if isinstance(atom, VectorTestFunction):
                        l_col = atom.shape[0]

                else:
                    raise ValueError('> Could not find {}'.format(atom))

            if isinstance(t, (Matrix, ImmutableDenseMatrix)):
                M_loc = M[i_row:i_row+l_row, i_col:i_col+l_col]
                if M_loc.shape == t.shape:
                    M_loc += t

                # TODO must check if a trial/test is used as test/trial
                elif M_loc.shape == t.shape[::-1]:
                    M_loc += t.transpose()

                else:
                    raise ValueError('Wrong sizes')

                M[i_row:i_row+l_row, i_col:i_col+l_col] += M_loc

            else:
                raise NotImplementedError('TODO')

        return M

    dim = a.ldim
    domain = a.domain
    tests = a.test_functions
    trials = a.trial_functions

    assert(len(a.test_spaces) == 1)
    assert(len(a.trial_spaces) == 1)
    assert(len(tests) == 1)
    assert(len(trials) == 1)

    V = a.test_spaces[0]
    U = a.trial_spaces[0]

    # the result of evaluate is a list of KernelExpression
    kernels = evaluate(a)
    kernels = [i.expr for i in kernels]

    expressions = []
    for kernel in kernels:
        if isinstance(kernel, (Matrix, ImmutableDenseMatrix)):

            n_rows, n_cols = kernel.shape

            # ... subs indexed test/trial functions by a new symbol
            tmp_tests = []
            tmp_trials = []
            for i_row in range(0, n_rows):
                for i_col in range(0, n_cols):
                    e = kernel[i_row,i_col]
                    indexed = e.atoms(IndexedTestTrial)
                    for ui in indexed:
                        i = ui.indices
                        if len(i) > 1:
                            raise ValueError('Expecting one index')
                        i = i[0]

                        space_name = 'VTmp{}'.format(i)
                        Vi = FunctionSpace(space_name, V.domain)
                        vi = TestFunction(Vi, '{test}{i}'.format(test=ui.base.name, i=i))
                        e = e.subs({ui: vi})

                        if ui.base in tests:
                            tmp_tests.append(vi)

                        elif ui.base in trials:
                            tmp_trials.append(vi)

                    kernel[i_row,i_col] = e
            # ...

            tmp_tests += [i for i in tests if not(i in tmp_tests)]
            tmp_trials += [i for i in trials if not(i in tmp_trials)]

            # ...
            lines = []
            for i_row in range(0, n_rows):
                line = []
                for i_col in range(0, n_cols):
                    e = kernel[i_row,i_col]

                    atoms = e.atoms(TestFunction)
                    _tests = [i for i in atoms if i in tmp_tests]
                    _trials = [i for i in atoms if i in tmp_trials]

                    eij = _tensorize_core(e, dim, _tests, _trials)

                    line.append(eij)

                lines.append(line)

            expr = Matrix(lines)
            # ...

        else:
            expr = _tensorize_core(kernel, dim, tests, trials)

        expressions.append(expr)
    # ...

    # TODO
    # looking for weighted atomic forms
    # this should be done if a flag is True
    # and used for LinearOperator Kron
#    expr = _tensorize_weights(expr)

    return expr

#==============================================================================
def subs_mul(expr):
    """substitute Mul with TensorProduct"""

    if isinstance(expr,(Add, Mul)):
        args = expr.args
        args = [subs_mul(arg) for arg in args]

        if isinstance(expr, Mul):
            args = expr.args
            forms = [i for i in args if isinstance(i, BilinearAtomicForm)]
            others = [i for i in args if not( i in forms )]

            t = TensorProduct(*forms)
            m = Mul(*others)
            return m*t

        elif isinstance(expr, Add):

            return Add(*args)
    else:

        return expr

#==============================================================================
# form here is a BilinearForm
def subs_form(form, newargs):
#    print('>>> subs_form : ', form)

    if isinstance( form, BilinearForm ):
        calls = form.expr.atoms(FormCall)
        if calls:
            return subs_form(form.expr, newargs)

        else:
            return _subs_bilinear_form_core(form, newargs)

    elif isinstance( form, LinearForm ):
        calls = form.expr.atoms(FormCall)
        if calls:
            return subs_form(form.expr, newargs)

        else:
            return _subs_linear_form_core(form, newargs)

    elif isinstance( form, FormCall ):
        return form

    elif isinstance( form, Add ):
        args = [subs_form(a, newargs) for a in form.args]
        return Add(*args)

    elif isinstance( form, Mul ):
        coeffs  = [i for i in form.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in form.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [subs_form(a, newargs) for a in vectors]
            j = Mul(*args)

        return Mul(i, j)

    elif isinstance( form, Pow ):

        b = subs_form(form.base, newargs)
        e = form.exp

        return Pow(b, e)

    else:
        return form


#==============================================================================
def _subs_bilinear_form_core(form, newargs):
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

    # ...
    d_tmp = {}
    for x in trial_functions:
        name = random_string( 6 )
        if isinstance(x, TestFunction):
            X = TestFunction(x.space, name=name)

        elif isinstance(x, VectorTestFunction):
            X = VectorTestFunction(x.space, name=name)

        else:
            raise TypeError('Only TestFunction and VectorTestFunction are available')

        d_tmp[X] = x
    # ...

    expr = form.expr

    # ... replacing trial functions by tmp symbols
    for k,v in zip(form.trial_functions, d_tmp):
        expr = expr.subs(k,v)
    # ...

    # ... replacing test functions
    for k,v in zip(form.test_functions, test_functions):
        expr = expr.subs(k,v)
    # ...

    # ... replacing trial functions from tmp symbols
    for k,v in d_tmp.items():
        expr = expr.subs(k,v)
    # ...

    # ...
    if len(test_functions)  == 1: test_functions  = test_functions[0]
    if len(trial_functions) == 1: trial_functions = trial_functions[0]

    test_trial = (test_functions, trial_functions)
    # ...

    return BilinearForm(test_trial, expr, name=form.name, check=False)


#==============================================================================
def _subs_linear_form_core(form, newargs):

    # ...
    test_functions = _sanitize_form_arguments(newargs, form, is_linear=True)
    # TODO is it ok to do this?
    test_functions = test_functions[0]

    if isinstance(test_functions, (TestFunction, VectorTestFunction)):
        test_functions = [test_functions]

    elif isinstance(test_functions, (tuple, list, Tuple)):
        test_functions = list(*test_functions)
    # ...

    expr = form.expr

    # ... replacing test functions
    for k,v in zip(form.test_functions, test_functions):
        expr = expr.subs(k,v)
    # ...

    if len(test_functions) == 1: test_functions = test_functions[0]

    return LinearForm(test_functions, expr, name=form.name, check=False)

#==============================================================================
def is_bilinear_form(expr, args):
    """checks if an expression is bilinear with respect to the given arguments."""
    # ...
    test_trial = _sanitize_form_arguments(args, expr, is_bilinear=True)

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

    # ...
    newargs = []
    d_newargs = {}
    for x in test_functions:
        tag    = random_string( 4 )
        c_left  = Constant('alpha_' + tag)
        c_right = Constant('beta_'  + tag)

        if isinstance(x, TestFunction):
            left  = TestFunction(x.space, name='l_' + tag)
            right = TestFunction(x.space, name='r_' + tag)

        elif isinstance(x, VectorTestFunction):
            left  = VectorTestFunction(x.space, name='l_' + tag)
            right = VectorTestFunction(x.space, name='r_' + tag)

        else:
            raise TypeError('Only TestFunction and VectorTestFunction are available')

        newargs += [c_left*left + c_right*right]
        d_newargs[x] = {'left':  (c_left, left),
                        'right': (c_right, right)}

    test_newargs   = newargs
    d_test_newargs = d_newargs
    # ...

    # ...
    newargs = []
    d_newargs = {}
    for x in trial_functions:
        tag    = random_string( 4 )
        c_left  = Constant('alpha_' + tag)
        c_right = Constant('beta_'  + tag)

        if isinstance(x, TestFunction):
            left  = TestFunction(x.space, name='l_' + tag)
            right = TestFunction(x.space, name='r_' + tag)

        elif isinstance(x, VectorTestFunction):
            left  = VectorTestFunction(x.space, name='l_' + tag)
            right = VectorTestFunction(x.space, name='r_' + tag)

        else:
            raise TypeError('Only TestFunction and VectorTestFunction are available')

        newargs += [c_left*left + c_right*right]
        d_newargs[x] = {'left':  (c_left, left),
                        'right': (c_right, right)}

    trial_newargs   = newargs
    d_trial_newargs = d_newargs
    # ...

    # ... replacing test functions for a1*u1+a2*u2
    test_newexpr = expr
    for k,v in zip(test_functions, test_newargs):
        test_newexpr = test_newexpr.subs(k,v)
    # ...

    # ... replacing test functions for a1*u1
    test_left_expr = expr
    for old in test_functions:
        c,u = d_test_newargs[old]['left']
        test_left_expr = test_left_expr.subs(old,c*u)
    # ...

    # ... replacing test functions for a2*u2
    test_right_expr = expr
    for old in test_functions:
        c,u = d_test_newargs[old]['right']
        test_right_expr = test_right_expr.subs(old,c*u)
    # ...

    # ... replacing trial functions for a1*u1+a2*u2
    trial_newexpr = expr
    for k,v in zip(trial_functions, trial_newargs):
        trial_newexpr = trial_newexpr.subs(k,v)
    # ...

    # ... replacing trial functions for a1*u1
    trial_left_expr = expr
    for old in trial_functions:
        c,u = d_trial_newargs[old]['left']
        trial_left_expr = trial_left_expr.subs(old,c*u)
    # ...

    # ... replacing trial functions for a2*u2
    trial_right_expr = expr
    for old in trial_functions:
        c,u = d_trial_newargs[old]['right']
        trial_right_expr = trial_right_expr.subs(old,c*u)
    # ...

    test_condition  = expand(test_newexpr) == test_left_expr + test_right_expr
    trial_condition = expand(trial_newexpr) == trial_left_expr + trial_right_expr
    if test_condition and trial_condition:
        return True

    return False

#==============================================================================
def is_linear_form(expr, args):
    """checks if an expression is linear with respect to the given arguments."""
    # ...
    test_functions = _sanitize_form_arguments(args, expr, is_linear=True)
    # TODO is it ok to do this?
    test_functions = test_functions[0]

    if isinstance(test_functions, (TestFunction, VectorTestFunction)):
        test_functions = [test_functions]

    elif isinstance(test_functions, (tuple, list, Tuple)):
        test_functions = list(*test_functions)
    # ...

    # ...
    newargs = []
    d_newargs = {}
    for x in test_functions:
        tag    = random_string( 4 )
        c_left  = Constant('alpha_' + tag)
        c_right = Constant('beta_'  + tag)

        if isinstance(x, TestFunction):
            left  = TestFunction(x.space, name='l_' + tag)
            right = TestFunction(x.space, name='r_' + tag)

        elif isinstance(x, VectorTestFunction):
            left  = VectorTestFunction(x.space, name='l_' + tag)
            right = VectorTestFunction(x.space, name='r_' + tag)

        else:
            raise TypeError('Only TestFunction and VectorTestFunction are available')

        newargs += [c_left*left + c_right*right]
        d_newargs[x] = {'left':  (c_left, left),
                        'right': (c_right, right)}
    # ...

    # ... replacing test functions for a1*u1+a2*u2
    newexpr = expr
    for k,v in zip(test_functions, newargs):
        newexpr = newexpr.subs(k,v)
    # ...

    # ... replacing test functions for a1*u1
    left_expr = expr
    for old in test_functions:
        c,u = d_newargs[old]['left']
        left_expr = left_expr.subs(old,c*u)
    # ...

    # ... replacing test functions for a2*u2
    right_expr = expr
    for old in test_functions:
        c,u = d_newargs[old]['right']
        right_expr = right_expr.subs(old,c*u)
    # ...

    if expand(newexpr) == left_expr + right_expr:
        return True

    return False


#==============================================================================
def linearize(form, fields, trials=None):
    """linearize a LinearForm around the fields."""
    # ...
    if not isinstance(form, LinearForm):
        raise TypeError('> Expecting a LinearForm')

    if not isinstance(fields, (list, tuple, Tuple)):
        fields = [fields]

    for f in fields:
        if not isinstance(f, (Field, VectorField)):
            raise TypeError('{} is not Field/VectorField'.format(f))

    if not(trials is None):
        if not isinstance(trials, (list, tuple, Tuple)):
            trials = [trials]

        assert( all([isinstance(i, str) for i in trials]) )
        assert( len(fields) == len(trials) )
    # ...

    expr           = form.expr
    test_functions = form.test_functions
    fields          = Tuple(*fields)

    # ...
    calls = list(expr.atoms(FormCall))
    if calls:
        raise NotImplementedError()
    # ...

    # ...
    trial_functions = []
    coeffs          = []
    newargs         = []
    for i,x in enumerate(fields):
        tag  = random_string( 4 )
        eps  = Constant('eps_' + tag)

        if trials is None:
            name = x.name + '_' + tag
        else:
            name = trials[i]

        if isinstance(x, Field):
            trial  = TestFunction(x.space, name=name)

        elif isinstance(x, VectorField):
            trial  = VectorTestFunction(x.space, name=name)

        else:
            raise TypeError('Only TestFunction and VectorTestFunction are available')

        newargs         += [x + eps*trial]
        trial_functions += [trial]
        coeffs          += [eps]
    # ...

    # ... replacing test functions for a1*u1+a2*u2
    newexpr = expr
    for k,v in zip(fields, newargs):
        newexpr = newexpr.subs(k,v)
    # ...

    newexpr = expand(newexpr)

    args = []
    for eps in coeffs:
        e = newexpr.series(eps,0,2)
        d = collect(e, eps, evaluate=False)
        args += [d[eps]]

    expr = Add(*args)

    test_trial = (trial_functions, test_functions)
    return BilinearForm(test_trial, expr, check=True)
