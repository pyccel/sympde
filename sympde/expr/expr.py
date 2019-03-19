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
from sympde.topology.space import ScalarTestFunction
from sympde.topology.space import VectorTestFunction
from sympde.topology.space import IndexedTestTrial
from sympde.topology.space import Unknown, VectorUnknown
from sympde.topology.space import Trace
from sympde.topology.space import ScalarField, VectorField, IndexedVectorField
from sympde.topology.measure import CanonicalMeasure
from sympde.topology.measure import CartesianMeasure
from sympde.topology.measure import Measure

from .errors import UnconsistentLinearExpressionError
from .basic import BasicForm
from .basic  import BasicExpr
from .basic  import is_linear_form, _sanitize_arguments



#==============================================================================
class LinearExpr(BasicExpr):
    is_linear = True

    def __new__(cls, arguments, expr, check=False):

        # ...
        if expr.atoms(BasicIntegral):
            raise TypeError('')
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

    def __call__(self, *args):
        args = _sanitize_arguments(args, is_linear=True)
        args = Tuple(*args)
        return self.expr.xreplace(dict(list(zip(self.variables, args))))

    def _eval_nseries(self, x, n, logx):
        return self.expr._eval_nseries(x, n, logx)

#==============================================================================
class BilinearExpr(BasicExpr):
    is_bilinear = True

    def __new__(cls, arguments, expr, check=False):

        # ...
        if expr.atoms(BasicIntegral):
            raise TypeError('')
        # ...

        # ...
        if not isinstance(arguments, (tuple, list, Tuple)):
            raise TypeError('(test, trial) must be a tuple, list or Tuple')

        if not(len(arguments) == 2):
            raise ValueError('Expecting a couple (test, trial)')
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

    def __call__(self, *args):
        args = _sanitize_arguments(args, is_bilinear=True)
        left,right = args
        if not is_sequence(left):
            left = [left]

        if not is_sequence(right):
            right = [right]

        args = Tuple(*left, *right)

        variables = Tuple(*self.variables[0], *self.variables[1])
        return self.expr.xreplace(dict(list(zip(self.variables, args))))


#==============================================================================
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

        if not isinstance(expr, (BasicExpr, Expr)):
            raise TypeError('')

        if isinstance(expr, BasicExpr):
            expr = expr.expr

        if isinstance(expr, Add):
            args = [cls.eval(a) for a in expr.args]
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

        return cls(expr, evaluate=False)

#==============================================================================
def _get_domain(a):
    # expr is an integral of BasicExpr or Add of Integral of BasicExpr
    if isinstance(a, BoundaryIntegral):
        domains = list(a.atoms(Boundary))
        if len(domains) == 1:
            return domains[0]

        else:
            return domains

    elif isinstance(a, DomainIntegral):
        expr = a._args[0]

    else:
        expr = a

    atoms  = []
    atoms += list(expr.atoms(ScalarTestFunction))
    atoms += list(expr.atoms(VectorTestFunction))
    atoms += list(expr.atoms(ScalarField))
    atoms += list(expr.atoms(VectorField))

    if len(atoms) == 0:
        raise ValueError('could not find any test function or field')

    space = atoms[0].space

    domains = list(expr.atoms(Boundary)) + [space.domain]
    if len(domains) == 1:
        return domains[0]

    else:
        return domains


#==============================================================================
class Functional(BasicForm):
    is_functional = True

    def __new__(cls, expr, domain):

        expr = BasicIntegral(expr)
        obj = Basic.__new__(cls, expr, domain)

        # compute dim from fields if available
        ls = list(expr.atoms((ScalarField, VectorField)))
        if ls:
            F = ls[0]
            space = F.space

        else:
            tag = random_string( 3 )
            space_name = 'space_{}'.format(tag)
            space = FunctionSpace(space_name, domain)
            # TODO vector case

        obj._ldim = domain.dim
        obj._space = space

        return obj

    @property
    def expr(self):
        return self._args[0]

    @property
    def domain(self):
        return self._args[1]

    @property
    def coordinates(self):
        return self.domain.coordinates

    @property
    def space(self):
        return self._space

    # TODO do we need it?
#    def _eval_nseries(self, x, n, logx):
#        return self.expr._eval_nseries(x, n, logx)


#==============================================================================
class LinearForm(BasicForm):
    is_linear = True

    def __new__(cls, arguments, expr):

        # ...
        integrals = expr.atoms(BasicIntegral)
        if integrals:
            for integral in integrals:
                expr = expr.subs(integral, integral._args[0])
        # ...

        # ...
        expr = expand(expr)
        # ...

        args = _sanitize_arguments(arguments, is_linear=True)
        expr = BasicIntegral(expr)
        obj = Basic.__new__(cls, args, expr)

        # ...
        domain = _get_domain(expr)
        obj._domain = domain
        # ...

        return obj

    @property
    def variables(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

    @property
    def body(self):
        if self._body is None:
            expr = self.expr
            integrals = expr.atoms(BasicIntegral)
            if integrals:
                for integral in integrals:
                    expr = expr.subs(integral, integral._args[0])

            self._body = expr

        return self._body

    @property
    def test_functions(self):
        return self.variables

    @property
    def test_spaces(self):
        return [u.space for u in self.test_functions]

    @property
    def coordinates(self):
        return self.test_spaces[0].coordinates

    @property
    def ldim(self):
        return self.test_spaces[0].ldim

    def __call__(self, *args, **kwargs):

        # ... use free variables if given and available
        expr = self._update_free_variables(**kwargs)
        # ...

        # ...
        args = Tuple(*args)

        variables = self.variables
        expr = expr.xreplace(dict(list(zip(variables, args))))
        # ...

        return expr


#==============================================================================
class BilinearForm(BasicForm):
    is_bilinear = True
    _is_symmetric = None

    def __new__(cls, arguments, expr):

        # ...
        integrals = expr.atoms(BasicIntegral)
        if integrals:
            for integral in integrals:
                expr = expr.subs(integral, integral._args[0])
        # ...

        # ...
        expr = expand(expr)
        # ...

        args = _sanitize_arguments(arguments, is_bilinear=True)
        expr = BasicIntegral(expr)
        obj = Basic.__new__(cls, args, expr)

        # ...
        domain = _get_domain(expr)
        obj._domain = domain
        # ...
        return obj

    @property
    def variables(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

    @property
    def body(self):
        if self._body is None:
            expr = self.expr
            integrals = expr.atoms(BasicIntegral)
            if integrals:
                for integral in integrals:
                    expr = expr.subs(integral, integral._args[0])

            self._body = expr

        return self._body

    @property
    def test_functions(self):
        return self.variables[1]

    @property
    def trial_functions(self):
        return self.variables[0]

    @property
    def test_spaces(self):
        return [u.space for u in self.test_functions]

    @property
    def trial_spaces(self):
        return [u.space for u in self.trial_functions]

    @property
    def coordinates(self):
        return self.test_spaces[0].coordinates

    @property
    def ldim(self):
        return self.test_spaces[0].ldim

    @property
    def is_symmetric(self):
        if self._is_symmetric is None:
            left, right = self.variables
            a1 = self(left, right)
            a2 = self(right, left)
            a1 = expand(a1)
            a2 = expand(a2)
#            print(a1)
#            print(a2)
            value = a1 == a2

            self._is_symmetric = value

        return self._is_symmetric

    def __call__(self, *args, **kwargs):

        # ... use free variables if given and available
        expr = self._update_free_variables(**kwargs)
        # ...

        # ...
        assert(len(args) == 2)
        
        new_args = []
        
        for arg in args:
        
            if is_sequence(arg):
                new_args += list(arg)
            else:
                new_args.append(arg) 
            
        args = Tuple(*new_args)

        variables = Tuple(*self.variables[0], *self.variables[1])
        expr = expr.xreplace(dict(list(zip(variables, args))))
        # ...

        return expr

#==============================================================================
class Norm(Functional):
    def __new__(cls, expr, domain, kind='l2'):
#        # ...
#        tests = expr.atoms((ScalarTestFunction, VectorTestFunction))
#        if tests:
#            msg = '> Expecting an Expression without test functions'
#            raise UnconsistentArgumentsError(msg)
#
#        if not isinstance(expr, (Expr, Matrix, ImmutableDenseMatrix)):
#            msg = '> Expecting Expr, Matrix, ImmutableDenseMatrix'
#            raise UnconsistentArgumentsError(msg)
#        # ...

        # ...
        if not(kind in ['l2', 'h1']):
            raise ValueError('> Only L2, H1 norms are available')
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

        obj = Functional.__new__(cls, expr, domain)
        obj._exponent = exponent

        return obj

    @property
    def exponent(self):
        return self._exponent


#==============================================================================
def linearize(form, fields, trials=None):
    """linearize a LinearForm around the fields."""
    # ...
    if not isinstance(form, (LinearExpr, LinearForm)):
        raise TypeError('> Expecting a LinearExpr or LinearForm')

    if not isinstance(fields, (list, tuple, Tuple)):
        fields = [fields]

    for f in fields:
        if not isinstance(f, (ScalarField, VectorField)):
            raise TypeError('{} is not ScalarField/VectorField'.format(f))

    if not(trials is None):
        if not isinstance(trials, (list, tuple, Tuple)):
            trials = [trials]

        assert( all([isinstance(i, (str, ScalarTestFunction, VectorTestFunction)) for i in trials]) )
        assert( len(fields) == len(trials) )

        newtrials = []
        for i in trials:
            if isinstance(i, (ScalarTestFunction, VectorTestFunction)):
                newtrials += [i.name]

            else:
                newtrials += [i]

        trials = newtrials
    # ...

    if isinstance(form, LinearForm):
        is_form = True
        expr = form.body

    else:
        is_form = False
        expr = form.expr

    test_functions = form.variables
    fields         = Tuple(*fields)

    # ...
    trial_functions = []
    newargs         = []
    eps  = Constant('eps_' + random_string( 4 ))
    for i,x in enumerate(fields):
        tag  = random_string( 4 )

        if trials is None:
            name = x.name + '_' + tag
        else:
            name = trials[i]

        if isinstance(x, ScalarField):
            trial  = ScalarTestFunction(x.space, name=name)

        elif isinstance(x, VectorField):
            trial  = VectorTestFunction(x.space, name=name)

        else:
            raise TypeError('Only ScalarTestFunction and VectorTestFunction are available')

        newargs         += [x + eps*trial]
        trial_functions += [trial]
    # ...

    # ...
    newexpr = expr
    for k,v in zip(fields, newargs):
        newexpr = newexpr.subs(k,v)
    # ...

    newexpr = expand(newexpr)

    e = newexpr.series(eps, 0, 2)
    d = collect(e, eps, evaluate=False)
    expr = d[eps]

#    print('> linearize = ', expr)
#    import sys; sys.exit(0)

    test_trial = (trial_functions, test_functions)

    if is_form:
        return BilinearForm(test_trial, expr)

    else:
        return BilinearExpr(test_trial, expr)
