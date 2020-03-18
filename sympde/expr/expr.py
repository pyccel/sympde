# coding: utf-8

from sympy import Dummy
from sympy import Matrix
from sympy.core import Basic
from sympy.core import Expr, Add, Mul
from sympy.core.expr import AtomicExpr
from sympy.core.numbers import Zero as sy_Zero
from sympy.core.containers import Tuple
from sympy.core.compatibility import is_sequence

from sympde.core.basic import CalculusFunction
from sympde.core.basic import Constant
from sympde.core.utils import random_string
from sympde.calculus import Dot, Inner
from sympde.calculus import Grad, Hessian
from sympde.topology import BasicDomain, Union
from sympde.topology import NormalVector
from sympde.topology import Boundary, Interface
from sympde.topology.space import ScalarFunctionSpace
from sympde.topology.space import ProductSpace
from sympde.topology.space import ScalarTestFunction
from sympde.topology.space import VectorTestFunction
from sympde.topology.space import Trace, trace_0, trace_1
from sympde.topology.space import ScalarField, VectorField

from .errors import UnconsistentLinearExpressionError
from .basic  import BasicForm
from .basic  import BasicExpr
from .basic  import _sanitize_arguments

#==============================================================================
def expand(expr):
    """
    Expand an expression by making sure that Mul objects are evaluated by using
    the * operator (which might be overloaded by the types in the expression).

    """
    from sympy import expand as _expand
    from operator  import mul
    from functools import reduce

    if isinstance(expr, Tuple):
        return expr

    coeff, args = _expand(expr).as_coeff_add()
    newargs = [c * reduce(mul, m) for a in args for c, m in [a.as_coeff_mul()]]

    return Add(coeff, *newargs)

#==============================================================================
def _get_domain(expr):
    # expr is an integral of BasicExpr or Add of Integral of BasicExpr
    if isinstance(expr, (DomainIntegral, BoundaryIntegral, InterfaceIntegral)):
        return expr.domain

    elif isinstance(expr, (Add,Mul)):
        domains = set()
        for a in expr.args:
            a = _get_domain(a)
            if isinstance(a, Union):
                domains = domains.union(a.args)
            elif isinstance(a, BasicDomain):
                domains = domains.union([a])
        if len(domains) == 1:
            return tuple(domains)[0]
        return Union(*domains)

#==============================================================================
class LinearExpr(BasicExpr):
    is_linear = True

    def __new__(cls, arguments, expr):

        if expr.atoms(DomainIntegral, BoundaryIntegral, InterfaceIntegral):
            raise TypeError('')

        args = _sanitize_arguments(arguments, is_linear=True)

        if not is_linear_expression(expr, args):
            msg = '> Expression is not linear'
            raise UnconsistentLinearExpressionError(msg)

        return Basic.__new__(cls, args, expr)

    @property
    def variables(self):
        return self._args[0]

    @property
    def expr(self):
        return self._args[1]

    def __call__(self, *args):
        args = _sanitize_arguments(args, is_linear=True)
        expr, _ =  self.expr._xreplace(dict(list(zip(self.variables, args))))
        return expr

    def _eval_nseries(self, x, n, logx):
        return self.expr._eval_nseries(x, n, logx)

#==============================================================================
class Integral(CalculusFunction):

    def __new__(cls, expr, domain, **options):
        # (Try to) sympify args first

        assert isinstance(domain, BasicDomain)
        return Integral.eval(expr, domain)

    @staticmethod
    def eval(expr, domain):
        """."""
        expr = expand(expr)

        if not isinstance(expr, Expr):
            raise TypeError('only Expr are accepted')

        if isinstance(expr, sy_Zero):
            return sy_Zero

        if isinstance(expr, Add):
            args = [Integral.eval(a, domain) for a in expr.args]
            return expr._new_rawargs(*args)

        if isinstance(domain, Union):
            expr = [Integral.eval(expr, d) for d in domain.args]
            return Add(*expr)

        if isinstance(domain, Boundary):
            return BoundaryIntegral(expr, domain)

        elif isinstance(domain, Interface):
            return InterfaceIntegral(expr, domain)

        else:
            return DomainIntegral(expr, domain)

#==============================================================================
class DomainIntegral(AtomicExpr):
    _op_priority = 20
    @property
    def expr(self):
        return self._args[0]
    @property
    def domain(self):
        return self._args[1]

    def __mul__(self, o):
        return DomainIntegral(self.expr*o, self.domain)

    def __rmul__(self, o):
        return DomainIntegral(self.expr*o, self.domain)

    def __eq__(self, a):
        if isinstance(a, DomainIntegral):
            eq = self.domain == a.domain
            eq = eq and self.expr == a.expr
            return eq
        return False

    def __hash__(self):
        return hash(self.expr) + hash(self.domain)

#==============================================================================
class BoundaryIntegral(AtomicExpr):
    _op_priority = 20

    def __new__(cls, expr, domain):

        atoms_1 = list(expr.atoms(Dot,Trace))

        for i in range(len(atoms_1)):
            a = atoms_1[i]
            if isinstance(a, Dot):
                if not isinstance(a.args[0], NormalVector):
                    if not isinstance(a.args[1], NormalVector):
                        atoms_1.remove(a)

        subs_1  = {a:Dummy() for a in atoms_1}
        expr, _ = expr._xreplace(subs_1)

        atoms_2 = expr.atoms(ScalarTestFunction, VectorTestFunction)
        subs_2  = {a:trace_0(a, domain) for a in atoms_2}
        expr, _ = expr._xreplace(subs_2)

        subs_3 = {}

        for key,val in subs_1.items():

            if isinstance(key, Dot):
                args = key.args
                if isinstance(args[0], NormalVector):
                    v = args[1]
                elif isinstance(args[1], NormalVector):
                    v = args[0]
                subs_3[val] = trace_1(v, domain)
            else:
                subs_3[val] = key

        expr, _ = expr._xreplace(subs_3)


        return Basic.__new__(cls, expr, domain)

    @property
    def expr(self):
        return self._args[0]
    @property
    def domain(self):
        return self._args[1]

    def __mul__(self, o):
        return BoundaryIntegral(self.expr*o, self.domain)

    def __rmul__(self, o):
        return BoundaryIntegral(self.expr*o, self.domain)

    def __eq__(self, a):
        if isinstance(a, BoundaryIntegral):
            eq = self.domain == a.domain
            eq = eq and self.expr == a.expr
            return eq
        return False

    def __hash__(self):
        return hash(self.expr) + hash(self.domain)

#==============================================================================
class InterfaceIntegral(AtomicExpr):
    _op_priority = 20

    @property
    def expr(self):
        return self._args[0]

    @property
    def domain(self):
        return self._args[1]

    def __mul__(self, o):
        return InterfaceIntegral(self.expr*o, self.domain)

    def __rmul__(self, o):
        return InterfaceIntegral(self.expr*o, self.domain)

    def __eq__(self, a):
        if isinstance(a, InterfaceIntegral):
            eq = self.domain == a.domain
            eq = eq and self.expr == a.expr
            return eq
        return False

    def __hash__(self):
        return hash(self.expr) + hash(self.domain)

#==============================================================================
class Functional(BasicForm):
    is_functional = True

    def __new__(cls, expr, domain, eval=True):

        if eval:
            expr = Integral(expr, domain)
        obj = Basic.__new__(cls, expr, domain)

        # compute dim from fields if available
        ls = tuple(expr.atoms((ScalarField, VectorField)))
        if ls:
            F = ls[0]
            space = F.space

        else:
            tag = random_string( 3 )
            space_name = 'space_{}'.format(tag)
            space = ScalarFunctionSpace(space_name, domain)
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

        # Trivial case: null expression
        if expr == 0:
            return sy_Zero

        # Check that integral expression is given
        integral_types = DomainIntegral, BoundaryIntegral, InterfaceIntegral
        if not expr.atoms(*integral_types):
            raise ValueError('Expecting integral Expression')

        # Expand integral expression and sanitize arguments
        # TODO: why do we 'sanitize' here?
        expr = expand(expr)
        args = _sanitize_arguments(arguments, is_linear=True)

        # Check linearity with respect to the given arguments
        if not is_linear_expression(expr, args):
            msg = 'Expression is not linear w.r.t [{}]'.format(args)
            raise UnconsistentLinearExpressionError(msg)

        # Create new object of type LinearForm
        obj = Basic.__new__(cls, args, expr)

        # Compute 'domain' property (scalar or tuple)
        # TODO: is this is useful?
        obj._domain = _get_domain(expr)

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
            integrals = expr.atoms(DomainIntegral)
            if integrals:
                for integral in integrals:
                    expr = expr.subs(integral, integral.expr)

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

    def __call__(self, *tests, **kwargs):

        # Use free variables if given and available
        expr = self._update_free_variables(**kwargs)

        # Make sure that 'values' is always a list
        if len(tests) == 1:
            values = tests[0]
            if not is_sequence(values):
                values = [values]
        else:
            values = tests

        # Substitute variables with given values in linear expression
        variables = self.variables
        subs      = dict(zip(variables, values))
        expr, _   = expr._xreplace(subs)

        return expr

#==============================================================================
class BilinearForm(BasicForm):
    is_bilinear = True
    _is_symmetric = None

    def __new__(cls, arguments, expr):

        # Trivial case: null expression
        if expr == 0:
            return sy_Zero

        # Check that integral expression is given
        integral_types = DomainIntegral, BoundaryIntegral, InterfaceIntegral
        if not expr.atoms(*integral_types):
            raise ValueError('Expecting integral Expression')

        # Expand integral expression and sanitize arguments
        # TODO: why do we 'sanitize' here?
        expr = expand(expr)
        args = _sanitize_arguments(arguments, is_bilinear=True)

        # Distinguish between trial and test functions
        trial_functions, test_functions = args

        # Check linearity with respect to trial functions
        if not is_linear_expression(expr, trial_functions):
            msg = ' Expression is not linear w.r.t trial functions {}'\
                    .format(trial_functions)
            raise UnconsistentLinearExpressionError(msg)

        # Check linearity with respect to test functions
        if not is_linear_expression(expr, test_functions):
            msg = ' Expression is not linear w.r.t test functions {}'\
                    .format(test_functions)
            raise UnconsistentLinearExpressionError(msg)

        # Create new object of type BilinearForm
        obj = Basic.__new__(cls, args, expr)

        # Compute 'domain' property (scalar or tuple)
        # TODO: is this is useful?
        obj._domain = _get_domain(expr)

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
            integrals = expr.atoms(DomainIntegral)
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
            a1 = expand(self(left, right))
            a2 = expand(self(right, left))
            self._is_symmetric = (a1 == a2)
        return self._is_symmetric

    def __call__(self, trials, tests, **kwargs):

        # Use free variables if given and available
        expr = self._update_free_variables(**kwargs)

        # If needed, convert positional arguments to lists
        if not is_sequence(trials): trials = [trials]
        if not is_sequence(tests ): tests  = [tests ]

        # Concatenate input values into single list
        values = [*trials, *tests]

        # Substitute variables with given values in bilinear expression
        variables = [*self.variables[0], *self.variables[1]]
        subs      = dict(zip(variables, values))
        expr, _   = expr._xreplace(subs)

        return expr

#==============================================================================
class Norm(Functional):
    is_norm = True

    def __new__(cls, expr, domain, kind='l2', eval=True):
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
        if not(kind in ['l2', 'h1', 'h2']):
            raise ValueError('> Only L2, H1, H2 norms are available')
        # ...

        # ...
        is_vector = isinstance(expr, (Matrix, Tuple, list, tuple))
        if is_vector:
            expr = Matrix(expr)
        # ...

        # ...
        exponent = None
        if kind == 'l2' and eval:
            exponent = 2

            if not is_vector:
                expr = expr*expr

            else:
                if not( expr.shape[1] == 1 ):
                    raise ValueError('Wrong expression for Matrix. must be a row')

                v = Tuple(*expr[:,0])
                expr = Dot(v, v)

        elif kind == 'h1'and eval :
            exponent = 2

            if not is_vector:
                expr = Dot(Grad(expr), Grad(expr))

            else:
                if not( expr.shape[1] == 1 ):
                    raise ValueError('Wrong expression for Matrix. must be a row')

                v = Tuple(*expr[:,0])
                expr = Inner(Grad(v), Grad(v))

        elif kind == 'h2'and eval :
            exponent = 2

            if not is_vector:
                expr = Dot(Hessian(expr), Hessian(expr))

            else:
                raise NotImplementedError('TODO')
        # ...

        obj = Functional.__new__(cls, expr, domain, eval=eval)
        obj._exponent = exponent

        return obj

    @property
    def exponent(self):
        return self._exponent

#==============================================================================
def linearize(form, fields, trials=None):
    """
    Parameters
    ----------
    form : LinearForm
        Linear form F(v; u) to be linearized. F is nonlinear in the parameter u

    fields : [Scalar|Vector]TestFunction or tuple
        Function u with respect to which F should be differentiated.
        Tuple if element of ProductSpace.

    trials : [Scalar|Vector]TestFunction or tuple
        Trial function to be used in resulting bilinear form.

    Results
    -------
    a : BilinearForm
        a(u, v) is bilinear form obtained from linearization of F(v; u).

    """
    if not isinstance(form, LinearForm):
        raise TypeError('> Expecting a LinearForm')

    # ...
    if not isinstance(fields, (list, tuple, Tuple)):
        fields = [fields]

    for field in fields:
        if not isinstance(field, (ScalarTestFunction, VectorTestFunction)):
            raise TypeError('> Expecting a TestFunction', field)

    # ...
    if trials:
        if not isinstance(trials, (list, tuple, Tuple)):
            trials = [trials]

        for trial in trials:
            if not isinstance(trial, (ScalarTestFunction, VectorTestFunction)):
                raise TypeError('> Expecting a TestFunction', trial)
    else:
        trials = [x.space.element(x.name + '_trial_' + random_string(4)) for x in fields]

    # ...

    eps        = Constant('eps_' + random_string(4))
    new_fields = [field + eps*trial for field, trial in zip(fields, trials)]

    integrals     = form.expr.args if isinstance(form.expr, Add) else [form.expr]
    new_integrals = []

    for I in integrals:

        g0 = I.expr
        g1 = expand(I.expr.subs(zip(fields, new_fields)))
        dg_du = (g1 - g0).as_leading_term(eps).as_coefficient(eps)

        if dg_du:
            new_I = integral(I.domain, dg_du)
            new_integrals.append(new_I)

    bilinear_expr = Add(*new_integrals)
    tests = form.variables

    return BilinearForm((trials, tests), bilinear_expr)

#==============================================================================
def is_linear_expression(expr, args, debug=True):
    """checks if an expression is linear with respect to the given arguments."""
    # ...
    left_args  = []
    right_args = []

    for arg in args:
        tag    = random_string( 4 )

        if isinstance(arg, ScalarTestFunction):
            left  = ScalarTestFunction(arg.space, name='l_' + tag)
            right = ScalarTestFunction(arg.space, name='r_' + tag)

        elif isinstance(arg, VectorTestFunction):
            left  = VectorTestFunction(arg.space, name='l_' + tag)
            right = VectorTestFunction(arg.space, name='r_' + tag)
        else:
            raise TypeError('argument must be a TestFunction')

        left_args  += [left]
        right_args += [right]
    # ...

    # ... check addition
    newexpr = expr
    for arg, left, right in zip(args, left_args, right_args):
        newarg  = left + right
        newexpr = newexpr.subs(arg, newarg)

    integrals = newexpr.atoms(DomainIntegral, BoundaryIntegral, InterfaceIntegral)

    for a in integrals:
        newexpr,_ = newexpr._xreplace({a:integral(a.domain, a.expr)})

    left_expr = expr
    for arg, left in zip(args, left_args):
        left_expr = left_expr.subs(arg, left)

    right_expr = expr
    for arg, right in zip(args, right_args):
        right_expr = right_expr.subs(arg, right)

    if not( expand(newexpr) == expand(left_expr + right_expr) ):
        # TODO use a warning or exception?
        if debug:
            print('Failed to assert addition property')
            print('{} != {}'.format(newexpr, left_expr + right_expr))

        return False

    # ...

    # ... check multiplication
    tag   = random_string( 4 )
    coeff = Constant('alpha_' + tag)

    newexpr = expr
    for arg, left in zip(args, left_args):
        newarg  = coeff * left
        newexpr = newexpr.subs(arg, newarg)

    left_expr = expr
    for arg, left in zip(args, left_args):
        left_expr = left_expr.subs(arg, left)

    left_expr = coeff * left_expr.subs(arg, left)

    if not( expand(newexpr) == expand(left_expr)):
        # TODO use a warning or exception?
        if debug:
            print('Failed to assert multiplication property')
            print('{} != {}'.format(newexpr, left_expr))

        return False
    # ...

    return True

#==============================================================================
def integral(domain, expr):
    """."""
    return Integral(expr, domain)
