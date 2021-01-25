# coding: utf-8

from sympy import Dummy
from sympy import Matrix, ImmutableDenseMatrix
from sympy.core import Basic, S
from sympy.core import Expr, Add, Mul
from sympy.core.expr import AtomicExpr
from sympy.core.numbers import Zero as sy_Zero
from sympy.core.containers import Tuple
from sympy.core.compatibility import is_sequence

from sympde.core.basic import CalculusFunction
from sympde.core.basic import Constant
from sympde.core.utils import random_string
from sympde.calculus import Dot, Inner, BasicOperator
from sympde.calculus import Grad, Hessian
from sympde.topology import BasicDomain, Union
from sympde.topology import NormalVector
from sympde.topology import Boundary, Interface, Domain, InteriorDomain
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
from sympy   import cacheit

from operator  import mul, add
from functools import reduce

class IntAdd(Add):
    _op_priority  = 20
    def __new__(cls, *args, **options):
        newargs = []
        for i in args:
            if isinstance(i, IntAdd):
                newargs += list(i.args)
            else:
                newargs += [i]
        newargs = [a for a in newargs if not a == 0]
        domains = list(set([i.domain for i in newargs]))
        groupes = [[i.expr for i in newargs if i.domain == d] for d in domains]
        args    = [Integral(reduce(add, g, S.Zero), d) for g,d in zip(groupes, domains)]
        return Add.__new__(cls, *args)

    def __mul__(self, o):
        return IntAdd(*[a*o for a in self.args])

    def __rmul__(self, o):
        return IntAdd(*[a*o for a in self.args])

    def __div__(self, o):
        return IntAdd(*[a/o for a in self.args])

    def __rdiv__(self, o):
        return IntAdd(*[a/o for a in self.args])

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def _sympystr(self, printer):
        args = [printer._print(i) for i in self.args]
        return 'IntAdd({})'.format(','.join(args))
#==============================================================================
def expand(expr):
    """
    Expand an expression by making sure that Mul objects are evaluated by using
    the * operator (which might be overloaded by the types in the expression).

    """
    from sympy import expand as _expand

    if isinstance(expr, Tuple):
        return expr
    coeff, args = _expand(expr, MatMul=True).as_coeff_add()
    newargs = [c * reduce(mul, m, 1) for a in args for c, m in [a.as_coeff_mul()]]
    expr    = coeff
    for e in newargs:
        expr += e
    return expr

#==============================================================================
def _get_domain(expr):
    # expr is an integral of BasicExpr or Add of Integral of BasicExpr
    if isinstance(expr, Integral):
        return expr.domain

    elif isinstance(expr, (Add, Mul)):
        domains = []
        for a in expr.args:
            a = _get_domain(a)
            if isinstance(a, Union):
                domains.extend(list(a.args))
            elif isinstance(a, BasicDomain):
                domains.append(a)
        return Union(*domains)

#==============================================================================
class LinearExpr(BasicExpr):
    is_linear = True

    def __new__(cls, arguments, expr, **options):

        if expr.atoms(Integral):
            raise TypeError('')

        args = _sanitize_arguments(arguments, is_linear=True)

        if not is_linear_expression(expr, args, integral=False):
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

    _op_priority          = 20
    is_domain_integral    = None
    is_boundary_integral  = None
    is_interface_integral = None

    def __new__(cls, expr, domain, **options):
        # (Try to) sympify args first
        if domain is None or expr == 0:
            return sy_Zero()

        assert isinstance(domain, BasicDomain)

        if not isinstance(expr, Expr):
            raise TypeError('only Expr are accepted')

        if isinstance(domain, Union):
            exprs = [cls(expr, d) for d in domain]
            return IntAdd(*exprs)
        elif isinstance(domain, Domain):
            interiors = domain.interior if isinstance(domain.interior, Union) else [domain.interior]
            exprs = [cls(expr, d) for d in interiors]
            return IntAdd(*exprs)
        elif isinstance(domain, Boundary):
            expr = cls.subs_boundary_expr(expr, domain)
            obj = CalculusFunction.__new__(cls, expr, domain)
            obj.is_boundary_integral = True

        elif isinstance(domain, Interface):
            obj = CalculusFunction.__new__(cls, expr, domain)
            obj.is_interface_integral = True

        elif isinstance(domain, InteriorDomain):
            obj = CalculusFunction.__new__(cls, expr, domain)
            obj.is_domain_integral = True
        else:
            raise TypeError(domain)

        return obj

    @property
    def expr(self):
        return self._args[0]

    @property
    def domain(self):
        return self._args[1]

    def __neg__(self):
        return Integral(-self.expr, self.domain)

    def __mul__(self, o):
        return Integral(self.expr*o, self.domain)

    def __rmul__(self, o):
        return Integral(self.expr*o, self.domain)

    def __div__(self, o):
        return Integral(self.expr/o, self.domain)

    def __rdiv__(self, o):
        return Integral(self.expr/o, self.domain)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __add__(self, o):
        if o == 0:
            return self
        else:
            return IntAdd(self, o)

    def __radd__(self, o):
        if o == 0:
            return self
        else:
            return IntAdd(self, o)

    def __eq__(self, a):
        if isinstance(a, Integral):
            eq = self.domain == a.domain
            eq = eq and ((self.expr == a.expr) or (self.expr - a.expr).expand() == 0)
            return eq
        return False

    def _hashable_content(self):
        return (self.expr , self.domain)

    def __hash__(self):
        return hash(self._hashable_content())

    @classmethod
    def subs_boundary_expr(cls, expr, domain):
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
        return expr

    def _sympystr(self, printer):
        domain = printer._print(self.domain)
        expr   = printer._print(self.expr)
        return 'Integral({}, {})'.format(domain, expr)

#==============================================================================
class Functional(BasicForm):
    is_functional = True

    def __new__(cls, expr, domain, evaluate=True, **options):

        # compute dim from fields if available
        ls = tuple(expr.atoms(ScalarField, VectorField, ScalarTestFunction, VectorTestFunction))
        if ls:
            F = ls[0]
            space = F.space
        else:
            space = None
            #TODO raise an error ?

        if isinstance(domain, Domain):
            domain = domain.interior

        if evaluate:
            expr = Integral(expr, domain)

        obj         = Basic.__new__(cls, expr)
        obj._ldim   = domain.dim
        obj._space  = space
        obj._domain = domain

        return obj

    @property
    def expr(self):
        return self._args[0]

    @property
    def domain(self):
        return self._domain

    @property
    def coordinates(self):
        return self.domain.coordinates

    @property
    def space(self):
        return self._space

#==============================================================================
class LinearForm(BasicForm):
    is_linear = True

    def __new__(cls, arguments, expr, **options):

        # Trivial case: null expression
        if expr == 0:
            return sy_Zero()

        # Check that integral expression is given
        if not expr.atoms(Integral):
            raise ValueError('Expecting an integral Expression')

        # TODO: why do we 'sanitize' here?
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
    def domain(self):
        return self._domain

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

    def __new__(cls, arguments, expr, **options):

        # Trivial case: null expression
        if expr == 0:
            return sy_Zero()

        # Check that integral expression is given

        if not expr.atoms(Integral):
            raise ValueError('Expecting integral Expression')

        # TODO: why do we 'sanitize' here?
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
    def domain(self):
        return self._domain

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

    def __new__(cls, expr, domain, kind='l2', evaluate=True, **options):
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
        kind = kind.lower()
        if not(kind in ['l2', 'h1', 'h2']):
            raise ValueError('> Only L2, H1, H2 norms are available')
        # ...

        # ...
        is_vector = isinstance(expr, (Matrix, ImmutableDenseMatrix, Tuple, list, tuple))
        if is_vector:
            expr = ImmutableDenseMatrix(expr)
        # ...

        # ...
        exponent = None
        if kind == 'l2' and evaluate:
            exponent = 2

            if not is_vector:
                expr = expr*expr

            else:
                if not( expr.shape[1] == 1 ):
                    raise ValueError('Wrong expression for Matrix. must be a row')

                v = Tuple(*expr[:,0])
                expr = Dot(v, v)

        elif kind == 'h1'and evaluate :
            exponent = 2

            if not is_vector:
                a    = Grad(expr)
                expr = Dot(a, a)

            else:
                if not( expr.shape[1] == 1 ):
                    raise ValueError('Wrong expression for Matrix. must be a row')

                v = Tuple(*expr[:,0])
                a = Grad(v)
                expr = Inner(a, a)

        elif kind == 'h2'and evaluate :
            exponent = 2

            if not is_vector:
                a    = Hessian(expr)
                expr = Dot(a, a)

            else:
                raise NotImplementedError('TODO')
        # ...

        obj = Functional.__new__(cls, expr, domain, evaluate=evaluate)
        obj._exponent = exponent
        obj._kind     = kind

        return obj

    @property
    def exponent(self):
        return self._exponent

    @property
    def kind(self):
        return self._kind

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
        g1 = I.expr.subs(zip(fields, new_fields)).expand()
        dg_du = ((g1-g0)/eps).series(eps, 0, 2).subs(eps, 0)

        if dg_du:
            new_I = integral(I.domain, dg_du)
            new_integrals.append(new_I)

    bilinear_expr = reduce(add, new_integrals)
    tests = form.variables

    return BilinearForm((trials, tests), bilinear_expr)

#==============================================================================
def is_linear_expression(expr, args, integral=True, debug=True):
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
    newargs = [left + right for left, right in zip(left_args, right_args)]

    newexpr    = expr.subs(zip(args, newargs))
    left_expr  = expr.subs(zip(args, left_args))
    right_expr = expr.subs(zip(args, right_args))

    a = newexpr
    b = left_expr + right_expr

    if not( (a-b).expand() == 0 ):
        # TODO use a warning or exception?
        if debug:
            print('Failed to assert addition property')
            print('{} != {}'.format(a.expand(), b.expand()))

        return False

    # ...

    # ... check multiplication
    tag   = random_string( 4 )
    coeff = Constant('alpha_' + tag)

    newexpr = expr
    for arg, left in zip(args, left_args):
        newarg  = coeff * left
        newexpr = newexpr.subs(arg, newarg)

    atoms     = list(newexpr.atoms(BasicOperator))
    subs      = [e.func(*e.args, evaluate=True) for e in atoms]
    newexpr   = newexpr.subs(zip(atoms, subs))

    left_expr = expr.subs(zip(args, left_args))
    left_expr = coeff * left_expr.subs(arg, left)

    if not( (newexpr-left_expr).expand() == 0):
        # TODO use a warning or exception?
        if debug:
            print('Failed to assert multiplication property')
            print('{} != {}'.format(newexpr, left_expr))

        return False
    # ...

    return True

#==============================================================================
def integral(domain, expr):
    return Integral(expr, domain)

def mul_integral(expr):
    args = expr.args
    return reduce(mul, args)

def mul_add(expr):
    args_int = [a for a in expr.args if isinstance(a, IntAdd)]
    coeff    = Mul(*[a for a in expr.args if not a in args_int])
    args_int[0] = [IntAdd(*[coeff*i for i in args_int[0].args])]
    return Mul(*args_int)

def add_int(expr):
    return IntAdd(*expr.args)

Basic._constructor_postprocessor_mapping[Integral] = {
    "Add": [add_int],
    "Mul": [mul_integral],
}

Basic._constructor_postprocessor_mapping[IntAdd] = {
    "Mul": [mul_add],
}

