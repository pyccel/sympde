# coding: utf-8
"""
The calculus subpackage provides different operators, as generic as possible,
but knowning properties between them. For instance, the following statement will
naturally give 0

>>> from sympde.calculus import grad, curl
>>> from sympde.topology import Domain
>>> from sympde.topology import ScalarFunctionSpace
>>> from sympde.topology import ScalarFunction
>>> from sympde.topology import element_of

>>> domain = Domain('Omega', dim=2)
>>> V = ScalarFunctionSpace('V', domain)
>>> u,u1,u2 = [ScalarFunction(V, name=i) for i in ['u', 'u1', 'u2']]

>>> curl(grad(u))
0


>>> domain = Domain('Omega', dim=2)

>>> V = ScalarFunctionSpace('V', domain)
>>> W = VectorFunctionSpace('W', domain)

>>> alpha, beta, gamma = [constant(i, dtype=float) for i in ['alpha','beta','gamma']]

>>> f,g,h = [element_of(V, name=i) for i in ['f','g','h']]
>>> F,G,H = [element_of(W, i) for i in ['F','G','H']]

Generic properties
******************

scalar gradient properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

>>> assert( grad(f+g) == grad(f) + grad(g) )
>>> assert( grad(alpha*f) == alpha*grad(f) )
>>> assert( grad(alpha*f + beta*g) == alpha*grad(f) + beta*grad(g)  )

>>> assert( grad(f*g) == f*grad(g) + g*grad(f) )
>>> assert( grad(f/g) == -f*grad(g)/g**2 + grad(f)/g )

>>> assert( expand(grad(f*g*h)) == f*g*grad(h) + f*h*grad(g) + g*h*grad(f) )

vector gradient properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

>>> assert( grad(F+G) == grad(F) + grad(G) )
>>> assert( grad(alpha*F) == alpha*grad(F) )
>>> assert( grad(alpha*F + beta*G) == alpha*grad(F) + beta*grad(G)  )

>>> assert( grad(dot(F,G)) == convect(F, G) + convect(G, F) + cross(F, curl(G)) - cross(curl(F), G) )

curl properties
^^^^^^^^^^^^^^^

>>> assert( curl(f+g) == curl(f) + curl(g) )
>>> assert( curl(alpha*f) == alpha*curl(f) )
>>> assert( curl(alpha*f + beta*g) == alpha*curl(f) + beta*curl(g)  )

divergence properties
^^^^^^^^^^^^^^^^^^^^^

>>> assert( div(F+G) == div(F) + div(G) )
>>> assert( div(alpha*F) == alpha*div(F) )
>>> assert( div(alpha*F + beta*G) == alpha*div(F) + beta*div(G)  )

>>> assert( div(cross(F,G)) == -dot(F, curl(G)) + dot(G, curl(F)) )

2D specific properties
**********************

rot properties
^^^^^^^^^^^^^^

>>> assert( rot(F+G) == rot(F) + rot(G) )
>>> assert( rot(alpha*F) == alpha*rot(F) )
>>> assert( rot(alpha*F + beta*G) == alpha*rot(F) + beta*rot(G)  )

3D specific properties
**********************


"""

from operator  import mul, add
from functools import reduce

from sympy                    import expand
from sympy                    import Indexed, sympify
from sympy                    import Matrix, ImmutableDenseMatrix
from sympy                    import cacheit
from sympy.core               import Basic
from sympy.core               import Add, Mul, Pow
from sympy.core.containers    import Tuple
from sympy.core.singleton     import S
from sympy.core.decorators    import call_highest_priority
from sympy                    import Symbol, Number, NumberSymbol

from sympde.old_sympy_utilities import is_sequence
from sympde.core.basic        import CalculusFunction
from sympde.core.basic        import _coeffs_registery
from sympde.core.basic        import ScalarConstant, VectorConstant, MatrixConstant
from sympde.topology.space    import ScalarFunction, VectorFunction
from sympde.topology.domain   import NormalVector, MinusNormalVector, PlusNormalVector
from sympde.topology.datatype import H1SpaceType, HcurlSpaceType
from sympde.topology.datatype import HdivSpaceType, UndefinedSpaceType

from .errors import ArgumentTypeError

__all__ = (
    'Average',
    'BasicOperator',
    'BasicOperatorAdd',
    'Convolution',
    'Cross',
    'Curl',
    'DiffOperator',
    'Div',
    'Dot',
    'Grad',
    'Inner',
    'Jump',
    'MinusInterfaceOperator',
    'NormalDerivative',
    'Outer',
    'PlusInterfaceOperator',
    'Rot',
#
    'add_basicop',
    'has',
    'is_constant',
    'is_scalar',
    'is_zero',
#
    'Dn',
    '_diff_ops',
    '_generic_ops',
    '_is_op_test_function',
    'avg',
    'bracket',
    'conv',
    'cross',
    'curl',
    'div',
    'dot',
    'grad',
    'hessian',
    'inner',
    'jump',
    'convect',
    'laplace',
    'D',
    'minus',
    'outer',
    'plus',
    'rot',
)

#==============================================================================
@cacheit
def has(obj, types):
    if hasattr(obj, 'args'):
        return isinstance(obj, types) or any(has(i, types) for i in obj.args)
    else:
        return isinstance(obj, types)

@cacheit
def is_zero(x):
    if isinstance(x, (Matrix, ImmutableDenseMatrix)):
        return all( i==0 for i in x[:])
    else:
        return x == 0

#==============================================================================
class BasicOperator(CalculusFunction):
    """
    Basic class for calculus operators.
    """

    _op_priority   = 10.005
    _unconsistent_types = ()
    _unevaluated_types = ()

    def __hash__(self):
        return hash(self.args)

    @call_highest_priority('__radd__')
    def __add__(self, o):
        return BasicOperatorAdd(self, o)

    @call_highest_priority('__add__')
    def __radd__(self, o):
        return BasicOperatorAdd(self, o)

class DiffOperator(CalculusFunction):
    """
    Basic class for calculus operators.
    """
    is_commutative = False
    _op_priority   = 10.005
    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

    def __hash__(self):
        return hash(self.args)

class BasicOperatorAdd(Add):
    _op_priority   = 10.005
    def __new__(cls, *args, **options):

        newargs = []
        for i in args:
            if isinstance(i, BasicOperatorAdd):
                newargs += list(i.args)
            elif not i == 0:
                newargs += [i]

        obj = Add.__new__(cls, *newargs)
        return obj

    @call_highest_priority('__rmul__')
    def __mul__(self, o):
        return BasicOperatorAdd(*[a*o for a in self.args])

    @call_highest_priority('__mul__')
    def __rmul__(self, o):
        return BasicOperatorAdd(*[a*o for a in self.args])

    @call_highest_priority('__radd__')
    def __add__(self, o):
        return BasicOperatorAdd(self, o)

    @call_highest_priority('__add__')
    def __radd__(self, o):
        return BasicOperatorAdd(self, o)

#==============================================================================
def is_constant(atom):
    """ Determine whether the given atom represents a constant number.
    """
    if isinstance(atom, _coeffs_registery):
        return True

    if (isinstance(atom, Pow) and
        isinstance(atom.base, _coeffs_registery) and
        isinstance(atom.exp , _coeffs_registery)):
        return True

    return False

#==============================================================================
def is_scalar(atom):
    """ Determine whether the given atom represents a scalar quantity.
    """
    return is_constant(atom) or isinstance(atom, ScalarFunction)

#==============================================================================
def _collect_mul_expr(expr, number=False, scalar=False, vector=False, matrix=False):
    """
    """
    dtypes = {}
    dtypes['number'] = (int, float, complex, Number, NumberSymbol)
    dtypes['scalar'] = (ScalarConstant, ScalarFunction)
    dtypes['vector'] = (VectorConstant, VectorFunction)
    dtypes['matrix'] = (MatrixConstant,)

    d = {}

    # ...
    if number:
        types = dtypes['number']
        args = []
        for i in expr.args:
            if isinstance(i, types):
                args.append(i)
            elif isinstance(i, Pow):
                flag = isinstance(i.base, types)
                flag = flag and isinstance(i.exp , types)
                if flag:
                    args.append(i)
        d['number'] = args
    # ...

    # ...
    if scalar:
        types = dtypes['scalar']
        args = []
        for i in expr.args:
            if isinstance(i, types):
                args.append(i)
            elif isinstance(i, Pow):
                flag = isinstance(i.base, types) or isinstance(i.base, dtypes['number'])
                flag = flag and ( isinstance(i.exp , types) or isinstance(i.exp, dtypes['number']) )
                if flag:
                    args.append(i)
        d['scalar'] = args
    # ...

    # ...
    if vector:
        types = dtypes['vector']
        d['vector'] = [i for i in expr.args if isinstance(i, types)]
    # ...

    # ...
    if matrix:
        types = dtypes['matrix']
        d['matrix'] = [i for i in expr.args if isinstance(i, types)]
    # ...

    # ...
    values = [i for values in d.values() for i in values]
    d['other'] = [i for i in expr.args if not i in values]
    # ...

    return d

#==============================================================================
def _split_mul_args(expr, number=False, scalar=False, constant=False):
    """
    splits a Mul expression with respect to specified arguments.
    this functions assumes that expr is Mul.
    Return:
        a dictionary
    """
    if scalar or constant:
        number = True

    dtypes = {}
    dtypes['number'] = (int, float, complex, Number, NumberSymbol)
    if constant:
        dtypes['scalar'] = (ScalarConstant,)
        dtypes['vector'] = (VectorConstant,)
    else:
        dtypes['scalar'] = (ScalarConstant, ScalarFunction)
        dtypes['vector'] = (VectorConstant, VectorFunction)
    # TODO must be treated as scalar & vector
    dtypes['matrix'] = (MatrixConstant,)

    d = {}
    if number:
        types = dtypes['number']
        args = []
        for i in expr.args:
            if isinstance(i, types):
                args.append(i)
            elif isinstance(i, Pow):
                flag = isinstance(i.base, types)
                flag = flag and isinstance(i.exp , types)
                if flag:
                    args.append(i)
        d['number'] = args

    if scalar:
        types = dtypes['scalar']

        args = []
        for i in expr.args:
            if isinstance(i, types):
                args.append(i)
            elif isinstance(i, Pow):
                flag = isinstance(i.base, types) or isinstance(i.base, dtypes['number'])
                flag = flag and ( isinstance(i.exp , types) or isinstance(i.exp, dtypes['number']) )
                if flag:
                    args.append(i)
        d['scalar'] = args

    values = [i for values in d.values() for i in values]
    d['other'] = [i for i in expr.args if not i in values]
    return d

#==============================================================================
def _split_mul_expr(expr, number=False, scalar=False, constant=False):
    """
    splits a Mul expression into two Mul expression: the first expression
    contains arguments that are filtered through the provided flags, while the
    second expression contains the other arguments.
    """
    d = _split_mul_args(expr, number=number, scalar=scalar, constant=constant)

    numbers = reduce(mul, d['number'], S.One)
    scalars = reduce(mul, d['scalar'], S.One)
    others  = reduce(mul, d['other'], S.One)
    c = numbers * scalars
    arg  = others

    return c, arg

#==============================================================================
# TODO add dot(u,u) +2*dot(u,v) + dot(v,v) = dot(u+v,u+v)
# now we only have dot(u,u) + dot(u,v)+ dot(v,u) + dot(v,v) = dot(u+v,u+v)
# add dot(u,v) = dot(v,u)
class Dot(BasicOperator):
    """
    Represents a generic Dot operator, without knowledge of the dimension.

    Examples

    >>> from sympde.calculus import Dot
    >>> from sympy import Tuple
    >>> from sympy.abc import x,y

    >>> a = Tuple(x,1)
    >>> b = Tuple(1,y)
    >>> dot(a,b)
    x + y

    This operator implements the properties of addition and multiplication

    >>> from sympde.core import constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = VectorFunctionSpace('V', domain)
    >>> u,u1,u2 = [VectorFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [VectorFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = constant('alpha', dtype=float)

    >>> dot(u1+u2,v1)
    Dot(u1, v1) + Dot(u2, v1)

    >>> dot(u1,alpha*v1)
    alpha*Dot(u1, v1)
    """

    is_scalar      = True
    is_commutative = True
    is_real        = False
    is_positive    = False
    _unconsistent_types = (int, float, complex, ScalarConstant, MatrixConstant)
    _unevaluated_types  = (VectorConstant, VectorFunction)

    def __new__(cls, arg1, arg2, **options):
        _eval = True

        # If one argument is the zero vector, return 0
        if is_zero(arg1) or is_zero(arg2):
            return S.Zero

        if isinstance(arg1, cls._unconsistent_types):
            raise TypeError('arg1 must be of a vector type, given {}'.format(type(arg1)))

        if isinstance(arg2, cls._unconsistent_types):
            raise TypeError('arg2 must be of a vector type, given {}'.format(type(arg2)))

        if isinstance(arg1, cls._unevaluated_types) and  isinstance(arg2, cls._unevaluated_types):
            _eval = False

        if isinstance(arg1, Add):
            args = [cls(i, arg2, evaluate=_eval) for i in arg1.args]
            return reduce(add, args, S.Zero)

        if isinstance(arg2, Add):
            args = [cls(arg1, i, evaluate=_eval) for i in arg2.args]
            return reduce(add, args, S.Zero)

        c1 = 1
        if isinstance(arg1, Mul):
            c1, arg1 = _split_mul_expr(arg1, scalar=True)

        c2 = 1
        if isinstance(arg2, Mul):
            c2, arg2 = _split_mul_expr(arg2, scalar=True)

        c = c1*c2

        if str(arg1) > str(arg2):
            arg1,arg2 = arg2,arg1

        obj = Basic.__new__(cls, arg1, arg2)

        return c*obj

#==============================================================================
class Cross(BasicOperator):
    """
    This operator represents the cross product between two expressions,
    regardless of the dimension.
    """
    is_scalar      = False
    is_commutative = False
    _unconsistent_types = (int, float, complex, ScalarConstant, MatrixConstant)
    _unevaluated_types  = (VectorConstant, VectorFunction)

    def __new__(cls, arg1, arg2, **options):
        _eval = True

        # Operator is anti-commutative, hence cross(u, u) = 0
        if arg1 == arg2:
            return S.Zero

        # If one argument is the zero vector, return 0
        if is_zero(arg1) or is_zero(arg2):
            return S.Zero

        if isinstance(arg1, cls._unconsistent_types):
            raise TypeError('arg1 must be of a vector type, given {}'.format(type(arg1)))

        if isinstance(arg2, cls._unconsistent_types):
            raise TypeError('arg2 must be of a vector type, given {}'.format(type(arg2)))

        if isinstance(arg1, cls._unevaluated_types) and  isinstance(arg2, cls._unevaluated_types):
            _eval = False

        if isinstance(arg1, Add):
            args = [cls(i, arg2, evaluate=_eval) for i in arg1.args]
            return reduce(add, args, S.Zero)

        if isinstance(arg2, Add):
            args = [cls(arg1, i, evaluate=_eval) for i in arg2.args]
            return reduce(add, args, S.Zero)

        c1 = 1
        if isinstance(arg1, Mul):
            c1, arg1 = _split_mul_expr(arg1, scalar=True)

        c2 = 1
        if isinstance(arg2, Mul):
            c2, arg2 = _split_mul_expr(arg2, scalar=True)

        c = c1*c2

        if str(arg1) > str(arg2):
            arg1,arg2 = arg2,arg1
            c   = -c

        obj = Basic.__new__(cls, arg1, arg2)

        return c*obj

#==============================================================================
# TODO add inner(u,u) +2*inner(u,v) + inner(v,v) = inner(u+v,u+v)
# now we only have inner(u,u) + inner(u,v)+ inner(v,u) + inner(v,v) = inner(u+v,u+v)
# add inner(u,v) = inner(v,u)
class Inner(BasicOperator):
    """
    Represents a generic Frobenius inner operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication

    Examples

    >>> from sympde.calculus import inner, grad
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = VectorFunctionSpace('V', domain)
    >>> u,u1,u2 = [VectorFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [VectorFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> inner(grad(u), grad(v))
    Inner(Grad(u), Grad(v))

    >>> inner(grad(u1+u2), grad(v))
    Inner(Grad(u1), Grad(v)) + Inner(Grad(u2), Grad(v))
    """
    is_scalar      = True
    is_commutative = True
    is_real        = False
    is_positive    = False
    _unconsistent_types = (int, float, complex, ScalarConstant, MatrixConstant)
    _unevaluated_types  = (VectorConstant, VectorFunction)

    def __new__(cls, arg1, arg2, **options):
        _eval = True

        # If one argument is the zero vector, return 0
        if is_zero(arg1) or is_zero(arg2):
            return S.Zero

        if isinstance(arg1, cls._unconsistent_types):
            raise TypeError('arg1 must be of a vector type, given {}'.format(type(arg1)))

        if isinstance(arg2, cls._unconsistent_types):
            raise TypeError('arg2 must be of a vector type, given {}'.format(type(arg2)))

        if isinstance(arg1, cls._unevaluated_types) and  isinstance(arg2, cls._unevaluated_types):
            _eval = False

        if isinstance(arg1, Add):
            args = [cls(i, arg2, evaluate=_eval) for i in arg1.args]
            return reduce(add, args, S.Zero)

        if isinstance(arg2, Add):
            args = [cls(arg1, i, evaluate=_eval) for i in arg2.args]
            return reduce(add, args, S.Zero)

        c1 = 1
        if isinstance(arg1, Mul):
            c1, arg1 = _split_mul_expr(arg1, scalar=True)

        c2 = 1
        if isinstance(arg2, Mul):
            c2, arg2 = _split_mul_expr(arg2, scalar=True)

        c = c1*c2

        if str(arg1) > str(arg2):
            arg1,arg2 = arg2,arg1

        obj = Basic.__new__(cls, arg1, arg2)

        return c*obj

#==============================================================================
class Outer(BasicOperator):
    """
    Represents a generic outer operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication


    """
    is_scalar      = False
    is_commutative = False
    _unconsistent_types = (int, float, complex, ScalarConstant, MatrixConstant)
    _unevaluated_types  = (VectorConstant, VectorFunction)

    def __new__(cls, arg1, arg2, **options):
        _eval = True

        # If one argument is the zero vector, return 0
        if is_zero(arg1) or is_zero(arg2):
            return S.Zero

        if isinstance(arg1, cls._unconsistent_types):
            raise TypeError('arg1 must be of a vector type, given {}'.format(type(arg1)))

        if isinstance(arg2, cls._unconsistent_types):
            raise TypeError('arg2 must be of a vector type, given {}'.format(type(arg2)))

        if isinstance(arg1, cls._unevaluated_types) and  isinstance(arg2, cls._unevaluated_types):
            _eval = False

        if isinstance(arg1, Add):
            args = [cls(i, arg2, evaluate=_eval) for i in arg1.args]
            return reduce(add, args, S.Zero)

        if isinstance(arg2, Add):
            args = [cls(arg1, i, evaluate=_eval) for i in arg2.args]
            return reduce(add, args, S.Zero)

        c1 = 1
        if isinstance(arg1, Mul):
            c1, arg1 = _split_mul_expr(arg1, scalar=True)

        c2 = 1
        if isinstance(arg2, Mul):
            c2, arg2 = _split_mul_expr(arg2, scalar=True)

        c = c1*c2

        obj = Basic.__new__(cls, arg1, arg2)

        return c*obj

#==============================================================================
# TODO add grad(a*A) = dot(grad(a),A) + a*grad(A) where A = (a1, a2, ...)
class Grad(DiffOperator):
    """
    Represents a generic Grad operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = ScalarFunctionSpace('V', domain)
    >>> u,u1,u2 = [ScalarFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [ScalarFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = constant('alpha', dtype=float)

    >>> grad(u1+u2,v1)
    Grad(u1, v1) + Grad(u2, v1)

    >>> grad(alpha*u1)
    alpha*Grad(u1)

    >>> grad(2)
    0
    """
    is_scalar      = False
    is_commutative = False
    _unevaluated_types  = (ScalarFunction, VectorFunction)

    def __new__(cls, expr, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(sympify(expr))
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, expr, **options)
        else:
            return r

    @classmethod
    def eval(cls, expr):
        if is_constant(expr):
            return S.Zero

        if isinstance(expr, cls._unevaluated_types):
            return cls(expr, evaluate=False)

        if isinstance(expr, Add):
            args = [cls(i, evaluate=True) for i in expr.args]
            return reduce(add, args, S.Zero)

        elif isinstance(expr, Mul):
            c, expr = _split_mul_expr(expr, scalar=True, constant=True)
            if isinstance(expr, Mul):
                a1 = expr.args[0]
                a2 = expr.func(*expr.args[1:])

                d_a1  = cls(expr.args[0], evaluate=True)
                d_a2  = cls(expr.func(*expr.args[1:]), evaluate=True)

                return expand(c * a1 * d_a2 + c * d_a1 * a2)
            else:
                return c * cls(expr, evaluate=True)

        elif isinstance(expr, Pow):  # TODO: fix this for the case where e is not a number
            b = expr.base
            e = expr.exp
            a = cls(b)
            expr = expr.func(b, e-1)
            if isinstance(a, Add):
                expr = reduce(add, [e*expr*i for i in a.args])
            else:
                expr = e*a*expr
            return expr

        # ... check consistency between space type and the operator
        if _is_sympde_atom(expr):
            if not isinstance(expr.space.kind, (UndefinedSpaceType, H1SpaceType)):
                msg = '> Wrong space kind, given {}'.format(expr.space.kind)
                raise ArgumentTypeError(msg)

        return cls(expr, evaluate=False)

#==============================================================================
class Curl(DiffOperator):
    """
    Represents a generic Curl operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = VectorFunctionSpace('V', domain)
    >>> u,u1,u2 = [VectorFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [VectorFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = constant('alpha', dtype=float)

    >>> curl(u1+u2,v1)
    Curl(u1, v1) + Curl(u2, v1)

    >>> curl(alpha*u1)
    alpha*Curl(u1)
    """

    is_scalar      = False
    is_commutative = False
    _unevaluated_types  = (ScalarFunction, VectorFunction)

    def __new__(cls, expr, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(sympify(expr))
        else:
            r = None

        if r is None:
            r = Basic.__new__(cls, expr, **options)

        if r and _is_sympde_atom(expr):
            r.is_scalar      = expr.space.domain.dim == 2
            r.is_commutative = expr.space.domain.dim == 2

        return r

    @classmethod
    def eval(cls, expr):

        if is_constant(expr):
            return S.Zero

        if isinstance(expr, Grad):
            return S.Zero

        if isinstance(expr, cls._unevaluated_types):
            return cls(expr, evaluate=False)

        if isinstance(expr, Add):
            args = [cls(i, evaluate=True) for i in expr.args]
            return reduce(add, args, S.Zero)

        elif isinstance(expr, Mul):
            d = _collect_mul_expr(expr, number=True, scalar=True, vector=False, matrix=False)

            number = reduce(mul, d['number'], S.One)
            scalar = reduce(mul, d['scalar'], S.One)
            other  = reduce(mul, d['other'], S.One)

            u = number*scalar
            F = other

            d_u  = Grad(u, evaluate=True)
            d_F  = cls(F, evaluate=True)

            return expand(u * d_F + Cross(d_u, F))

        elif isinstance(expr, Pow):
            raise NotImplementedError('{}'.format(type(expr)))

        # ... check consistency between space type and the operator
        if _is_sympde_atom(expr):
            if not isinstance(expr.space.kind, (UndefinedSpaceType,
                                                HcurlSpaceType,
                                                H1SpaceType)):
                msg = '> Wrong space kind, given {}'.format(expr.space.kind)
                raise ArgumentTypeError(msg)

        return cls(expr, evaluate=False)

#==============================================================================
class Rot(DiffOperator):
    """
    Represents a generic 2D rotational operator.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = ScalarFunctionSpace('V', domain)
    >>> u,u1,u2 = [ScalarFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [ScalarFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = constant('alpha', dtype=float)

    >>> rot(u1+u2,v1)
    Rot(u1, v1) + Rot(u2, v1)

    >>> rot(alpha*u1)
    alpha*Rot(u1)
    """

    is_scalar      = False
    is_commutative = False
    _unevaluated_types  = (ScalarFunction, VectorFunction)

    def __new__(cls, expr, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(sympify(expr))
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, expr, **options)
        else:
            return r

    @classmethod
    def eval(cls, expr):

        if is_constant(expr):
            return S.Zero

        if isinstance(expr, Grad):
            return S.Zero

        if isinstance(expr, cls._unevaluated_types):
            return cls(expr, evaluate=False)

        if isinstance(expr, Add):
            args = [cls(i, evaluate=True) for i in expr.args]
            return reduce(add, args, S.Zero)

        elif isinstance(expr, Mul):
            d = _collect_mul_expr(expr, number=True, scalar=True, vector=False, matrix=False)

            number = reduce(mul, d['number'], S.One)
            scalar = reduce(mul, d['scalar'], S.One)
            other  = reduce(mul, d['other'], S.One)

            u = number*scalar
            F = other

            d_u  = Grad(u, evaluate=True)
            d_F  = cls(F, evaluate=True)

            return expand(u * d_F + Cross(d_u, F))

        elif isinstance(expr, Pow):
            raise NotImplementedError('{}'.format(type(expr)))

        # ... check consistency between space type and the operator
        # TODO add appropriate space types
        if _is_sympde_atom(expr):
            if not isinstance(expr.space.kind, (UndefinedSpaceType, H1SpaceType)):
                msg = '> Wrong space kind, given {}'.format(expr.space.kind)
                raise ArgumentTypeError(msg)
        # ...

        return cls(expr, evaluate=False)

#==============================================================================
class Div(DiffOperator):
    """
    Represents a generic Div operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = VectorFunctionSpace('V', domain)
    >>> u,u1,u2 = [VectorFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [VectorFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = constant('alpha', dtype=float)

    >>> div(u1+u2,v1)
    Div(u1, v1) + Div(u2, v1)

    >>> div(alpha*u1)
    alpha*Div(u1)
    """
    is_commutative = True
    is_scalar      = True
    _unevaluated_types  = (ScalarFunction, VectorFunction)

    def __new__(cls, expr, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(sympify(expr))
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, expr, **options)
        else:
            return r

    @classmethod
    def eval(cls, expr):

        if is_constant(expr):
            return S.Zero

        if isinstance(expr, Curl):
            return S.Zero

        if isinstance(expr, Cross):
            a,b = expr._args
            return Dot(b, Curl(a)) - Dot(a, Curl(b))

        if isinstance(expr, cls._unevaluated_types):
            return cls(expr, evaluate=False)

        if isinstance(expr, Add):
            args = [cls(i, evaluate=True) for i in expr.args]
            return reduce(add, args, S.Zero)

        elif isinstance(expr, Mul):
            d = _collect_mul_expr(expr, number=True, scalar=True, vector=False, matrix=False)

            number = reduce(mul, d['number'], S.One)
            scalar = reduce(mul, d['scalar'], S.One)
            other  = reduce(mul, d['other'], S.One)

            u = number*scalar
            F = other

            d_u  = Grad(u, evaluate=True)
            d_F  = cls(F, evaluate=True)

            return expand(u * d_F + Dot(d_u, F))

        elif isinstance(expr, Pow):
            raise NotImplementedError('{}'.format(type(expr)))

        # ... check consistency between space type and the operator
        if _is_sympde_atom(expr):
            if not isinstance(expr.space.kind, (UndefinedSpaceType,
                                                HdivSpaceType,
                                                H1SpaceType)):
                msg = '> Wrong space kind, given {}'.format(expr.space.kind)
                raise ArgumentTypeError(msg)
        # ...

        return cls(expr, evaluate=False)

#==============================================================================
class Convolution(BasicOperator):
    """
    Represents a generic Convolution operator, without knowledge of the dimension.

    """

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

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        if not len(_args) == 2:
            raise ValueError('Expecting two arguments')

        left,right = _args
        if (left == 0) or (right == 0):
            return 0

        # ...
        if isinstance(left, (list, tuple, Tuple)) and isinstance(right, (list, tuple, Tuple)):
            assert( len(left) == len(right) )
            n = len(left)
            args = [left[i]*right[i] for i in range(0,n)]
            return Add(*args)
        # ...

        # ...
        if isinstance(left, Add):
            args = [cls.eval(i, right) for i in left.args]
            return Add(*args)
        # ...

        # ...
        if isinstance(right, Add):
            args = [cls.eval(left, i) for i in right.args]
            return Add(*args)
        # ...

        # ... from now on, we construct left and right with some coeffs
        #     return is done at the end
        alpha = S.One
        if isinstance(left, Mul):
            coeffs  = [a for a in left.args if isinstance(a, _coeffs_registery)]
            for a in left.args:
                if ( isinstance(a, Pow) and
                     isinstance(a.base, _coeffs_registery) and
                     isinstance(a.exp, _coeffs_registery) ):
                    coeffs += [a]

            vectors = [a for a in left.args if not(a in coeffs)]

            a = S.One
            if coeffs:
                a = Mul(*coeffs)

            b = S.One
            if vectors:
                b = Mul(*vectors)

            alpha *= a
            left   = b

        if isinstance(right, Mul):
            coeffs  = [a for a in right.args if isinstance(a, _coeffs_registery)]
            for a in right.args:
                if ( isinstance(a, Pow) and
                     isinstance(a.base, _coeffs_registery) and
                     isinstance(a.exp, _coeffs_registery) ):
                    coeffs += [a]

            vectors = [a for a in right.args if not(a in coeffs)]

            a = S.One
            if coeffs:
                a = Mul(*coeffs)

            b = S.One
            if vectors:
                b = Mul(*vectors)

            alpha *= a
            right  = b
        # ...

        # ... this is a hack to ensure commutativity
        #     TODO to be improved
        try:
            if str(right) < str(left):
                return alpha*cls(right, left, evaluate=False)

        except:
            pass
        # ...

        return alpha*cls(left, right, evaluate=False)

#==============================================================================
class NormalDerivative(DiffOperator):
    """
    Represents the normal derivative.

    This operator implements the properties of addition and multiplication

    Examples

    """

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

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        if isinstance(expr, Add):
            args = [cls.eval(a) for a in expr.args]
            return Add(*args)

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if isinstance(a, _coeffs_registery)]
            vectors = [a for a in expr.args if not(a in coeffs)]

            a = S.One
            if coeffs:
                a = Mul(*coeffs)

            b = S.One
            if vectors:
                try:
                    if len(vectors) == 1:
                        f = vectors[0]
                        b = cls(f)

                    elif len(vectors) == 2:
                        f,g = vectors
                        b = f*cls(g) + g*cls(f)

                    else:
                        left = vectors[0]
                        right = Mul(*vectors[1:])

                        f_left  = cls(left, evaluate=True)
                        f_right = cls(right, evaluate=True)

                        b = left * f_right + f_left * right

                except:
                    b = cls(Mul(*vectors), evaluate=False)

            return Mul(a, b)

        return cls(expr, evaluate=False)

#==============================================================================
class Jump(BasicOperator):
    """
    Represents the jump of an expression at the interface of two subdomains.

    This operator implements the properties of addition and multiplication

    Examples

    """

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

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        if isinstance(expr, Add):
            args = [cls.eval(a) for a in expr.args]
            return Add(*args)

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if isinstance(a, _coeffs_registery)]
            vectors = [a for a in expr.args if not(a in coeffs)]

            a = S.One
            if coeffs:
                a = Mul(*coeffs)

            b = S.One
            if vectors:
                try:
                    if len(vectors) == 1:
                        f = vectors[0]
                        b = cls(f)

                    elif len(vectors) == 2:
                        f,g = vectors
                        b = f*cls(g) + g*cls(f)

                    else:
                        left = vectors[0]
                        right = Mul(*vectors[1:])

                        f_left  = cls(left, evaluate=True)
                        f_right = cls(right, evaluate=True)

                        b = left * f_right + f_left * right

                except:
                    b = cls(Mul(*vectors), evaluate=False)

            return Mul(a, b)

        return cls(expr, evaluate=False)

#==============================================================================
class Average(BasicOperator):
    """
    Represents the average of an expression at the interface of two subdomains.

    This operator implements the properties of addition and multiplication

    Examples

    """

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

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        if isinstance(expr, Add):
            args = [cls.eval(a) for a in expr.args]
            return Add(*args)

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if isinstance(a, _coeffs_registery)]
            vectors = [a for a in expr.args if not(a in coeffs)]

            a = S.One
            if coeffs:
                a = Mul(*coeffs)

            b = S.One
            if vectors:
                try:
                    if len(vectors) == 1:
                        f = vectors[0]
                        b = cls(f)

                    elif len(vectors) == 2:
                        f,g = vectors
                        b = f*cls(g) + g*cls(f)

                    else:
                        left = vectors[0]
                        right = Mul(*vectors[1:])

                        f_left  = cls(left, evaluate=True)
                        f_right = cls(right, evaluate=True)

                        b = left * f_right + f_left * right

                except:
                    b = cls(Mul(*vectors), evaluate=False)

            return Mul(a, b)

        return cls(expr, evaluate=False)

#==============================================================================
class MinusInterfaceOperator(BasicOperator):
    """
    The minus operator represents the value of an expression on the first side
    of an interface.
    """

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

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        if isinstance(expr, Add):
            args = [cls.eval(a) for a in expr.args]
            return Add(*args)

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if isinstance(a, _coeffs_registery)]
            vectors = [a for a in expr.args if not(a in coeffs)]

            a = S.One
            if coeffs:
                a = Mul(*coeffs)

            b = S.One
            if vectors:
                try:
                    if len(vectors) == 1:
                        f = vectors[0]
                        b = cls(f)

                    elif len(vectors) == 2:
                        f,g = vectors
                        b = f*cls(g) + g*cls(f)

                    else:
                        left = vectors[0]
                        right = Mul(*vectors[1:])

                        f_left  = cls(left, evaluate=True)
                        f_right = cls(right, evaluate=True)

                        b = left * f_right + f_left * right

                except:
                    b = cls(Mul(*vectors), evaluate=False)

            return Mul(a, b)

        elif isinstance(expr, NormalDerivative):
            n = NormalVector('n')
            u = expr._args[0]

            return Dot(Grad(cls(u)), cls(n))

        elif isinstance(expr, NormalVector):
            return MinusNormalVector('n')

        elif expr.is_zero:
            return S.Zero

        elif isinstance(expr, (Matrix, ImmutableDenseMatrix)):
            newexpr = Matrix.zeros(*expr.shape)
            for i in range(expr.shape[0]):
                for j in range(expr.shape[1]):
                    newexpr[i,j] = cls(expr[i,j])
            return type(expr)(newexpr)

        return cls(expr, evaluate=False)

    @property
    def space(self):
        return self.args[0].space

    def __getitem__(self, key):
        return MinusInterfaceOperator(self.args[0][key])

#==============================================================================
class PlusInterfaceOperator(BasicOperator):
    """
    The plus operator represents the value of an expression on the second side
    of an interface.
    """

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

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        if isinstance(expr, Add):
            args = [cls.eval(a) for a in expr.args]
            return Add(*args)

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if isinstance(a, _coeffs_registery)]
            vectors = [a for a in expr.args if not(a in coeffs)]

            a = S.One
            if coeffs:
                a = Mul(*coeffs)

            b = S.One
            if vectors:
                try:
                    if len(vectors) == 1:
                        f = vectors[0]
                        b = cls(f)

                    elif len(vectors) == 2:
                        f,g = vectors
                        b = f*cls(g) + g*cls(f)

                    else:
                        left = vectors[0]
                        right = Mul(*vectors[1:])

                        f_left  = cls(left, evaluate=True)
                        f_right = cls(right, evaluate=True)

                        b = left * f_right + f_left * right

                except:
                    b = cls(Mul(*vectors), evaluate=False)

            return Mul(a, b)

        elif isinstance(expr, NormalDerivative):
            n = NormalVector('n')
            u = expr._args[0]

            return Dot(Grad(cls(u)), cls(n))

        elif isinstance(expr, NormalVector):
            return PlusNormalVector('n')

        elif expr.is_zero:
            return S.Zero

        elif isinstance(expr, (Matrix, ImmutableDenseMatrix)):
            newexpr = Matrix.zeros(*expr.shape)
            for i in range(expr.shape[0]):
                for j in range(expr.shape[1]):
                    newexpr[i,j] = cls(expr[i,j])
            return type(expr)(newexpr)

        return cls(expr, evaluate=False)

    @property
    def space(self):
        return self.args[0].space

    def __getitem__(self, key):
        return PlusInterfaceOperator(self.args[0][key])

#==============================================================================

_generic_ops  = (Dot, Cross, Inner, Outer, Convolution)

_diff_ops  = (Grad, Curl, Rot, Div)

# ... alias for ufl compatibility
# user friendly functions
cross = Cross
dot   = Dot
inner = Inner
outer = Outer

grad = Grad
curl = Curl
rot = Rot
div = Div

convect = lambda a,u: dot(a, grad(u))
laplace = lambda _: div(grad(_))
D       = lambda w: (grad(w) + Transpose(grad(w))) / 2
hessian = None # TODO to be defined
bracket = lambda u,v: dot(grad(u), rot(v))
conv    = Convolution

jump  = Jump
avg   = Average
Dn    = NormalDerivative
minus = MinusInterfaceOperator
plus  = PlusInterfaceOperator
# ...

_is_op_test_function = lambda op: (isinstance(op, (Grad, Curl, Div)) and
                                   isinstance(op._args[0], (ScalarFunction, VectorFunction)))

_is_sympde_atom   = lambda a: isinstance(a, (ScalarFunction, VectorFunction, minus, plus))

def add_basicop(expr):
    return BasicOperatorAdd(*expr.args)

Basic._constructor_postprocessor_mapping[BasicOperator] = {
    "Add": [add_basicop],
}
