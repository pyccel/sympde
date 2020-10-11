# coding: utf-8
"""
The calculus subpackage provides different operators, as generic as possible,
but knowning properties between them. For instance, the following statement will
naturally give 0

>>> from sympde.calculus import grad, curl
>>> from sympde.topology import Domain
>>> from sympde.topology import ScalarFunctionSpace
>>> from sympde.topology import ScalarTestFunction

>>> domain = Domain('Omega', dim=2)
>>> V = ScalarFunctionSpace('V', domain)
>>> u,u1,u2 = [ScalarTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]

>>> curl(grad(u))
0


>>> domain = Domain('Omega', dim=2)

>>> V = ScalarFunctionSpace('V', domain)
>>> W = VectorFunctionSpace('W', domain)

>>> alpha, beta, gamma = [Constant(i) for i in ['alpha','beta','gamma']]

>>> f,g,h = [ScalarField(V, name=i) for i in ['f','g','h']]
>>> F,G,H = [VectorField(W, i) for i in ['F','G','H']]

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

laplace properties
^^^^^^^^^^^^^^^^^^

>>> assert( laplace(f+g) == laplace(f) + laplace(g) )
>>> assert( laplace(alpha*f) == alpha*laplace(f) )
>>> assert( laplace(alpha*f + beta*g) == alpha*laplace(f) + beta*laplace(g)  )

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

from sympy.core.compatibility import is_sequence
from sympy.core import Basic
from sympy      import Indexed, IndexedBase, sympify
from sympy      import Matrix, ImmutableDenseMatrix
from sympy.core import Add, Mul, Pow, Symbol
from sympy.core.containers import Tuple
from sympy.core.singleton  import S


from sympde.core.basic import CalculusFunction
from sympde.core.basic import _coeffs_registery

from sympde.topology.space import ScalarTestFunction, VectorTestFunction, IndexedTestTrial
from sympde.topology.space import ScalarField, VectorField, IndexedVectorField
from sympde.topology.space import _is_sympde_atom
from sympde.topology.domain import NormalVector, MinusNormalVector, PlusNormalVector
from sympde.topology.datatype import H1SpaceType, HcurlSpaceType
from sympde.topology.datatype import HdivSpaceType, L2SpaceType, UndefinedSpaceType

from .matrices import MatrixSymbolicExpr
from .errors import ArgumentTypeError
from sympy.core.decorators import call_highest_priority
from operator  import mul,add
from functools import reduce
from sympy import cacheit

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
    return is_constant(atom) or isinstance(atom, (ScalarField, ScalarTestFunction))

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

    >>> from sympde.core import Constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = VectorFunctionSpace('V', domain)
    >>> u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = Constant('alpha', is_real=True)

    >>> dot(u1+u2,v1)
    Dot(u1, v1) + Dot(u2, v1)

    >>> dot(u1,alpha*v1)
    alpha*Dot(u1, v1)
    """

    is_scalar      = True
    is_commutative = True
    is_real        = False
    is_positive    = False
    def __new__(cls, arg1, arg2, **options):
        # (Try to) sympify args first

        # If one argument is the zero vector, return 0
        if is_zero(arg1) or is_zero(arg2):
            return S.Zero

        types = (VectorTestFunction, ScalarTestFunction)
        if isinstance(arg1, Add):
            a = [i for i in arg1.args if has(i, types)]
            b = [i for i in arg1.args if i not in a]

            a = [cls(i, arg2) for i in a]
            b = cls(arg1.func(*b), arg2)
            a = reduce(add, a, S.Zero)
            return a+b

        if isinstance(arg2, Add):
            a = [i for i in arg2.args if has(i, types)]
            b = [i for i in arg2.args if i not in a]

            a = [cls(arg1, i) for i in a]
            b = cls(arg1, arg2.func(*b))
            a = reduce(add, a, S.Zero)
            return a+b

        if isinstance(arg1, Mul):
            a = arg1.args
        else:
            a = [arg1]

        if isinstance(arg2, Mul):
            b = arg2.args
        else:
            b = [arg2]

        args_1 = [i for i in a if not i.is_commutative]
        c1     = [i for i in a if not i in args_1]
        args_2 = [i for i in b if not i.is_commutative]
        c2     = [i for i in b if not i in args_2]

        a = reduce(mul, args_1)
        b = reduce(mul, args_2)
        c = Mul(*c1)*Mul(*c2)

        if str(a) > str(b):
            a,b = b,a

        obj = Basic.__new__(cls, a, b)

        if a == b:
            obj.is_real     = True
            obj.is_positive = True

        return c*obj


#==============================================================================
#TODO add cross(u,v) + cross(v,u) = 0
class Cross(BasicOperator):
    """
    This operator represents the cross product between two expressions,
    regardless of the dimension.
    """
    is_scalar      = False
    is_commutative = False

    def __new__(cls, arg1, arg2, **options):

        # Operator is anti-commutative, hence cross(u, u) = 0
        if arg1 == arg2:
            return S.Zero

        # If one argument is the zero vector, return 0
        if is_zero(arg1) or is_zero(arg2):
            return S.Zero

        types = (VectorTestFunction, ScalarTestFunction)
        if isinstance(arg1, Add):
            a = [i for i in arg1.args if has(i, types)]
            b = [i for i in arg1.args if i not in a]

            a = [cls(i, arg2) for i in a]
            b = cls(arg1.func(*b), arg2)
            a = reduce(add, a, S.Zero)
            return a+b

        if isinstance(arg2, Add):
            a = [i for i in arg2.args if has(i, types)]
            b = [i for i in arg2.args if i not in a]

            a = [cls(arg1, i) for i in a]
            b = cls(arg1, arg2.func(*b))
            a = reduce(add, a, S.Zero)
            return a+b

        if isinstance(arg1, Mul):
            a = arg1.args
        else:
            a = [arg1]

        if isinstance(arg2, Mul):
            b = arg2.args
        else:
            b = [arg2]

        args_1 = [i for i in a if not i.is_commutative]
        c1     = [i for i in a if not i in args_1]
        args_2 = [i for i in b if not i.is_commutative]
        c2     = [i for i in b if not i in args_2]

        a = reduce(mul, args_1)
        b = reduce(mul, args_2)
        c = Mul(*c1)*Mul(*c2)

        if str(a) > str(b):
            a,b = b,a
            c   = -c

        obj = Basic.__new__(cls, a, b)

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
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = VectorFunctionSpace('V', domain)
    >>> u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> inner(grad(u), grad(v))
    Inner(Grad(u), Grad(v))

    >>> inner(grad(u1+u2), grad(v))
    Inner(Grad(u1), Grad(v)) + Inner(Grad(u2), Grad(v))
    """
    is_scalar      = True
    is_commutative = True
    is_real        = False
    is_positive    = False
    def __new__(cls, arg1, arg2, **options):
        # (Try to) sympify args first

        # If one argument is the zero vector, return 0
        if is_zero(arg1) or is_zero(arg2):
            return S.Zero

        types = (VectorTestFunction, ScalarTestFunction)
        if isinstance(arg1, Add):
            a = [i for i in arg1.args if has(i, types)]
            b = [i for i in arg1.args if i not in a]

            a = [cls(i, arg2) for i in a]
            b = cls(arg1.func(*b), arg2)
            a = reduce(add, a, S.Zero)
            return a+b

        if isinstance(arg2, Add):
            a = [i for i in arg2.args if has(i, types)]
            b = [i for i in arg2.args if i not in a]

            a = [cls(arg1, i) for i in a]
            b = cls(arg1, arg2.func(*b))
            a = reduce(add, a, S.Zero)
            return a+b

        if isinstance(arg1, Mul):
            a = arg1.args
        else:
            a = [arg1]

        if isinstance(arg2, Mul):
            b = arg2.args
        else:
            b = [arg2]

        args_1 = [i for i in a if not i.is_commutative]
        c1     = [i for i in a if not i in args_1]
        args_2 = [i for i in b if not i.is_commutative]
        c2     = [i for i in b if not i in args_2]

        a = reduce(mul, args_1)
        b = reduce(mul, args_2)
        c = Mul(*c1)*Mul(*c2)

        if str(a) > str(b):
            a,b = b,a

        obj = Basic.__new__(cls, a, b)

        if a == b:
            obj.is_real     = True
            obj.is_positive = True

        return c*obj

#==============================================================================
class Outer(BasicOperator):
    """
    Represents a generic outer operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication


    """

    is_scalar      = False
    is_commutative = False

    def __new__(cls, arg1, arg2, **options):
        # (Try to) sympify args first
        # If one argument is the zero vector, return 0
        if is_zero(arg1) or is_zero(arg2):
            return S.Zero

        types = (VectorTestFunction, ScalarTestFunction)
        if isinstance(arg1, Add):
            a = [i for i in arg1.args if has(i, types)]
            b = [i for i in arg1.args if i not in a]

            a = [cls(i, arg2) for i in a]
            b = cls(arg1.func(*b), arg2)
            a = reduce(add, a, S.Zero)
            return a+b

        if isinstance(arg2, Add):
            a = [i for i in arg2.args if has(i, types)]
            b = [i for i in arg2.args if i not in a]

            a = [cls(arg1, i) for i in a]
            b = cls(arg1, arg2.func(*b))
            a = reduce(add, a, S.Zero)
            return a+b

        if isinstance(arg1, Mul):
            a = arg1.args
        else:
            a = [arg1]

        if isinstance(arg2, Mul):
            b = arg2.args
        else:
            b = [arg2]

        args_1 = [i for i in a if not i.is_commutative]
        c1     = [i for i in a if not i in args_1]
        args_2 = [i for i in b if not i.is_commutative]
        c2     = [i for i in b if not i in args_2]

        a = reduce(mul, args_1)
        b = reduce(mul, args_2)
        c = Mul(*c1)*Mul(*c2)

        obj = Basic.__new__(cls, a, b)

        return c*obj

#==============================================================================
# TODO add it to evaluation
# Convect(F, G) = dot(F, nabla) G
class Convect(BasicOperator):
    r"""
    This operator represents the convection operator defined as
    :math:`convect(F, G) := (F \cdot \\nabla) G`.

    This operator implements the properties of addition and multiplication

    Examples

    >>> domain = Domain('Omega', dim=2)
    >>> V = ScalarFunctionSpace('V', domain)
    >>> W = VectorFunctionSpace('W', domain)
    >>> alpha, beta, gamma = [Constant(i) for i in ['alpha','beta','gamma']]
    >>> f,g,h = [ScalarField(V, name=i) for i in ['f','g','h']]
    >>> F,G,H = [VectorField(W, i) for i in ['F','G','H']]

    >>> convect(F+G, H)
    convect(F,H) + convect(G,H)

    >>> convect(alpha*F,H)
    alpha*convect(F,H)

    >>> convect(F,alpha*H)
    alpha*convect(F,H)
    """
    is_scalar      = False
    is_commutative = False

    def __new__(cls, arg1, arg2, **options):
        # (Try to) sympify args first

        arg1 , arg2 = sympify(arg1), sympify(arg2)
        # If one argument is the zero vector, return 0
        if is_zero(arg1) or arg2.is_number:
            return S.Zero

        types = (VectorTestFunction, ScalarTestFunction)
        if isinstance(arg1, Add):
            a = [i for i in arg1.args if has(i, types)]
            b = [i for i in arg1.args if i not in a]

            a = [cls(i, arg2) for i in a]
            b = cls(arg1.func(*b), arg2)
            a = reduce(add, a, S.Zero)
            return a+b

        if isinstance(arg2, Add):
            a = [i for i in arg2.args if has(i, types)]
            b = [i for i in arg2.args if i not in a]

            a = [cls(arg1, i) for i in a]
            b = cls(arg1, arg2.func(*b))
            a = reduce(add, a, S.Zero)
            return a+b

        if isinstance(arg1, Mul):
            a = arg1.args
        else:
            a = [arg1]

        if isinstance(arg2, Mul):
            b = arg2.args
        else:
            b = [arg2]

        args_1 = [i for i in a if not i.is_commutative]
        c1     = [i for i in a if not i in args_1]
        args_2 = [i for i in b if not i.is_commutative]
        c2     = [i for i in b if not i in args_2]

        a = reduce(mul, args_1)
        b = reduce(mul, args_2)
        c = Mul(*c1)*Mul(*c2)

        obj = Basic.__new__(cls, a, b)

        return c*obj

#==============================================================================
# TODO add grad(a*A) = dot(grad(a),A) + a*grad(A) where A = (a1, a2, ...)
class Grad(DiffOperator):
    """
    Represents a generic Grad operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import Constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = ScalarFunctionSpace('V', domain)
    >>> u,u1,u2 = [ScalarTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [ScalarTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = Constant('alpha', is_real=True)

    >>> grad(u1+u2,v1)
    Grad(u1, v1) + Grad(u2, v1)

    >>> grad(alpha*u1)
    alpha*Grad(u1)

    >>> grad(2)
    0
    """
    is_scalar      = False
    is_commutative = False
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
        types = (VectorTestFunction, ScalarTestFunction)
        if not has(expr, types):
            if expr.is_number:
                return S.Zero
            return cls(expr, evaluate=False)

        if isinstance(expr, Add):
            a = [i for i in expr.args if has(i, types)]
            b = [i for i in expr.args if i not in a]
            a = [cls(i) for i in a]
            b = cls(expr.func(*b))
            return reduce(add, a) + b

        elif isinstance(expr, Mul):

            commutative           = [a for a in expr.args if a.is_commutative]
            non_commutative       = [a for a in expr.args if not a in commutative]

            coeffs                = [a for a in commutative if a.is_number]
            commutative_free_expr = [a for a in commutative if not a in coeffs and not has(a, types)]
            commutative           = [a for a in commutative if not a in commutative_free_expr+coeffs]


            a  = reduce(mul, coeffs, S.One)
            b1 = reduce(mul, commutative_free_expr, S.One)
            b2 = reduce(mul, commutative, S.One)
            b3 = reduce(mul, non_commutative, S.One)

            if not b3 == 1:
                b = cls( b1*b2*b3, evaluate=False)
                return a*b
            elif not b1 == 1:

                d_b1  = cls(b1, evaluate=False)
                d_b2  = cls(b2, evaluate=True)

                return a * b1 * d_b2 + a * d_b1 * b2
            elif not  b2 == 1:
                if not isinstance(b2, Mul):
                    return a * cls(b2, evaluate=True)
                else:
                    arg1 = b2.args[0]
                    arg2 = b2.func(*b2.args[1:])

                    d_arg1  = cls(b2.args[0], evaluate=True)
                    d_arg2  = cls(b2.func(*b2.args[1:]), evaluate=True)

                return a * arg1 * d_arg2 + a * d_arg1 * arg2
            else:
                return S.Zero

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


    >>> from sympde.core import Constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = VectorFunctionSpace('V', domain)
    >>> u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = Constant('alpha', is_real=True)

    >>> curl(u1+u2,v1)
    Curl(u1, v1) + Curl(u2, v1)

    >>> curl(alpha*u1)
    alpha*Curl(u1)
    """

    is_scalar      = False
    is_commutative = False

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

        types = (VectorTestFunction, ScalarTestFunction)
        if not has(expr, types):
            if expr.is_number:
                return S.Zero
            return cls(expr, evaluate=False)

        if isinstance(expr, Add):
            a = [i for i in expr.args if has(i, types)]
            b = [i for i in expr.args if i not in a]
            a = [cls(i) for i in a]
            b = cls(expr.func(*b))
            return reduce(add, a) + b

        elif isinstance(expr, Mul):

            coeffs  = [a for a in expr.args if a.is_number]
            vectors = [a for a in expr.args if not(a in coeffs)]

            a = Mul(*coeffs)
            b = cls(expr.func(*vectors), evaluate=False)

            return a*b

        elif isinstance(expr, Grad):
            return S.Zero

        # ... check consistency between space type and the operator
        if _is_sympde_atom(expr):
            if not isinstance(expr.space.kind, (UndefinedSpaceType,
                                                HcurlSpaceType,
                                                H1SpaceType)):
                msg = '> Wrong space kind, given {}'.format(expr.space.kind)
                raise ArgumentTypeError(msg)
        # ...
        return cls(expr, evaluate=False)

#==============================================================================
class Rot(DiffOperator):
    """
    Represents a generic 2D rotational operator.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import Constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = ScalarFunctionSpace('V', domain)
    >>> u,u1,u2 = [ScalarTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [ScalarTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = Constant('alpha', is_real=True)

    >>> rot(u1+u2,v1)
    Rot(u1, v1) + Rot(u2, v1)

    >>> rot(alpha*u1)
    alpha*Rot(u1)
    """

    is_scalar      = False
    is_commutative = False

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

        types = (VectorTestFunction, ScalarTestFunction)
        if not has(expr, types):
            if expr.is_number:
                return S.Zero
            return cls(expr, evaluate=False)

        if isinstance(expr, Add):
            a = [i for i in expr.args if has(i, types)]
            b = [i for i in expr.args if i not in a]
            a = [cls(i) for i in a]
            b = cls(expr.func(*b))
            return reduce(add, a) + b

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if a.is_number]
            vectors = [a for a in expr.args if not(a in coeffs)]

            a = Mul(*coeffs)
            b = cls(expr.func(*vectors), evaluate=False)

            return a*b

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


    >>> from sympde.core import Constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = VectorFunctionSpace('V', domain)
    >>> u,u1,u2 = [VectorTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [VectorTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = Constant('alpha', is_real=True)

    >>> div(u1+u2,v1)
    Div(u1, v1) + Div(u2, v1)

    >>> div(alpha*u1)
    alpha*Div(u1)
    """
    is_commutative = True
    is_scalar      = True

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


        types = (VectorTestFunction, ScalarTestFunction)
        if not has(expr, types):
            if expr.is_number:
                return S.Zero
            return cls(expr, evaluate=False)

        if isinstance(expr, Add):
            a = [i for i in expr.args if has(i, types)]
            b = [i for i in expr.args if i not in a]
            a = [cls(i) for i in a]
            b = cls(expr.func(*b))
            return reduce(add, a) + b

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if a.is_number]
            vectors = [a for a in expr.args if not(a in coeffs)]

            a = Mul(*coeffs)

            b = S.One
            if vectors:
                if len(vectors) == 2:
                    a,b = vectors
                    # TODO remove try/except using regularity from space
                    try:
                        if isinstance(a, (Tuple, VectorTestFunction, VectorField)):
                            f = b ; F = a
                            return f*Div(F) + Dot(F, grad(f))

                        elif isinstance(b, (Tuple, VectorTestFunction, VectorField)):
                            f = a ; F = b
                            return f*Div(F) + Dot(F, grad(f))

                    except:
                        return cls(expr.func(*vectors), evaluate=False)

                b = cls(expr.func(*vectors), evaluate=False)

            return a*b

        elif isinstance(expr, Cross):
            a,b = expr._args
            return Dot(b, Curl(a)) - Dot(a, Curl(b))

        elif isinstance(expr, Curl):
            return S.Zero

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
class Laplace(DiffOperator):
    """
    Represents a generic Laplace operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import Constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = ScalarFunctionSpace('V', domain)
    >>> u,u1,u2 = [ScalarTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [ScalarTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = Constant('alpha', is_real=True)

    >>> laplace(u1+u2,v1)
    Laplace(u1, v1) + Laplace(u2, v1)

    >>> laplace(alpha*u1)
    alpha*Laplace(u1)
    """

    is_scalar      = True
    is_commutative = True
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

        types = (VectorTestFunction, ScalarTestFunction)
        if not has(expr, types):
            if expr.is_number:
                return S.Zero
            return cls(expr, evaluate=False)

        if isinstance(expr, Add):
            a = [i for i in expr.args if has(i, types)]
            b = [i for i in expr.args if i not in a]
            a = [cls(i) for i in a]
            b = cls(expr.func(*b))
            return reduce(add, a) + b

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if a.is_number]
            vectors = [a for a in expr.args if not(a in coeffs)]

            a = Mul(*coeffs)

            if len(vectors) == 2:
                f,g = vectors
                b = f*cls(g) + g*cls(f) + 2 * Dot(Grad(f), Grad(g))

            else:
                b = cls(Mul(*vectors), evaluate=False)

            return a*b

        # ... check consistency between space type and the operator
        # TODO add appropriate space types
        if _is_sympde_atom(expr):
            if not isinstance(expr.space.kind, UndefinedSpaceType):
                msg = '> Wrong space kind, given {}'.format(expr.space.kind)
                raise ArgumentTypeError(msg)
        # ...

        return cls(expr, evaluate=False)

#==============================================================================
class Hessian(DiffOperator):
    """
    Represents a generic Hessian operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import Constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = ScalarFunctionSpace('V', domain)
    >>> u,u1,u2 = [ScalarTestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [ScalarTestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = Constant('alpha', is_real=True)

    >>> hessian(u1+u2,v1)
    Hessian(u1, v1) + Hessian(u2, v1)

    >>> hessian(alpha*u1)
    alpha*Hessian(u1)
    """

    is_scalar      = False
    is_commutative = False

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

        types = (VectorTestFunction, ScalarTestFunction)
        if not has(expr, types):
            if expr.is_number:
                return S.Zero
            return cls(expr, evaluate=False)

        if isinstance(expr, Add):
            a = [i for i in expr.args if has(i, types)]
            b = [i for i in expr.args if i not in a]
            a = [cls(i) for i in a]
            b = cls(expr.func(*b))
            return reduce(add, a) + b

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if a.is_number]
            vectors = [a for a in expr.args if not(a in coeffs)]

            a = Mul(*coeffs)
            b = cls(Mul(*vectors), evaluate=False)

            return a*b

        # ... check consistency between space type and the operator
        # TODO add appropriate space types
        if _is_sympde_atom(expr):
            if not isinstance(expr.space.kind, UndefinedSpaceType):
                msg = '> Wrong space kind, given {}'.format(expr.space.kind)
                raise ArgumentTypeError(msg)
        # ...

        return cls(expr, evaluate=False)

#==============================================================================
class Bracket(DiffOperator):
    """
    This operator represents the Poisson bracket between two expressions.
    """

    is_scalar      = True
    is_commutative = True

    def __new__(cls, arg1, arg2, **options):


        if options.pop('evaluate', True):
            return cls.eval(sympify(arg1), sympify(arg2))
        else:
            return Basic.__new__(cls, arg1, arg2, **options)

    @classmethod
    def eval(cls, arg1, arg2):

        # Operator is anti-commutative, hence [u, u] = 0

        if arg1.is_number or arg2.is_number:
            return S.Zero

        if arg1 == arg2:
            return S.Zero

        # Recursive application of differentiation rules to both arguments
        for expr, args in (arg1, lambda a: (a, arg2)), \
                          (arg2, lambda a: (arg1, a)):

            # Derivative of sum: d(a+b+c) = da + db + dc
            if isinstance(expr, Add):
                return Add(*[cls.eval(*args(a)) for a in expr.args])

            # Derivative of product: d(a*b*c) = da * (b*c) + a * d(b*c)
            elif isinstance(expr, Mul):
                coeffs = [a for a in expr.args if isinstance(a, _coeffs_registery)]
                fields = [a for a in expr.args if a not in coeffs]

                terms  = []
                for i in range(len(fields)):
                    factors = [(cls.eval(*args(f)) if i == j else f)
                               for j, f in enumerate(fields)]
                    terms.append(Mul(*factors))

                a = Mul(*coeffs)
                b = Add(*terms )

                return Mul(a, b)

            # Derivative of constant: d(const) = 0
            elif isinstance(expr, _coeffs_registery):
                return S.Zero

        # Automatic evaluation to canonical form: reorder arguments by using
        # anti-commutativity property [v, u] = -[u, v] and stop recursion.
        if str(arg1) > str(arg2):
            return -cls(arg2, arg1, evaluate=False)

        # Stop recursion
        return cls(arg1, arg2, evaluate=False)

#==============================================================================
# TODO improve
class StrainTensor(Grad):
    """
    This operator represents the strain tensor (grad(u) + transpose(grad(u)))/2
    """

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{D}({arg})'.format(D=sstr('D'), arg=sstr(self.args[0]))

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

        return cls(expr, evaluate=False)

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

        return cls(expr, evaluate=False)

#==============================================================================

_generic_ops  = (Dot, Cross, Inner, Outer, Convect, Convolution)

_diff_ops  = (Grad, Curl, Rot, Div,
              Bracket, Laplace, Hessian)

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

convect = Convect
bracket = Bracket
laplace = Laplace
hessian = Hessian
D       = StrainTensor
conv    = Convolution

jump  = Jump
avg   = Average
Dn    = NormalDerivative
minus = MinusInterfaceOperator
plus  = PlusInterfaceOperator
# ...

_is_op_test_function = lambda op: (isinstance(op, (Grad, Curl, Div)) and
                                   isinstance(op._args[0], (ScalarTestFunction, VectorTestFunction)))

_is_op_field         = lambda op: (isinstance(op, (Grad, Curl, Div)) and
                                   isinstance(op._args[0], (ScalarField, VectorField)))

def add_basicop(expr):
    return BasicOperatorAdd(*expr.args)

Basic._constructor_postprocessor_mapping[BasicOperator] = {
    "Add": [add_basicop],
}

