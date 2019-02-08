# coding: utf-8
"""
The calculus subpackage provides different operators, as generic as possible,
but knowning properties between them. For instance, the following statement will
naturally give 0

>>> from sympde.calculus import grad, curl
>>> from sympde.topology import Domain
>>> from sympde.topology import FunctionSpace
>>> from sympde.topology import TestFunction

>>> domain = Domain('Omega', dim=2)
>>> V = FunctionSpace('V', domain)
>>> u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]

>>> curl(grad(u))
0


>>> domain = Domain('Omega', dim=2)

>>> V = FunctionSpace('V', domain)
>>> W = VectorFunctionSpace('W', domain)

>>> alpha, beta, gamma = [Constant(i) for i in ['alpha','beta','gamma']]

>>> f,g,h = [Field(V, name=i) for i in ['f','g','h']]
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
from sympy import Indexed, IndexedBase
from sympy.core import Add, Mul, Pow
from sympy.core.containers import Tuple
from sympy.core.singleton import S

from sympde.core.basic import CalculusFunction
from sympde.core.basic import _coeffs_registery

from sympde.topology.space import TestFunction, VectorTestFunction, IndexedTestTrial
from sympde.topology.space import Field, VectorField, IndexedVectorField


#==============================================================================
class BasicOperator(CalculusFunction):
    """
    Basic class for calculus operators.
    """

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

#==============================================================================
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
        """."""

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

                elif isinstance(a, (Field, TestFunction)):
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

                elif isinstance(a, (Field, TestFunction)):
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
# TODO add properties
class Cross(BasicOperator):
    """
    This operator represents the cross product between two expressions,
    regardless of the dimension.
    """
    pass

#==============================================================================
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
        """."""

        if not _args:
            return

        if not len(_args) == 2:
            raise ValueError('Expecting two arguments')

        left,right = _args
        if (left == 0) or (right == 0):
            return 0

        if isinstance(left, Add):
            args = [cls.eval(i, right) for i in left.args]
            return Add(*args)

        if isinstance(right, Add):
            args = [cls.eval(left, i) for i in right.args]
            return Add(*args)

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

                elif isinstance(a, (Field, TestFunction)):
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

                elif isinstance(a, (Field, TestFunction)):
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
# TODO add it to evaluation
# Convect(F, G) = dot(F, nabla) G
class Convect(BasicOperator):
    """
    This operator represents the convection operator defined as
    :math:`convect(F, G) := (F \cdot \\nabla) G`.

    This operator implements the properties of addition and multiplication

    Examples

    >>> domain = Domain('Omega', dim=2)
    >>> V = FunctionSpace('V', domain)
    >>> W = VectorFunctionSpace('W', domain)
    >>> alpha, beta, gamma = [Constant(i) for i in ['alpha','beta','gamma']]
    >>> f,g,h = [Field(V, name=i) for i in ['f','g','h']]
    >>> F,G,H = [VectorField(W, i) for i in ['F','G','H']]

    >>> convect(F+G, H)
    convect(F,H) + convect(G,H)

    >>> convect(alpha*F,H)
    alpha*convect(F,H)

    >>> convect(F,alpha*H)
    alpha*convect(F,H)
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
        """."""

        if not _args:
            return

        if not len(_args) == 2:
            raise ValueError('Expecting two arguments')

        left,right = _args
        if (left == 0) or (right == 0):
            return 0

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

        return alpha*cls(left, right, evaluate=False)

#==============================================================================
class Grad(BasicOperator):
    """
    Represents a generic Grad operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import Constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = FunctionSpace('V', domain)
    >>> u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = Constant('alpha', is_real=True)

    >>> grad(u1+u2,v1)
    Grad(u1, v1) + Grad(u2, v1)

    >>> grad(alpha*u1)
    alpha*Grad(u1)

    >>> grad(2)
    0
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
        """."""

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

            return Mul(a, b)

        elif isinstance(expr, Pow):
            b = expr.base
            e = expr.exp
            return e*cls(b)*Pow(b, e-1)

        elif isinstance(expr, Dot):
            a,b = expr._args
            return Cross(a, Curl(b)) - Cross(Curl(a), b) + Convect(a,b) + Convect(b,a)

        return cls(expr, evaluate=False)

#==============================================================================
class Curl(BasicOperator):
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
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        if isinstance(expr, Add):
            args = expr.args
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
                if len(vectors) == 2:
                    a,b = vectors
                    if isinstance(a, (Tuple, VectorTestFunction, VectorField,
                                      Grad)):
                        f = b ; F = a
                        return f*Curl(F) + Cross(grad(f), F)

                    elif isinstance(b, (Tuple, VectorTestFunction, VectorField,
                                       Grad)):
                        f = a ; F = b
                        return f*Curl(F) + Cross(grad(f), F)

                b = cls(Mul(*vectors), evaluate=False)

            return Mul(a, b)

        elif isinstance(expr, Cross):
            a,b = expr._args
            return a * Div(b) - b*Div(a) + Convect(b, a) - Convect(a, b)

        elif isinstance(expr, Grad):
            return 0

        elif isinstance(expr, Curl):
            f = expr._args[0]
            return Grad(Div(f)) - Laplace(f)

        return cls(expr, evaluate=False)

#==============================================================================
class Rot(BasicOperator):
    """
    Represents a generic 2D rotational operator.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import Constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = FunctionSpace('V', domain)
    >>> u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = Constant('alpha', is_real=True)

    >>> rot(u1+u2,v1)
    Rot(u1, v1) + Rot(u2, v1)

    >>> rot(alpha*u1)
    alpha*Rot(u1)
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
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        if isinstance(expr, Add):
            args = expr.args
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
                b = cls(Mul(*vectors), evaluate=False)

            return Mul(a, b)

        return cls(expr, evaluate=False)

#==============================================================================
class Div(BasicOperator):
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
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        if isinstance(expr, Add):
            args = expr.args
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
                if len(vectors) == 2:
                    a,b = vectors
                    if isinstance(a, (Tuple, VectorTestFunction, VectorField)):
                        f = b ; F = a
                        return f*Div(F) + Dot(F, grad(f))

                    elif isinstance(b, (Tuple, VectorTestFunction, VectorField)):
                        f = a ; F = b
                        return f*Div(F) + Dot(F, grad(f))

                b = cls(Mul(*vectors), evaluate=False)

            return Mul(a, b)

        elif isinstance(expr, Cross):
            a,b = expr._args
            return Dot(b, Curl(a)) - Dot(a, Curl(b))

        elif isinstance(expr, Curl):
            return 0

        return cls(expr, evaluate=False)

#==============================================================================
class Laplace(BasicOperator):
    """
    Represents a generic Laplace operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import Constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = FunctionSpace('V', domain)
    >>> u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = Constant('alpha', is_real=True)

    >>> laplace(u1+u2,v1)
    Laplace(u1, v1) + Laplace(u2, v1)

    >>> laplace(alpha*u1)
    alpha*Laplace(u1)
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
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        if isinstance(expr, Add):
            args = expr.args
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
                if len(vectors) == 1:
                    f = vectors[0]
                    b = cls(f)

                elif len(vectors) == 2:
                    f,g = vectors
                    b = f*cls(g) + g*cls(f) + 2 * Dot(Grad(f), Grad(g))

                else:
                    b = cls(Mul(*vectors), evaluate=False)

            return Mul(a, b)

        return cls(expr, evaluate=False)

#==============================================================================
class Hessian(BasicOperator):
    """
    Represents a generic Hessian operator, without knowledge of the dimension.

    This operator implements the properties of addition and multiplication

    Examples


    >>> from sympde.core import Constant
    >>> from sympde.topology import Domain
    >>> from sympde.topology import VectorFunctionSpace
    >>> from sympde.topology import VectorTestFunction

    >>> domain = Domain('Omega', dim=2)
    >>> V = FunctionSpace('V', domain)
    >>> u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    >>> v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    >>> alpha = Constant('alpha', is_real=True)

    >>> hessian(u1+u2,v1)
    Hessian(u1, v1) + Hessian(u2, v1)

    >>> hessian(alpha*u1)
    alpha*Hessian(u1)
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
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        if isinstance(expr, Add):
            args = expr.args
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
                b = cls(Mul(*vectors), evaluate=False)

            return Mul(a, b)

        return cls(expr, evaluate=False)

#==============================================================================
# TODO add properties
class Bracket(BasicOperator):
    """
    This operator represents the Poisson bracket between two expressions.
    """
    pass


Outer = Dot # TODO add the Outer class Function

# ...
_generic_ops  = (Dot, Cross, Inner, Outer,
                 Grad, Curl, Rot, Div, Convect,
                 Bracket, Laplace, Hessian)
# ...

# ... alias for ufl compatibility
cross = Cross
dot = Dot
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
# ...
