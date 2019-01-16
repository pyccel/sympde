# coding: utf-8

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

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

#==============================================================================
class Dot(BasicOperator):

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

        a,b = _args
        if (a == 0) or (b == 0):
            return 0

        if isinstance(a, (list, tuple, Tuple)) and isinstance(b, (list, tuple, Tuple)):
            assert( len(a) == len(b) )
            n = len(a)
            args = [a[i]*b[i] for i in range(0,n)]
            return Add(*args)

        if isinstance(a, Add):
            args = [cls.eval(i, b) for i in a.args]
            return Add(*args)

        if isinstance(b, Add):
            args = [cls.eval(a, i) for i in b.args]
            return Add(*args)

        return cls(a, b, evaluate=False)

#==============================================================================
class Cross(BasicOperator):
    pass

#==============================================================================
class Inner(BasicOperator):

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

        if isinstance(left, Mul):
            coeffs  = [i for i in left.args if isinstance(i, _coeffs_registery)]
            vectors = [i for i in left.args if not(i in coeffs)]

            a = S.One
            if coeffs:
                a = Mul(*coeffs)

            b = S.One
            if vectors:
                b = Mul(*vectors)

            return Mul(a, cls.eval(b, right))

        if isinstance(right, Mul):
            coeffs  = [i for i in right.args if isinstance(i, _coeffs_registery)]
            vectors = [i for i in right.args if not(i in coeffs)]

            a = S.One
            if coeffs:
                a = Mul(*coeffs)

            b = S.One
            if vectors:
                b = Mul(*vectors)

            return Mul(a, cls.eval(left, b))

        return cls(left, right, evaluate=False)

#==============================================================================
# TODO add it to evaluation
# Convect(F, G) = dot(F, nabla) G
class Convect(BasicOperator):
    pass

#==============================================================================
class Grad(BasicOperator):

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
