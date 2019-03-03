# coding: utf-8

# TODO - use BasicOperator instead of LinearOperator

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Integer, Float
from sympy.core.singleton import Singleton
from sympy.core.compatibility import with_metaclass
from sympy.core import Add, Mul
from sympy.core.singleton import S

from sympde.core.basic import _coeffs_registery
from sympde.core import LinearOperator
from sympde.core.basic import CalculusFunction

from .form import DifferentialForm


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
class ExteriorDerivative(LinearOperator):

    nargs = None

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

        if not( len(_args) == 1):
            raise ValueError('Expecting one argument')

        expr = _args[0]

        if isinstance(expr, ExteriorDerivative):
            return 0

        if isinstance(expr, _coeffs_registery):
            return 0

        if isinstance(expr, DifferentialForm):
            if expr.index.index == expr.dim:
                return 0

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

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{d}({arg})'.format(d=sstr('d'), arg=sstr(self.args[0]))

#==============================================================================
class ExteriorProduct(LinearOperator):

    nargs = None

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

        if not( len(_args) == 2):
            raise ValueError('Expecting two arguments')

        left = _args[0]
        right = _args[1]
        # TODO add properties in the spirit of ExteriorDerivative

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

    @property
    def math_symbol(self):
        math_str = '/\\'
        return math_str

    def _sympystr(self, printer):
        sstr = printer.doprint
        left = self.args[0]
        right = self.args[1]
        return '{left} {math} {right}'.format( math=sstr(self.math_symbol),
                                                left=sstr(left),
                                                right=sstr(right) )

#==============================================================================
class InteriorProduct(LinearOperator):

    nargs = None

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

        if not( len(_args) == 2):
            raise ValueError('Expecting two arguments')

        left = _args[0]
        right = _args[1]
        # TODO add properties?

        return cls(left, right, evaluate=False)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return 'ip({left}, {right})'.format( left=sstr(self.args[0]),
                                             right=sstr(self.args[1]) )

#==============================================================================
class LieDerivative(LinearOperator):

    nargs = None

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

        if not( len(_args) == 2):
            raise ValueError('Expecting two arguments')

        left = _args[0]
        right = _args[1]
        # TODO add properties?

        return cls(left, right, evaluate=False)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return 'ld({left}, {right})'.format( left=sstr(self.args[0]),
                                             right=sstr(self.args[1]) )


#==============================================================================
# TODO improve and test
class PullBack(LinearOperator):

    nargs = None
    _name = None

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        name = options.pop('name', None)

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            obj = Basic.__new__(cls, *args, **options)
            obj._name = name
            return obj
        else:
            r._name = name
            return r

    @property
    def name(self):
        return self._name

    @classmethod
    def eval(cls, *_args):
        """."""

        if not _args:
            return

        if not( len(_args) == 1):
            raise ValueError('Expecting one argument')

        expr = _args[0]
        # TODO add properties?
        if isinstance(expr, ExteriorProduct):
            left = expr.args[0]
            right = expr.args[1]
            return ExteriorProduct(cls(left, evaluate=False),
                                   cls(right, evaluate=False))

        return cls(expr, evaluate=False)

#    def __call__(self, *args):
#        return PullBack(*args, name=self.name)

    def _sympystr(self, printer):
        sstr = printer.doprint
        name = self.name
        if not name:
            name = 'PullBack'
        arg = self.args[0]
        return '{name}({arg})'.format( name=sstr(name), arg=sstr(arg) )


#==============================================================================
class Adjoint(LinearOperator):

    nargs = None

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

        if not( len(_args) == 1):
            raise ValueError('Expecting one argument')

        expr = _args[0]

        return cls(expr, evaluate=False)

#==============================================================================
# TODO add properties
class AdjointExteriorDerivative(LinearOperator):

    nargs = None

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

        if not( len(_args) == 1):
            raise ValueError('Expecting one argument')

        expr = _args[0]

        if isinstance(expr, AdjointExteriorDerivative):
            return 0

        if isinstance(expr, _coeffs_registery):
            return 0

        if isinstance(expr, DifferentialForm):
            if expr.index.index == 0:
                return 0

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

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{d}({arg})'.format(d=sstr('delta'), arg=sstr(self.args[0]))

#==============================================================================
# TODO add properties
class AdjointInteriorProduct(LinearOperator):

    nargs = None

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

        if not( len(_args) == 2):
            raise ValueError('Expecting two arguments')

        left = _args[0]
        right = _args[1]
        # TODO add properties?

        return cls(left, right, evaluate=False)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return 'jp({left}, {right})'.format( left=sstr(self.args[0]),
                                             right=sstr(self.args[1]) )

#==============================================================================
# TODO add properties
class AdjointLieDerivative(LinearOperator):

    nargs = None

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

        if not( len(_args) == 2):
            raise ValueError('Expecting two arguments')

        left = _args[0]
        right = _args[1]
        # TODO add properties?

        return cls(left, right, evaluate=False)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return 'Ld({left}, {right})'.format( left=sstr(self.args[0]),
                                             right=sstr(self.args[1]) )

#==============================================================================
# TODO even/odd dim
class Hodge(LinearOperator):

    nargs = None

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

        if not( len(_args) == 1):
            raise ValueError('Expecting one argument')

        expr = _args[0]

        if isinstance(expr, Hodge):
            arg = expr._args[0]
            if isinstance(arg, DifferentialForm):
                k = arg.index.index
                n = arg.dim
                c = (-1)**(k*(n-k))
                return c*arg

        if isinstance(expr, _coeffs_registery):
            return 0

        elif isinstance(expr, Add):
            args = expr.args
            args = [cls.eval(a) for a in expr.args]
            return Add(*args)

        # TODO improve
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

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{hodge}({arg})'.format(hodge=sstr('hodge'), arg=sstr(self.args[0]))

#==============================================================================
#      user friendly names
d = ExteriorDerivative
wedge = ExteriorProduct
ip = InteriorProduct
ld = LieDerivative

delta = AdjointExteriorDerivative
jp = AdjointInteriorProduct
Ld = AdjointLieDerivative

hodge = Hodge
# ...
