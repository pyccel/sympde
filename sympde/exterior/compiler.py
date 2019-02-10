# coding: utf-8

# TODO - assert the space type is not Undefined

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
from sympde.core.basic import CalculusFunction
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import H1Space, HcurlSpace, HdivSpace, L2Space, UndefinedSpace
from sympde.topology import TestFunction, ScalarTestFunction, VectorTestFunction
from sympde.topology import Field, ScalarField, VectorField
from sympde.topology.space import _is_sympde_atom
from sympde.calculus import Grad, Curl, Div
from sympde.calculus import Dot, Cross
#from sympde.calculus import grad, dot, inner, cross, rot, curl, div

from .form import DifferentialForm
from .calculus import d, wedge, ip
from .calculus import delta, jp



#==============================================================================
class ExteriorCalculusExpr(CalculusFunction):

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args, **options)
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
    def eval(cls, *_args, **kwargs):
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        index = kwargs.pop('index', None)

        if isinstance(expr, Grad):
            arg = expr._args[0]
            name = arg.name

            if index is None: index = 0

            if index == 0:
                return d(DifferentialForm(name, index=index))

            elif index == 3:
                return -delta(DifferentialForm(name, index=index))

            else:
                raise ValueError('')

        if isinstance(expr, Curl):
            arg = expr._args[0]
            name = arg.name

            if index is None: index = 1

            if index == 1:
                return d(DifferentialForm(name, index=index))

            elif index == 2:
                return delta(DifferentialForm(name, index=index))

            else:
                raise ValueError('')

        if isinstance(expr, Div):
            arg = expr._args[0]
            name = arg.name

            if index is None: index = 2

            if index == 2:
                return d(DifferentialForm(name, index=index))

            elif index == 1:
                return -delta(DifferentialForm(name, index=index))

            else:
                raise ValueError('')

        if isinstance(expr, Dot):
            # TODO ORDER OF LEFT AND RIGHT DEPEND ON THE STR!!
            #      TO BE IMPROVED
            left, right = expr._args[:]

            if _is_sympde_atom(right):
                if index is None: index = 1

                right = DifferentialForm(right.name, index=index)

                if index == 1:

                    return ip(left, right)

                elif index == 2:

                    return jp(left, right)

                else:
                    raise NotImplementedError('')

            else:
                raise NotImplementedError('')

        if isinstance(expr, Cross):
            # TODO ORDER OF LEFT AND RIGHT DEPEND ON THE STR!!
            #      TO BE IMPROVED
            left, right = expr._args[:]

            if _is_sympde_atom(right):
                if index is None: index = 2

                right = DifferentialForm(right.name, index=index)

                if index == 2:

                    return ip(left, right)

                elif index == 1:

                    return -jp(left, right)

                else:
                    raise NotImplementedError('')

            else:
                raise NotImplementedError('')

        # TODO improve
        if isinstance(expr, Mul):
            if len(expr.args) == 2:

                convert = False
                if isinstance(expr.args[0], ScalarTestFunction):
                    left  = expr.args[1]
                    right = expr.args[0]
                    convert = True

                elif isinstance(expr.args[1], ScalarTestFunction):
                    left  = expr.args[0]
                    right = expr.args[1]
                    convert = True

                if convert:
                    if index is None: index = 3

                    right = DifferentialForm(right.name, index=index)

                    if index == 3:

                        return ip(left, right)

                    elif index == 0:

                        return jp(left, right)

                    else:
                        raise NotImplementedError('')

            else:
                raise NotImplementedError('')

        return cls(expr, evaluate=False)
