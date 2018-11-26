# coding: utf-8

from itertools import chain

from sympy.core import Symbol
from sympy import Tuple

from pyccel.codegen.printing.pycode import PythonCodePrinter as PyccelPythonCodePrinter

from sympde.core.generic import Dot, Inner, Cross
from sympde.core.generic import Grad, Rot, Curl, Div
from sympde.core.geometry import Line, Square, Cube
from sympde.core.derivatives import _partial_derivatives
from sympde.core.derivatives import print_expression


class PythonCodePrinter(PyccelPythonCodePrinter):

    def _print_dx(self, expr):
        arg = expr.args[0]
        if isinstance(arg, _partial_derivatives):
            arg = print_expression(arg, mapping_name=False)

        else:
            arg = self._print(arg) + '_'

        return arg + 'x'

    def _print_dy(self, expr):
        arg = expr.args[0]
        if isinstance(arg, _partial_derivatives):
            arg = print_expression(arg, mapping_name=False)

        else:
            arg = self._print(arg) + '_'

        return arg + 'y'

    def _print_dz(self, expr):
        arg = expr.args[0]
        if isinstance(arg, _partial_derivatives):
            arg = print_expression(arg, mapping_name=False)

        else:
            arg = self._print(arg) + '_'

        return arg + 'z'

    def _print_IndexedTestTrial(self, expr):
        base = self._print(expr.base)
        index = self._print(expr.indices[0])
        return  '{base}_{i}'.format(base=base, i=index)

    def _print_IndexedVectorField(self, expr):
        base = self._print(expr.base)
        index = self._print(expr.indices[0])
        return  '{base}_{i}'.format(base=base, i=index)


def pycode(expr, **settings):
    """ Converts an expr to a string of Python code
    Parameters
    ==========
    expr : Expr
        A SymPy expression.
    fully_qualified_modules : bool
        Whether or not to write out full module names of functions
        (``math.sin`` vs. ``sin``). default: ``True``.
    Examples
    ========
    >>> from sympy import tan, Symbol
    >>> from sympy.printing.pycode import pycode
    >>> pycode(tan(Symbol('x')) + 1)
    'math.tan(x) + 1'
    """
    return PythonCodePrinter(settings).doprint(expr)
