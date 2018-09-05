# -*- coding: utf-8 -*-
#

from sympy.core import Symbol
from sympy import Mul
from sympy.printing.latex import LatexPrinter as LatexPrinterSympy

from sympde.core.expr import BilinearForm, LinearForm, FunctionForm
from sympde.core.generic import Dot, Inner, Cross
from sympde.core.generic import Grad, Rot, Curl, Div
from sympde.core.geometry import Line, Square, Cube
from sympde.core.derivatives import sort_partial_derivatives
from sympde.core.derivatives import get_index_derivatives
from sympde.core.derivatives import get_atom_derivatives

def print_integral(domain, expr, measure=None):
    measure_str = ''
    if not(measure is None):
        measure_str = latex(measure)

    expr_str = latex(expr)

    if isinstance(domain, Line):
        xmin, xmax = domain.bounds

        xmin = latex(xmin)
        xmax = latex(xmax)

        int_str = r'\int_{' + xmin + '}^{' + xmax + '}'

    elif isinstance(domain, (Square, Cube)):
        xmins, xmaxs = domain.bounds

        int_str = ''
        for xmin, xmax in zip(xmins, xmaxs):
            xmin = latex(xmin)
            xmax = latex(xmax)

            int_str += r'\int_{' + xmin + '}^{' + xmax + '}'

    return '{integral} {expr} {measure}'.format(integral=int_str,
                                                expr=expr_str,
                                                measure=measure_str)

class LatexPrinter(LatexPrinterSympy):

    # ... differential operators
    def _print_Grad(self, expr):
        return r'\nabla{' + self._print(expr.args[0]) + '}'

    def _print_Curl(self, expr):
        return r'\nabla \times ' + self._print(expr.args[0])

    def _print_Div(self, expr):
        return r'\nabla \cdot ' + self._print(expr.args[0])

    def _print_Rot(self, expr):
        return r'\mathrm{rot} ' + self._print(expr.args[0])

    def _print_DifferentialOperator(self, expr):
        coord = self._print(expr.coordinate)
        arg = self._print(expr.args[0])
        return r'\partial_{' + coord + '}' +  arg
    # ...

    # ...
    def _print_VectorTestFunction(self, expr):
        return r'\mathbf{' + self._print(expr.name) + '}'
    # ...

    # ... algebraic operators
    def _print_Dot(self, expr):
        left = expr.args[0]
        left = self._print(left)

        right = expr.args[1]
        right = self._print(right)

        return r'{left} \cdot {right}'.format(left=left, right=right)

    def _print_Inner(self, expr):
        left = expr.args[0]
        left = self._print(left)

        right = expr.args[1]
        right = self._print(right)

        return r'{left} : {right}'.format(left=left, right=right)

    def _print_Cross(self, expr):
        left = expr.args[0]
        left = self._print(left)

        right = expr.args[1]
        right = self._print(right)

        return r'{left} \times {right}'.format(left=left, right=right)
    # ...

    # ... forms
    def _print_BilinearForm(self, expr):
        txt = print_integral(expr.domain, expr.expr,
                             measure=expr.measure)
        return txt

    def _print_LinearForm(self, expr):
        txt = print_integral(expr.domain, expr.expr,
                             measure=expr.measure)
        return txt

    def _print_FunctionForm(self, expr):
        txt = print_integral(expr.domain, expr.expr,
                             measure=expr.measure)
        return txt

    def _print_Kron(self, expr):
        raise NotImplementedError('TODO')

    # TODO this implementation works only for the 1D case
    def _print_BilinearAtomicForm(self, expr):
        # TODO move this to __new__ of BilinearAtomicForm
        if not isinstance(expr.expr, Mul):
            raise TypeError('> BilinearAtomicForm must always be of instance Mul')

        coord = expr.coordinates
        args = expr.expr.args
        code = ''
        for i in args:
            atom = get_atom_derivatives(i)
            d = get_index_derivatives(i)
            n_deriv = d[coord.name]
            deriv = '\prime'*n_deriv

            if deriv:
                code = '{code} {atom}^{deriv}'.format(code=code, atom=self._print(atom),
                                                      deriv=deriv)
            else:
                code = '{code} {atom}'.format(code=code, atom=self._print(atom))

        txt = print_integral(expr.domain, code,
                             measure=expr.measure)
        return txt
    # ...

    def _print_Measure(self, expr):
        txt = ''.join(j for j in ['d{}'.format(self._print(i)) for i in expr.args])
        return txt

def latex(expr, **settings):

    return LatexPrinter(settings).doprint(expr)
