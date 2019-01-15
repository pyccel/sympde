# -*- coding: utf-8 -*-
#

from sympy.core import Symbol
from sympy import Mul, Tuple
from sympy.printing.latex import LatexPrinter as LatexPrinterSympy

from pyccel.ast.core import Nil

from sympde.expr.expr import BilinearForm, LinearForm, Integral, FormCall
from sympde.calculus import Dot, Inner, Cross
from sympde.calculus import Grad, Rot, Curl, Div
from sympde.topology import Line, Square, Cube, Domain
from sympde.topology.derivatives import sort_partial_derivatives
from sympde.topology.derivatives import get_index_derivatives
from sympde.topology.derivatives import get_atom_derivatives
from sympde.topology.space import ProductSpace

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

    def _print_Bracket(self, expr):
        u,v = [self._print(i) for i in expr.args]
        return r'[' + u + ',' + v + ']'

    def _print_DifferentialOperator(self, expr):
        coord = self._print(expr.coordinate)
        arg = self._print(expr.args[0])
        return r'\partial_{' + coord + '}' +  arg
    # ...

    # ...
    def _print_VectorTestFunction(self, expr):
        return r'\mathbf{' + self._print(expr.name) + '}'
    # ...

    # ...
    def _print_VectorField(self, expr):
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
        calls = expr.atoms(FormCall)
        if not calls:
            domain = expr.domain
            if expr.boundary:
                domain = expr.boundary
            return self._write_integral_expr(domain, expr.expr,
                                             measure=expr.measure)
        else:
            return self._print(expr.expr)

    def _print_LinearForm(self, expr):
        calls = expr.atoms(FormCall)
        if not calls:
            domain = expr.domain
            if expr.boundary:
                domain = expr.boundary
            return self._write_integral_expr(domain, expr.expr,
                                             measure=expr.measure)
        else:
            return self._print(expr.expr)

    def _print_Integral(self, expr):
        calls = expr.atoms(FormCall)
        if not calls:
            domain = expr.domain
            if expr.boundary:
                domain = expr.boundary
            return self._write_integral_expr(domain, expr.expr,
                                             measure=expr.measure)
        else:
            return self._print(expr.expr)

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

        txt = self._write_integral_expr(expr.domain, code, measure=expr.measure)
        return txt
    # ...

    def _print_Nil(self, expr):
        return '\ldots'

    def _print_Domain(self, expr):
        return '{}'.format(self._print(expr.name))

    def _print_Interval(self, expr):
        return '{}'.format(self._print(expr.name))

    def _print_ProductSpace(self, expr):
        spaces = [self._print(i) for i in expr.spaces]
        return r' \times '.join(spaces)

    def _print_FormCall(self, expr):
        form = expr.expr
        name = form.name
        if name is None:
            raise ValueError('> undefined name for a form')

        args = []

        # ...
        test = [self._print(i) for i in form.test_functions]
        test_str = ','.join(i for i in test)

        if len(test) == 1:
            test = test_str

        else:
            test = '({})'.format(test_str)

        args += [test]
        # ...

        # ...
        if isinstance(form, BilinearForm):
            # ...
            trial = [self._print(i) for i in form.trial_functions]
            trial_str = ','.join(i for i in trial)

            if len(trial) == 1:
                trial = trial_str

            else:
                trial = '({})'.format(trial_str)

            args += [trial]
            # ...
        # ...

        args = ', '.join(args)
        name = self._print(name)
        return '{name}({args})'.format(name=name, args=args)

    def _print_Equation(self, expr):
        lhs = self._print(expr.lhs)
        rhs = self._print(expr.rhs)

        prelude  = ''
        epilogue = ''

        tests, trials = expr.lhs.arguments

        # ... we print tests function here, since sympy printing does not look
        #     nice for a tuple (quad)
        if isinstance(tests, (list, tuple, Tuple)):
            if len(tests) == 1:
                tests = tests[0]
                test_spaces = tests.space

                tests = self._print(tests)

            else:
                test_spaces = [u.space for u in tests]
                test_spaces = ProductSpace(*test_spaces)

                tests = ','.join([self._print(i) for i in tests])
                tests = '({})'.format(tests)

        else:
            test_spaces = tests.space
            tests = self._print(tests)
        # ...

        # ...
        if isinstance(trials, (list, tuple, Tuple)):
            if len(trials) == 1:
                trials = trials[0]
                trial_spaces = trials.space

                trials = self._print(trials)

            else:
                trial_spaces = [u.space for u in trials]
                trial_spaces = ProductSpace(*trial_spaces)

                trials = ','.join([self._print(i) for i in trials])
                trials = '({})'.format(trials)

        else:
            trial_spaces = trials.space
            trials = self._print(trials)
        # ...

        test_spaces = self._print(test_spaces)
        trial_spaces = self._print(trial_spaces)

        prelude = 'find $' + trials + ' \in ' + trial_spaces + '$' + ' such that'
        epilogue = r',\quad \forall ' + tests + ' \in ' + test_spaces

        if expr.is_undefined:
            body = ''
        else:
            body = '{lhs} &= {rhs}'.format(lhs=lhs, rhs=rhs)

        code = (prelude + '\n' +
                r'\begin{align*}' +
                body +
                epilogue +
                r'\end{align*}')
        return code

    def _print_Measure(self, expr):
        txt = ''.join(j for j in ['d{}'.format(self._print(i)) for i in expr.args])
        return txt

    def _print_dx(self, expr):
        arg = self._print(expr.args[0])
        return r'\partial_x {}'.format(arg)

    def _print_dy(self, expr):
        arg = self._print(expr.args[0])
        return r'\partial_y {}'.format(arg)

    def _print_dz(self, expr):
        arg = self._print(expr.args[0])
        return r'\partial_z {}'.format(arg)

    def _print_Boundary(self, expr):
        return self._print(expr.name)

    def _print_BoundaryVector(self, expr):
        name = self._print(expr.name)
        return r'\mathbf{' + name + '}'

    def _write_integral_expr(self, domain, expr, measure=None):
        measure_str = ''
        if not(measure is None):
            measure_str = latex(measure)

        expr_str = latex(expr)

        # TODO improve this using latex for ProductDomain
        int_str = r'\int_{' + self._print(domain) + '}'
        measure_str = ''

        return '{integral} {expr} {measure}'.format(integral=int_str,
                                                    expr=expr_str,
                                                    measure=measure_str)

def latex(expr, **settings):

    return LatexPrinter(settings).doprint(expr)
