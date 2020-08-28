# -*- coding: utf-8 -*-
#

from sympy.core import Symbol
from sympy import Mul, Tuple
from sympy.printing.latex import LatexPrinter as LatexPrinterSympy
from sympy.printing.latex import translate
from sympy import Indexed, IndexedBase, Matrix, ImmutableDenseMatrix

from sympde.topology import NormalVector, TangentVector
from sympde.topology import Line, Square, Cube, Domain
from sympde.topology.derivatives import sort_partial_derivatives
from sympde.topology.derivatives import get_index_derivatives
from sympde.topology.derivatives import get_atom_derivatives
from sympde.topology.space import ProductSpace
from sympde.topology.space import ScalarTestFunction, VectorTestFunction

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

    def _print_Convect(self, expr):
        raise NotImplementedError('TODO')

    def _print_StrainTensor(self, expr):
        return r'\mathrm{D}\left( ' + self._print(expr.args[0]) + r' \right)'

    def _print_Convolution(self, expr):
        left = self._print(expr.args[0])
        right = self._print(expr.args[1])
        return left + r' \ast ' + right

    def _print_Laplace(self, expr):
        return r'\nabla^2' + self._print(expr.args[0])

    def _print_Hessian(self, expr):
        return r'\mathrm{hess}( ' + self._print(expr.args[0]) + ' )'

    def _print_DifferentialOperator(self, expr):
        coord = self._print(expr.coordinate)
        arg = self._print(expr.args[0])
        return r'\partial_{' + coord + '}' +  arg

    def _print_Trace(self, expr):
        return self._print(expr.expr)

    def _print_MinusInterfaceOperator(self, expr):
        arg = expr.args[0]
        if isinstance(arg, (ScalarTestFunction, VectorTestFunction,
                            Symbol, IndexedBase, Indexed)):
            return self._print(arg) + '_{-}'

        else:
            return r'\left( ' + self._print(arg) + r' \right)_{-}'

    def _print_PlusInterfaceOperator(self, expr):
        arg = expr.args[0]
        if isinstance(arg, (ScalarTestFunction, VectorTestFunction,
                            Symbol, IndexedBase, Indexed)):
            return self._print(arg) + '_{+}'

        else:
            return r'\left( ' + self._print(arg) + r' \right)_{+}'

    def _print_NormalDerivative(self, expr):
        nn = self._print(NormalVector('n'))
        arg = expr.args[0]
        if isinstance(arg, (ScalarTestFunction, VectorTestFunction)):
            arg = r'\nabla ' + self._print(arg)

        else:
            arg = r'\nabla \left( ' + self._print(arg) + r' \right)'

        return arg + r' \cdot ' + nn
    # ...

    # ...
    def _print_Jump(self, expr):
        return r'[\![ ' + self._print(expr.args[0]) + r' ]\!]'
    # ...

    # ...
    def _print_Average(self, expr):
        return r'\{ ' + self._print(expr.args[0]) + r' \}'
    # ...

    # ...
    def _print_NormalVector(self, expr):
        return r'\mathbf{n}'

    def _print_MinusNormalVector(self, expr):
        return r'\mathbf{n}_{-}'

    def _print_PlusNormalVector(self, expr):
        return r'\mathbf{n}_{+}'

    def _print_TangentVector(self, expr):
        return r'\mathbf{t}'
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

    def _print_Outer(self, expr):
        raise NotImplementedError('TODO')
    # ...

    # ...
    def _print_Integral(self, expr):

        if expr.is_domain_integral:
            dim = expr.domain.dim
            if dim == 1:
                dx = r'~dx'

            else:
                dx = r'~d\mathbf{x}'
        else:
            dx = r'~ds'

        domain = self._print(expr.domain)
        expr   = self._print(expr.expr)

        integral = r'\int_{' + domain + '}'

        return '{integral} {expr} {dx}'.format(integral=integral,
                                               expr=expr,
                                               dx=dx)

    # ... forms
    def _print_BilinearForm(self, expr):
        trials, tests = expr.variables
        if len(trials) == 1: trials = trials[0]
        if len(tests) == 1: tests = tests[0]
        trials = self._print(trials)
        tests  = self._print(tests)
        expr = self._print(expr.expr)
        pattern = r'\left( {trials}, {tests} \right) \mapsto {expr}'
        return pattern.format(trials=trials, tests=tests, expr=expr)

    def _print_LinearForm(self, expr):
        tests = expr.variables
        if len(tests) == 1: tests = tests[0]
        tests  = self._print(tests)
        expr = self._print(expr.expr)
        pattern = r'{tests} \mapsto {expr}'
        return pattern.format(tests=tests, expr=expr)

    def _print_Functional(self, expr):
        raise NotImplementedError('TODO')

    def _print_Kron(self, expr):
        raise NotImplementedError('TODO')
    # ...

    def _print_Nil(self, expr):
        return r'\ldots'

    def _print_Domain(self, expr):
        return translate(expr.name)

    def _print_Interval(self, expr):
        return '{}'.format(self._print(expr.name))

    def _print_ProductSpace(self, expr):
        spaces = [self._print(i) for i in expr.spaces]
        return r' \times '.join(spaces)

    def _print_Equation(self, expr):
        lhs = expr.lhs
        rhs = expr.rhs
        trials, tests = lhs.variables

        lhs = self._print(lhs(trials, tests))
        rhs = self._print(rhs(tests))

        prelude  = ''
        epilogue = ''

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

        prelude = ( r'\mbox{find} ~'
                   + trials + r' \in ' + trial_spaces
                   + r' ~\mbox{such that}\\')
        epilogue = r',\quad \forall~ ' + tests + r' \in ' + test_spaces

        body = '{lhs} = {rhs}'.format(lhs=lhs, rhs=rhs)

        code = prelude + body + epilogue
        return code

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

    # ........................................
    #            EXTERIOR CALCULUS
    # ........................................
    def _print_ExteriorDerivative(self, expr):
        arg = self._print(expr._args[0])
        return r'\mathrm{d} ' + arg

    def _print_AdjointExteriorDerivative(self, expr):
        arg = self._print(expr._args[0])
        # TODO which one to take?
#        return r'\mathrm{\delta} ' + arg
        return r'\star\mathrm{d} ' + arg

    def _print_ExteriorProduct(self, expr):
        left, right = expr._args[:]
        left = self._print(left)
        right = self._print(right)
        return left + r' \wedge ' + right

    def _print_InteriorProduct(self, expr):
        left, right = expr._args[:]
        left = self._print(left)
        right = self._print(right)
        return left + r' \lrcorner ' + right

    def _print_AdjointInteriorProduct(self, expr):
        left, right = expr._args[:]
        left = self._print(left)
        right = self._print(right)
        return r'\star \left(' + left + r' \lrcorner ' + right + r' \right)'

    def _print_Hodge(self, expr):
        arg = self._print(expr._args[0])
        return r'\star ' + arg

    def _print_PullBack(self, expr):
        raise NotImplementedError('')
    # ........................................

    def _print_KernelExpression(self, expr):
        target = self._print(expr.target)
        expr   = self._print(expr.expr)
        return r'\textbf{target}~' + target + r'\\' + expr

def latex(expr, **settings):

    return LatexPrinter(settings).doprint(expr)
