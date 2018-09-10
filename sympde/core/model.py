# coding: utf-8

from numpy import unique
from collections import OrderedDict

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Integer, Float
from sympy.printing.preview import preview as sympy_preview

from pyccel.ast.core import Assign
from pyccel.ast.core import Nil

from .expr import BasicForm, BilinearForm, LinearForm, FunctionForm


class Model(Basic):
    """
    Represents a mathematical model.

    Examples

    """
    def __new__(cls, **kwargs):
        # TODO this does not work if a form is multiplied by a coefficient

        # ...
        nil = Nil()
        equations = []
        d_forms = {}
        for name, i in list(kwargs.items()):
            if isinstance(i, BasicForm):
                # ...
                signature = None
                if isinstance(i, BilinearForm):
                    test = i.test_functions
                    if len(test) == 1:
                        test = test[0]

                    else:
                        test = Tuple(*test)

                    trial = i.trial_functions
                    if len(trial) == 1:
                        trial = trial[0]

                    else:
                        trial = Tuple(*trial)

                    signature = Tuple(test, trial)

                elif isinstance(i, LinearForm):
                    test = i.test_functions
                    if len(test) == 1:
                        test = test[0]

                    else:
                        test = Tuple(*test)

                    signature = test
                # ...

                # ...
                if signature:
                    lhs = Function(name)(*signature)
                else:
                    lhs = Symbol(name)
                # ...

                stmt = Assign(lhs, i)
                equations.append(stmt)
                d_forms[name] = i

            elif isinstance(i, (tuple, list, Tuple)):
                if not(len(i) == 2):
                    raise ValueError('> Expecting two elements in tuple/list/Tuple')

                lhs = i[0] ; rhs = i[1]
                if not isinstance(lhs, BilinearForm):
                    raise TypeError('> Expecting a BilinearForm as a lhs')

                if not isinstance(rhs, LinearForm):
                    raise TypeError('> Expecting a LinearForm as a rhs')

                stmt = Assign(lhs, rhs)
                equations.append(stmt)
                # TODO add Equation class, that is like a dictionary
                d_forms[name] = lhs
                raise NotImplementedError('TODO')
        # ...

        name = kwargs.pop('name', None)

        obj = Basic.__new__(cls, equations)

        # ...
        if not name:
            name = obj.__class__.__name__
        obj._name = name
        # ...

        # ...
        obj._forms = OrderedDict(sorted(d_forms.items()))
        # ...

        return obj

    @property
    def name(self):
        return self._name

    @property
    def equations(self):
        return self._args[0]

    @property
    def forms(self):
        return self._forms

    def preview(self, euler=False, packages=None,
                output='dvi', outputTexFile=None):
        # ...
        default_packages = ('amsmath', 'amsfonts')
        if packages:
            packages = tuple(packages)
        else:
            packages = ()

        includes = '\n'.join(['\\usepackage{%s}' % p for p in default_packages + packages])

        preamble = r"\documentclass[a4paper,9pt]{article}"
        preamble = preamble + '\n' + includes
        preamble = preamble + '\n' + r'\begin{document}'
        # ...

        from sympde.printing.latex import latex
        sympy_preview(latex(self), euler=euler, preamble=preamble,
                      output=output, outputTexFile=outputTexFile)
