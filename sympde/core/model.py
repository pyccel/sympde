# coding: utf-8

# TODO problem if a form is multiplied by a coefficient
#      must change _print_form_call once FormCall is created


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

from .expr import BasicForm, BilinearForm, LinearForm, Integral

class Equation(Assign):
    pass


class Model(Basic):
    """
    Represents a mathematical model.

    Examples

    """
    _name = None
    _forms = None
    _equations = None

    def __new__(cls, **kwargs):

        obj = Basic.__new__(cls)

        # ...
        forms = kwargs.pop('forms', None)
        equations = kwargs.pop('equations', [])
        # ...

        # ...
        if forms is None:
            raise ValueError('> forms must be provided')
        # ...

        # ... set form name
        d_forms = OrderedDict(sorted(forms.items()))
        for name, form in list(d_forms.items()):
            form.set_name(name)
        # ...

        # ...
        name = kwargs.pop('name', None)
        if not name:
            name = obj.__class__.__name__
        # ...

        # ...
        obj._name = name
        obj._forms = d_forms
        obj._equations = equations
        # ...

        return obj

    @property
    def forms(self):
        return self._forms

    @property
    def equations(self):
        return self._equations

    @property
    def name(self):
        return self._name

    def add_equations(self, equations):
        # TODO do we need to check that lhs/rhs are Bilinear/Linear forms?

        if isinstance(equations, Equation):
            equations = [equations]

        elif not isinstance(i, (tuple, list, Tuple)):
            raise TypeError('> Expecting an iterable')

        self._equations += list(equations)

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
