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

from .expr import BasicForm, BilinearForm, LinearForm, Integral, FormCall

class Equation(Basic):
    def __new__(cls, lhs, rhs):
        # ...
        if lhs is None or isinstance(lhs, Nil):
            lhs = Nil()

        elif not(isinstance(lhs, FormCall) and isinstance(lhs.expr, BilinearForm)):
            raise TypeError('> wrong type for lhs')
        # ...

        # ...
        if rhs is None or isinstance(rhs, Nil):
            rhs = Nil()

        elif not(isinstance(rhs, FormCall) and isinstance(rhs.expr, LinearForm)):
            raise TypeError('> wrong type for rhs')
        # ...


        return Basic.__new__(cls, lhs, rhs)

    @property
    def lhs(self):
        return self._args[0]

    @property
    def rhs(self):
        return self._args[1]

    @property
    def is_undefined(self):
        return isinstance(self.lhs, Nil) or isinstance(self.rhs, Nil)


class Model(Basic):
    """
    Represents a mathematical model.

    Examples

    """
    _name = None
    _domain = None
    _forms = None
    _equation = None

    def __new__(cls, **kwargs):

        obj = Basic.__new__(cls)

        # ...
        forms = kwargs.pop('forms', None)
        equation = kwargs.pop('equation', None)
        # ...

        # ...
        if forms is None:
            raise ValueError('> forms must be provided')

        elif not isinstance(forms, (tuple, list, Tuple)):
            raise TypeError('> Expecting an iterable')
        # ...

        # ... set form name
        d_forms = {}
        for form in forms:
            name = form.name
            if name is None:
                raise ValueError('> Bilinear and Linear forms must be assigned a name')

            d_forms[name] = form

        d_forms = OrderedDict(sorted(d_forms.items()))
        # ...

        # ...
        name = kwargs.pop('name', None)
        if not name:
            name = obj.__class__.__name__
        # ...

        # ...
        obj._name = name
        obj._forms = d_forms
        obj._equation = equation
        # ...

        return obj

    @property
    def forms(self):
        return self._forms

    @property
    def equation(self):
        return self._equation

    @property
    def name(self):
        return self._name

    @property
    def domain(self):
        return self._domain

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
