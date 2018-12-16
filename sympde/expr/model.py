# coding: utf-8

from collections import OrderedDict

from sympy.core import Basic
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy.printing.preview import preview as sympy_preview

from pyccel.ast.core import Nil

from sympde.core.generic import Dot, Inner, Cross
from sympde.core.generic import Grad, Rot, Curl, Div
from sympde.topology.basic import Boundary
from sympde.topology.space import TestFunction, VectorTestFunction
from sympde.topology.space import FunctionSpace

from sympde.expr.errors import UnconsistentLhsError, UnconsistentRhsError
from sympde.expr.errors import UnconsistentBCError
from sympde.expr.errors import UnconsistentArgumentsError
from sympde.expr.expr import BasicForm, BilinearForm, LinearForm, Integral, FormCall
from sympde.expr.expr import is_sum_of_form_calls
from sympde.expr.utils import random_string
from sympde.expr.equation import Equation, BasicBoundaryCondition, DirichletBC


#==============================================================================
class Model(Basic):
    """
    Represents a mathematical model.

    Examples

    """
    _name = None
    _domain = None
    _forms = None
    _equation = None

    def __new__(cls, domain, **kwargs):

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
        obj._domain = domain
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

    def set_equation(self, *args, **kwargs):
        self._equation = Equation(*args, **kwargs)

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
