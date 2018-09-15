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
from .expr import is_sum_of_form_calls
from .geometry import Boundary
from .errors import UnconsistentLhsError, UnconsistentRhsError
from .errors import UnconsistentBCError

class BasicBoundaryCondition(Basic):

    def __new__(cls, boundary, value=None):
        return Basic.__new__(cls, boundary, value)

    @property
    def boundary(self):
        return self._args[0]

class DirichletBC(BasicBoundaryCondition):
    pass

# TODO add check on test/trial functions between lhs/rhs
class Equation(Basic):
    def __new__(cls, lhs, rhs, bc=None):
        # ...
        if lhs is None or isinstance(lhs, Nil):
            lhs = Nil()

        elif not isinstance(lhs, FormCall) and not is_sum_of_form_calls (lhs):
            raise UnconsistentLhsError('> lhs must be a call or sum of calls')
        # ...

        # ...
        if rhs is None or isinstance(rhs, Nil):
            rhs = Nil()

        elif not isinstance(rhs, FormCall) and not is_sum_of_form_calls (rhs):
            raise UnconsistentRhsError('> rhs must be a call or sum of calls')
        # ...

        # ...
        if bc:
            if isinstance(bc, BasicBoundaryCondition):
                bc = [bc]

            elif isinstance(bc, (list, tuple, Tuple)):
                for i in bc:
                    if not isinstance(i, BasicBoundaryCondition):
                        msg = '> Expecting a list of BasicBoundaryCondition'
                        raise TypeError(msg)

            else:
                raise TypeError('> Wrong type for bc')

            bc = Tuple(*bc)

            # ... check that the same boundary is not used in the weak
            #     formulation and strong condition
            lhs_bnd = lhs.atoms(Boundary)
            rhs_bnd = rhs.atoms(Boundary)
            bc_bnd  = bc.atoms(Boundary)

            if lhs_bnd & bc_bnd:
                msg = '> {} used for lhs and Dirichlet'.format(lhs_bnd & bc_bnd)
                raise UnconsistentBCError(msg)

            if rhs_bnd & bc_bnd:
                msg = '> {} used for rhs and Dirichlet'.format(rhs_bnd & bc_bnd)
                raise UnconsistentBCError(msg)
            # ...
        # ...

        return Basic.__new__(cls, lhs, rhs, bc)

    @property
    def lhs(self):
        return self._args[0]

    @property
    def rhs(self):
        return self._args[1]

    @property
    def bc(self):
        return self._args[2]

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
