# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Integer, Float

from .expr import BasicForm, BilinearForm, LinearForm, FunctionForm

class Model(Basic):
    """
    Represents a mathematical model.

    Examples

    """
    def __new__(cls, *args, **kwargs):
        forms = [i for i in args if isinstance(i, BasicForm)]
        name = kwargs.pop('name', None)

        obj = Basic.__new__(cls, forms)

        # ...
        if not name:
            name = obj.__class__.__name__
        obj._name = name
        # ...

        return obj

    @property
    def name(self):
        return self._name

    @property
    def forms(self):
        return self._args[0]

    @property
    def bilinear_forms(self):
        return [i for i in self.forms if isinstance(i, BilinearForm)]

    @property
    def linear_forms(self):
        return [i for i in self.forms if isinstance(i, LinearForm)]

    @property
    def function_forms(self):
        return [i for i in self.forms if isinstance(i, FunctionForm)]
