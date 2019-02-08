# coding: utf-8
# TODO: - Unknown is not used here (mlhipy) remove it?

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Integer, Float
from sympy.core.singleton import Singleton
from sympy.core.compatibility import with_metaclass
from sympy.core import Add, Mul
from sympy.core.singleton import S

from sympde.core.basic import _coeffs_registery
from sympde.core import LinearOperator


# ...
class FormType(with_metaclass(Singleton, Basic)):
    """Base class representing differential form types"""
    _index = None

    @property
    def index(self):
        return self._index

    @property
    def name(self):
        if self.index is None:
            return None

        else:
            return 'V_{}'.format(self.index)

    def __str__(self):
        return str(self.name)

class TypeZeroForm(FormType):
    _index = 0
    pass

class TypeOneForm(FormType):
    _index = 1
    pass

class TypeTwoForm(FormType):
    _index = 2
    pass

class TypeThreeForm(FormType):
    _index = 3
    pass

class TypeFourForm(FormType):
    _index = 4
    pass

class TypeFiveForm(FormType):
    _index = 5
    pass

class TypeSixForm(FormType):
    _index = 6
    pass

dtype_registry = {0: TypeZeroForm(),
                  1: TypeOneForm(),
                  2: TypeTwoForm(),
                  3: TypeThreeForm(),
                  4: TypeFourForm(),
                  5: TypeFiveForm(),
                  6: TypeSixForm()}

def get_index_form(dtype):
    if isinstance(dtype, FormType):
        return dtype

    elif isinstance(dtype, int):
        if not(dtype in list(dtype_registry.keys())):
            raise ValueError('Only forms from 0 to 6 are available.')

        return dtype_registry[dtype]

    else:
        raise TypeError('> Expecting an integer or a datatype')
# ...
