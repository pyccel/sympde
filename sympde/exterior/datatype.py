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


#==============================================================================
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

    def __lt__(self, other):
        assert(isinstance(other, FormType))
        return self.index < other.index

    def __le__(self, other):
        assert(isinstance(other, FormType))
        return self.index <= other.index

    def __gt__(self, other):
        assert(isinstance(other, FormType))
        return self.index > other.index

    def __ge__(self, other):
        assert(isinstance(other, FormType))
        return self.index >= other.index

class ZeroFormType(FormType):
    _index = 0
    pass

class OneFormType(FormType):
    _index = 1
    pass

class TwoFormType(FormType):
    _index = 2
    pass

class ThreeFormType(FormType):
    _index = 3
    pass

class FourFormType(FormType):
    _index = 4
    pass

class FiveFormType(FormType):
    _index = 5
    pass

class SixFormType(FormType):
    _index = 6
    pass


#==============================================================================
# user friendly

ZeroForm  = ZeroFormType()
OneForm   = OneFormType()
TwoForm   = TwoFormType()
ThreeForm = ThreeFormType()
FourForm  = FourFormType()
FiveForm  = FiveFormType()
SixForm   = SixFormType()

dtype_registry = {0: ZeroForm,
                  1: OneForm,
                  2: TwoForm,
                  3: ThreeForm,
                  4: FourForm,
                  5: FiveForm,
                  6: SixForm}


#==============================================================================
def FormTypeFactory(argnames=["_index"], index=None):

    name = 'Dynamic{}FormType'.format(index)

    # ...
    def __init__(self, **kwargs):
        for key, value in list(kwargs.items()):
            # here, the argnames variable is the one passed to the
            # DataTypeFactory call
            if key not in argnames:
                raise TypeError("Argument %s not valid for %s"
                    % (key, self.__class__.__name__))
            setattr(self, key, value)
        FormType.__init__(self)
    # ...

    if index is None:
        raise ValueError('index must be given')

    elif isinstance(index, int) and index <= 6:
        raise ValueError('use dtype_registry for index <= 6')

    newclass = type(name, (FormType,),
                    {"__init__":          __init__,
                     "_index":            index})
    return newclass


#==============================================================================
def get_index_form(dtype):
    if isinstance(dtype, FormType):
        return dtype

    elif isinstance(dtype, int):
        if not(dtype in list(dtype_registry.keys())):
            raise ValueError('Only forms from 0 to 6 are available.')

        return dtype_registry[dtype]

    elif isinstance(dtype, Symbol):
        dtype = FormTypeFactory(index=dtype)
        return dtype()

    else:
        raise TypeError('> Expecting an integer or a datatype')
