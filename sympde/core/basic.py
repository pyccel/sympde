# coding: utf-8

from sympy        import Atom, Expr
from sympy        import Function
from sympy        import Number
from sympy        import NumberSymbol
from sympy        import Tuple
from sympy.core   import Basic
from sympy.core   import Symbol
from sympy.tensor import IndexedBase

#==============================================================================
class IndexedElement(Expr):
    is_number = True

    def __new__(cls, base, indices, **options):
        return Expr.__new__(cls, base, indices)

    @property
    def base(self):
        return self._args[0]

    @property
    def indices(self):
        return self._args[1]

    def _sympystr(self, printer):
        sstr = printer.doprint
        var = self.args[0]
        key = self.args[1]
        if isinstance(key, (tuple, list, Tuple)):
            key = ','.join(sstr(_) for _ in key)
            return '{}[{}]'.format(sstr(var),key)
        else:
            return '{}[{}]'.format(sstr(var),sstr(key))

#==============================================================================
class BasicConstant(Expr):
    is_commutative = False
    _op_priority   = 20.0
    _shape = None

    @property
    def name(self):
        return self.args[0]

    @property
    def value(self):
        return self._value

    @property
    def dtype(self):
        return self._dtype

    @property
    def shape(self):
        return self._shape

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

#==============================================================================
class ScalarConstant(BasicConstant):
    is_number = True

    def __new__(cls, *args, **options):
        assert len(args) == 1

        name = args[0]
        if not isinstance(name, str):
            raise TypeError('Keyword-only argument `name` must be a string.')

        value = options.pop('value', None)
        dtype = options.pop('dtype', None)

        if value is None:
            if dtype is None:
                raise ValueError('dtype must be provided, since value is not  available.')
        else:
            if not dtype is None:
                value = dtype(value)
            else:
                dtype = type(value)

        assert dtype in (int, float, complex)

        obj = Expr.__new__(cls, name)
        obj._value = value
        obj._dtype = dtype
        obj._shape = 1
        obj.is_integer = dtype == int
        obj.is_real    = dtype == float
        obj.is_complex = dtype == complex
        obj.is_commutative = True

        return obj

#==============================================================================
class VectorConstant(BasicConstant):

    def __new__(cls, *args, **options):
        assert len(args) == 1

        name = args[0]
        if not isinstance(name, str):
            raise TypeError('Keyword-only argument `name` must be a string.')

        value = options.pop('value', None)
        dtype = options.pop('dtype', None)
        shape = options.pop('shape', None)

        if value is None:
            if dtype is None:
                raise ValueError('dtype must be provided, since value is not  available.')
            if shape is None:
                raise ValueError('shape must be provided, since value is not  available.')

            assert isinstance(shape, int)
            assert dtype in (int, float, complex)
        else:
            assert isinstance(value, (tuple, list))
            if not len({type(_) for _ in value}) == 1:
                raise ValueError('Elements of value must be of the same type.')

            if not dtype is None:
                value = [dtype(_) for _ in value]
            else:
                dtype = type(value[0])

            if not shape is None:
                assert shape == len(value)
            else:
                shape = len(value)

            value = Tuple(*value)

        obj = Expr.__new__(cls, name)
        obj._value = value
        obj._dtype = dtype
        obj._shape = shape
        obj.is_integer = dtype == int
        obj.is_real    = dtype == float
        obj.is_complex = dtype == complex

        return obj

    def __getitem__(self, key):
        assert isinstance(key, (int, Number, Symbol))
        return IndexedElement(self, key)

#==============================================================================
class MatrixConstant(BasicConstant):

    def __new__(cls, *args, **options):
        assert len(args) == 1

        name = args[0]
        if not isinstance(name, str):
            raise TypeError('Keyword-only argument `name` must be a string.')

        value = options.pop('value', None)
        dtype = options.pop('dtype', None)
        shape = options.pop('shape', None)

        if value is None:
            if dtype is None:
                raise ValueError('dtype must be provided, since value is not  available.')
            if shape is None:
                raise ValueError('shape must be provided, since value is not  available.')

            assert isinstance(shape, (tuple, list))
            assert dtype in (int, float, complex)
        else:
            assert isinstance(value, (tuple, list))
            for row in value:
                if not isinstance(row, (tuple, list)):
                    raise TypeError('Each row of `mat` should be a list or a tuple.')

            row_sizes = {len(row) for row in value}
            if len(row_sizes) != 1:
                raise ValueError('Each row of `value` should have the same length.')

            if not dtype is None:
                value = [[dtype(_) for _ in row] for row in value]
            else:
                dtype = type(value[0][0])

            _shape = (len(value), len(value[0]))
            if not shape is None:
                shape = tuple(shape)
                assert shape == _shape
            else:
                shape = _shape

            value = Tuple(*[Tuple(*row) for row in value])

        obj = Expr.__new__(cls, name)
        obj._value = value
        obj._dtype = dtype
        obj._shape = shape
        obj.is_integer = dtype == int
        obj.is_real    = dtype == float
        obj.is_complex = dtype == complex

        return obj

    def __getitem__(self, key):
        assert isinstance(key, (tuple, list))
        for x in key:
            assert isinstance(x, (int, Number, Symbol))
        return IndexedElement(self, key)

#==============================================================================
def constant(name, shape=None, dtype=None, value=None):
    if value is None:
        assert dtype in (int, float, complex)

        if (shape == 0) or (shape is None):
            return ScalarConstant(name, dtype=dtype)
        elif isinstance(shape, int):
            return VectorConstant(name, shape=shape, dtype=dtype)
        elif isinstance(shape, (tuple, list)):
            return MatrixConstant(name, shape=shape, dtype=dtype)
        else:
            raise ValueError('Wrong arguments')
    else:
        if isinstance(value, (int, float, complex)):
            assert (shape == 0) or (shape is None)
            return ScalarConstant(name, dtype=dtype, value=value)
        elif isinstance(value, (tuple, list)):
            if isinstance(shape, int):
                return VectorConstant(name, shape=shape, dtype=dtype, value=value)
            elif isinstance(shape, (tuple, list)):
                return MatrixConstant(name, shape=shape, dtype=dtype, value=value)
            else:
                raise ValueError('Wrong arguments')
        else:
            raise ValueError('Wrong arguments')

#==============================================================================
# TODO to be removed
class Constant(Symbol):
    """
    Represents a constant symbol.

    Examples

    """
    _label = ''
    is_number = True
    def __new__(cls, *args, **kwargs):
        label = kwargs.pop('label', '')

        obj = Symbol.__new__(cls, *args, **kwargs)
        obj._label = label
        return obj

    @property
    def label(self):
        return self._label

#==============================================================================
class CalculusFunction(Function):
    """this class is needed to distinguish between functions and calculus
    functions when manipulating our expressions"""
    pass

#==============================================================================
class BasicMapping(IndexedBase):
    """
    Represents a basic class for mapping.
    """
    pass

#==============================================================================
class BasicDerivable(Basic):
    pass

#==============================================================================
_coeffs_registery = (int, float, complex, Number, NumberSymbol,
                     Constant,
                     ScalarConstant, VectorConstant, MatrixConstant)
