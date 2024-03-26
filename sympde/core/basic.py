# coding: utf-8

from sympy        import cacheit
from sympy        import Atom, Expr, AtomicExpr
from sympy        import Function
from sympy        import Number
from sympy        import NumberSymbol
from sympy        import Tuple
from sympy.core   import Basic
from sympy.core   import Symbol
from sympy.tensor import IndexedBase

#==============================================================================
class NoneValue(Symbol):
    """
    Represents a None symbol.
    """
    pass

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
#class BasicConstant(AtomicExpr, Boolean):
#class BasicConstant(AtomicExpr):
class BasicConstant(Symbol):
#    _op_priority    = 10.0

    __slots__ = ('name', 'shape', 'value')

    def __new__(cls, name, shape, value, **assumptions):
        obj = Symbol.__new__(cls, name, **assumptions)
        obj.shape = shape
        obj.value = value

        return obj

    def __getnewargs_ex__(self):
        return ((self.name, self.shape, self.value), self.assumptions0)

    def _hashable_content(self):
        # Note: user-specified assumptions not hashed, just derived ones
        return (self.name, self.shape, self.value) + tuple(sorted(self.assumptions0.items()))

    @property
    def _diff_wrt(self):
        return True

    @property
    def free_symbols(self):
        return {self}

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

#==============================================================================
class ScalarConstant(BasicConstant):
    pass

#==============================================================================
class VectorConstant(BasicConstant):

    def __getitem__(self, key):
        assert isinstance(key, (int, Number, Symbol))
        return IndexedElement(self, key)

#==============================================================================
class MatrixConstant(BasicConstant):

    def __getitem__(self, key):
        assert isinstance(key, (tuple, list))
        for x in key:
            assert isinstance(x, (int, Number, Symbol))
        return IndexedElement(self, key)

#==============================================================================
# TODO ARA must check consistency between shape and value
def constant(name, shape=None, dtype=None, value=None):
    if not isinstance(name, str):
        raise TypeError('Keyword-only argument `name` must be a string.')

    if value is None and dtype is None:
        raise ValueError('dtype must be provided, since value is not  available.')

    # ... check validity and set default arguments values
    if value is None:
        # scalar case
        if shape is None:
            shape = 0

        value = NoneValue('None')

        if not isinstance(shape, (int, tuple, list)):
            raise TypeError('shape must be an integer, tuple or list, but was given {}'.format(type(shape)))
    else:
        if isinstance(value, (int, float, complex)):
            # scalar case
            if dtype is None:
                dtype = type(value)
            elif not (dtype in (int, float, complex)):
                raise TypeError('dtype must be (int, float, complex), but {} was given'.format(dtype))

            # cast type
            value = dtype(value)

            # set default value for shape
            shape = 0

        elif isinstance(value, (tuple, list)):
            if not len({type(_) for _ in value}) == 1:
                raise ValueError('Elements of value must be of the same type.')

            if isinstance(value[0], (int, float, complex)):
                # vector case
                if dtype is None:
                    dtype = type(value[0])
                elif not (dtype in (int, float, complex)):
                    raise TypeError('dtype must be (int, float, complex), but {} was given'.format(dtype))

                # cast type
                value = [dtype(_) for _ in value]
                value = Tuple(*value)

                # set default value for shape
                shape = len(value)

            elif isinstance(value[0], (tuple, list)):
                # matrix case
                for row in value:
                    if not isinstance(row, (tuple, list)):
                        raise TypeError('Each row of `mat` should be a list or a tuple.')

                row_sizes = {len(row) for row in value}
                if len(row_sizes) != 1:
                    raise ValueError('Each row of `value` should have the same length.')

                if dtype is None:
                    dtype = type(value[0][0])
                elif not (dtype in (int, float, complex)):
                    raise TypeError('dtype must be (int, float, complex), but {} was given'.format(dtype))

                # cast type
                value = [[dtype(_) for _ in row] for row in value]
                value = Tuple(*[Tuple(*row) for row in value])

                # set default value for shape
                shape = (len(value), len(value[0]))

            else:
                raise TypeError('Elements must be (int, float, complex) or (tuple, list).')
    # ...

    assumptions = {}
    if dtype == int:
        assumptions['integer'] = True
    elif dtype == float:
        assumptions['real'] = True
    elif dtype == complex:
        assumptions['complex'] = True
    else:
        raise NotImplementedError('')

    if isinstance(shape, int) and shape == 0:
        assumptions['number'] = True
        assumptions['commutative'] = True
        return ScalarConstant(name, shape, value, **assumptions)
    elif isinstance(shape, int) and shape > 0:
        return VectorConstant(name, shape, value, **assumptions)
    elif isinstance(shape, (tuple, list)) and len(shape) == 2:
        return MatrixConstant(name, shape, value, **assumptions)
    else:
        raise ValueError('Wrong arguments, shape = {}'.format(shape))

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
                     ScalarConstant, VectorConstant, MatrixConstant)
