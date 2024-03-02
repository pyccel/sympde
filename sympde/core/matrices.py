from sympy                 import Expr, S
from sympy                 import Add, Mul, Pow
from sympy                 import Tuple
from sympy                 import sympify
from sympy.core.decorators import call_highest_priority
from sympde.core.basic     import _coeffs_registery, Basic

from sympy import ImmutableDenseMatrix

import functools


class MatrixSymbolicExpr(Expr):
    is_commutative = False
    _op_priority   = 20.0
    is_MatrixSymbolicExpr = True
    is_Matrix     = True
    is_ZeroMatrix = False
    is_Identity   = False

    def __new__(cls, *args, **options):
        return Expr.__new__(cls, *args)

    def inv(self):
        return Inverse(self)

    def inverse(self):
        return Inverse(self)

    def transpose(self):
        return Transpose(self)

    @property
    def T(self):
        return self.transpose()

    def det(self):
        return SymbolicDeterminant(self)

    def trace(self):
        return SymbolicTrace(self)

    # The following is adapted from the core Expr object
    def __neg__(self):
        return MatSymbolicMul(S.NegativeOne, self)

    def __abs__(self):
        return MatAbs(self)

    def __getitem__(self, key):
        return MatrixElement(self, key)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        return MatSymbolicAdd(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return MatSymbolicAdd(other, self)

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return MatSymbolicAdd(self, -other)

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return MatSymbolicAdd(other, -self)

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return MatSymbolicMul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return MatSymbolicMul(other, self)

    @call_highest_priority('__pow__')
    def __pow__(self, other):
        return MatSymbolicPow(self, other)

    @call_highest_priority('__div__')
    def __div__(self, other):
        return MatSymbolicMul(self , other**S.NegativeOne)

    @call_highest_priority('__rdiv__')
    def __rdiv__(self, other):
        return MatSymbolicMul(other, Pow(self, S.NegativeOne))

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

class Inverse(MatrixSymbolicExpr):
    is_Matrix     = False
    is_commutative = False
    def __new__(cls, *args, **options):
        assert len(args) == 1
        if isinstance(args[0], Inverse):
            return args[0].arg
        return Expr.__new__(cls, *args)

    @property
    def arg(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        arg  = self.arg
        if not arg.is_Matrix:
            return '({})**(-1)'.format(sstr(arg))
        return '{}**(-1)'.format(sstr(arg))

class Transpose(MatrixSymbolicExpr):
    is_commutative = False
    is_Matrix      = False
    def __new__(cls, *args, **options):
        assert len(args) == 1
        if isinstance(args[0], Transpose):
            return args[0].arg
        elif isinstance(args[0], Add):
            return Add(*[Transpose(a) for a in args[0].args])
        elif isinstance(args[0], Mul):
            coeffs = [a for a in args[0].args if a.is_commutative]
            args   = [a for a in args[0].args if not a.is_commutative]
            return Mul(*coeffs)*Expr.__new__(cls, Mul(*args))
        return Expr.__new__(cls, *args)

    @property
    def arg(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        arg  = self.arg
        if not arg.is_Matrix:
            return '({}).T'.format(sstr(arg))
        return '{}.T'.format(sstr(arg))

class MatSymbolicMul(MatrixSymbolicExpr, Mul):
    is_Matrix     = False
    def __new__(cls, *args, **options):
        args = [sympify(a) for a in args if a != 1]

        for i,a in enumerate(args):
            if isinstance(a, Add):
                newargs = []
                for e in a.args:
                    t = args[:i] + [e] + args[i+1:]
                    newargs.append(MatSymbolicMul(*t))
                return type(a)(*newargs)

        newargs = []
        for a in args:
            if isinstance(a, Mul):
                newargs += list(a.args)
            else:
                newargs.append(a)

        args   = [a for a in newargs if not a.is_commutative]
        coeffs = [a for a in newargs if a.is_commutative]

        if coeffs:
            c = Mul(*coeffs)
            if not args:
                return c
            elif isinstance(c, Mul):
                args = list(c.args) + args
            elif c != 1:
                args = [c] + args

        if len(args) == 0:
            return S.One
        elif len(args) == 1:
            return args[0]

        return Expr.__new__(cls, *args)

    def _sympystr(self, printer):
        sstr = printer.doprint
        code = sstr(self.args[0])
        for e in self.args[1:]:
            if isinstance(e, Add):
                code += ' * ({})'.format(sstr(e))
            else:
                code += ' * {}'.format(sstr(e))

        return code

    @staticmethod
    def _expandsums(sums):
        """
        Helper function for _eval_expand_mul.
        sums must be a list of instances of Basic.
        """

        L = len(sums)
        if L == 1:
            return sums[0].args
        terms = []
        left = MatSymbolicMul._expandsums(sums[:L//2])
        right = MatSymbolicMul._expandsums(sums[L//2:])
        terms = [MatSymbolicMul(a, b) for a in left for b in right]
        added = MatSymbolicAdd(*terms)
        return Add.make_args(added) # it may have collapsed down to one term

class MatSymbolicAdd(MatrixSymbolicExpr, Add):
    is_Matrix     = False
    def __new__(cls, *args, **options):
        args = [sympify(a) for a in args if a != 0]
        newargs = []
        for i in args:
            if isinstance(i, (MatSymbolicAdd, Add)):
                newargs += list(i.args)
            else:
                newargs.append(i)

        args = newargs

        if len(args) == 0:
            return S.Zero
        elif len(args) == 1:
            return args[0]
        args = sorted(args, key=str)
        return Expr.__new__(cls, *args)

    def _sympystr(self, printer):
        sstr = printer.doprint
        #return ' + '.join((sstr(i) for i in self.args))
        return 'MatSymbolicAdd({})'.format(','.join((sstr(i) for i in self.args)))

class MatSymbolicPow(MatrixSymbolicExpr, Pow):
    is_Matrix     = False
    def _sympystr(self, printer):
        sstr = printer.doprint
        #return '**'.join((sstr(i) for i in self.args))
        return 'MatSymbolicPow({})'.format(','.join((sstr(i) for i in self.args)))

class MatSymbolicAbs(MatrixSymbolicExpr):
    def _sympystr(self, printer):
        sstr = printer.doprint
        return 'Abs({})'.format(sstr(self.args[0]))

class SymbolicDeterminant(Expr):
    is_commutative = True
    def __new__(cls, *args, **options):
        assert len(args) == 1
        return Basic.__new__(cls, *args)

    @property
    def arg(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        arg = sstr(self.arg)
        return 'det({})'.format(arg)

class SymbolicTrace(Expr):
    is_commutative = True
    def __new__(cls, *args, **options):
        assert len(args) == 1
        arg = args[0]
        if isinstance(arg, Add):
            return Add(*[SymbolicTrace(a) for a in arg.args])

        elif isinstance(arg, Mul):
            coeffs = [a for a in arg.args if isinstance(a, _coeffs_registery)]
            mats   = [a for a in arg.args if not a in coeffs]
            if coeffs:
                return Mul(*coeffs)*SymbolicTrace(arg.func(*mats))

        return Basic.__new__(cls, *args)

    @property
    def arg(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        arg = sstr(self.arg)
        return 'tr({})'.format(arg)

class MatrixElement(Expr):
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
        return '{}[{}]'.format(sstr(self.args[0]),sstr(self.args[1]))

class Matrix(MatrixSymbolicExpr):
    def __new__(cls, mat, *, name):
        if not isinstance(mat, list):
            raise TypeError('Positional argument `mat` should be a list of lists.')

        for row in mat:
            if not isinstance(row, list):
                raise TypeError('Each row of `mat` should be a list.')

        row_sizes = {len(row) for row in mat}
        if len(row_sizes) != 1:
            raise ValueError('Each row of `mat` should have the same length.')

        if not isinstance(name, str):
            raise TypeError('Keyword-only argument `name` must be a string.')

        # Call constructor of base class Expr with immutable arguments
        new_mat = Tuple(*[Tuple(*row) for row in mat])
        obj = Expr.__new__(cls, new_mat, name)
        return obj

    @property
    def name(self):
        return self.args[1]

    def __getitem__(self, key):
        i,j = key
        return self.args[0][i][j]

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

    def __hash__(self):
        return hash((self.name, self.args))

    def __eq__(self, other):
        if not isinstance(other, Matrix):
            return False

        if not self.name == other.name:
            return False

        result = all(x == y for x, y in zip(self.args[0], other.args[0]))
        return result

    def to_sympy(self):
        return ImmutableDenseMatrix(self.args[0])


class Vector(MatrixSymbolicExpr):

    def __new__(cls, vec, *, name):
        if not isinstance(vec, list):
            raise TypeError('Positional argument `vec` should be a list.')

        if not isinstance(name, str):
            raise TypeError('Keyword-only argument `name` must be a string.')

        # Call constructor of base class Expr with immutable arguments
        new_vec = Tuple(*vec)
        obj = Expr.__new__(cls, new_vec, name)
        return obj

    @property
    def name(self):
        return self.args[1]

    def __getitem__(self, key):
        i = key
        return self.args[0][i]

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}'.format(sstr(self.name))

    def __hash__(self):
        return hash((self.name, self.args))

    def __eq__(self, other):
        # TODO BUG
        # We should check that other is a Vector
        # right now, there is a problem when using Cross
        # see the linearity test of
        # f = lambda u,v: cross(b*u, b*v)
        # where b is a Vector
        if not isinstance(other, Vector):
            return False

        if not self.name == other.name:
            return False

        result = self.args[0] == other.args[0]
        return result

    def to_sympy(self):
        args = [[a] for a in self.args[0]]
        return ImmutableDenseMatrix(args)


Basic._constructor_postprocessor_mapping[MatrixSymbolicExpr] = {
    "Mul": [lambda x: MatSymbolicMul(*x.args)],
    "Add": [lambda x: MatSymbolicAdd(*x.args)]
}

