from sympy                 import Expr, S
from sympy                 import Add, Mul, Pow
from sympy                 import sympify
from sympy.core.decorators import call_highest_priority
from sympde.core.basic     import _coeffs_registery, Basic

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
        return '{}**(-1)'.format(sstr(self.arg))

class Transpose(MatrixSymbolicExpr):
    is_commutative = False
    def __new__(cls, *args, **options):
        assert len(args) == 1
        if isinstance(args[0], Transpose):
            return args[0].arg
        return Expr.__new__(cls, *args)

    @property
    def arg(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{}.T'.format(sstr(self.arg))

class MatSymbolicMul(MatrixSymbolicExpr, Mul):

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

    def __new__(cls, *args, **options):
        args = [sympify(a) for a in args if a != 0]
        newargs = []
        for i in args:
            if isinstance(i, MatSymbolicAdd):
                newargs += list(i.args)
            else:
                newargs.append(i)

        args = newargs

        if len(args) == 0:
            return S.Zero
        elif len(args) == 1:
            return args[0]

        return Expr.__new__(cls, *args)

    def _sympystr(self, printer):
        sstr = printer.doprint
        #return ' + '.join((sstr(i) for i in self.args))
        return 'MatSymbolicAdd({})'.format(','.join((sstr(i) for i in self.args)))

class MatSymbolicPow(MatrixSymbolicExpr, Pow):
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

def f1(*x):
    return MatSymbolicMul(*x)
def f2(*x):
    return MatSymbolicAdd(*x)

Basic._constructor_postprocessor_mapping[MatrixSymbolicExpr] = {
    "Mul": [f1],
    "Add": [f2]
}

