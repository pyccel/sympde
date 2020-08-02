from sympy                 import Expr, Abs, S
from sympy                 import Add, Mul, Pow
from sympy.core.expr       import AtomicExpr
from sympy.core.decorators import call_highest_priority
from sympy import Basic


class MatrixSymbolicExpr(Expr):
    is_commutative = False
    _op_priority   = 20.0
    is_MatrixSymbolicExpr = True

    def __new__(cls, *args):
        return Expr.__new__(cls, *args)

    def inv(self):
        return Inverse(self)

    def transpose(self):
        return Transpose(self)

    @property
    def T(self):
        return self.transpose()

    def det(self):
        return SymbolicDeterminant(self)

    # The following is adapted from the core Expr object
    def __neg__(self):
        return MatMul(S.NegativeOne, self)

    def __abs__(self):
        return MatAbs(self)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        return MatAdd(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return MatAdd(other, self)

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return MatAdd(self, -other)

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return MatAdd(other, -self)

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return MatMul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return MatMul(other, self)

    @call_highest_priority('__pow__')
    def __pow__(self, other):
        return MatPow(self, other)

    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return self * other**S.NegativeOne

    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        return MatMul(other, Pow(self, S.NegativeOne))

class Inverse(MatrixSymbolicExpr):
    is_commutative = False
    def __new__(cls, arg):
        if isinstance(arg, Inverse):
            return arg.arg
        return Expr.__new__(cls, arg)

    @property
    def arg(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        return 'Inverse({})'.format(sstr(self.arg))
        
class Transpose(MatrixSymbolicExpr):
    is_commutative = False
    def __new__(cls, arg):
        if isinstance(arg, Transpose):
            return arg.arg
        return Expr.__new__(cls, arg)

    @property
    def arg(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        return 'Transpose({})'.format(sstr(self.arg))

class MatMul(MatrixSymbolicExpr, Mul):
    is_MatMul = True
    def __new__(cls, *args):
        if len(args) == 0:
            return S.One
        elif len(args) == 1:
            return args[0]
        newargs = []
        for a in args:
            if isinstance(a, Mul):
                newargs += list(a.args)
            else:
                newargs.append(a)
        mats   = [a for a in newargs if not a.is_commutative]
        coeffs = [a for a in newargs if a.is_commutative]
        newargs = coeffs + mats
        return Expr.__new__(cls, *newargs)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return 'MatMul({})'.format(','.join((sstr(i) for i in self.args)))

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
        left = MatMul._expandsums(sums[:L//2])
        right = MatMul._expandsums(sums[L//2:])
        
        terms = [MatMul(a, b) for a in left for b in right]
        added = MatAdd(*terms)
        return Add.make_args(added) # it may have collapsed down to one term


class MatAdd(MatrixSymbolicExpr, Add):
    is_MatAdd = True
    def __new__(cls, *args):
        if len(args) == 0:
            return S.Zero
        elif len(args) == 1:
            return args[0]

        return Expr.__new__(cls, *args)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return 'MatAdd({})'.format(','.join((sstr(i) for i in self.args)))

class MatPow(MatrixSymbolicExpr, Pow):
    def _sympystr(self, printer):
        sstr = printer.doprint
        return 'MatPow({})'.format(','.join((sstr(i) for i in self.args)))

class MatAbs(MatrixSymbolicExpr, Mul):
    def _sympystr(self, printer):
        sstr = printer.doprint
        return 'Abs({})'.format(sstr(self.args[0]))

class SymbolicDeterminant(Expr):

    def __new__(cls, arg):
        assert(isinstance(arg, MatrixSymbolicExpr))
        return Basic.__new__(cls, arg)

    @property
    def arg(self):
        return self._args[0]       
       
    def _sympystr(self, printer):
        sstr = printer.doprint
        arg = sstr(self.arg)
        return 'det({})'.format(arg)


