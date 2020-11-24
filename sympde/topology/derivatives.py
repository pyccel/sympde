# coding: utf-8

# TODO add action of diff operators on sympy known functions

import numpy as np
from itertools   import groupby
from collections import OrderedDict

from sympy import Basic
from sympy import Symbol
from sympy import Expr
from sympy import Tuple
from sympy import Matrix, ImmutableDenseMatrix
from sympy import Mul, Add, Pow
from sympy import Derivative
from sympy import S
from sympy import Indexed, IndexedBase
from sympy import diff
from sympy import log
from sympy import preorder_traversal
from sympy import cacheit

from sympy.core.function      import AppliedUndef
from sympy.core.function      import UndefinedFunction
from sympy.core.compatibility import is_sequence

from sympde.core.basic    import CalculusFunction
from sympde.core.basic    import _coeffs_registery
from sympde.core.basic    import BasicMapping
from sympde.core.algebra  import LinearOperator
from sympde.calculus.core import minus, plus
from sympde.calculus.core import has

from .space   import ScalarTestFunction, VectorTestFunction, IndexedTestTrial
from .space   import ScalarField, VectorField, IndexedVectorField

#==============================================================================
class DifferentialOperator(LinearOperator):
    """
    This class is a linear operator that applies the Leibniz formula

    Parameters
    ----------
    expr : Expr
        expr represents a Sympy expression

    """
    coordinate = None
    logical    = None

    # this is a hack since sometimes we use Matrix in the eval of abstract
    # operators
    def __getitem__(self, indices, **kw_args):
        return self

    @classmethod
    @cacheit
    def eval(cls, expr):

        types = (VectorTestFunction, ScalarTestFunction,
                DifferentialOperator, ScalarField, VectorField)

        if isinstance(expr, _logical_partial_derivatives):
            atom    = get_atom_logical_derivatives(expr)
            indices = get_index_logical_derivatives(expr)

            if cls in _logical_partial_derivatives:
                indices[cls.coordinate] += 1

            for i,n in enumerate(list(indices.values())[::-1]):
                d = _logical_partial_derivatives[-i-1]
                for _ in range(n):
                    atom = d(atom, evaluate=False)

            if cls in _partial_derivatives:
                atom = cls(atom, evaluate=False)
            elif cls not in _logical_partial_derivatives:
                raise NotImplementedError('TODO')

            return atom

        if isinstance(expr, (VectorTestFunction, VectorField)):
            n = expr.shape[0]
            args = [cls(expr[i], evaluate=False) for i in range(0, n)]
            args = Tuple(*args)
            return Matrix([args])
        elif isinstance(expr, (list, tuple, Tuple, Matrix, ImmutableDenseMatrix)):
            args = [cls(i, evaluate=True) for i in expr]
            args = Tuple(*args)
            return Matrix([args])
        elif isinstance(expr, (IndexedTestTrial, IndexedVectorField, DifferentialOperator)):
            return cls(expr, evaluate=False)

        elif isinstance(expr, (ScalarField, ScalarTestFunction)):
            return cls(expr, evaluate=False)

        elif isinstance(expr, (minus, plus)):
            return cls(expr, evaluate=False)

        elif isinstance(expr, Indexed) and isinstance(expr.base, BasicMapping):
            return cls(expr, evaluate=False)
        elif not has(expr, types):
            if expr.is_number:
                return S.Zero

            elif isinstance(expr, Expr):
                x = Symbol(cls.coordinate)
                if cls.logical:
                    M = expr.atoms(Mapping)
                    if len(M)>0:
                        M = list(M)[0]
                        expr_primes = [diff(expr, M[i]) for i in range(M.pdim)]
                        Jj = Jacobian(M)[:,cls.grad_index]
                        expr_prime = sum([ei*Jji for ei,Jji in zip(expr_primes, Jj)])
                        return expr_prime + diff(expr, x)
                return diff(expr, x)


        if isinstance(expr, Add):
            args = [cls(a, evaluate=True) for a in expr.args]
            v    = Add(*args)
            return v

        elif isinstance(expr, Mul):
            coeffs  = [a for a in expr.args if isinstance(a, _coeffs_registery)]
            vectors = [a for a in expr.args if not(a in coeffs)]

            c = S.One
            if coeffs:
                c = Mul(*coeffs)

            V = S.Zero
            if vectors:
                if len(vectors) == 1:
                    # do we need to use Mul?
                    V = cls(vectors[0], evaluate=True)

                elif len(vectors) == 2:
                    a = vectors[0]
                    b = vectors[1]

                    fa = cls(a, evaluate=True)
                    fb = cls(b, evaluate=True)

                    V = a * fb + fa * b

                else:
                    a = vectors[0]
                    b = Mul(*vectors[1:])

                    fa = cls(a, evaluate=True)
                    fb = cls(b, evaluate=True)

                    V = a * fb + fa * b

            v = Mul(c, V)
            return v

        elif isinstance(expr, Pow):
            b = expr.base
            e = expr.exp
            v = (log(b)*cls(e, evaluate=True) + e*cls(b, evaluate=True)/b) * b**e
            return v

        else:
            msg = '{expr} of type {type}'.format(expr=expr, type=type(expr))
            raise NotImplementedError(msg)

    def _sympystr(self, printer):
        sstr = printer.doprint
        args = ','.join(sstr(i) for i in self.args)
        return 'd{}({})'.format(sstr(self.coordinate), args)
#==============================================================================
class dx(DifferentialOperator):
    coordinate = 'x'
    grad_index = 0 # index in grad


class dy(DifferentialOperator):
    coordinate = 'y'
    grad_index = 1 # index in grad


class dz(DifferentialOperator):
    coordinate = 'z'
    grad_index = 2 # index in grad


_partial_derivatives = (dx, dy, dz)

#==============================================================================
class dx1(DifferentialOperator):
    coordinate = 'x1'
    grad_index = 0 # index in grad
    logical    = True

class dx2(DifferentialOperator):
    coordinate = 'x2'
    grad_index = 1 # index in grad
    logical    = True

class dx3(DifferentialOperator):
    coordinate = 'x3'
    grad_index = 2 # index in grad
    logical    = True

_logical_partial_derivatives = (dx1, dx2, dx3)

#==============================================================================
@cacheit
def find_partial_derivatives(expr):
    """
    returns all partial derivative expressions
    """
    if isinstance(expr, (Add, Mul)):
        return find_partial_derivatives(expr.args)

    elif isinstance(expr, Pow):
        return find_partial_derivatives(expr.base)

    elif isinstance(expr, (list, tuple, Tuple)):
        args = []
        for a in expr:
            args += find_partial_derivatives(a)
        return args

    elif isinstance(expr, _partial_derivatives):
        return (expr,)

    elif isinstance(expr, _logical_partial_derivatives):
        return (expr,)

    return ()

#==============================================================================
@cacheit
def get_number_derivatives(expr):
    """
    returns the number of partial derivatives in expr.
    this is still an experimental version, and it assumes that expr is of the
    form d(a) where a is a single atom.
    """
    n = 0
    if isinstance(expr, _partial_derivatives):
        assert(len(expr.args) == 1)

        n += 1 + get_number_derivatives(expr.args[0])
    return n

#==============================================================================
@cacheit
def sort_partial_derivatives(expr):
    """returns the partial derivatives of an expression, sorted.
    """
    ls = []

    args = find_partial_derivatives(expr)

#    # ... Note
#    #     group by is given the wrong answer for expr =mu * u + dx(u) + dx(dx(u))
#    for key, group in groupby(args, lambda x: get_number_derivatives(x)):
#        g = [a for a in group]
#        for a in group:
#            ls.append(a)
#    # ...

    # ...
    d = {}
    for a in args:
        n = get_number_derivatives(a)
        if n in d.keys():
            d[n] += [a]
        else:
            d[n] = [a]
    # ...

    # ...
    if not d:
        return []
    # ...

    # ... sort keys from high to low
    keys = np.asarray(list(d.keys()))
    keys.sort()
    keys = keys[::-1]
    # ...

    # ... construct a list of partial derivatives from high to low order
    ls = []
    for k in keys:
        ls += d[k]
    # ...

    return tuple(ls)

#==============================================================================
@cacheit
def get_index_derivatives(expr):
    """
    """
    coord = ['x','y','z']

    d = OrderedDict()
    for c in coord:
        d[c] = 0

    ops = [a for a in preorder_traversal(expr) if isinstance(a, _partial_derivatives)]
    for i in ops:
        op = type(i)

        if isinstance(i, dx):
            d['x'] += 1

        elif isinstance(i, dy):
            d['y'] += 1

        elif isinstance(i, dz):
            d['z'] += 1

    return d

#==============================================================================
@cacheit
def get_atom_derivatives(expr):
    """
    """

    if isinstance(expr, _partial_derivatives):
        assert(len(expr.args) == 1)

        return get_atom_derivatives(expr.args[0])

    else:
        return expr

#==============================================================================
def get_index_logical_derivatives(expr):
    """
    """
    coord = ['x1','x2','x3']

    d = OrderedDict()
    for c in coord:
        d[c] = 0

    ops = [a for a in preorder_traversal(expr) if isinstance(a, _logical_partial_derivatives)]
    for i in ops:
        op = type(i)

        if isinstance(i, dx1):
            d['x1'] += 1

        elif isinstance(i, dx2):
            d['x2'] += 1

        elif isinstance(i, dx3):
            d['x3'] += 1

    return d

#==============================================================================
@cacheit
def get_atom_logical_derivatives(expr):
    """
    """

    if isinstance(expr, _logical_partial_derivatives):
        assert(len(expr.args) == 1)

        return get_atom_logical_derivatives(expr.args[0])

    else:
        return expr

#==============================================================================
class DotBasic(CalculusFunction):

    nargs = None
    name = 'Dot'

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

class Dot_1d(DotBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        if not( len(_args) == 2):
            raise ValueError('Expecting two arguments')

        u = _args[0]
        v = _args[1]

        return u * v

class Dot_2d(DotBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        if not( len(_args) == 2):
            raise ValueError('Expecting two arguments')

        u = _args[0]
        v = _args[1]

        if isinstance(u, (Matrix, ImmutableDenseMatrix)):
            if isinstance(v, (Matrix, ImmutableDenseMatrix)):
                raise NotImplementedError('TODO')

            else:
                return Tuple(u[0,0]*v[0] + u[0,1]*v[1], u[1,0]*v[0] + u[1,1]*v[1])

        else:
            if isinstance(v, (Matrix, ImmutableDenseMatrix)):
                raise NotImplementedError('TODO')

            else:
                return u[0]*v[0] + u[1]*v[1]

class Dot_3d(DotBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        if not( len(_args) == 2):
            raise ValueError('Expecting two arguments')

        u = _args[0]
        v = _args[1]

        if isinstance(u, (Matrix, ImmutableDenseMatrix)):
            if isinstance(v, (Matrix, ImmutableDenseMatrix)):
                raise NotImplementedError('TODO')

            else:
                return Tuple(u[0,0]*v[0] + u[0,1]*v[1] + u[0,2]*v[2],
                             u[1,0]*v[0] + u[1,1]*v[1] + u[1,2]*v[2],
                             u[2,0]*v[0] + u[2,1]*v[1] + u[2,2]*v[2])

        else:
            if isinstance(v, (Matrix, ImmutableDenseMatrix)):
                raise NotImplementedError('TODO')

            else:
                return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

#==============================================================================
class CrossBasic(CalculusFunction):

    nargs = None
    name = 'Cross'

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

class Cross_2d(CrossBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        v = _args[1]

        return u[0]*v[1] - u[1]*v[0]

class Cross_3d(CrossBasic):

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        v = _args[1]

        return Tuple(u[1]*v[2] - u[2]*v[1],
                     u[2]*v[0] - u[0]*v[2],
                     u[0]*v[1] - u[1]*v[0])

#==============================================================================
class GradBasic(CalculusFunction):

    nargs = None
    name = 'Grad'

    @cacheit
    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

class Grad_1d(GradBasic):

    @classmethod
    @cacheit
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return dx(u)

class Grad_2d(GradBasic):

    @classmethod
    @cacheit
    def eval(cls, *_args):

        if not _args:
            return

        u  = _args[0]
        du = (dx(u), dy(u))

        if isinstance(du[0], (Tuple, Matrix, ImmutableDenseMatrix)):
            lines = [list(d[:]) for d in du]
        else:
            lines = [[d] for d in du]

        v = ImmutableDenseMatrix(lines)
        return v

class Grad_3d(GradBasic):

    @classmethod
    @cacheit
    def eval(cls, *_args):

        if not _args:
            return

        u  = _args[0]
        du = (dx(u), dy(u), dz(u))

        if isinstance(du[0], (Tuple, Matrix, ImmutableDenseMatrix)):
            lines = [list(d[:]) for d in du]
        else:
            lines = [[d] for d in du]

        v = ImmutableDenseMatrix(lines)

        return v

#==============================================================================
class CurlBasic(CalculusFunction):

    nargs = None
    name = 'Curl'

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

class Curl_2d(CurlBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return dx(u[1])-dy(u[0])

class Curl_3d(CurlBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return ImmutableDenseMatrix([[dy(u[2]) - dz(u[1])],
                                     [dz(u[0]) - dx(u[2])],
                                     [dx(u[1]) - dy(u[0])]])

#==============================================================================
class Rot_2d(CalculusFunction):

    nargs = None
    name = 'Grad'

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return ImmutableDenseMatrix([[dy(u)],[-dx(u)]])

#==============================================================================
class DivBasic(CalculusFunction):

    nargs = None
    name = 'Div'

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

class Div_1d(DivBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        return dx(u[0])

class Div_2d(DivBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return dx(u[0]) + dy(u[1])

class Div_3d(DivBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return dx(u[0]) + dy(u[1]) + dz(u[2])

#==============================================================================
class LaplaceBasic(CalculusFunction):

    nargs = None
    name = 'Laplace'

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

class Laplace_1d(LaplaceBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return dx(dx(u))

class Laplace_2d(LaplaceBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return dx(dx(u)) + dy(dy(u))

class Laplace_3d(LaplaceBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return dx(dx(u)) + dy(dy(u)) + dz(dz(u))

#==============================================================================
class HessianBasic(CalculusFunction):

    nargs = None
    name = 'Hessian'

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

class Hessian_1d(HessianBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return dx(dx(u))

class Hessian_2d(HessianBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return ImmutableDenseMatrix([[dx(dx(u)), dx(dy(u))],
                       [dx(dy(u)), dy(dy(u))]])

class Hessian_3d(HessianBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return ImmutableDenseMatrix([[dx(dx(u)), dx(dy(u)), dx(dz(u))],
                       [dx(dy(u)), dy(dy(u)), dy(dz(u))],
                       [dx(dz(u)), dy(dz(u)), dz(dz(u))]])

#==============================================================================
class BracketBasic(CalculusFunction):
    pass

class Bracket_2d(BracketBasic):

    nargs = None
    name = 'Bracket'

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        v = _args[1]

        return dx(u)*dy(v) - dy(u)*dx(v)

#==============================================================================
class LogicalGrad_1d(GradBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return dx1(u)

class LogicalGrad_2d(GradBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u  = _args[0]
        du = (dx1(u), dx2(u))

        if isinstance(du[0], (Tuple, Matrix, ImmutableDenseMatrix)):
            lines = [list(d[:]) for d in du]
        else:
            lines = [[d] for d in du]

        v = ImmutableDenseMatrix(lines)
        return v

class LogicalGrad_3d(GradBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u  = _args[0]
        du = (dx1(u), dx2(u), dx3(u))

        if isinstance(du[0], (Tuple, Matrix, ImmutableDenseMatrix)):
            lines = [list(d[:]) for d in du]
        else:
            lines = [[d] for d in du]

        v = ImmutableDenseMatrix(lines)
        return v

#==============================================================================
class LogicalCurl_2d(CurlBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return dx1(u[1])-dx2(u[0])

class LogicalCurl_3d(CurlBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return ImmutableDenseMatrix([[dx2(u[2]) - dx3(u[1])],
                     [dx3(u[0]) - dx1(u[2])],
                     [dx1(u[1]) - dx2(u[0])]])

#==============================================================================
class LogicalRot_2d(CalculusFunction):

    nargs = None
    name = 'Grad'

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return ImmutableDenseMatrix([[dx2(u)],[-dx1(u)]])

#==============================================================================
class LogicalDiv_1d(DivBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return dx1(u[0])

class LogicalDiv_2d(DivBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return dx1(u[0]) + dx2(u[1])

class LogicalDiv_3d(DivBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]

        return dx1(u[0]) + dx2(u[1]) + dx3(u[2])

#==============================================================================
class LogicalLaplace_1d(LaplaceBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return dx1(dx1(u))

class LogicalLaplace_2d(LaplaceBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return dx1(dx1(u)) + dx2(dx2(u))

class LogicalLaplace_3d(LaplaceBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return dx1(dx1(u)) + dx2(dx2(u)) + dx3(dx3(u))

#==============================================================================
class LogicalHessian_1d(HessianBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return dx1(dx1(u))

class LogicalHessian_2d(HessianBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return ImmutableDenseMatrix([[dx1(dx1(u)), dx1(dx2(u))],
                                     [dx1(dx2(u)), dx2(dx2(u))]])

class LogicalHessian_3d(HessianBasic):

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        if isinstance(u, (VectorTestFunction, VectorField)):
            raise NotImplementedError('TODO')

        return ImmutableDenseMatrix([[dx1(dx1(u)), dx1(dx2(u)), dx1(dx3(u))],
                                     [dx1(dx2(u)), dx2(dx2(u)), dx2(dx3(u))],
                                     [dx1(dx3(u)), dx2(dx3(u)), dx3(dx3(u))]])

#==============================================================================

class LogicalBracket_2d(BracketBasic):

    nargs = None
    name = 'Bracket'

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    @classmethod
    def eval(cls, *_args):

        if not _args:
            return

        u = _args[0]
        v = _args[1]

        return dx1(u)*dx2(v) - dx2(u)*dx1(v)

#==============================================================================
def get_index_derivatives_atom(expr, atom, verbose=False):
    """This function return a dictionary of partial derivative indices for
    a given atom.
    it must be called after atomizing the expression.
    """
    ops = sort_partial_derivatives(expr)
    if verbose:
        print('> ops = ', ops)

    indices = []
    for i in ops:
        a = get_atom_derivatives(i)
        if a == atom:
            index = get_index_derivatives(i)
            indices.append(index)

    return tuple(indices)

#==============================================================================
def get_index_logical_derivatives_atom(expr, atom, verbose=False):
    """This function return a dictionary of partial derivative indices for
    a given atom.
    it must be called after atomizing the expression.
    """
    ops = sort_partial_derivatives(expr)
    if verbose:
        print('> ops = ', ops)

    indices = []
    for i in ops:
        a = get_atom_logical_derivatives(i)
        if a == atom:
            index = get_index_logical_derivatives(i)
            indices.append(index)

    return tuple(indices)

#==============================================================================
def get_max_partial_derivatives(expr, F=None):
    if F is None:
        Fs = (list(expr.atoms(ScalarTestFunction)) +
              list(expr.atoms(VectorTestFunction)) +
              list(expr.atoms(IndexedTestTrial)) +
              list(expr.atoms(VectorField)) +
              list(expr.atoms(IndexedVectorField)) +
              list(expr.atoms(ScalarField)))

        indices = []
        for F in Fs:
            indices += get_index_derivatives_atom(expr, F)
    else:
        indices = get_index_derivatives_atom(expr, F)

    d = {'x':0, 'y':0, 'z':0}
    for dd in indices:
        for k,v in dd.items():
            if v > d[k]: d[k] = v
    return d

#==============================================================================
def get_max_logical_partial_derivatives(expr, F=None):
    if F is None:
        Fs = (list(expr.atoms(ScalarTestFunction)) +
              list(expr.atoms(VectorTestFunction)) +
              list(expr.atoms(IndexedTestTrial)) +
              list(expr.atoms(VectorField)) +
              list(expr.atoms(IndexedVectorField)) +
              list(expr.atoms(ScalarField)))

        indices = []
        for F in Fs:
            indices += get_index_logical_derivatives_atom(expr, F)
    else:
        indices = get_index_logical_derivatives_atom(expr, F)

    d = {'x1':0, 'x2':0, 'x3':0}
    for dd in indices:
        for k,v in dd.items():
            if v > d[k]: d[k] = v
    return d

from .mapping import Mapping, Jacobian

