# coding: utf-8

# TODO - transpose of BilinearForm
#      - add unknown status if only one space is given to the BilinearForm
#        => can not be gelatized except if it is called on a couple test/trial
#      - check that a BilinearForm is bilinear (using Constant)
#      - check that a LinearForm is linear
#      - add is_symmetric property for BilinearForm
#      - treat Function atom in atomize, normalize, matricize
#      - treat Field atom in atomize, normalize, matricize
#      - shall we matricize a FunctionForm or not?

from numpy import zeros

from sympy.core import Basic
from sympy.core import Symbol
from sympy.core import Function
from sympy.core import Expr, Add, Mul, Pow
from sympy import S
from sympy.core.containers import Tuple
from sympy import preorder_traversal
from sympy import Indexed, IndexedBase, Matrix, ImmutableDenseMatrix
from sympy.physics.quantum import TensorProduct
from sympy import expand
from sympy import Integer, Float
from sympy.core.expr import AtomicExpr
from sympy.physics.quantum import TensorProduct

from .measure import CanonicalMeasure
from .measure import CartesianMeasure
from .measure import Measure
from .geometry import Line, Square, Cube
from .derivatives import _partial_derivatives
from .derivatives import partial_derivative_as_symbol
from .derivatives import sort_partial_derivatives
from .derivatives import get_atom_derivatives
from .derivatives import dx, dy, dz
from .derivatives import (Grad_1d, Div_1d,
                          Grad_2d, Curl_2d, Rot_2d, Div_2d,
                          Grad_3d, Curl_3d, Div_3d)

from .basic import _coeffs_registery
from .basic import CalculusFunction
from .basic import Field, Constant

from .generic import Dot, Inner, Cross
from .generic import Grad, Rot, Curl, Div
from .generic import _generic_ops

from .algebra import (Dot_1d,
                      Dot_2d, Inner_2d, Cross_2d,
                      Dot_3d, Inner_3d, Cross_3d)

from .space import BasicSobolevSpace
from .space import ProductSpace
from .space import TestFunction
from .space import VectorTestFunction
from spl.linalg.stencil import StencilVector

class BasicForm(Expr):

    def _init_domain(domain, ldim):
        if not( domain is None ):
            return domain

        if ldim == 1:
            return Line()

        elif ldim == 2:
            return Square()

        elif ldim == 3:
            return Cube()

        else:
            raise NotImplementedError('')

    def _init_measure(measure, coordinates):
        if not( measure is None ):
            return measure

        if not( coordinates is None ):
            return Measure(coordinates)

        else:
            raise NotImplementedError('')


    @property
    def fields(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Field)]
        # no redanduncy
        return sorted(list(set(ls)))

    @property
    def constants(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Constant)]
        # no redanduncy
        return list(set(ls))

    @property
    def domain(self):
        return self._domain

    @property
    def measure(self):
        return self._measure

    @property
    def mapping(self):
        return self._mapping


# TODO we should check that the only free symbols are fields, constants or coordinates
class FunctionForm(BasicForm):
    """

    Examples

    """
    _ldim = None
    _coordinates = None
    def __new__(cls, expr, coordinates=None, domain=None, measure=None, mapping=None):

        # ... check that there are no test functions in the expression
        ls = [a for a in expr.free_symbols if isinstance(a, (TestFunction, VectorTestFunction))]
        if not(len(ls) == 0):
            raise TypeError('Cannot use test functions in FunctionForm')
        # ...

        # ... compute dim from fields if available
        ls = [a for a in expr.free_symbols if isinstance(a, Field)]
        if ls:
            F = ls[0]
            ldim = F.space.ldim

            if coordinates is None:
                coordinates = F.space.coordinates

        else:
            if coordinates is None:
                raise ValueError('> Coordinates must be provided if the expression has no fields')

            ldim = len(coordinates)
            if ldim == 1:
                coordinates = coordinates[0]
        # ...

        domain = BasicForm._init_domain(domain, ldim)
        measure = BasicForm._init_measure(measure, coordinates)

        obj = Basic.__new__(cls, expr)
        obj._ldim = ldim
        obj._coordinates = coordinates
        obj._domain = domain
        obj._measure = measure
        obj._mapping = mapping

        return obj

    @property
    def expr(self):
        return self._args[0]

    @property
    def ldim(self):
        return self._ldim

    @property
    def coordinates(self):
        return self._coordinates

    def _sympystr(self, printer):
        sstr = printer.doprint
        expr = self.expr
        return sstr(expr)

    # TODO how to implement this?
    def __call__(self, *args):
        raise NotImplementedError('')


class LinearForm(BasicForm):
    """

    Examples

    """
    def __new__(cls, test_functions, expr, domain=None, measure=None, mapping=None):
        # ...
        if isinstance(test_functions, (TestFunction, VectorTestFunction)):
            test_functions = [test_functions]
            test_functions = Tuple(*test_functions)

        elif isinstance(test_functions, (tuple, list, Tuple)):
            spaces = [i.space for i in test_functions]
            V = ProductSpace(*spaces)
            name = ''.join(i.name for i in test_functions)
            v = VectorTestFunction(V, name=name)

            i = 0
            for w in test_functions:
                if isinstance(w, TestFunction):
                    expr = expr.subs({w: v[i]})

                elif isinstance(w, VectorTestFunction):
                    for j in range(0, w.shape[0]):
                        expr = expr.subs({w[j]: v[i+j]})

                i += w.space.shape

            test_functions = [v]

        else:
            raise TypeError('Wrong type for test function(s)')

        test_functions = Tuple(*test_functions)
        # ...

        obj = Basic.__new__(cls, expr, test_functions)

        domain = BasicForm._init_domain(domain, obj.ldim)
        measure = BasicForm._init_measure(measure, obj.coordinates)

        obj._domain = domain
        obj._measure = measure
        obj._mapping = mapping
        return obj

    @property
    def expr(self):
        return self._args[0]

    @property
    def test_functions(self):
        return self._args[1]

    @property
    def ldim(self):
        return self.test_spaces[0].ldim

    @property
    def coordinates(self):
        return self.test_spaces[0].coordinates

    @property
    def test_spaces(self):
        return [u.space for u in self.test_functions]

    def _sympystr(self, printer):
        sstr = printer.doprint
        expr = self.expr
        return sstr(expr)

    # TODO: when subs args there is a possible bug when a new var exists in the
    #       old args in a different order
    def __call__(self, *args):
        # ...
        if not(len(args) == 1):
            raise ValueError('Expecting one argument')
        # ...

        # ...
        test_functions = args[0]
        if isinstance(test_functions, (TestFunction, VectorTestFunction)):
            test_functions = [test_functions]

        elif isinstance(test_functions, (tuple, list, Tuple)):
            test_functions = list(*test_functions)
        # ...

        # ... we need to atomize the expression here, since Grad, Div etc are
        #     not defined yet.
        #     TODO: this may be improved, if we fuse 1d, 2d, 3d definitions of Generic
        #     operators into one single Function expression
        #     in this case, operators like Cross, Curl must be able to act on
        #     slices
        expr = atomize(self.expr)
        # ...

        # ... replacing test functions
        d = {}
        for k,v in zip(self.test_functions, test_functions):
            d[k] = v
        expr = expr.subs(d)
        # ...

        return expr


class BilinearForm(BasicForm):
    """

    Examples

    """
    def __new__(cls, test_trial, expr, domain=None, measure=None, mapping=None):
        # ...
        if not isinstance(test_trial, (tuple, list, Tuple)):
            raise TypeError('(test, trial) must be a tuple, list or Tuple')

        if not(len(test_trial) == 2):
            raise ValueError('Expecting a couple (test, trial)')

        # ...
        test_functions = test_trial[0]
        if isinstance(test_functions, (TestFunction, VectorTestFunction)):
            test_functions = [test_functions]

        elif isinstance(test_functions, (tuple, list, Tuple)):
            spaces = [i.space for i in test_functions]
            V = ProductSpace(*spaces)
            name = ''.join(i.name for i in test_functions)
            v = VectorTestFunction(V, name=name)

            i = 0
            for w in test_functions:
                if isinstance(w, TestFunction):
                    expr = expr.subs({w: v[i]})

                elif isinstance(w, VectorTestFunction):
                    for j in range(0, w.shape[0]):
                        expr = expr.subs({w[j]: v[i+j]})

                i += w.space.shape

            test_functions = [v]

        else:
            raise TypeError('Wrong type for test function(s)')

        test_functions = Tuple(*test_functions)
        # ...

        # ...
        trial_functions = test_trial[1]
        if isinstance(trial_functions, (TestFunction, VectorTestFunction)):
            trial_functions = [trial_functions]

        elif isinstance(trial_functions, (tuple, list, Tuple)):
            spaces = [i.space for i in trial_functions]
            V = ProductSpace(*spaces)
            name = ''.join(i.name for i in trial_functions)
            v = VectorTestFunction(V, name=name)

            i = 0
            for w in trial_functions:
                if isinstance(w, TestFunction):
                    expr = expr.subs({w: v[i]})

                elif isinstance(w, VectorTestFunction):
                    for j in range(0, w.shape[0]):
                        expr = expr.subs({w[j]: v[i+j]})

                i += w.space.shape

            trial_functions = [v]

        else:
            raise TypeError('Wrong type for trial function(s)')

        trial_functions = Tuple(*trial_functions)
        # ...

        obj = Basic.__new__(cls, expr, test_functions, trial_functions)

        domain = BasicForm._init_domain(domain, obj.ldim)
        measure = BasicForm._init_measure(measure, obj.coordinates)

        obj._domain = domain
        obj._measure = measure
        obj._mapping = mapping
        return obj

    @property
    def expr(self):
        return self._args[0]

    @property
    def test_functions(self):
        return self._args[1]

    @property
    def trial_functions(self):
        return self._args[2]

    @property
    def ldim(self):
        return self.test_spaces[0].ldim

    @property
    def coordinates(self):
        return self.test_spaces[0].coordinates

    @property
    def test_spaces(self):
        return [u.space for u in self.test_functions]

    @property
    def trial_spaces(self):
        return [u.space for u in self.trial_functions]

    def _sympystr(self, printer):
        sstr = printer.doprint
        expr = self.expr
        return sstr(expr)

    # TODO: when subs args there is a possible bug when a new var exists in the
    #       old args in a different order
    def __call__(self, *args):
        # ...
        test_trial = args
        if not isinstance(test_trial, (tuple, list, Tuple)):
            raise TypeError('(test, trial) must be a tuple, list or Tuple')

        if not(len(test_trial) == 2):
            raise ValueError('Expecting a couple (test, trial)')
        # ...

        # ...
        test_functions = test_trial[0]
        if isinstance(test_functions, (TestFunction, VectorTestFunction)):
            test_functions = [test_functions]

        elif isinstance(test_functions, (tuple, list, Tuple)):
            test_functions = list(*test_functions)
        # ...

        # ...
        trial_functions = test_trial[1]
        if isinstance(trial_functions, (TestFunction, VectorTestFunction)):
            trial_functions = [trial_functions]

        elif isinstance(trial_functions, (tuple, list, Tuple)):
            trial_functions = list(*trial_functions)
        # ...

        # ... we need to atomize the expression here, since Grad, Div etc are
        #     not defined yet.
        #     TODO: this may be improved, if we fuse 1d, 2d, 3d definitions of Generic
        #     operators into one single Function expression
        #     in this case, operators like Cross, Curl must be able to act on
        #     slices
        expr = atomize(self.expr)
        # ...

        # in order to avoid problems when swapping indices, we need to create
        # temp symbols

        # ...
        d_tmp = {}
        for x in trial_functions:
            name = '{name}_{hash}'.format(name=x.name, hash=abs(hash(x)))
            X = x.duplicate(name)

            d_tmp[X] = x
        # ...

        # ... replacing trial functions by tmp symbols
        d = {}
        for k,v in zip(self.trial_functions, d_tmp):
            d[k] = v
        expr = expr.subs(d)
        # ...

        # ... replacing test functions
        d = {}
        for k,v in zip(self.test_functions, test_functions):
            d[k] = v
        expr = expr.subs(d)
        # ...

        # ... replacing trial functions from tmp symbols
        expr = expr.subs(d_tmp)
        # ...

        return expr


class BilinearAtomicForm(BilinearForm, AtomicExpr):
    """

    Examples

    """
    _name = None

    @property
    def name(self):
        return self._name

    def _sympystr(self, printer):
        sstr = printer.doprint
        name = sstr(self.name)

        test = [sstr(i) for i in self.test_functions]
        test = ','.join(i for i in test)

        trial = [sstr(i) for i in self.trial_functions]
        trial = ','.join(i for i in trial)

        return '{name}({test},{trial})'.format(name=name, trial=trial, test=test)

class Mass(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'Mass'
    def __new__(cls, test, trial):

        test_trial = [test, trial]
        expr = test * trial

        return BilinearForm.__new__(cls, test_trial, expr)

class Stiffness(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'Stiffness'
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = d(test) * d(trial)

        return BilinearForm.__new__(cls, test_trial, expr)

class Advection(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'Advection'
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = test * d(trial)

        return BilinearForm.__new__(cls, test_trial, expr)

class AdvectionT(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'AdvectionT'
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = d(test) * trial

        return BilinearForm.__new__(cls, test_trial, expr)


class Kron(BilinearAtomicForm):

    """."""

    def __new__(cls, ls):
        
        ln = len(ls)
        dim = len(ls[0])
        obj = Basic.__new__(cls, ls)
        if dim >= 3:
            raise NotImplementedError('TODO')
         
        args1 = ['A%s'%i for i in range(ln)]
        args2 = ['B%s'%i for i in range(ln)]
        try:
            import importlib
            package = importlib.import_module("symfe.codegen.templates.kron")
        except:
            raise ImportError('could not import kron_dot')
        name = 'kron_dot_2d'
        template = getattr(package,name)
        body = [template['body'].format(MAT1=args1[i],MAT2=args2[i]) for i in range(ln)]
        body = '\n'.join(i for i in body)
        args = ','.join(arg for arg in args1+args2)
        function = template['function'].format(___MAT_ARGS___=args,__BODY__=body)
        args_types = ','.join('double[:,:]' for i in range(2*ln))
        header = template['header'].format(__ARGS_TYPES__=args_types)
        dot = compile(function,'','single')
        dic = {}
        eval(dot,dic)
        _dot = dic[name]
        setattr(obj, '_dot',_dot)
        return obj

    @property
    def args(self):
        return self._args[0]
 
    def dot(self, x):
        space = x.space
        args = list(zip(*self.args))
        args1 = args[0]
        args2 = args[1]
        args1 = [arg._data for arg in args1]
        args2 = [arg._data for arg in args2]
        #args1 = [arg._data.T for arg in args1]
        #args2 = [arg._data.T for arg in args2]
        starts = space.starts
        ends   = space.ends
        pads   = space.pads
        Y      = StencilVector(space)
        X_tmp  = StencilVector(space)
        #self._dot(starts,ends,pads,x._data.T,Y._data.T,X_tmp._data.T,*args1,*args2)
        self._dot(starts,ends,pads,x._data,Y._data,X_tmp._data,*args1,*args2)
        return Y
        
  
    def __str__(self):
        return 'Kron'
  
    def _sympystr(self, printer):
        return 'Kron'

# ...
def atomize(expr, dim=None):
    """
    """
#    if not isinstance(expr, (Add, Mul,
#                             _partial_derivatives, _generic_ops,
#                             TestFunction, VectorTestFunction, Indexed,
#                             Field, Constant, Symbol, Function,
#                             Integer, Float, Matrix, ImmutableDenseMatrix,
#                             list, tuple, Tuple)):
#        msg = ('> Wrong input type.')
#
#        raise TypeError(msg, ', given ', expr, type(expr))

    if not isinstance(expr, (Expr,
                             _partial_derivatives, _generic_ops,
                             TestFunction, VectorTestFunction, Indexed,
                             Field, Constant, Symbol, Function,
                             Integer, Float, Matrix, ImmutableDenseMatrix,
                             list, tuple, Tuple)):
        msg = ('> Wrong input type.')

        raise TypeError(msg, ', given ', expr, type(expr))

#    print('> expr [atomize] = ', expr, type(expr))

    # ... compute dim if None
    if dim is None:
        ls = [i for i in expr.free_symbols if isinstance(i, (TestFunction,
                                                             VectorTestFunction,
                                                             Field))]

        if ls:
            atom = ls[0]
            if atom.space is None:
                raise ValueError('Expecting atom to be associated to a space')

            dim = atom.space.ldim
    # ...

    if isinstance(expr, (list, tuple, Tuple)):
        args = [atomize(i, dim=dim) for i in expr]
        return Tuple(*args)

    elif isinstance(expr, Add):
        args = [atomize(i, dim=dim) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in expr.args if not(i in coeffs)]
#        print('expr    = ', expr)
#        print('vectors = ', vectors)
#        print('coeffs  = ', coeffs )

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [atomize(i, dim=dim) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)

    elif isinstance(expr, Pow):

        b = atomize(expr.base, dim=dim)
        e = expr.exp

        return Pow(b, e)

    elif isinstance(expr, (Dot, Inner, Cross, Grad, Rot, Curl, Div)):
        # if i = Dot(...) then type(i) is Grad
        op = type(expr)
        new  = eval('{0}_{1}d'.format(op, dim))

        args = [atomize(i, dim=dim) for i in expr.args]
        return new(*args)

    elif isinstance(expr, Matrix):
        n,m = expr.shape
        lines = []
        for i in range(0, n):
            line = []
            for j in range(0, m):
                line.append(atomize(expr[i,j], dim=dim))
            lines.append(line)
        return Matrix(lines)

    return expr
# ...

# ...
def normalize(expr, basis=None, enable_fields=False):
    """
    must be applied after calling atomize

    basis: dict
        for every space we give the name of the basis function symbol
    """
    # ... compute dim
    ls = [i for i in expr.free_symbols if isinstance(i, (TestFunction,
                                                         VectorTestFunction,
                                                         Field))]

    if ls:
        atom = ls[0]
        if atom.space is None:
            raise ValueError('Expecting atom to be associated to a space')

        dim = atom.space.ldim
    # ...

    # ...
    expr = atomize(expr)
    # ...

#    print('> expr [normalize] = ', expr, type(expr))

    if isinstance(expr, (list, tuple, Tuple)):
        args = [normalize(i, basis=basis, enable_fields=enable_fields) for i in expr]
        return Tuple(*args)

    elif isinstance(expr, Add):
        args = [normalize(i, basis=basis, enable_fields=enable_fields) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [normalize(i, basis=basis, enable_fields=enable_fields) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)

    elif isinstance(expr, Pow):

        b = normalize(expr.base, basis=basis, enable_fields=enable_fields)
        e = expr.exp

        return Pow(b, e)

    elif isinstance(expr, _partial_derivatives):
        ops = sort_partial_derivatives(expr)

        # ...
        for i in ops:

            if not(len(i.args) == 1):
                raise ValueError('expecting only one argument for partial derivatives')

            arg = i.args[0]

            name = None
            shape = dim
            if not(basis is None):
                atom = get_atom_derivatives(i)
                if isinstance(atom, (TestFunction, VectorTestFunction)):
                    if atom in list(basis.keys()):
                        name = basis[atom]

                    if isinstance(atom, VectorTestFunction):
                        shape = atom.shape[0]

                elif isinstance(atom, Indexed):
                    base = atom.base
                    shape = base.shape[0]
                    if base in list(basis.keys()):
                        name = basis[base]

            # terms like dx(..)
            if  enable_fields or not isinstance(arg, Field):
                new = partial_derivative_as_symbol(i, name=name, dim=shape)
                expr = expr.subs(i, new)
        # ...

        return expr

    elif isinstance(expr, (TestFunction, VectorTestFunction)):
        if not(basis is None):
            if expr in list(basis.keys()):
                name = basis[expr]
                return Symbol(name)

    elif isinstance(expr, Indexed):
        if not(basis is None):
            base = expr.base
            if base in list(basis.keys()):
                name = basis[base]
                indices = expr.indices
                shape = base.shape[0]
                return IndexedBase(name, shape=shape)[indices]

    return expr
# ...

# ...
def matricize(expr):
    """
    must be applied after calling normalize
    """

    # ... we need first to expand the expression
    expr = expand(expr)
    # ...

#    print('> expr = ', expr, type(expr))

    if isinstance(expr, (list, tuple, Tuple)):
        args = [matricize(i) for i in expr]
        return Tuple(*args)

    elif isinstance(expr, Add):
        args = [matricize(i) for i in expr.args]
        # we cannot return Add(*args)
        # since it gives the error:
        # TypeError: cannot add <class 'sympy.matrices.immutable.ImmutableDenseMatrix'> and <class 'sympy.core.numbers.Zero'>
        # when args are of type Matrix
        r = args[0]
        for a in args[1:]:
            r += a
        return r

    elif isinstance(expr, Mul):
        # a coeff can be a symbol, otherwise the expression rot(v) * rot(u) + c * div(v) * div(u)
        # raises an error
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery) or isinstance(i, Symbol)]
        others = [i for i in expr.args if not(i in coeffs)]
        for i in others:
            if isinstance(i, Function) and not(isinstance(i, CalculusFunction)):
                coeffs.append(i)
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [matricize(i) for i in vectors]

            if len(args) == 1:
                j = args[0]

            elif len(args) == 2:
                # TODO how to be sure about who is left/right? test/trial?
                left = args[0]
                right = args[1]

                if isinstance(left, Matrix) and isinstance(right, Matrix):
                    j = TensorProduct(left.transpose(), right)
                else:
                    j = Mul(*args)

            else:
                raise ValueError('Expecting one or two arguments')


        # we cannot return Mul(i, j)
        # since it gives the error:
        # TypeError: cannot add <class 'sympy.matrices.immutable.ImmutableDenseMatrix'> and <class 'sympy.core.mul.Mul'>
        # when args are of type Matrix
        return i * j

    elif isinstance(expr, Indexed):
        base = expr.base
        if not(len(expr.indices) == 1):
            raise ValueError('Expecting exactly one index')

        index = expr.indices[0]
        name = '{}'.format(base)

        dim = base.shape[0]

        M = Matrix(zeros(dim))
        M[index] = Symbol(name)

        return M

    return expr
# ...

# ... TODO compute basis if not given
def gelatize(a, basis=None, verbose=False):

    if not isinstance(a, (BasicForm, Add, Mul)):
        raise TypeError('Expecting a BasicForm, Add or Mul')

    if isinstance(a, Add):
        args = [gelatize(i, basis=basis, verbose=verbose) for i in a.args]
        return Add(*args)

    elif isinstance(a, Mul):
        # a coeff can be a symbol, otherwise the expression c1 * a
        # raises an error
        coeffs  = [i for i in a.args if isinstance(i, _coeffs_registery) or isinstance(i, Symbol)]
        vectors = [i for i in a.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [gelatize(i, basis=basis, verbose=verbose) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)

    dim = a.ldim
    expr = a.expr

    expr = atomize(expr, dim=dim)
    if verbose:
        print('> atomized   >>> {0}'.format(expr))

    expr = normalize(expr, basis=basis)
    if verbose:
        print('> normalized >>> {0}'.format(expr))

    # TODO is it ok to keep this?
    if isinstance(a, FunctionForm):
        return expr

    expr = matricize(expr)
    if verbose:
        print('> matricized >>> {0}'.format(expr))

    return expr
# ...

# TODO - get dim from atoms
#      - check coefficinets/functions
def _tensorize_core(expr, dim, tests, trials):

    
    

    if isinstance(expr, Add):
        args = [_tensorize_core(i, dim, tests, trials) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        args = [i for i in expr.args if not(i in coeffs)]

        d_atoms = {}
        _coordinates = ['x', 'y', 'z']
        test_trial = list(tests) + list(trials)
        for a in test_trial:
            d_atoms[a] = []

            new = S.One
            for i in range(0, dim):
                coord = _coordinates[i]
                Vi = BasicSobolevSpace('V_{}'.format(i),
                                       ldim=1,
                                       coordinates=[coord])

                ai = TestFunction(Vi, '{test}{i}'.format(test=a.name, i=i))
                d_atoms[a].append(ai)

                new *= ai
            expr = expr.subs({a: new})

        # make sure we have sum of products
        expr = expand(expr)

        # ...
        # TODO - improve this later
        #      - must distinguish between test/trial
        assert(len(tests) == 1)
        assert(len(trials) == 1)

        v = tests[0]
        u = trials[0]

        ops = {'x': dx, 'y': dy, 'z': dz}
        
        for ui,vi in zip(d_atoms[u], d_atoms[v]):
            coord = ui.space.coordinates.name
            d = ops[coord]

            # ... Mass
            old = vi*ui
            new = Mass(vi,ui)

            expr = expr.subs({old: new})
            # ...

            # ... Stiffness
            old = d(vi)*d(ui)
            new = Stiffness(vi,ui)

            expr = expr.subs({old: new})
            # ...

            # ... Advection
            old = vi*d(ui)
            new = Advection(vi,ui)

            expr = expr.subs({old: new})
            # ...

            # ... Transpose of Advection
            old = d(vi)*ui
            new = AdvectionT(vi,ui)

            expr = expr.subs({old: new})
            # ...
        
        
        expr = subs_mul(expr)
        # ...
    

    return expr

# ... TODO not working yet for high order derivatives
def inv_normalize(expr, ns):
    """Returns an expression using partial differential operators from an
    expression based on Symbols like u_x to define the derivative of u w.r.t x.
    a Space for scalar functions will be created.

    ns: dict
        similar to the one in subs method of Expr
    """
    if not isinstance(ns, dict):
        raise TypeError('> Expecting a dictionary')

    # ...
    names = []
    tests = []
    for i,u in list(ns.items()):
        if isinstance(i, Symbol):
            n = str(i.name)

        elif isinstance(i, str):
            n = i

        else:
            raise TypeError('> Found a wrong type for test function name')

        names.append(n)
        tests.append(u)
    # ...

    # ... replace symbols by test functions
    for n,u in zip(names, tests):
        expr = expr.subs({n: u})
    # ...

    # ... all order derivatives
    d_ops = {'x': dx, 'y': dy, 'z': dz}

    def _op_from_str(u, suffixes):
        atom = u
        for coord in suffixes[::-1]:
            d = d_ops[coord]
            atom = d(atom)
        return atom

    free_symbols = expr.free_symbols
    for n,u in zip(names, tests):
        suffixes = [i.name.split('_')[-1] for i in free_symbols if i.name.startswith('{}_'.format(n))]
        for i in suffixes:
            new = _op_from_str(u, i)
            old = '{name}_{suffix}'.format(name=n, suffix=i)
            expr = expr.subs({old: new})
    # ...

    return expr

# ...

def tensorize(a):

    if not isinstance(a, BilinearForm):
        raise TypeError('Expecting a BilinearForm')

    dim = a.ldim
    expr = a.expr
    tests = a.test_functions
    trials = a.trial_functions

    # ... TODO this is not the best thing to do
    #          we should modify matricize and call it here
    e = gelatize(a)
    if isinstance(e, (Matrix, ImmutableDenseMatrix)):

        assert(len(a.test_spaces) == 1)
        assert(len(a.trial_spaces) == 1)

        assert(len(tests) == 1)
        assert(len(trials) == 1)

        V = a.test_spaces[0]
        U = a.trial_spaces[0]

        v = tests[0]
        u = trials[0]

        from symfe.core import H1Space
        V = H1Space(V.name, ldim=V.ldim, coordinates=[i.name for i in V.coordinates])
        U = H1Space(U.name, ldim=U.ldim, coordinates=[i.name for i in U.coordinates])

        vv = TestFunction(V, name=v.name*2)
        uu = TestFunction(U, name=u.name*2)

        tests = [TestFunction(V, name=i.name*2) for i in tests]
        trials = [TestFunction(U, name=i.name*2) for i in trials]

        print(e)
        expr = inv_normalize(e, {v: vv, u: uu})
        print(expr)
        n_rows, n_cols = expr.shape
        lines = []
        for i in range(0, n_rows):
            line = []
            for j in range(0, n_cols):
                eij = _tensorize_core(expr[i,j], dim, tests, trials)
                line.append(eij)
            lines.append(line)
        expr = Matrix(lines)

    else:

        expr = atomize(expr)
        expr = _tensorize_core(expr, dim, tests, trials)
    # ...

    return expr

def subs_mul(expr):
    """substitute Mul with TensorProduct"""

    if isinstance(expr,(Add, Mul)):
        args = expr.args
        args = [subs_mul(arg) for arg in args]

        if isinstance(expr, Mul):

            return TensorProduct(*args)
        elif isinstance(expr, Add):

            return Add(*args)
    else:

        return expr
 

