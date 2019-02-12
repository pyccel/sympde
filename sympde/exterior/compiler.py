# coding: utf-8

# TODO - assert the space type is not Undefined

from numpy import unique
from collections import OrderedDict

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
from sympde.core.basic import CalculusFunction
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import H1SpaceType, HcurlSpaceType, HdivSpaceType
from sympde.topology import L2SpaceType, UndefinedSpaceType
from sympde.topology import TestFunction, ScalarTestFunction, VectorTestFunction
from sympde.topology import Field, ScalarField, VectorField
from sympde.topology.space import _is_sympde_atom
from sympde.topology.space import _is_test_function
from sympde.calculus.core import _is_op_test_function
from sympde.calculus import Grad, Curl, Div
from sympde.calculus import Dot, Inner, Cross
#from sympde.calculus import grad, dot, inner, cross, rot, curl, div

from .form import DifferentialForm
from .calculus import d, wedge, ip
from .calculus import delta, jp, hodge
from .calculus import AdjointExteriorDerivative



#==============================================================================
class ExteriorCalculusExpr(CalculusFunction):

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args, **options)
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
    def eval(cls, *_args, **kwargs):
        """."""

        if not _args:
            return

        if not len(_args) == 1:
            raise ValueError('Expecting one argument')

        expr = _args[0]
        tests = kwargs.pop('tests', [])
        atoms = kwargs.pop('atoms', {})

        assert(isinstance(atoms, (dict, OrderedDict)))

        if isinstance(expr, (ScalarField, VectorField)):
            return expr

        if isinstance(expr, (tuple, list, Tuple)):
            return expr

        if isinstance(expr, _coeffs_registery):
            return expr

        if _is_test_function(expr):
            name = expr.name
            dim  = expr.space.ldim
            kind = expr.space.kind

            if expr in atoms.keys():
                return atoms[expr]

            elif isinstance(kind, H1SpaceType):
                if not(expr in tests):
                    return DifferentialForm(name, index=0, dim=dim)

                else:
                    return DifferentialForm(name, index=3, dim=dim)

            elif isinstance(kind, HcurlSpaceType):
                if not(expr in tests):
                    return DifferentialForm(name, index=1, dim=dim)

                else:
                    return DifferentialForm(name, index=2, dim=dim)

            elif isinstance(kind, HdivSpaceType):
                if not(expr in tests):
                    return DifferentialForm(name, index=2, dim=dim)

                else:
                    return DifferentialForm(name, index=1, dim=dim)

            elif isinstance(kind, L2SpaceType):
                if not(expr in tests):
                    return DifferentialForm(name, index=3, dim=dim)

                else:
                    return DifferentialForm(name, index=0, dim=dim)

            else:
                msg = 'Cannot convert {} to differential form'.format(expr)
                raise TypeError(msg)

        if isinstance(expr, Grad):
            arg = expr._args[0]
            newarg = cls.eval(arg, tests=tests, atoms=atoms)

            if not(arg in tests):
                return d(newarg)

            else:
                return -delta(newarg)

        if isinstance(expr, Curl):
            arg = expr._args[0]
            newarg = cls.eval(arg, tests=tests, atoms=atoms)

            if not(arg in tests):
                return d(newarg)

            else:
                return delta(newarg)

        if isinstance(expr, Div):
            arg = expr._args[0]
            newarg = cls.eval(arg, tests=tests, atoms=atoms)

            if not(arg in tests):
                return d(newarg)

            else:
                return -delta(newarg)

        if isinstance(expr, Dot):
            left, right = expr._args[:]

            newleft = cls.eval(left, tests=tests, atoms=atoms)
            newright = cls.eval(right, tests=tests, atoms=atoms)

            if _is_test_function(left) and _is_test_function(right):
                if left in tests:
                    return wedge(newright, hodge(newleft))

                elif right in tests:
                    return wedge(newleft, hodge(newright))

                else:
                    raise ValueError('argument not appears as a test function')

            elif _is_op_test_function(left) and _is_test_function(right):
                if left._args[0] in tests:
                    return wedge(newright, hodge(newleft))

                elif right in tests:
                    return wedge(newleft, hodge(newright))

                else:
                    raise ValueError('argument not appears as a test function')

            elif _is_test_function(left) and _is_op_test_function(right):
                return wedge(newleft, hodge(newright))

            elif _is_test_function(right):

                if not(right in tests):
                    return ip(newleft, newright)

                else:
                    return jp(newleft, newright)

            elif _is_test_function(left):
                raise NotImplementedError('')

            else:
                raise NotImplementedError('')

        if isinstance(expr, Cross):
            # TODO ORDER OF LEFT AND RIGHT DEPEND ON THE STR!!
            #      TO BE IMPROVED
            left, right = expr._args[:]

            if _is_test_function(right):
                dim = right.space.ldim

                if not(right in tests):
                    right = DifferentialForm(right.name, index=2, dim=dim)
                    return ip(left, right)

                else:
                    right = DifferentialForm(right.name, index=1, dim=dim)
                    return -jp(left, right)

            else:
                raise NotImplementedError('')

        if isinstance(expr, Inner):
            left, right = expr._args[:]

            if isinstance(left, Grad) and isinstance(right, Grad):
                if (_is_test_function(left._args[0]) and
                    _is_test_function(right._args[0])):

                    left  = left._args[0]
                    right = right._args[0]

                    newleft = cls.eval(left, tests=tests, atoms=atoms)
                    newright = cls.eval(right, tests=tests, atoms=atoms)

                    if left in tests:
                        expr  = wedge(d(newright), hodge(d(newleft)))
                        expr += wedge(d(delta(newright)), hodge(newleft))

                    elif right in tests:
                        expr  = wedge(d(newleft), hodge(d(newright)))
                        expr += wedge(d(delta(newleft)), hodge(newright))

                    else:
                        raise ValueError('either left/right must be a test function')

                    return expr

                else:
                    raise NotImplementedError('')

            else:
                raise NotImplementedError('')

        if isinstance(expr, Add):
            args = [cls.eval(a, tests=tests, atoms=atoms) for a in expr.args]
            return Add(*args)

        if isinstance(expr, Mul):

            coeffs  = [a for a in expr.args if isinstance(a, _coeffs_registery)]
            vectors = [a for a in expr.args if not(a in coeffs)]

            alpha = S.One
            if coeffs:
                alpha = Mul(*coeffs)

            if len(vectors) == 2:
                left, right = vectors

                newleft = cls.eval(left, tests=tests, atoms=atoms)
                newright = cls.eval(right, tests=tests, atoms=atoms)
                if _is_test_function(left) and _is_test_function(right):

                    if left in tests:
                        return alpha*wedge(newright, hodge(newleft))

                    elif right in tests:
                        return alpha*wedge(newleft, hodge(newright))

                    else:
                        raise ValueError('argument not appears as a test function')

                elif _is_op_test_function(left) and _is_test_function(right):
                    if left._args[0] in tests:
                        return alpha*wedge(newright, hodge(newleft))

                    elif right in tests:
                        return alpha*wedge(newleft, hodge(newright))

                    else:
                        raise ValueError('argument not appears as a test function')

                elif _is_test_function(left) and _is_op_test_function(right):
                    if right._args[0] in tests:
                        return alpha*wedge(newleft, hodge(newright))

                    elif left in tests:
                        return alpha*wedge(newright, hodge(newleft))

                    else:
                        raise ValueError('argument not appears as a test function')

                else:
                    convert = False
                    if isinstance(vectors[0], ScalarTestFunction):
                        left  = vectors[1]
                        right = vectors[0]
                        convert = True

                    elif isinstance(vectors[1], ScalarTestFunction):
                        left  = vectors[0]
                        right = vectors[1]
                        convert = True

                    if convert:
                        dim = right.space.ldim

                        if not(right in tests):
                            right = DifferentialForm(right.name, index=3, dim=dim)
                            return alpha*ip(left, right)

                        else:
                            right = DifferentialForm(right.name, index=0, dim=dim)
                            return alpha*jp(left, right)

            else:

                raise NotImplementedError('')

        return cls(expr, evaluate=False)


#==============================================================================
# TODO improve
def augmented_expression(expr, tests, atoms):
    trials = [atoms[i] for i in set(atoms.keys()) - set(tests) ]
    tests  = [atoms[i] for i in tests]

    constraints = {}

    deltas = list(expr.atoms(AdjointExteriorDerivative))
    trials_delta = [u for u in trials if delta(u) in deltas]

    # ... replace delta(u) for all u in trials
    for u in trials_delta:
        name  = 'aux_{}'.format(u.name)
        dim   = u.dim
        index = u.index.index-1

        aux = DifferentialForm(name, index, dim)

        expr = expr.subs(delta(u), aux)

        constraints[delta(u)] = aux
    # ...

    if constraints:
        constraints = [k-v for k,v in constraints.items()]

    return expr, constraints
