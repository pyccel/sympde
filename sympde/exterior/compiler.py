# coding: utf-8

# TODO - assert the space type is not Undefined

from collections import OrderedDict

from sympy import Basic
from sympy import Add, Mul
from sympy import S
from sympy import Indexed
from sympy import Tuple

from sympde.core.basic import _coeffs_registery
from sympde.core.basic import CalculusFunction
from sympde.topology   import H1SpaceType, HcurlSpaceType, HdivSpaceType
from sympde.topology   import L2SpaceType, UndefinedSpaceType
from sympde.topology   import ScalarTestFunction, VectorTestFunction
from sympde.calculus   import Grad, Curl, Div
from sympde.calculus   import Dot, Inner, Cross
#from sympde.calculus import grad, dot, inner, cross, rot, curl, div

from .form     import DifferentialForm
from .calculus import d, wedge, ip, ld
from .calculus import delta, jp, hodge
from .calculus import AdjointExteriorDerivative

#==============================================================================

_is_function      = lambda u: isinstance(u, (ScalarTestFunction, VectorTestFunction))
_is_test_function = lambda u: _is_function(u) and not isinstance(u.space.kind, UndefinedSpaceType)
_is_field         = lambda u: _is_function(u) and     isinstance(u.space.kind, UndefinedSpaceType)
_is_proxy         = lambda u: _is_field(u) or isinstance(u, Tuple)

_is_op_test_function = lambda op: (isinstance(op, (Grad, Curl, Div)) and
                                   _is_test_function(op.args[0]))

#==============================================================================
def _decompose_lie_derivative(*args):

    #---------------------------------------
    # Try to extract Grad and Cross terms
    #---------------------------------------

    # ... 1 form
    grad_term  = [i for i in args if (isinstance(i, Grad) and
                                      isinstance(i._args[0], Dot))]

    # TODO improve
    cross_term = [i for i in args if (isinstance(i, Cross) and
                                      (isinstance(i._args[0], Curl) or
                                       isinstance(i._args[1], Curl)))]

    if grad_term and cross_term:
        others = [i for i in args if i not in (grad_term + cross_term)]
        return grad_term, cross_term, others

    #---------------------------------------
    # Try to extract Curl and Mul terms
    #---------------------------------------
    curl_term = []
    mul_term  = []
    others    = []

    for i in args:

        # ... 1 form
        if isinstance(i, Curl) and isinstance(i.args[0], Cross):
            curl_term.append(i)

        # Special case to handle anti-commutativity of Curl operator:
        # Mul(-1, Curl(Cross(a, b))) --> Curl(Cross(b, a))
        # TODO avoid this hack
        elif (isinstance(i, Mul) and len(i.args) == 2 and i.args[0] == -1 and
              isinstance(i.args[1], Curl) and
              isinstance(i.args[1].args[0], Cross)
              ):
            a, b = i.args[1].args[0].args
            new_cross = Cross(b, a, evaluate=False)
            new_curl  = Curl(new_cross)
            curl_term.append(new_curl)

        # TODO improve
        elif (isinstance(i, Mul) and (isinstance(i.args[0], Div) or
                                      isinstance(i.args[1], Div)) ):
            mul_term.append(i)

        else:
            others.append(i)

    if curl_term and mul_term:
        return curl_term, mul_term, others

    #---------------------------------------
    # Unsuccessful: return input arguments
    #---------------------------------------

    return [], [], args

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

        if _is_field(expr):
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

            if isinstance(arg, Mul):
                if len(arg.args) == 2:
                    left, right = arg.args[:]

                    newleft = cls.eval(left, tests=tests, atoms=atoms)
                    newright = cls.eval(right, tests=tests, atoms=atoms)

                    if  _is_test_function(right) and _is_proxy(left):
                        return ld(newleft, newright)

                    elif  _is_test_function(left) and _is_proxy(right):
                        return ld(newright, newleft)

                else:
                    raise NotImplementedError('')

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

            elif  _is_op_test_function(right) and _is_proxy(left):
                raise NotImplementedError('')

            elif  _is_op_test_function(left) and _is_proxy(right):
                if isinstance(left, Grad):
                    if left._args[0] in tests:
                        raise NotImplementedError('')

                    else:
                        return ld(newright, newleft._args[0])

                else:
                    raise NotImplementedError('')

            else:
                print(type(left))
                print(type(right))
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
                if (_is_function(left._args[0]) and
                    _is_function(right._args[0])):

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
            # looking for Lie derivative terms
            a, b, c = _decompose_lie_derivative(*expr.args)

            newargs = []

            for ai in a:
                for bi in b:

                    ai_tests = [i for i in ai.atoms() if _is_test_function(i)]
                    bi_tests = [i for i in bi.atoms() if _is_test_function(i)]

                    test = set(ai_tests) & set(bi_tests)
                    test = list(test)
                    if len(test) == 1:
                        test = test[0]

                        # to distinguish between lie deriv of 1 and 2 forms
                        if isinstance(ai, Grad) and isinstance(bi, Cross):
                            if isinstance(ai._args[0], Dot):
                                # TODO implement anti commutativity for
                                #      Cross
                                if (isinstance(bi._args[0], Curl) and
                                    (_is_field(bi._args[1]) or
                                     isinstance(bi._args[1], Tuple))):
                                    bi_test  = bi._args[0]._args[0]
                                    bi_proxy = bi._args[1]
                                    # grad(dot(,)) args
                                    ai_left,ai_right = ai._args[0]._args[:]
                                    if _is_test_function(ai_left):
                                        ai_proxy = ai_right
                                        ai_test  = ai_left

                                    elif _is_test_function(ai_right):
                                        ai_proxy = ai_left
                                        ai_test  = ai_right

                                    if ((ai_test == bi_test) and
                                        (ai_proxy == bi_proxy)):

                                        if not ai_test == test:
                                            raise NotImplementedError('')

                                        ai_proxy = cls.eval(ai_proxy, tests=tests, atoms=atoms)
                                        ai_test  = cls.eval(ai_test, tests=tests, atoms=atoms)

                                        newargs = [ld(ai_proxy, ai_test)]

                                else:
                                    raise NotImplementedError('')

                        # to distinguish between lie deriv of 1 and 2 forms
                        elif isinstance(ai, Curl) and isinstance(bi, Mul):
                            if isinstance(ai._args[0], Cross) and len(bi.args) == 2:
                                # TODO implement anti commutativity for
                                #      Cross
                                if (_is_test_function(bi._args[0]._args[0]) and
                                    (_is_field(bi._args[1]) or
                                     isinstance(bi._args[1], Tuple))):

                                    bi_test  = bi._args[0]._args[0]
                                    bi_proxy = bi._args[1]
                                    # grad(dot(,)) args
                                    ai_left,ai_right = ai._args[0]._args[:]
                                    if _is_test_function(ai_left):
                                        ai_proxy = ai_right
                                        ai_test  = ai_left

                                    elif _is_test_function(ai_right):
                                        ai_proxy = ai_left
                                        ai_test  = ai_right

                                    if ((ai_test == bi_test) and
                                        (ai_proxy == bi_proxy)):

                                        if not ai_test == test:
                                            raise NotImplementedError('')

                                        ai_proxy = cls.eval(ai_proxy, tests=tests, atoms=atoms)
                                        ai_test  = cls.eval(ai_test, tests=tests, atoms=atoms)

                                        newargs = [ld(ai_proxy, ai_test)]

                                else:
                                    raise NotImplementedError('')

                    elif len(test) > 1:
                        raise ValueError('wrong expression')

            args = [cls.eval(i, tests=tests, atoms=atoms) for i in c] + newargs
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
def augmented_expression(expr, tests, atoms, weak=True):
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

    if not weak:
        return expr, constraints

    else:
        raise NotImplementedError('')
