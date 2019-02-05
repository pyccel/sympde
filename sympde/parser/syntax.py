# coding: utf-8

from sympy import Tuple
from sympy import Function
from sympy import pi
from sympy import Pow

#from vale.utilities import (grad, d_var, inner, outer, cross, dot, \
#                           replace_symbol_derivatives)

DEBUG = False
#DEBUG = True

from pyccel.ast.utilities import math_functions as _known_functions_math

from sympde.topology import Domain              as sym_Domain
from sympde.topology import Boundary            as sym_Boundary
from sympde.topology import NormalVector        as sym_NormalVector
from sympde.topology import TangentVector       as sym_TangentVector
from sympde.topology import FunctionSpace       as sym_FunctionSpace
from sympde.topology import VectorFunctionSpace as sym_VectorFunctionSpace
from sympde.topology import ProductSpace        as sym_ProductSpace
from sympde.topology import Field               as sym_Field
from sympde.topology import VectorField         as sym_VectorField
from sympde.topology import TestFunction        as sym_TestFunction
from sympde.topology import VectorTestFunction  as sym_VectorTestFunction
from sympde.topology import Constant            as sym_Constant

from sympde.expr     import LinearForm          as sym_LinearForm
from sympde.expr     import BilinearForm        as sym_BilinearForm
from sympde.expr     import Equation            as sym_Equation
from sympde.expr     import EssentialBC         as sym_EssentialBC

from sympde.calculus import grad, dot, inner, cross, rot, curl, div
from sympde.calculus import laplace, hessian, bracket
from sympde.topology import (dx, dy, dz)

_known_operators = {
    'dx':      dx,
    'dy':      dy,
    'dz':      dz,
    'grad':    grad,
    'dot':     dot,
    'inner':   inner,
    'cross':   cross,
    'rot':     rot,
    'curl':    curl,
    'div':     div,
    'laplace': laplace,
    'hessian': hessian,
    'bracket': bracket,
}

_known_constants_math = {'pi':pi}

# Global variable namespace
namespace = {}
stack     = {}
settings  = {}

#======================================================================
def insert_namespace(key, value):
    if key in namespace.keys():
        raise ValueError('{} already defined'.format(key))

    namespace[key] = value

#======================================================================
class BasicPDE(object):

    def __init__(self, **kwargs):
        self.namespace = namespace

#======================================================================
class PDE(BasicPDE):
    """Class for PDE syntax."""
    def __init__(self, **kwargs):
        self.declarations = kwargs.pop('declarations')
        BasicPDE.__init__(self, **kwargs)

#======================================================================
class Domain(BasicPDE):
    """Class representing a Domain."""
    def __init__(self, **kwargs):
        name = kwargs.pop('name')
        dim  = kwargs.pop('dim', None)
        filename = kwargs.pop('filename', None)

        expr = self
        if not( dim is None ):
            atom = sym_Domain(name, dim=dim)
            insert_namespace(name, atom)

            # insert coordinates
            for x in atom.coordinates:
                insert_namespace(x.name, x)

            # insert normal and tangent vectors
            insert_namespace('nn', sym_NormalVector('nn'))
            insert_namespace('tt', sym_TangentVector('tt'))

        elif not( filename is None ):
            raise NotImplementedError('')

        self.name = name
        BasicPDE.__init__(self, **kwargs)

#======================================================================
# TODO kind is not used yet
class FunctionSpace(BasicPDE):
    """Class representing a Finite Element FunctionSpace."""
    def __init__(self, **kwargs):
        name   = kwargs.pop('name')
        domain = kwargs.pop('domain')
        kind   = kwargs.pop('kind', 'h1')

        domain = namespace[domain]
        V = sym_FunctionSpace(name, domain)

        insert_namespace(name, V)

        self.name = name
        BasicPDE.__init__(self, **kwargs)

#======================================================================
# TODO kind is not used yet
class VectorFunctionSpace(BasicPDE):
    """Class representing a Finite Element VectorFunctionSpace."""
    def __init__(self, **kwargs):
        name   = kwargs.pop('name')
        domain = kwargs.pop('domain')
        kind   = kwargs.pop('kind', 'h1')

        domain = namespace[domain]
        V = sym_VectorFunctionSpace(name, domain)

        namespace[name] = V

        self.name = name
        BasicPDE.__init__(self, **kwargs)


#======================================================================
class Field(BasicPDE):
    """Class representing a Field."""
    def __init__(self, **kwargs):
        name  = kwargs.pop('name')
        space = kwargs.pop('space')

        space = namespace[space]
        if isinstance(space, sym_FunctionSpace):
            v = sym_Field(space, name=name)

        elif isinstance(space, sym_VectorFunctionSpace):
            v = sym_VectorField(name, space=space)

        insert_namespace(name, v)

        self.name = name
        BasicPDE.__init__(self, **kwargs)

#======================================================================
class Alias(BasicPDE):
    """Class representing an Alias."""
    def __init__(self, **kwargs):
        name  = kwargs.pop('name')
        rhs = kwargs.pop('rhs')

        insert_namespace(name, rhs.expr)

        self.name = name
        BasicPDE.__init__(self, **kwargs)

#======================================================================
# TODO - add bc
#      - add name/label
class Equation(BasicPDE):
    """Class representing an Equation."""
    def __init__(self, **kwargs):
        name  = kwargs.pop('name', 'equation')

        tests  = kwargs.pop('tests')
        trials = kwargs.pop('trials')
        lhs    = kwargs.pop('lhs')
        rhs    = kwargs.pop('rhs')
        bc     = kwargs.pop('bc', None)

        # ... create test functions
        args = tests
        space = namespace[args.space]

        if isinstance(space, sym_ProductSpace):
            space = space.spaces
        else:
            space = [space]

        assert(len(space) == len(args.functions))

        functions = []
        for i,V in zip(args.functions, space):
            if isinstance(V, sym_FunctionSpace):
                v = sym_TestFunction(V, name=i.name)

            elif isinstance(V, sym_VectorFunctionSpace):
                v = sym_VectorTestFunction(V, name=i.name)

            functions.append(v)

        for v in functions:
            insert_namespace(v.name, v)

        if len(functions) == 1:
            functions = functions[0]

        test_functions = functions
        # ...

        # ... create trial functions
        args = trials
        space = namespace[args.space]

        if isinstance(space, sym_ProductSpace):
            space = space.spaces
        else:
            space = [space]

        assert(len(space) == len(args.functions))
        functions = []
        for i,V in zip(args.functions, space):
            if isinstance(V, sym_FunctionSpace):
                v = sym_TestFunction(V, name=i.name)

            elif isinstance(V, sym_VectorFunctionSpace):
                v = sym_VectorTestFunction(V, name=i.name)

            functions.append(v)

        for v in functions:
            insert_namespace(v.name, v)

        if len(functions) == 1:
            functions = functions[0]

        trial_functions = functions
        # ...

        # ... prepare boundary conditions
        if bc:
            # TODO get domain from space
            domain = [i for i in namespace.values() if isinstance(i, sym_Domain)]
            domain = domain[0]

            _bc = []
            for b in bc:
                bnd     = b.boundary
                bnd     = sym_Boundary(bnd, domain)

                bnd_lhs = b.lhs.expr
                bnd_rhs = b.rhs.expr

                sym_bc = sym_EssentialBC(bnd_lhs, bnd_rhs, bnd)
                _bc.append(sym_bc)

            bc = _bc
        # ...

        # ... define sympde Equation
        rhs = rhs.expr
        lhs = lhs.expr

        atom = sym_Equation(lhs, rhs, bc=bc)
        insert_namespace(name, atom)
        # ...

        # ... clean namespace
        if not isinstance(test_functions, (tuple, list, Tuple)):
            test_functions = [test_functions]

        if not isinstance(trial_functions, (tuple, list, Tuple)):
            trial_functions = [trial_functions]

        for v in test_functions + trial_functions:
            namespace.pop(v.name)
        # ...

        self.name = name
        BasicPDE.__init__(self, **kwargs)

#======================================================================
class Real(BasicPDE):
    """Class representing a Real number."""
    def __init__(self, **kwargs):
        name  = kwargs.pop('name')

        atom = sym_Constant(name, real=True)
        insert_namespace(name, atom)

        self.name = name
        BasicPDE.__init__(self, **kwargs)

#======================================================================
class Complex(BasicPDE):
    """Class representing a Complex number."""
    def __init__(self, **kwargs):
        name  = kwargs.pop('name')

        atom = sym_Constant(name, complex=True)
        insert_namespace(name, atom)

        self.name = name
        BasicPDE.__init__(self, **kwargs)

# TODO
class Function(object):
    """Class representing a Function."""
    def __init__(self, **kwargs):
        """
        A Function has the following attributs

        * name
        * parameters

        .. note::
            The grammar rule to define a Function is

            Function:
            "Function" LPAREN parameters*=ID[','] RPAREN  DEF name=ID
            ;
        """
        self.name       = kwargs.pop('name')
        self.parameters = kwargs.pop('parameters', {})

        namespace[self.name] = self


#======================================================================
class LinearForm(BasicPDE):
    """Class representing a Linear Form."""

    def __init__(self, **kwargs):
        name = kwargs.pop('name')
        args = kwargs.pop('args')
        body = kwargs.pop('body')

        # ... create test functions
        space = namespace[args.space]

        if isinstance(space, sym_ProductSpace):
            space = space.spaces
        else:
            space = [space]

        assert(len(space) == len(args.functions))

        functions = []
        for i,V in zip(args.functions, space):
            if isinstance(V, sym_FunctionSpace):
                v = sym_TestFunction(V, name=i.name)

            elif isinstance(V, sym_VectorFunctionSpace):
                v = sym_VectorTestFunction(V, name=i.name)

            functions.append(v)

        for v in functions:
            insert_namespace(v.name, v)

        if len(functions) == 1:
            functions = functions[0]
        # ...

        # ...
        if isinstance(body, SimpleBodyForm):
            expression = body.expression.expr

        elif isinstance(body, Expression):
            expression = body.expr
        # ...

        atom = sym_LinearForm(functions, expression)
        insert_namespace(name, atom)

        # ... clean namespace
        if not isinstance(functions, (tuple, list, Tuple)):
            functions = [functions]

        for v in functions:
            namespace.pop(v.name)
        # ...

        self.name = name
        BasicPDE.__init__(self, **kwargs)

#======================================================================
class BilinearForm(BasicPDE):
    """Class representing a Bilinear Form."""

    def __init__(self, **kwargs):
        name       = kwargs.pop('name')
        args_test  = kwargs.pop('args_test')
        args_trial = kwargs.pop('args_trial')
        body       = kwargs.pop('body')

        # ... create test functions
        args = args_test
        space = namespace[args.space]

        if isinstance(space, sym_ProductSpace):
            space = space.spaces
        else:
            space = [space]

        assert(len(space) == len(args.functions))

        functions = []
        for i,V in zip(args.functions, space):
            if isinstance(V, sym_FunctionSpace):
                v = sym_TestFunction(V, name=i.name)

            elif isinstance(V, sym_VectorFunctionSpace):
                v = sym_VectorTestFunction(V, name=i.name)

            functions.append(v)

        for v in functions:
            insert_namespace(v.name, v)

        if len(functions) == 1:
            functions = functions[0]

        test_functions = functions
        # ...

        # ... create trial functions
        args = args_trial
        space = namespace[args.space]

        if isinstance(space, sym_ProductSpace):
            space = space.spaces
        else:
            space = [space]

        assert(len(space) == len(args.functions))
        functions = []
        for i,V in zip(args.functions, space):
            if isinstance(V, sym_FunctionSpace):
                v = sym_TestFunction(V, name=i.name)

            elif isinstance(V, sym_VectorFunctionSpace):
                v = sym_VectorTestFunction(V, name=i.name)

            functions.append(v)

        for v in functions:
            insert_namespace(v.name, v)

        if len(functions) == 1:
            functions = functions[0]

        trial_functions = functions
        # ...

        # ...
        if isinstance(body, SimpleBodyForm):
            expression = body.expression.expr

        elif isinstance(body, Expression):
            expression = body.expr
        # ...

        args = (test_functions, trial_functions)
        atom = sym_BilinearForm(args, expression)
        insert_namespace(name, atom)

        # ... clean namespace
        if not isinstance(test_functions, (tuple, list, Tuple)):
            test_functions = [test_functions]

        if not isinstance(trial_functions, (tuple, list, Tuple)):
            trial_functions = [trial_functions]

        for v in test_functions + trial_functions:
            namespace.pop(v.name)
        # ...

        self.name = name
        BasicPDE.__init__(self, **kwargs)


#======================================================================
class BodyForm(object):
    """Class representing the body of a linear/bilinear form."""
    def __init__(self, **kwargs):

        super(BodyForm, self).__init__()


#======================================================================
class SimpleBodyForm(BodyForm):
    """Class representing the body of a simple linear/bilinear form."""
    def __init__(self, **kwargs):

        self.expression = kwargs.pop('expression')
        self.domain = kwargs.pop('domain', None)

        super(SimpleBodyForm, self).__init__()


#======================================================================
class Trailer(object):
    """Class representing a trailer."""
    def __init__(self, **kwargs):
        self.args = kwargs.pop('args', [])

    @property
    def expr(self):
        return [i.expr for i in self.args]

#======================================================================
class Power(object):
    """Class representing a power."""
    def __init__(self, **kwargs):
        self.arg = kwargs.pop('arg', None)
        self.op  = kwargs.pop('op',  None)

    @property
    def expr(self):
        return Pow(self.op.expr, self.arg.expr, evaluate=False)

#======================================================================
class ExpressionElement(object):
    """Class representing an element of an expression."""
    def __init__(self, **kwargs):

        # textX will pass in parent attribute used for parent-child
        # relationships. We can use it if we want to.
        self.parent = kwargs.get('parent', None)

        # We have 'op' attribute in all grammar rules
        self.op = kwargs['op']

        super(ExpressionElement, self).__init__()


#======================================================================
class Factor(ExpressionElement):
    """Class representing a signed factor."""
    def __init__(self, **kwargs):
        self.sign = kwargs.pop('sign', '+')
        self.trailer = kwargs.pop('trailer', [])
        super(Factor, self).__init__(**kwargs)

    @property
    def expr(self):
        expr = self.op.expr
        trailer = self.trailer
        if trailer:
            expr = expr(*trailer.expr)

        return -expr if self.sign == '-' else expr


#======================================================================
class Term(ExpressionElement):
    @property
    def expr(self):
        ret = self.op[0].expr
        for operation, operand in zip(self.op[1::2], self.op[2::2]):
            if operation == '*':
                ret *= operand.expr
            else:
                ret /= operand.expr
        return ret


#======================================================================
class Expression(ExpressionElement):
    @property
    def expr(self):
        ret = self.op[0].expr
        for operation, operand in zip(self.op[1::2], self.op[2::2]):
            if operation == '+':
                ret += operand.expr
            else:
                ret -= operand.expr
        return ret

#======================================================================
class Operand(ExpressionElement):
    @property
    def expr(self):
#        if DEBUG:
#        if True:
#            print("> Operand ")
#            print(self.op, type(self.op))

        op = self.op
        if type(op) in {int, float}:
            return op

        elif isinstance(op, str):
            if op in _known_operators.keys():
                return _known_operators[op]

            if op in _known_functions_math.keys():
                return _known_functions_math[op]

            if op in _known_constants_math.keys():
                return _known_constants_math[op]

            if not(op in namespace.keys()):
                raise ValueError('{} not found'.format(op))

            return namespace[op]

#        elif type(op) == list:
#            # op is a list
#            for O in op:
#                if O in namespace:
#                    # TODO use isinstance
#                    if type(namespace[O]) in [Field, Function, Real]:
#                        return namespace[O].expr
#                    else:
#                        return namespace[O]
#                elif O in stack:
#                    if DEBUG:
#                        print((">>> found local variables: " + O))
#                    return Symbol(O)
#                else:
#                    raise Exception('Unknown variable "{}" at position {}'
#                                    .format(O, self._tx_position))
#        elif isinstance(op, ExpressionElement):
#            return op.expr
#        elif op in stack:
#            if DEBUG:
#                print((">>> found local variables: " + op))
#            return Symbol(op)
#        elif op in namespace:
#            # TODO use isinstance
#            if type(namespace[op]) in [Field, Function, Real]:
#                return namespace[op].expr
#            else:
#                return namespace[op]
        else:
            raise Exception('Unknown variable "{}" at position {}'
                            .format(op, self._tx_position))
