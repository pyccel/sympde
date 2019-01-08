# coding: utf-8

from sympy import Symbol, sympify
from sympy import Function

#from vale.utilities import (grad, d_var, inner, outer, cross, dot, \
#                           replace_symbol_derivatives)

DEBUG = False
#DEBUG = True

# Global variable namespace
namespace = {}
stack     = {}
settings  = {}

from sympde.topology import Domain              as sym_Domain
from sympde.topology import FunctionSpace       as sym_FunctionSpace
from sympde.topology import VectorFunctionSpace as sym_VectorFunctionSpace
from sympde.topology import Field               as sym_Field
from sympde.topology import VectorField         as sym_VectorField
from sympde.topology import TestFunction        as sym_TestFunction
from sympde.topology import VectorTestFunction  as sym_VectorTestFunction
from sympde.topology import Constant            as sym_Constant

from sympde.expr     import LinearForm          as sym_LinearForm
from sympde.expr     import BilinearForm        as sym_BilinearForm


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
            namespace[name] = sym_Domain(name, dim=dim)

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

        namespace[name] = V

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
class TestFunction(BasicPDE):
    """Class representing a test function."""
    def __init__(self, **kwargs):
        """
        A Field has the following attributs

        * name
        """
        name  = kwargs.pop('name')

        # the appropriate object will created later in the Linear/Bilinear form
        namespace[name] = Symbol(name)

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
            v = sym_Field(name, space=space)

        elif isinstance(space, sym_VectorFunctionSpace):
            v = sym_VectorField(name, space=space)

        namespace[name] = v

        self.name = name
        BasicPDE.__init__(self, **kwargs)

#======================================================================
class Real(BasicPDE):
    """Class representing a Real number."""
    def __init__(self, **kwargs):
        name  = kwargs.pop('name')

        namespace[name] = sym_Constant(name, real=True)

        self.name = name
        BasicPDE.__init__(self, **kwargs)

#======================================================================
class Complex(BasicPDE):
    """Class representing a Complex number."""
    def __init__(self, **kwargs):
        name  = kwargs.pop('name')

        namespace[name] = sym_Constant(name, complex=True)

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

#    @property
#    def expr(self):
#        return Symbol(self.name)
#

#======================================================================
class LinearForm(BasicPDE):
    """Class representing a Linear Form."""

    def __init__(self, **kwargs):
        name = kwargs.pop('name')
        args = kwargs.pop('args')
        body = kwargs.pop('body')

        # ... create test functions
        space = namespace[args.space]
        functions = []
        if isinstance(space, sym_FunctionSpace):
            for i in args.functions:
                v = sym_TestFunction(space, name=i.name)
                functions.append(v)

        elif isinstance(space, sym_VectorFunctionSpace):
            for i in args.functions:
                v = sym_VectorTestFunction(space, name=i.name)
                functions.append(v)

        for v in functions:
            namespace[v.name] = v

        if len(functions) == 1:
            functions = functions[0]
        # ...

        # ...
        if isinstance(body, SimpleBodyForm):
            expression = body.expression.expr

        elif isinstance(body, ExpressionBodyForm):
            raise NotImplementedError('')
        # ...

        namespace[name] = sym_LinearForm(functions, expression,
                                         name=name)

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
        functions = []
        if isinstance(space, sym_FunctionSpace):
            for i in args.functions:
                v = sym_TestFunction(space, name=i.name)
                functions.append(v)

        elif isinstance(space, sym_VectorFunctionSpace):
            for i in args.functions:
                v = sym_VectorTestFunction(space, name=i.name)
                functions.append(v)

        for v in functions:
            namespace[v.name] = v

        if len(functions) == 1:
            functions = functions[0]

        test_functions = functions
        # ...

        # ... create trial functions
        args = args_trial
        space = namespace[args.space]
        functions = []
        if isinstance(space, sym_FunctionSpace):
            for i in args.functions:
                v = sym_TestFunction(space, name=i.name)
                functions.append(v)

        elif isinstance(space, sym_VectorFunctionSpace):
            for i in args.functions:
                v = sym_VectorTestFunction(space, name=i.name)
                functions.append(v)

        for v in functions:
            namespace[v.name] = v

        if len(functions) == 1:
            functions = functions[0]

        trial_functions = functions
        # ...

        # ...
        if isinstance(body, SimpleBodyForm):
            expression = body.expression.expr

        elif isinstance(body, ExpressionBodyForm):
            raise NotImplementedError('')
        # ...

        args = (test_functions, trial_functions)
        namespace[name] = sym_BilinearForm(args, expression,
                                           name=name)

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
        if DEBUG:
            print("> Trailer ")
        return [i.expr.expr for i in self.args]

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
class FactorSigned(ExpressionElement):
    """Class representing a signed factor."""
    def __init__(self, **kwargs):
        self.sign = kwargs.pop('sign', '+')
        self.trailer = kwargs.pop('trailer', [])
        super(FactorSigned, self).__init__(**kwargs)

    @property
    def expr(self):
        if DEBUG:
            print("> FactorSigned ")
        trailer = self.trailer
        if trailer:
            print(trailer.expr)
            print('PAR ICI')
            import sys; sys.exit(0)
        expr = self.op.expr
        return -expr if self.sign == '-' else expr


#======================================================================
class Term(ExpressionElement):
    @property
    def expr(self):
        if DEBUG:
            print("> Term ")
        ret = self.op[0].expr
        for operation, operand in zip(self.op[1::2], self.op[2::2]):
            if operation == '*':
                ret *= sympify(operand.expr)
            else:
                ret /= sympify(operand.expr)
        return ret


#======================================================================
class Expression(ExpressionElement):
    @property
    def expr(self):
        if DEBUG:
            print("> Expression ")
        ret = self.op[0].expr
        for operation, operand in zip(self.op[1::2], self.op[2::2]):
            if operation == '+':
                ret += sympify(operand.expr)
            else:
                ret -= sympify(operand.expr)
        return ret


#======================================================================
class Operand(ExpressionElement):
    @property
    def expr(self):
        if DEBUG:
            print("> Operand ")
            print(("> stack : ", stack))
            print((self.op))
#        op = self.op[0]
        op = self.op
        if type(op) in {int, float}:
            return op
        elif type(op) == list:
            # op is a list
            for O in op:
                if O in namespace:
                    # TODO use isinstance
                    if type(namespace[O]) in [Field, Function, Real]:
                        return namespace[O].expr
                    else:
                        return namespace[O]
                elif O in stack:
                    if DEBUG:
                        print((">>> found local variables: " + O))
                    return Symbol(O)
                else:
                    raise Exception('Unknown variable "{}" at position {}'
                                    .format(O, self._tx_position))
        elif isinstance(op, ExpressionElement):
            return op.expr
        elif op in stack:
            if DEBUG:
                print((">>> found local variables: " + op))
            return Symbol(op)
        elif op in namespace:
            # TODO use isinstance
            if type(namespace[op]) in [Field, Function, Real]:
                return namespace[op].expr
            else:
                return namespace[op]
        else:
            raise Exception('Unknown variable "{}" at position {}'
                            .format(op, self._tx_position))


#======================================================================
class ExpressionBodyForm(ExpressionElement):
    @property
    def expr(self):
        if DEBUG:
            print("> ExpressionBodyForm ")
        ret = self.op[0].expr
        for operation, operand in zip(self.op[1::2], self.op[2::2]):
            if operation in ['+', '-']:
                ret += operation + ' ' + operand.expr
            else:
                raise Exception('Unknown operation "{}" at position {}'
                                .format(operation, self._tx_position))
        return ret


#======================================================================
class TermForm(ExpressionElement):
    @property
    def expr(self):
        if DEBUG:
            print("> TermForm ")
        ret = self.op[0].expr
        for operation, operand in zip(self.op[1::2], self.op[2::2]):
            if operation == '*':
                ret += '*' + ' ' + operand.expr
            else:
                raise Exception('Unknown operation "{}" at position {}'
                                .format(operation, self._tx_position))
        return ret


#======================================================================
class CallForm(ExpressionBodyForm):
    """Class representing the call to a linear/bilinear form."""
    def __init__(self, **kwargs):

        self.name = kwargs.pop('name')
        self.args = kwargs.pop('args')

#        super(ExpressionBodyForm, self).__init__()

    @property
    def expr(self):
        if DEBUG:
            print("> CallForm")
        if self.name in namespace:
            b  = namespace[stack["parent"]]

            if isinstance(b, LinearForm):
                if len(self.args) != 1:
                    raise Exception('Expecting exactly one argument at position {}.'
                                    .format(self._tx_position))

                i_row = b.args.functions.index(self.args[0])
                bi = namespace[self.name]
                b.blocks[i_row] = bi
            elif isinstance(b, BilinearForm):
                if len(self.args) != 2:
                    raise Exception('Expecting exactly two argument at position {}.'
                                    .format(self._tx_position))

                i_row = b.args_test.functions.index(self.args[0])
                i_col = b.args_trial.functions.index(self.args[1])
                bi = namespace[self.name]
                b.blocks[i_row, i_col] = bi
            else:
                raise ValueError("Expecting a Linear or Bilinear form.")

        return self.name

