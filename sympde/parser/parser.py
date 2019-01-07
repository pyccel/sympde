# coding: utf-8

import os
from sympy import Symbol, sympify

#from .utilities import grad, d_var, inner, outer, cross, dot
from sympde.parser.syntax import (PDE,
                     Expression, Term, Operand,
                     FactorSigned, FactorUnary, FactorBinary,
                     LinearForm, BilinearForm,
                     BodyForm, SimpleBodyForm, ExpressionBodyForm,
                     TermForm, CallForm,
                     Domain, FunctionSpace, VectorFunctionSpace, Field, Function, Real)

from textx.metamodel import metamodel_from_str

# ...
def get_by_name(ast, name):
    """
    Returns an object from the AST by giving its name.
    """
    for token in ast.declarations:
        if token.name == name:
            return token
    return None
# ...

# ...
def ast_to_dict(ast):
    """
    Returns an object from the AST by giving its name.
    """
    tokens = {}
    for token in ast.declarations:
        tokens[token.name] = token
    return tokens
# ...

class BasicParser(object):
    """ Class for a Parser using TextX.

    A parser can be created from a grammar (str) or a filename. It is preferable
    to specify the list classes to have more control over the abstract grammar;
    for example, to use a namespace, and to do some specific anotation.

    >>> parser = Parser(filename="gammar.tx")

    Once the parser is created, you can parse a given set of instructions by
    calling

    >>> parser.parse(["Field(V) :: u"])

    or by providing a file to parse

    >>> parser.parse_from_file("tests/inputs/1d/poisson.vl")
    """
    def __init__(self, grammar=None, filename=None, \
                 classes=None):
        """Parser constructor.

        grammar : str
            abstract grammar describing the DSL.

        filename: str
            name of the file containing the abstract grammar.

        classes : list
            a list of Python classes to be used to describe the grammar. Take a
            look at TextX documentation for more details.
        """

        _grammar = grammar

        # ... read the grammar from a file
        if not (filename is None):
            dir_path = os.path.dirname(os.path.realpath(__file__))
            filename = os.path.join(dir_path, filename)

            f = open(filename)
            _grammar = f.read()
            _grammar.replace("\n", "")
            f.close()
        # ...

        # ...
        self.grammar = _grammar
        # ...

        # ...
        if classes is None:
            self.model = metamodel_from_str(_grammar)
        else:
            self.model = metamodel_from_str(_grammar, classes=classes)
        # ...

    def parse(self, instructions):
        """Parse a set of instructions with respect to the grammar.

        instructions: list
            list of instructions to parse.
        """
        # ... parse the DSL code
        return self.model.model_from_str(instructions)
        # ...

    def parse_from_file(self, filename):
        """Parse a set of instructions with respect to the grammar.

        filename: str
            a file containing the instructions to parse.
        """
        # ... read a DSL code
        f = open(filename)
        instructions = f.read()
        instructions.replace("\n", "")
        f.close()
        # ...

        # ... parse the DSL code
        return self.parse(instructions)
        # ...

# User friendly parser

class Parser(BasicParser):
    """A Class for SymPDE parser.

    This is an extension of the Parser class. Additional treatment is done for
    Linear and Bilinear Forms to define their dependencies: user_fields,
    user_functions and user_constants.

    """
    def __init__(self, **kwargs):
        """parser constructor.

        It takes the same arguments as the Parser class.
        """
        classes = [PDE,
                   Expression, Term, Operand,
                   FactorSigned, FactorUnary, FactorBinary,
                   LinearForm, BilinearForm,
                   BodyForm, SimpleBodyForm, ExpressionBodyForm,
                   TermForm, CallForm,
                   Domain, FunctionSpace, VectorFunctionSpace, Field, Function, Real
                   ]

        try:
            filename = kwargs["filename"]
        except:
            filename = "grammar.tx"

        super(Parser, self).__init__(filename = filename,
                                     classes=classes)

    def parse_from_file(self, filename):
        """Parse a set of instructions with respect to the grammar and returns
        the AST.

        filename: str
            a file containing the instructions to parse.
        """
        ast = super(Parser, self).parse_from_file(filename)
        print(ast)

        # ... annotating the AST
        for token in ast.declarations:
            print(token)
#            annotate_form(token, ast)
        # ...

        return ast

## ...
#def annotate_form(token, ast):
#    """Annotates a linear/bilinear form using the ast.
#
#    token: sympde.syntax.LinearForm or sympde.syntax.BilinearForm
#        a linear/bilinear form object to annotate
#
#    ast: list
#        abstract tree
#    """
#    if isinstance(token, LinearForm):
#        user_fields    = []
#        user_functions = []
#        user_constants = []
#
#        space  = get_by_name(ast, token.args.space)
#        domain = get_by_name(ast, space.domain)
#
#        expr = token.to_sympy()
#        if type(expr) == dict:
#            for key, form in list(expr.items()):
#                free_symbols = form.free_symbols
#                for symbol in free_symbols:
#                    var = get_by_name(ast, str(symbol))
#                    if isinstance(var, Field):
#                        user_fields.append(var.name)
#                    elif isinstance(var, Function):
#                        user_functions.append(var.name)
#                    elif isinstance(var, Real):
#                        user_constants.append(var.name)
#        else:
#            free_symbols = expr.free_symbols
#            for symbol in free_symbols:
#                var = get_by_name(ast, str(symbol))
#                if isinstance(var, Field):
#                    user_fields.append(var.name)
#                elif isinstance(var, Function):
#                    user_functions.append(var.name)
#                elif isinstance(var, Real):
#                    user_constants.append(var.name)
#
#        token.set("dim", domain.dim)
#        token.set("user_fields", user_fields)
#        token.set("user_functions", user_functions)
#        token.set("user_constants", user_constants)
#
#    elif isinstance(token, BilinearForm):
#        user_fields    = []
#        user_functions = []
#        user_constants = []
#
#        space_test  = get_by_name(ast, token.args_test.space)
#        space_trial = get_by_name(ast, token.args_trial.space)
#        domain      = get_by_name(ast, space_test.domain)
#
#        expr = token.to_sympy()
#        if type(expr) == dict:
#            for key, form in list(expr.items()):
#                free_symbols = form.free_symbols
#                for symbol in free_symbols:
#                    var = get_by_name(ast, str(symbol))
#                    if isinstance(var, Field):
#                        user_fields.append(var.name)
#                    elif isinstance(var, Function):
#                        user_functions.append(var.name)
#                    elif isinstance(var, Real):
#                        user_constants.append(var.name)
#        else:
#            free_symbols = expr.free_symbols
#            for symbol in free_symbols:
#                var = get_by_name(ast, str(symbol))
#                if isinstance(var, Field):
#                    user_fields.append(var.name)
#                elif isinstance(var, Function):
#                    user_functions.append(var.name)
#                elif isinstance(var, Real):
#                    user_constants.append(var.name)
#
#        token.set("dim", domain.dim)
#        token.set("user_fields", user_fields)
#        token.set("user_functions", user_functions)
#        token.set("user_constants", user_constants)
#
#    return token
## ...
