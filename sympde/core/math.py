# coding: utf-8

from sympy import Function
from sympy import preorder_traversal
from sympy import NumberSymbol
from sympy.printing.pycode import _known_functions_math
from sympy.printing.pycode import _known_constants_math

def math_atoms_as_str(expr):
    math_functions = [str(type(i)) for i in preorder_traversal(expr) if isinstance(i, Function)]
    math_functions = [i for i in math_functions if i in _known_functions_math.values()]
    math_functions = list(set(math_functions)) # remove redundancies

    math_constants = [str(i) for i in preorder_traversal(expr) if isinstance(i, NumberSymbol)]
    math_constants = [i for i in math_constants if i in _known_constants_math.values()]
    math_constants = list(set(math_constants)) # remove redundancies

    return math_functions + math_constants
