# coding: utf-8

from sympy.core.compatibility import is_sequence
from sympy import Indexed, IndexedBase

from .basic import CalculusFunction

# ... generic operators
class GenericFunction(CalculusFunction):

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            return Indexed(self, *indices, **kw_args)
        else:
            return Indexed(self, indices, **kw_args)

class Dot(GenericFunction):
    pass

class Cross(GenericFunction):
    pass

class Grad(GenericFunction):
    pass

class Curl(GenericFunction):
    pass

class Rot(GenericFunction):
    pass

class Div(GenericFunction):
    pass

class Inner(GenericFunction):
    pass

class Bracket(GenericFunction):
    pass

Outer = Dot # TODO add the Outer class Function

# ...
_generic_ops  = (Dot, Cross, Inner, Outer,
                 Grad, Curl, Rot, Div,
                 Bracket)
# ...

# ... alias for ufl compatibility
cross = Cross
dot = Dot
inner = Inner
outer = Outer

grad = Grad
curl = Curl
rot = Rot
div = Div

bracket = Bracket
# ...
