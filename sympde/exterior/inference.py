# coding: utf-8


from numpy import unique

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
from sympde.core import LinearOperator
from sympde.core.basic import CalculusFunction

from .form import DifferentialForm
from .datatype import get_index_form
from .calculus import ExteriorDerivative, AdjointExteriorDerivative
from .calculus import ExteriorProduct
from .calculus import Hodge


#==============================================================================
def _get_dim(expr):
    if isinstance(expr, DifferentialForm):
        return expr.dim

    else:
        atoms = list(expr.atoms(DifferentialForm))
        if not atoms:
            raise ValueError('Cannot compute dim')

        return atoms[0].dim

#==============================================================================
def infere_type(expr):
    if isinstance(expr, DifferentialForm):
        return expr.index

    elif isinstance(expr, ExteriorDerivative):
        arg = expr.args[0]
        k = infere_type(arg)

        return get_index_form(k.index + 1)

    elif isinstance(expr, AdjointExteriorDerivative):
        arg = expr.args[0]
        k = infere_type(arg)

        return get_index_form(k.index - 1)

    elif isinstance(expr, ExteriorProduct):
        left = expr.args[0]
        right = expr.args[1]
        k_left = infere_type(left)
        k_right = infere_type(right)

        return get_index_form(k_left.index + k_right.index)

    elif isinstance(expr, Hodge):
        arg = expr.args[0]
        k = infere_type(arg)
        n = _get_dim(arg)

        return get_index_form(n-k.index)

    elif isinstance(expr, Add):
        indices = set([infere_type(i) for i in expr.args])
        indices = list(indices)
        if not( len(indices) == 1 ):
            raise ValueError('> Incompatible types. Found {}'.format(indices))

        return indices[0]

    return None
