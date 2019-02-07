# coding: utf-8

import pytest

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import IndexedBase
from sympy import Matrix
from sympy import Function
from sympy import pi, cos, sin, exp
from sympy import srepr
from sympy.physics.quantum import TensorProduct

from sympde.core import Constant
from sympde.calculus import grad, dot, inner, cross, rot, curl, div
from sympde.calculus import laplace, hessian, bracket
from sympde.topology import (dx, dy, dz)
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import Field, VectorField
from sympde.topology import ProductSpace
from sympde.topology import TestFunction
from sympde.topology import VectorTestFunction
from sympde.topology import Unknown
from sympde.topology import Domain, Boundary, NormalVector, TangentVector
from sympde.topology import Trace, trace_0, trace_1
from sympde.topology import Square
from sympde.topology import ElementDomain
from sympde.topology import Area

from sympde.expr import BilinearForm, LinearForm
from sympde.expr import Projection
from sympde.expr import TerminalExpr

from sympde.expr.errors import UnconsistentLhsError
from sympde.expr.errors import UnconsistentRhsError
from sympde.expr.errors import UnconsistentBCError

from sympde.expr import Equation, EssentialBC, NewtonIteration

DIM = 2
domain = Domain('Omega', dim=DIM)

#==============================================================================
def test_newton_2d_1():

    # ... abstract model
    B1 = Boundary(r'\Gamma_1', domain)

    V = FunctionSpace('V', domain)

    x,y = domain.coordinates

    Un = Field(V, name='Un')

    v = TestFunction(V, name='v')

    f  = -4.*exp(-Un)
    l = LinearForm(v, dot(grad(v), grad(Un)) - f*v )

    eq = NewtonIteration(l, Un, trials='u')
    # ...

    u = TestFunction(V, name='u')

    # ...
    expected =  -4.0*u*v*exp(-Un) + dx(u)*dx(v) + dy(u)*dy(v)

    expr = TerminalExpr(eq.lhs)[0]
    assert(expr.expr == expected)
    # ...

    # ...
    expected = -4.0*v*exp(-Un) - dx(Un)*dx(v) - dy(Un)*dy(v)

    expr = TerminalExpr(eq.rhs)[0]
    assert(expr.expr == expected)
    # ...

    # ...
    bc = EssentialBC(u, 0, B1)
    eq = NewtonIteration(l, Un, bc=bc, trials=u)
    # ...


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
