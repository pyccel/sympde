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
from sympde.topology import ScalarField, VectorField
from sympde.topology import ProductSpace
from sympde.topology import ScalarTestFunction
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
from sympde.expr import find

from sympde.expr.errors import UnconsistentLhsError
from sympde.expr.errors import UnconsistentRhsError
from sympde.expr.errors import UnconsistentBCError

from sympde.expr import Equation, EssentialBC, NewtonIteration

DIM = 2
domain = Domain('Omega', dim=DIM)

#==============================================================================
def test_find_2d_1():

    V = FunctionSpace('V', domain)
    U = FunctionSpace('U', domain)

    v = ScalarTestFunction(V, name='v')
    u = ScalarTestFunction(U, name='u')

    x,y = domain.coordinates

    alpha = Constant('alpha')
    kappa = Constant('kappa', real=True)
    eps   = Constant('eps', real=True)

    B1 = Boundary(r'\Gamma_1', domain)
    B2 = Boundary(r'\Gamma_2', domain)
    B3 = Boundary(r'\Gamma_3', domain)

    # ... bilinear/linear forms
    expr = dot(grad(v), grad(u))
    a1 = BilinearForm((v,u), expr)

    expr = v*u
    a2 = BilinearForm((v,u), expr)

    expr = v*trace_0(u, B1)
    a_B1 = BilinearForm((v, u), expr)

    expr = x*y*v
    l1 = LinearForm(v, expr)

    expr = cos(x+y)*v
    l2 = LinearForm(v, expr)

    expr = x*y*trace_0(v, B2)
    l_B2 = LinearForm(v, expr)
    # ...

    # ...
    equation = find(u, forall=v, lhs=a1(u,v), rhs=l1(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + alpha*a2(v,u))
    equation = find(u, forall=v, lhs=a(u,v), rhs=l1(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    equation = find(u, forall=v, lhs=a(u,v), rhs=l1(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    l = LinearForm(v, l1(v) + l2(v))
    equation = find(u, forall=v, lhs=a(u,v), rhs=l(v))
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    l = LinearForm(v, l1(v) + alpha*l_B2(v))
    equation = find(u, forall=v, lhs=a(u,v), rhs=l(v))
    # ...

    # ... using bc
    equation = find(u, forall=v, lhs=a(u,v), rhs=l(v), bc=EssentialBC(u, 0, B1))
    # ...

##==============================================================================
# TODO shall we have 'find' with Newton?
#def test_find_newton_2d_1():
#
#    # ... abstract model
#    B1 = Boundary(r'\Gamma_1', domain)
#
#    V = FunctionSpace('V', domain)
#
#    x,y = domain.coordinates
#
#    Un = ScalarField(V, name='Un')
#
#    v = ScalarTestFunction(V, name='v')
#
#    f  = -4.*exp(-Un)
#    l = LinearForm(v, dot(grad(v), grad(Un)) - f*v )
#
#    eq = NewtonIteration(l, Un, trials='u')
#    # ...
#
#    u = ScalarTestFunction(V, name='u')
#
#    # ...
#    expected =  -4.0*u*v*exp(-Un) + dx(u)*dx(v) + dy(u)*dy(v)
#
#    expr = TerminalExpr(eq.lhs)[0]
#    assert(expr.expr == expected)
#    # ...
#
#    # ...
#    expected = -4.0*v*exp(-Un) - dx(Un)*dx(v) - dy(Un)*dy(v)
#
#    expr = TerminalExpr(eq.rhs)[0]
#    assert(expr.expr == expected)
#    # ...
#
#    # ...
#    bc = EssentialBC(u, 0, B1)
#
#    # can be called with rhs = 0
#    equation = find(u, forall=v, lhs=l, rhs=0, fields=Un, bc=bc)
#    print(equation)
#
#    # or without rhs
#    equation = find(u, forall=v, lhs=l, fields=Un, bc=bc)
#    print(equation)
#    # ...

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
