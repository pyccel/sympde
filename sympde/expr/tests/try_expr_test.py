# coding: utf-8

# TODO: - add assert to every test

import pytest

from sympy.core.containers import Tuple
from sympy import Function
from sympy import pi, cos, sin, exp
from sympy import ImmutableDenseMatrix as Matrix

from sympde.core     import Constant
from sympde.calculus import grad, dot, inner, rot, div
from sympde.calculus import laplace, bracket, convect
from sympde.calculus import jump, avg, Dn, minus, plus

from sympde.topology import dx1, dx2, dx3
from sympde.topology import dx, dy, dz
from sympde.topology import Mapping
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import element_of, elements_of
from sympde.topology import InteriorDomain, Union
from sympde.topology import Boundary, NormalVector
from sympde.topology import Domain
#from sympde.topology import trace_1  # TODO [YG, 27.01.2021]: fix trace
from sympde.topology import Square
from sympde.topology import ElementDomain
from sympde.topology import Area

from sympde.expr.expr import LinearExpr
from sympde.expr.expr import LinearForm, BilinearForm
from sympde.expr.expr import integral
from sympde.expr.expr import Functional, Norm
from sympde.expr.expr import linearize
from sympde.expr.evaluation import TerminalExpr

def try_terminal_expressions_for_navier_stokes():

    print("try this in a different file")

    domain = Square()
    x, y   = domain.coordinates

    mu = 1
    ux = cos(y*pi)
    uy = x*(x-1)
    ue = Matrix([[ux], [uy]])
    pe = sin(pi*y)
    # ...

    # ... Compute right-hand side
    a = TerminalExpr(-mu*laplace(ue), domain)
    b = TerminalExpr(    grad(ue), domain)
    c = TerminalExpr(    grad(pe), domain)

    # Verify that div(u) = 0
    assert (ux.diff(x) + uy.diff(y)).simplify() == 0

    d = TerminalExpr(div(ue), domain)
    assert d.simplify() == 0

    f = (a + b.T*ue + c).simplify()

    fx = -mu*(ux.diff(x, 2) + ux.diff(y, 2)) + ux*ux.diff(x) + uy*ux.diff(y) + pe.diff(x)
    fy = -mu*(uy.diff(x, 2) + uy.diff(y, 2)) + ux*uy.diff(x) + uy*uy.diff(y) + pe.diff(y)

    # [MCP 18.07.2024] for now, this test fails here because f is essentially 0: this should be fixed
    print(f'a = {a}')
    print(f'b = {b}')
    print(f'f = {f}')
    assert (f[0]-fx).simplify() == 0
    assert (f[1]-fy).simplify() == 0
    

if __name__ == '__main__':
    
    try_terminal_expressions_for_navier_stokes()