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
def test_equation_2d_1():

    V = FunctionSpace('V', domain)
    U = FunctionSpace('U', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

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

#    # ...
#    with pytest.raises(UnconsistentLhsError):
#        equation = Equation(a1, l1(v))
#
#    with pytest.raises(UnconsistentLhsError):
#        equation = Equation(l1(v), l1(v))
#
#    with pytest.raises(UnconsistentLhsError):
#        equation = Equation(a1(v,u) + alpha*a2(v,u), l1(v))
#    # ...
#
#    # ...
#    with pytest.raises(UnconsistentRhsError):
#        equation = Equation(a1(v,u), l1)
#
#    with pytest.raises(UnconsistentRhsError):
#        equation = Equation(a1(v,u), a1(v,u))
#
#    with pytest.raises(UnconsistentRhsError):
#        equation = Equation(a1(v,u), l1(v) + l2(v))
#    # ...

    # ...
    equation = Equation(a1, l1, tests=v, trials=u)
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + alpha*a2(v,u))
    equation = Equation(a, l1, tests=v, trials=u)
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    equation = Equation(a, l1, tests=v, trials=u)
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    l = LinearForm(v, l1(v) + l2(v))
    equation = Equation(a, l, tests=v, trials=u)
    # ...

    # ...
    a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
    l = LinearForm(v, l1(v) + alpha*l_B2(v))
    equation = Equation(a, l, tests=v, trials=u)
    # ...

    # ... using bc
    equation = Equation(a1, l1, tests=v, trials=u, bc=EssentialBC(u, 0, B1))
    # ...

    # ... Poisson with Nitsch method
    g = cos(pi*x)*cos(pi*y)
    a0 = BilinearForm((u,v), dot(grad(u),grad(v)))
    a_B1 = BilinearForm((u,v), - kappa * u*trace_1(grad(v), B1)
                               - v*trace_1(grad(u), B1)
                               + trace_0(u, B1) * trace_0(v, B1) / eps)
    a = BilinearForm((u,v), a0(u,v) + a_B1(u,v))
    l = LinearForm(v, g*v/eps - kappa * trace_1(v, B1) * g)
    equation = Equation(a, l, tests=v, trials=u)
    # ...

#    # ... TODO FIX THIS, NOT RAISED ANYMORE!
#    with pytest.raises(UnconsistentBCError):
#        a = BilinearForm((v,u), a1(v,u) + a_B1(v,u))
#        equation = Equation(a(v,u), l1(v), bc=EssentialBC(u,0,B1))
#    # ...

#    # ... TODO FIX THIS, NOT RAISED ANYMORE!
#    # ...
#    with pytest.raises(UnconsistentBCError):
#        l = LinearForm(v, l1(v) + alpha*l_B2(v))
#        equation = Equation(a1(v,u), l(v), bc=EssentialBC(u,0,B2))
#    # ...


##==============================================================================
## TODO
#def test_projection_2d():
#
#    V = FunctionSpace('V', domain)
#    x,y = domain.coordinates
#
#    alpha = Constant('alpha')
#
#    u = Projection(x**2+alpha*y, V)


#==============================================================================
def test_equation_2d_2():

    domain = Square()

    V = FunctionSpace('V', domain)

    x,y = domain.coordinates

    pn, wn = [Field(V, name=i) for i in ['pn', 'wn']]

    dp    = TestFunction(V, name='dp')
    dw    = TestFunction(V, name='dw')
    tau   = TestFunction(V, name='tau')
    sigma = TestFunction(V, name='sigma')

    Re    = Constant('Re', real=True)
    dt    = Constant('dt', real=True)
    alpha = Constant('alpha', real=True)

    s  = BilinearForm((tau,sigma), dot(grad(tau), grad(sigma)))
    m  = BilinearForm((tau,sigma), tau*sigma)
    b1 = BilinearForm((tau,dw), bracket(pn, dw) * tau)
    b2 = BilinearForm((tau,dp), bracket(dp, wn) * tau)

    l1 = LinearForm(tau, bracket(pn, wn)*tau - 1./Re * dot(grad(tau), grad(wn)))

    expr =  m(tau,dw) - alpha*dt*b1(tau,dw) - dt*b2(tau,dp) - (alpha*dt/Re)*s(tau,dw)
    a = BilinearForm(((tau, sigma),(dp,dw)), expr)

    l = LinearForm((tau, sigma), dt*l1(tau))

    bc  = [EssentialBC(dp, 0, domain.boundary)]
    bc += [EssentialBC(dw, 0, domain.boundary)]
    equation = Equation(a, l, tests=[tau, sigma], trials=[dp,dw], bc=bc)

#    # TODO not working yet!! gives the wrong result => result must be a vector
#    # and not a scalar
#    print(evaluate(l, verbose=True))
#    print(evaluate(equation.rhs.expr, verbose=True))

#==============================================================================
def test_equation_2d_3():

    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    x,y = domain.coordinates

    B1 = Boundary(r'\Gamma_1', domain)

    # ... bilinear/linear forms
    a1 = BilinearForm((v,u), dot(grad(v), grad(u)))
    a2 = BilinearForm((v,u), v*u)

    l1 = LinearForm(v, x*y*v)
    l2 = LinearForm(v, cos(x+y)*v)
    # ...

    # ...
    bc = EssentialBC(u, 0, B1)
    eq = Equation(a1, l1, tests=v, trials=u, bc=bc)
    # ...

    # ...
    nn = NormalVector('nn')
    bc = EssentialBC(dot(grad(u), nn), 0, B1)
    eq = Equation(a1, l1, tests=v, trials=u, bc=bc)
    # ...

#==============================================================================
def test_equation_2d_4():

    V = VectorFunctionSpace('V', domain)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')

    x,y = domain.coordinates

    B1 = Boundary(r'\Gamma_1', domain)

    # ... bilinear/linear forms
    a1 = BilinearForm((v,u), inner(grad(v), grad(u)))

    f = Tuple(x*y, sin(pi*x)*sin(pi*y))
    l1 = LinearForm(v, dot(f,v))
    # ...

    # ...
    bc = EssentialBC(u, 0, B1)
    eq = Equation(a1, l1, tests=v, trials=u, bc=bc)
    # ...

    # ...
    bc = EssentialBC(u[0], 0, B1)
    eq = Equation(a1, l1, tests=v, trials=u, bc=bc)
    # ...

    # ...
    nn = NormalVector('nn')
    bc = EssentialBC(dot(u, nn), 0, B1)
    eq = Equation(a1, l1, tests=v, trials=u, bc=bc)
    # ...


#==============================================================================
def test_equation_2d_5():
    domain = Square()
    x,y = domain.coordinates

    f0 = Tuple(2*pi**2*sin(pi*x)*sin(pi*y),
              2*pi**2*sin(pi*x)*sin(pi*y))

    f1 = cos(pi*x)*cos(pi*y)

    W = VectorFunctionSpace('W', domain)
    V = FunctionSpace('V', domain)
    X = ProductSpace(W, V)

    F = VectorField(W, name='F')
    G = Field(V, name='G')

    u,v = [VectorTestFunction(W, name=i) for i in ['u', 'v']]
    p,q = [      TestFunction(V, name=i) for i in ['p', 'q']]

    a0 = BilinearForm((v,u), inner(grad(v), grad(u)))
    print('     a0 done.')
    a1 = BilinearForm((q,p), p*q)
    print('     a1 done.')
    a  = BilinearForm(((v,q),(u,p)), a0(v,u) + a1(q,p))
    print('     a  done.')

    l0 = LinearForm(v, dot(f0, v))
    l1 = LinearForm(q, f1*q)
    l  = LinearForm((v,q), l0(v) + l1(q))

    print('****************************')
    bc = EssentialBC(u, 0, domain.boundary)
    equation = Equation(a, l, tests=[v,q], trials=[u,p], bc=bc)

    # ...
    print('=======')
    print(equation.lhs.expr)
    print('')
    # ...

    # ...
    print('=======')
    print(equation.rhs.expr)
    print('')
    # ...

#==============================================================================
def test_equation_2d_6():

    domain = Square()
    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)

    b1 = 1.
    b2 = 0.
    b = Tuple(b1, b2)

    # right hand side
    f = x*y

    e = ElementDomain()
    area = Area(e)

    V = FunctionSpace('V', domain)

    u,v = [TestFunction(V, name=i) for i in ['u', 'v']]

    # ...
    expr = kappa * dot(grad(u), grad(v)) + dot(b, grad(u)) * v
    a = BilinearForm((v,u), expr)
    # ...

    # ...
    expr = f * v
    l0 = LinearForm(v, expr)
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * dot(b, grad(v))
    s1 = BilinearForm((v,u), expr)

    expr = - f * dot(b, grad(v))
    l1 = LinearForm(v, expr)
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * ( dot(b, grad(v)) - kappa * laplace(v))
    s2 = BilinearForm((v,u), expr)

    expr = - f * ( dot(b, grad(v)) - kappa * laplace(v))
    l2 = LinearForm(v, expr)
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * ( dot(b, grad(v)) + kappa * laplace(v))
    s3 = BilinearForm((v,u), expr)

    expr = - f * ( dot(b, grad(v)) + kappa * laplace(v))
    l3 = LinearForm(v, expr)
    # ...

    # ...
    expr = a(v,u) + mu*area*s1(v,u)
    a1 = BilinearForm((v,u), expr)
    # ...

    # ...
    expr = a(v,u) + mu*area*s2(v,u)
    a2 = BilinearForm((v,u), expr)
    # ...

    # ...
    expr = a(v,u) + mu*area*s3(v,u)
    a3 = BilinearForm((v,u), expr)
    # ...

    bc = EssentialBC(u, 0, domain.boundary)

    # ...
    l = LinearForm(v, l0(v) + mu*area*l1(v))

    eq_1 = Equation(a1, l, tests=v, trials=u, bc=bc)
    # ...

    # ...
    l = LinearForm(v, l0(v) + mu*area*l2(v))

    eq_2 = Equation(a2, l, tests=v, trials=u, bc=bc)
    # ...

    # ...
    l = LinearForm(v, l0(v) + mu*area*l3(v))

    eq_3 = Equation(a3, l, tests=v, trials=u, bc=bc)
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
