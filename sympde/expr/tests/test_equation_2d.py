# coding: utf-8

#import pytest

from sympy.core.containers import Tuple
from sympy import pi, cos, sin
from sympy import ImmutableDenseMatrix as Matrix

from sympde.core        import Constant
from sympde.calculus    import grad, dot, inner
from sympde.calculus    import laplace, bracket
from sympde.topology    import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology    import element_of
from sympde.topology    import ProductSpace
from sympde.topology    import Domain, Boundary, NormalVector
from sympde.topology    import Square
from sympde.topology    import ElementDomain
from sympde.topology    import Area
from sympde.expr        import BilinearForm, LinearForm, integral
from sympde.expr        import Equation, EssentialBC
from sympde.expr.errors import UnconsistentLhsError
from sympde.expr.errors import UnconsistentRhsError
from sympde.expr.errors import UnconsistentBCError

DIM = 2
domain = Domain('Omega', dim=DIM)

#==============================================================================
def test_equation_2d_1():

    V = ScalarFunctionSpace('V', domain)
    U = ScalarFunctionSpace('U', domain)

    v = element_of(V, name='v')
    u = element_of(U, name='u')

    x,y = domain.coordinates

    alpha = Constant('alpha')
    kappa = Constant('kappa', real=True)
    eps   = Constant('eps', real=True)
    nn    = NormalVector('nn')

    B1 = Boundary(r'\Gamma_1', domain)
    B2 = Boundary(r'\Gamma_2', domain)
    B3 = Boundary(r'\Gamma_3', domain)

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)
    int_2 = lambda expr: integral(B2, expr)
    
    # ... bilinear/linear forms
    expr = dot(grad(v), grad(u))
    a1 = BilinearForm((v,u), int_0(expr))

    expr = v*u
    a2 = BilinearForm((v,u),int_0(expr))

    expr = v*u
    a_B1 = BilinearForm((v, u), int_1(expr))

    expr = x*y*v
    l1 = LinearForm(v, int_0(expr))

    expr = cos(x+y)*v
    l2 = LinearForm(v, int_0(expr))

    expr = x*y*v
    l_B2 = LinearForm(v, int_2(expr))
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
    a0 = BilinearForm((u,v), int_0(dot(grad(u),grad(v))))
    a_B1 = BilinearForm((u,v), - int_1(kappa * u*dot(grad(v), nn)
                               - v*dot(grad(u), nn)
                               + u*v/eps))


    a = BilinearForm((u,v), a0(u,v) + a_B1(u,v))
    l = LinearForm(v, int_0(g*v/eps) - int_1(kappa*v*g))
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


#==============================================================================
def test_equation_2d_2():

    domain = Square()

    V = ScalarFunctionSpace('V', domain)

    x,y = domain.coordinates

    pn, wn = [element_of(V, name=i) for i in ['pn', 'wn']]

    dp    = element_of(V, name='dp')
    dw    = element_of(V, name='dw')
    tau   = element_of(V, name='tau')
    sigma = element_of(V, name='sigma')

    Re    = Constant('Re', real=True)
    dt    = Constant('dt', real=True)
    alpha = Constant('alpha', real=True)

    int_0 = lambda expr: integral(domain , expr)
    
    s  = BilinearForm((tau,sigma), int_0(dot(grad(tau), grad(sigma))))
    m  = BilinearForm((tau,sigma), int_0(tau*sigma))
    b1 = BilinearForm((tau,dw), int_0(bracket(pn, dw) * tau))
    b2 = BilinearForm((tau,dp), int_0(bracket(dp, wn) * tau))

    l1 = LinearForm(tau, int_0(bracket(pn, wn)*tau - 1./Re * dot(grad(tau), grad(wn))))

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

    V = ScalarFunctionSpace('V', domain)

    v = element_of(V, name='v')
    u = element_of(V, name='u')

    x,y = domain.coordinates

    B1 = Boundary(r'\Gamma_1', domain)
    
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    # ... bilinear/linear forms
    a1 = BilinearForm((v,u), int_0(dot(grad(v), grad(u))))
    a2 = BilinearForm((v,u), int_0(v*u))

    l1 = LinearForm(v, int_0(x*y*v))
    l2 = LinearForm(v, int_0(cos(x+y)*v))
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

    v = element_of(V, name='v')
    u = element_of(V, name='u')
    x,y = domain.coordinates

    B1 = Boundary(r'\Gamma_1', domain)
    
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)

    # ... bilinear/linear forms
    a1 = BilinearForm((v,u), int_0(inner(grad(v), grad(u))))

    f = Matrix([x*y, sin(pi*x)*sin(pi*y)])
    l1 = LinearForm(v, int_0(dot(f,v)))
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

    f0 = Matrix([2*pi**2*sin(pi*x)*sin(pi*y),
                 2*pi**2*sin(pi*x)*sin(pi*y)])

    f1 = cos(pi*x)*cos(pi*y)

    W = VectorFunctionSpace('W', domain)
    V = ScalarFunctionSpace('V', domain)
    X = ProductSpace(W, V)

    F = element_of(W, name='F')
    G = element_of(V, name='G')

    u,v = [element_of(W, name=i) for i in ['u', 'v']]
    p,q = [element_of(V, name=i) for i in ['p', 'q']]

    int_0 = lambda expr: integral(domain , expr)
    
    a0 = BilinearForm((v,u), int_0(inner(grad(v), grad(u))))
    print('     a0 done.')
    a1 = BilinearForm((q,p), int_0(p*q))
    print('     a1 done.')
    a  = BilinearForm(((v,q),(u,p)), a0(v,u) + a1(q,p))
    print('     a  done.')

    l0 = LinearForm(v, int_0(dot(f0, v)))
    l1 = LinearForm(q, int_0(f1*q))
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
    b = Matrix([b1, b2])

    # right hand side
    f = x*y

    e = ElementDomain()
    area = Area(e)

    V = ScalarFunctionSpace('V', domain)

    u,v = [element_of(V, name=i) for i in ['u', 'v']]

    int_0 = lambda expr: integral(domain , expr)
    
    # ...
    expr = kappa * dot(grad(u), grad(v)) + dot(b, grad(u)) * v
    a = BilinearForm((v,u), int_0(expr))
    # ...
    # ...
    expr = f * v
    l0 = LinearForm(v, int_0(expr))
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * dot(b, grad(v))
    s1 = BilinearForm((v,u), int_0(expr))

    expr = - f * dot(b, grad(v))
    l1 = LinearForm(v, int_0(expr))
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * ( dot(b, grad(v)) - kappa * laplace(v))
    s2 = BilinearForm((v,u), int_0(expr))

    expr = - f * ( dot(b, grad(v)) - kappa * laplace(v))
    l2 = LinearForm(v, int_0(expr))
    # ...

    # ...
    expr = (- kappa * laplace(u) + dot(b, grad(u))) * ( dot(b, grad(v)) + kappa * laplace(v))
    s3 = BilinearForm((v,u), int_0(expr))

    expr = - f * ( dot(b, grad(v)) + kappa * laplace(v))
    l3 = LinearForm(v, int_0(expr))
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
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()

