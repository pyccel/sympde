# coding: utf-8

import pytest

from sympy import Function
from sympy import Integer, Float, Rational
from sympy import expand

from sympde.core     import Constant
from sympde.calculus import grad, dot, inner, outer, cross, rot, curl, div
from sympde.calculus import laplace, hessian, bracket, convect, D, conv
from sympde.calculus import ArgumentTypeError
from sympde.calculus import jump, avg, Dn, minus, plus
from sympde.topology import Domain
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import ProductSpace
from sympde.topology import element_of, elements_of

#==============================================================================
@pytest.mark.parametrize('dim', [1, 2, 3])
def test_Dot(dim):

    domain  = Domain('Omega', dim=dim)
    W       = VectorFunctionSpace('W', domain)
    a, b, c = elements_of(W, names='a, b, c')
    r       = Constant('r')

    # Commutativity (symmetry)
    assert dot(a, b) == dot(b, a)

    # Bilinearity: vector addition
    assert dot(a, b + c) == dot(a, b) + dot(a, c)
    assert dot(a + b, c) == dot(a, c) + dot(b, c)

    # Bilinearity: scalar multiplication
    assert dot(a, r * b) == r * dot(a, b)
    assert dot(r * a, b) == r * dot(a, b)

    # Special case: null vector
    assert dot(a, 0) == 0
    assert dot(0, a) == 0

    # Special case: two arguments are the same
    assert dot(a, a).is_real
    assert dot(a, a).is_positive

    # TODO: check exceptions

#==============================================================================
@pytest.mark.parametrize('dim', [1, 2, 3])
def test_Cross(dim):

    domain  = Domain('Omega', dim=dim)
    W       = VectorFunctionSpace('W', domain)
    a, b, c = elements_of(W, names='a, b, c')
    r       = Constant('r')

    # Anti-commutativity (anti-symmetry)
    assert cross(a, b) == -cross(b, a)

    # Bilinearity: vector addition
    assert cross(a, b + c) == cross(a, b) + cross(a, c)
    assert cross(a + b, c) == cross(a, c) + cross(b, c)

    # Bilinearity: scalar multiplication
    assert cross(a, r * b) == r * cross(a, b)
    assert cross(r * a, b) == r * cross(a, b)

    # Special case: null vector
    assert cross(a, 0) == 0
    assert cross(0, a) == 0

    # Special case: two arguments are the same
    assert cross(a, a) == 0

    # TODO: check exceptions

#==============================================================================
@pytest.mark.parametrize('dim', [1, 2, 3])
def test_Inner(dim):

    domain  = Domain('Omega', dim=dim)
    W       = VectorFunctionSpace('W', domain)
    a, b, c = elements_of(W, names='a, b, c')
    r       = Constant('r')

    # Commutativity
    assert inner(a, b) == inner(b, a)

    # Bilinearity: vector addition
    assert inner(a, b + c) == inner(a, b) + inner(a, c)
    assert inner(a + b, c) == inner(a, c) + inner(b, c)

    # Bilinearity: scalar multiplication
    assert inner(a, r * b) == r * inner(a, b)
    assert inner(r * a, b) == r * inner(a, b)

    # Special case: null vector
    assert inner(a, 0) == 0
    assert inner(0, a) == 0

    # Special case: two arguments are the same
    assert inner(a, a).is_real
    assert inner(a, a).is_positive

    # TODO: check exceptions

#==============================================================================
@pytest.mark.parametrize('dim', [1, 2, 3])
def test_Outer(dim):

    domain  = Domain('Omega', dim=dim)
    W       = VectorFunctionSpace('W', domain)
    a, b, c = elements_of(W, names='a, b, c')
    r       = Constant('r')

    # Not commutative
    assert outer(a, b) != outer(b, a)

    # Bilinearity: vector addition
    assert outer(a, b + c) == outer(a, b) + outer(a, c)
    assert outer(a + b, c) == outer(a, c) + outer(b, c)

    # Bilinearity: scalar multiplication
    assert outer(a, r * b) == r * outer(a, b)
    assert outer(r * a, b) == r * outer(a, b)

    # Special case: null vector
    assert outer(a, 0) == 0
    assert outer(0, a) == 0

    # TODO: check exceptions

#==============================================================================
def test_zero_derivative():

    assert grad(1)              == 0  # native int
    assert grad(2.3)            == 0  # native float
    assert grad(4+5j)           == 0  # native complex
    assert grad(Integer(1))     == 0  # sympy Integer
    assert grad(Float(2.3))     == 0  # sympy Float
    assert grad(Rational(6, 7)) == 0  # sympy Rational
    assert grad(Constant('a'))  == 0  # sympde Constant

    assert laplace(1)              == 0  # native int
    assert laplace(2.3)            == 0  # native float
    assert laplace(4+5j)           == 0  # native complex
    assert laplace(Integer(1))     == 0  # sympy Integer
    assert laplace(Float(2.3))     == 0  # sympy Float
    assert laplace(Rational(6, 7)) == 0  # sympy Rational
    assert laplace(Constant('a'))  == 0  # sympde Constant

    assert hessian(1)              == 0  # native int
    assert hessian(2.3)            == 0  # native float
    assert hessian(4+5j)           == 0  # native complex
    assert hessian(Integer(1))     == 0  # sympy Integer
    assert hessian(Float(2.3))     == 0  # sympy Float
    assert hessian(Rational(6, 7)) == 0  # sympy Rational
    assert hessian(Constant('a'))  == 0  # sympde Constant

    # 2D convection of constant scalar field
    domain = Domain('Omega', dim=2)
    W = VectorFunctionSpace('W', domain)
    F = element_of(W, name='F')

    assert convect(F, 1)              == 0  # native int
    assert convect(F, 2.3)            == 0  # native float
    assert convect(F, 4+5j)           == 0  # native complex
    assert convect(F, Integer(1))     == 0  # sympy Integer
    assert convect(F, Float(2.3))     == 0  # sympy Float
    assert convect(F, Rational(6, 7)) == 0  # sympy Rational
    assert convect(F, Constant('a'))  == 0  # sympde Constant

    # 3D convection of constant scalar field
    domain = Domain('Omega', dim=3)
    Z = VectorFunctionSpace('Z', domain)
    G = element_of(Z, name='G')

    assert convect(G, 1)              == 0  # native int
    assert convect(G, 2.3)            == 0  # native float
    assert convect(G, 4+5j)           == 0  # native complex
    assert convect(G, Integer(1))     == 0  # sympy Integer
    assert convect(G, Float(2.3))     == 0  # sympy Float
    assert convect(G, Rational(6, 7)) == 0  # sympy Rational
    assert convect(G, Constant('a'))  == 0  # sympde Constant

    # Poisson's bracket in 2D
    domain = Domain('Omega', dim=2)
    V = ScalarFunctionSpace('V', domain)
    u = element_of(V, name='V')

    assert bracket(u, 1)              == 0  # native int
    assert bracket(u, 2.3)            == 0  # native float
    assert bracket(u, 4+5j)           == 0  # native complex
    assert bracket(u, Integer(1))     == 0  # sympy Integer
    assert bracket(u, Float(2.3))     == 0  # sympy Float
    assert bracket(u, Rational(6, 7)) == 0  # sympy Rational
    assert bracket(u, Constant('a'))  == 0  # sympde Constant

    assert bracket(1             , u) == 0  # native int
    assert bracket(2.3           , u) == 0  # native float
    assert bracket(4+5j          , u) == 0  # native complex
    assert bracket(Integer(1)    , u) == 0  # sympy Integer
    assert bracket(Float(2.3)    , u) == 0  # sympy Float
    assert bracket(Rational(6, 7), u) == 0  # sympy Rational
    assert bracket(Constant('a') , u) == 0  # sympde Constant

#==============================================================================
def test_calculus_2d_1():
    domain = Domain('Omega', dim=2)

    V = ScalarFunctionSpace('V', domain)
    W = VectorFunctionSpace('W', domain)

    alpha, beta, gamma = [Constant(i) for i in ['alpha','beta','gamma']]

    f, g, h = elements_of(V, names='f, g, h')
    F, G, H = elements_of(W, names='F, G, H')

    # ... scalar gradient properties
    assert( grad(f+g) == grad(f) + grad(g) )
    assert( grad(alpha*h) == alpha*grad(h) )
    assert( grad(alpha*f + beta*g) == alpha*grad(f) + beta*grad(g)  )

    assert( grad(f*g) == f*grad(g) + g*grad(f) )
    assert( grad(f/g) == -f*grad(g)/g**2 + grad(f)/g )

    assert( expand(grad(f*g*h)) == f*g*grad(h) + f*h*grad(g) + g*h*grad(f) )
    # ...

    # ... vector gradient properties
    assert( grad(F+G) == grad(F) + grad(G) )
    assert( grad(alpha*H) == alpha*grad(H) )
    assert( grad(alpha*F + beta*G) == alpha*grad(F) + beta*grad(G)  )

    #assert( grad(dot(F,G)) == convect(F, G) + convect(G, F) + cross(F, curl(G)) - cross(curl(F), G) )
    # ...

    # ... curl properties
    assert( curl(f+g) == curl(f) + curl(g) )
    assert( curl(alpha*h) == alpha*curl(h) )
    assert( curl(alpha*f + beta*g) == alpha*curl(f) + beta*curl(g)  )
    # ...

    # ... laplace properties
    assert( laplace(f+g) == laplace(f) + laplace(g) )
    assert( laplace(alpha*h) == alpha*laplace(h) )
    assert( laplace(alpha*f + beta*g) == alpha*laplace(f) + beta*laplace(g)  )
    # ...

    # ... divergence properties
    assert( div(F+G) == div(F) + div(G) )
    assert( div(alpha*H) == alpha*div(H) )
    assert( div(alpha*F + beta*G) == alpha*div(F) + beta*div(G)  )

    #assert( div(cross(F,G)) == -dot(F, curl(G)) + dot(G, curl(F)) )
    # ...

    # ... rot properties
    assert( rot(F+G) == rot(F) + rot(G) )
    assert( rot(alpha*H) == alpha*rot(H) )
    assert( rot(alpha*F + beta*G) == alpha*rot(F) + beta*rot(G)  )
    # ...

    # ... Poisson's bracket properties
    assert bracket(alpha * f, g) == alpha * bracket(f, g)
    assert bracket(f, alpha * g) == alpha * bracket(f, g)
    assert bracket(f + h, g) == bracket(f, g) + bracket(h, g)
    assert bracket(f, g + h) == bracket(f, g) + bracket(f, h)
    assert bracket(f, f) == 0
    assert bracket(f, g) == -bracket(g, f)
    assert bracket(f, g * h) == g * bracket(f, h) + bracket(f, g) * h
    # ...

#==============================================================================
def test_calculus_2d_2():
    domain = Domain('Omega', dim=2)

    W = VectorFunctionSpace('W', domain)

    alpha, beta, gamma = [Constant(i) for i in ['alpha','beta','gamma']]

    F, G, H = elements_of(W, names='F, G, H')

    # ... vector gradient properties
    assert( convect(F+G, H) == convect(F,H) + convect(G,H) )
    assert( convect(alpha*F,H) == alpha*convect(F,H) )
    assert( convect(F,alpha*H) == alpha*convect(F,H) )
    # ...

#==============================================================================
def test_calculus_2d_3():
    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    V = ScalarFunctionSpace('V', domain)

    alpha, beta, gamma = [Constant(i) for i in ['alpha','beta','gamma']]

    f, g = elements_of(V, names='f, g')
    K = Function('K')(x,y)
    S = Function('S')(x,y)

    # ... scalar gradient properties
    assert(conv(K, alpha*f+beta*g) == alpha*conv(K,f) + beta*conv(K,g))
    assert(conv(alpha*K+beta*S, f) == alpha*conv(K,f) + beta*conv(S,f))
    # ...

#==============================================================================
def test_calculus_3d():
    domain = Domain('Omega', dim=3)

    V = ScalarFunctionSpace('V', domain)
    W = VectorFunctionSpace('W', domain)

    alpha, beta, gamma = [Constant(i) for i in ['alpha','beta','gamma']]

    f, g, h = elements_of(V, names='f, g, h')
    F, G, H = elements_of(W, names='F, G, H')

    # ... gradient properties
    assert( grad(f+g) == grad(f) + grad(g) )
    assert( grad(alpha*h) == alpha*grad(h) )
    assert( grad(alpha*f + beta*g) == alpha*grad(f) + beta*grad(g)  )

    assert( grad(f*g) == f*grad(g) + g*grad(f) )
    assert( grad(f/g) == -f*grad(g)/g**2 + grad(f)/g )

    assert( expand(grad(f*g*h)) == f*g*grad(h) + f*h*grad(g) + g*h*grad(f) )
    # ...

    # ... curl properties
    assert( curl(F+G) == curl(F) + curl(G) )
    assert( curl(alpha*H) == alpha*curl(H) )
    assert( curl(alpha*F + beta*G) == alpha*curl(F) + beta*curl(G)  )

    #assert( curl(cross(F,G)) == F*div(G) - G*div(F) - convect(F, G) + convect(G, F) )
    #assert( curl(f*F) == f*curl(F) + cross(grad(f), F) )
    # ...

    # ... laplace properties
    assert( laplace(f+g) == laplace(f) + laplace(g) )
    assert( laplace(alpha*h) == alpha*laplace(h) )
    assert( laplace(alpha*f + beta*g) == alpha*laplace(f) + beta*laplace(g)  )
    #assert( laplace(f*g) == f*laplace(g) + g*laplace(f) + 2*dot(grad(f), grad(g)) )
    # ...

    # ... divergence properties
    assert( div(F+G) == div(F) + div(G) )
    assert( div(alpha*H) == alpha*div(H) )
    assert( div(alpha*F + beta*G) == alpha*div(F) + beta*div(G)  )

    #assert( div(cross(F,G)) == -dot(F, curl(G)) + dot(G, curl(F)) )
    #assert( div(f*F) == f*div(F) + dot(F, grad(f)))
    # ...

    # ...
    assert( curl(grad(h)) == 0 )
    assert( div(curl(H)) == 0 )
    #assert( div(cross(grad(f), grad(g))) == 0 )
    #assert( curl(curl(F)) == grad(div(F)) - laplace(F))
    #assert( curl(f*grad(g)) == cross(grad(f), grad(g)) )
    # ...


#==============================================================================
def test_calculus_3d_3():
    domain = Domain('Omega', dim=3)

    H1    = ScalarFunctionSpace('V0', domain, kind='H1')
    Hcurl = VectorFunctionSpace('V1', domain, kind='Hcurl')
    Hdiv  = VectorFunctionSpace('V2', domain, kind='Hdiv')
    L2    = ScalarFunctionSpace('V3', domain, kind='L2')
    V     = ScalarFunctionSpace('V', domain, kind=None)
    W     = VectorFunctionSpace('W', domain, kind=None)

    X = ProductSpace(H1, Hcurl, Hdiv, L2)
    v0, v1, v2, v3 = element_of(X, ['v0', 'v1', 'v2', 'v3'])

    v = element_of(V, 'v')
    w = element_of(W, 'w')

    # ... consistency of grad operator
    # we can apply grad on an undefined space type or H1
    expr = grad(v0)
    expr = grad(v)
    expr = grad(w)

    # wa cannot apply grad to a function in Hcurl
    with pytest.raises(ArgumentTypeError): expr = grad(v1)

    # wa cannot apply grad to a function in Hdiv
    with pytest.raises(ArgumentTypeError): expr = grad(v2)

    # wa cannot apply grad to a function in L2
    with pytest.raises(ArgumentTypeError): expr = grad(v3)
    # ...

    # ... consistency of curl operator
    # we can apply curl on an undefined space type, H1 or Hcurl
    expr = curl(v0)
    expr = curl(v1)
    expr = curl(v)
    expr = curl(w)

    # wa cannot apply curl to a function in Hdiv
    with pytest.raises(ArgumentTypeError): expr = curl(v2)

    # wa cannot apply curl to a function in L2
    with pytest.raises(ArgumentTypeError): expr = curl(v3)
    # ...

    # ... consistency of div operator
    # we can apply div on an undefined space type, H1 or Hdiv
    expr = div(v0)
    expr = div(v2)
    expr = div(v)
    expr = div(w)

    # wa cannot apply div to a function in Hcurl
    with pytest.raises(ArgumentTypeError): expr = div(v1)

    # wa cannot apply div to a function in L2
    with pytest.raises(ArgumentTypeError): expr = div(v3)
    # ...

#==============================================================================
def test_calculus_3d_4():
    domain = Domain('Omega', dim=3)

    W = VectorFunctionSpace('W', domain)

    alpha, beta, gamma = [Constant(i) for i in ['alpha','beta','gamma']]

    F, G, H = elements_of(W, names='F, G, H')

    # ...
    expected = alpha*inner(D(F), D(G)) + beta*inner(D(F), D(H))
    assert(inner(D(F), D(alpha*G+beta*H)) == expected)
    # ...

#==============================================================================
def test_calculus_3d_5():
    domain = Domain('Omega', dim=3)

    W = VectorFunctionSpace('W', domain)

    alpha, beta, gamma = [Constant(i) for i in ['alpha','beta','gamma']]

    F, G, H = elements_of(W, names='F, G, H')

    # ...
    expected = alpha*outer(F, G) + beta*outer(F, H)
    assert(outer(F, alpha*G+beta*H) == expected)
    # ...

#    # ... TODO
#    expr = outer(E*F, G*H)
#    assert(outer(E*F, G*H) == outer(E,G)*outer(F,H))
#    # ...

#==============================================================================
def test_calculus_2d_4():

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain, kind=None)

    u, v = elements_of(V, names='u, v')

    a = Constant('a', is_real=True)

    # ... jump operator
    assert(jump(u+v) == jump(u) + jump(v))
    assert(jump(a*u) == a*jump(u))
    # ...

    # ... avg operator
    assert(avg(u+v) == avg(u) + avg(v))
    assert(avg(a*u) == a*avg(u))
    # ...

    # ... Dn operator
    assert(Dn(u+v) == Dn(u) + Dn(v))
    assert(Dn(a*u) == a*Dn(u))
    # ...

    # ... minus operator
    assert(minus(u+v) == minus(u) + minus(v))
    assert(minus(a*u) == a*minus(u))
    # ...

    # ... plus operator
    assert(plus(u+v) == plus(u) + plus(v))
    assert(plus(a*u) == a*plus(u))
    # ...

def test_calculus_2d_5():

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    V = VectorFunctionSpace('V', domain, kind=None)

    u, v = elements_of(V, names='u, v')

    a = Constant('a', is_real=True)

    # ... jump operator
    assert(jump(u+v) == jump(u) + jump(v))
    assert(jump(a*u) == a*jump(u))
    # ...

    # ... avg operator
    assert(avg(u+v) == avg(u) + avg(v))
    assert(avg(a*u) == a*avg(u))
    # ...

    # ... Dn operator
    assert(Dn(u+v) == Dn(u) + Dn(v))
    assert(Dn(a*u) == a*Dn(u))
    # ...

    # ... minus operator
    assert(minus(u+v)      == minus(u) + minus(v))
    assert(minus(a*u)      == a*minus(u))
    assert(div(minus(u+v)) == div(minus(u)) + div(minus(v)))
    # ...

    # ... plus operator
    assert(plus(u+v)      == plus(u) + plus(v))
    assert(plus(a*u)      == a*plus(u))
    assert(div(plus(u+v)) == div(plus(u)) + div(plus(v)))
#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()

#test_calculus_2d_1()
#test_calculus_2d_2()
#test_calculus_2d_3()
#test_calculus_2d_4()
#test_calculus_3d()

#test_calculus_3d_3()
#test_calculus_3d_4()
#test_calculus_3d_5()
