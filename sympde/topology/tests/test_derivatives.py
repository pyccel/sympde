# coding: utf-8

import pytest

from sympy import Tuple
from sympy import Matrix
from sympy import srepr
from sympy import Integer, Float, Rational, symbols
from sympy import expand, S, diff
from sympy import pi, sqrt
from sympy import sin, cos, tan
from sympy import cot, sinc, asin, acos, atan, acot, sec, csc, asec, acsc

from sympde.core import constant
from sympde.calculus import grad, dot, inner, outer, cross, rot, curl, div
from sympde.calculus import jump, avg, Dn, minus, plus
from sympde.topology import Domain, element_of
from sympde.topology import get_index_derivatives_atom
from sympde.topology import get_max_partial_derivatives
from sympde.topology import ScalarFunctionSpace
from sympde.topology import Mapping
from sympde.topology import (dx, dy, dz)
from sympde.topology import (dx1, dx2, dx3)
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import ProductSpace
from sympde.topology import element_of, elements_of

#==============================================================================
def indices_as_str(a):
    a = dict(sorted(a.items()))
    code = ''
    for k,n in list(a.items()):
        code += k*n
    return code

#==============================================================================
def test_partial_derivatives_1():

    dim = 3
    domain  = Domain('Omega', dim=dim)
    x,y,z = symbols('x,y,z')

    V  = ScalarFunctionSpace('V', domain)
    u1 = element_of(V, name='u1')
    u2 = element_of(V, name='u2')
    u3 = element_of(V, name='u3')

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    # ...
    for d in (dx,dy,dz,dx1,dx2,dx3):
        # ... constants
        assert d(1)              == 0  # native int
        assert d(2.3)            == 0  # native float
        # TODO
    #    assert d(4+5j)           == 0  # native complex
        assert d(Integer(1))     == 0  # sympy Integer
        assert d(Float(2.3))     == 0  # sympy Float
        assert d(Rational(6, 7)) == 0  # sympy Rational

        assert( d(c1) == 0 )
        assert( d(1/c1) == 0 )
        # ...

        assert d(u1+u2) == d(u1) + d(u2)
        assert d(c1*u1) == c1*d(u1)
        assert d(u1/c1) == d(u1)/c1
        assert d(c1*u1 + c2*u2) == c1*d(u1) + c2*d(u2)
        assert d(u1/c1 + u2/c2) == d(u1)/c1 + d(u2)/c2

        assert d(u1*u2) == u1*d(u2) + u2*d(u1)
        assert d(u1/u2) == -u1*d(u2)/u2**2 + d(u1)/u2

        assert d(u1*u2*u3) == u1*u2*d(u3) + u1*u3*d(u2) + u2*u3*d(u1)

        # TODO add more tests
        assert d(sin(u1)) == cos(u1)*d(u1)
        assert d(cos(u1)) == -sin(u1)*d(u1)
        assert d(tan(u1)) == (tan(u1)**2 + 1)*d(u1)

        assert d(cot(u1))  == (-1 - cot(u1)**2)*d(u1)
        assert d(sinc(u1)) == (cos(u1)/u1 - sin(u1)/u1**2)*d(u1)
        assert d(asin(u1)) == (1/sqrt(1 - u1**2))*d(u1)
        assert d(acos(u1)) == (-1/sqrt(1 - u1**2))*d(u1)
        assert d(atan(u1)) == (1/(1 + u1**2))*d(u1)
        assert d(acot(u1)) == (-1/(1 + u1**2))*d(u1)
        assert d(sec(u1))  == (tan(u1)*sec(u1))*d(u1)
        assert d(csc(u1))  == (-cot(u1)*csc(u1))*d(u1)
        assert d(asec(u1)) == (1/(u1**2*sqrt(1 - 1/u1**2)))*d(u1)
        assert d(acsc(u1)) == (-1/(u1**2*sqrt(1 - 1/u1**2)))*d(u1)
    # ...

#==============================================================================
def test_partial_derivatives_2():

    # ...
    x1, x2, x3 = symbols('x1, x2, x3')
    x, y, z = symbols('x, y, z')

    constants = ['c1', 'c2', 'c3', 'D', 'k', 'k1', 'k2', 'eps', 'b', 'R0', 'rmax', 'rmin']
    c1, c2, c3, D, k, k1, k2, eps, b, R0, rmax, rmin = [constant(i, dtype=float) for i in constants]

    expressions = []

    d = {}
    d['expr'] = c1 + (rmax*x1 + rmin*(-x1 + 1))*cos(x2)
    d['dx1'] = (rmax - rmin)*cos(x2)
    d['dx2'] = -(rmax*x1 + rmin*(-x1 + 1))*sin(x2)
    d['dx3'] = S.Zero
    expressions.append(d)

    d = {}
    d['expr'] = c2 + (rmax*x1 + rmin*(-x1 + 1))*sin(x2)
    d['dx1'] = (rmax - rmin)*sin(x2)
    d['dx2'] = (rmax*x1 + rmin*(-x1 + 1))*cos(x2)
    d['dx3'] = S.Zero
    expressions.append(d)

    d = {}
    d['expr'] = -D*x1**2 + c1 + x1*(-k + 1)*cos(x2)
    d['dx1'] = -2*D*x1 + (-k + 1)*cos(x2)
    d['dx2'] = -x1*(-k + 1)*sin(x2)
    d['dx3'] = S.Zero
    expressions.append(d)

    d = {}
    d['expr'] = c2 + x1*(k + 1)*sin(x2)
    d['dx1'] = (k + 1)*sin(x2)
    d['dx2'] = x1*(k + 1)*cos(x2)
    d['dx3'] = S.Zero
    expressions.append(d)

    d = {}
    d['expr'] = (-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 1)/eps
    d['dx1'] = -cos(x2)/sqrt(eps*(eps + 2*x1*cos(x2)) + 1)
    d['dx2'] = x1*sin(x2)/sqrt(eps*(eps + 2*x1*cos(x2)) + 1)
    d['dx3'] = S.Zero
    expressions.append(d)

    d = {}
    d['expr'] = b*x1*sin(x2)/(sqrt(-eps**2/4 + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)) + c2
    d['dx1'] = b*(eps*x1*sin(x2)*cos(x2)/(sqrt(-eps**2/4 + 1)*sqrt(eps*(eps + 2*x1*cos(x2)) + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)**2) + sin(x2)/(sqrt(-eps**2/4 + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)))
    d['dx2'] = b*x1*(-eps*x1*sin(x2)**2/(sqrt(eps*(eps + 2*x1*cos(x2)) + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)**2) + cos(x2)/(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2))/sqrt(-eps**2/4 + 1)
    d['dx3'] = S.Zero
    expressions.append(d)

    d = {}
    d['expr'] = 2.0*eps*sin(2.0*k1*pi*x1)*sin(2.0*k2*pi*x2) + 2.0*x1 - 1.0
    d['dx1'] = 4.0*eps*k1*pi*sin(2.0*k2*pi*x2)*cos(2.0*k1*pi*x1) + 2.0
    d['dx2'] = 4.0*eps*k2*pi*sin(2.0*k1*pi*x1)*cos(2.0*k2*pi*x2)
    d['dx3'] = S.Zero
    expressions.append(d)

    d = {}
    d['expr'] = 2.0*eps*sin(2.0*k1*pi*x1)*sin(2.0*k2*pi*x2) + 2.0*x2 - 1.0
    d['dx1'] = 4.0*eps*k1*pi*sin(2.0*k2*pi*x2)*cos(2.0*k1*pi*x1)
    d['dx2'] = 4.0*eps*k2*pi*sin(2.0*k1*pi*x1)*cos(2.0*k2*pi*x2) + 2.0
    d['dx3'] = S.Zero
    expressions.append(d)

    d = {}
    d['expr'] = (R0 + x1*cos(x2))*cos(x3)
    d['dx1'] = cos(x2)*cos(x3)
    d['dx2'] = -x1*sin(x2)*cos(x3)
    d['dx3'] = -(R0 + x1*cos(x2))*sin(x3)
    expressions.append(d)

    d = {}
    d['expr'] = (R0 + x1*cos(x2))*sin(x3)
    d['dx1'] = sin(x3)*cos(x2)
    d['dx2'] = -x1*sin(x2)*sin(x3)
    d['dx3'] = (R0 + x1*cos(x2))*cos(x3)
    expressions.append(d)

    d = {}
    d['expr'] = x1*sin(x2)
    d['dx1'] = sin(x2)
    d['dx2'] = x1*cos(x2)
    d['dx3'] = S.Zero
    expressions.append(d)

    d = {}
    d['expr'] = -D*x1**2 + c1 + x1*(-k + 1)*cos(x2)
    d['dx1'] = -2*D*x1 + (-k + 1)*cos(x2)
    d['dx2'] = -x1*(-k + 1)*sin(x2)
    d['dx3'] = S.Zero
    expressions.append(d)

    d = {}
    d['expr'] = c2 + x1*(k + 1)*sin(x2)
    d['dx1'] = (k + 1)*sin(x2)
    d['dx2'] = x1*(k + 1)*cos(x2)
    d['dx3'] = S.Zero
    expressions.append(d)

    d = {}
    d['expr'] = c3 + x1**2*x3*sin(2*x2)
    d['dx1'] = 2*x1*x3*sin(2*x2)
    d['dx2'] = 2*x1**2*x3*cos(2*x2)
    d['dx3'] = x1**2*sin(2*x2)
    expressions.append(d)
    # ...

    # ...
    for d in expressions:
        assert expand(dx1(d['expr'])) == expand(d['dx1'])
        assert expand(dx2(d['expr'])) == expand(d['dx2'])
        assert expand(dx3(d['expr'])) == expand(d['dx3'])

        assert expand(dx(d['expr'].subs(x1,x))) == expand(d['dx1'].subs(x1,x))
        assert expand(dy(d['expr'].subs(x2,y))) == expand(d['dx2'].subs(x2,y))
        assert expand(dz(d['expr'].subs(x3,z))) == expand(d['dx3'].subs(x3,z))
    # ...

#==============================================================================
def test_partial_derivatives_3():
    domain = Domain('Omega', dim=2)
    M      = Mapping('M', dim=2)

    mapped_domain = M(domain)

    x,y = mapped_domain.coordinates

    V = ScalarFunctionSpace('V', mapped_domain)

    F,u,v,w = [element_of(V, name=i) for i in ['F', 'u', 'v', 'w']]
    uvw = Tuple(u,v,w)

    alpha = constant('alpha', dtype=float)
    beta = constant('beta', dtype=float)

    assert(dx(x**2) == 2*x)
    assert(dy(x**2) == 0)
    assert(dz(x**2) == 0)

    assert(dx(y**2) == 0)
    assert(dy(y**2) == 2*y)
    assert(dz(y**2) == 0)

    assert(dx(x*F) == F + x*dx(F))
    assert(dx(uvw) == Matrix([[dx(u), dx(v), dx(w)]]))
    assert(dx(uvw) + dy(uvw) == Matrix([[dx(u) + dy(u),
                                         dx(v) + dy(v),
                                         dx(w) + dy(w)]]))

#    expected = Matrix([[alpha*dx(u) + beta*dy(u),
#                        alpha*dx(v) + beta*dy(v),
#                        alpha*dx(w) + beta*dy(w)]])
#    assert(alpha * dx(uvw) + beta * dy(uvw) == expected)
    # ...

#    expr = alpha * dx(uvw) + beta * dy(uvw)
#    print(expr)

#    print('> ', srepr(expr))
#    print('')

#==============================================================================
def test_partial_derivatives_4():
    domain = Domain('Omega', dim=2)
    M      = Mapping('M', dim=2)

    mapped_domain = M(domain)

    V = ScalarFunctionSpace('V', mapped_domain)
    F = element_of(V, name='F')

    alpha = constant('alpha', dtype=float)
    beta = constant('beta', dtype=float)

    # ...
    expr = alpha * dx(F)

    indices = get_index_derivatives_atom(expr, F)[0]
    assert(indices_as_str(indices) == 'x')
    # ...

    # ...
    expr = dy(dx(F))

    indices = get_index_derivatives_atom(expr, F)[0]
    assert(indices_as_str(indices) == 'xy')
    # ...

    # ...
    expr = alpha * dx(dy(dx(F)))

    indices = get_index_derivatives_atom(expr, F)[0]
    assert(indices_as_str(indices) == 'xxy')
    # ...

    # ...
    expr = alpha * dx(dx(F)) + beta * dy(F) + dx(dy(F))

    indices = get_index_derivatives_atom(expr, F)
    indices = [indices_as_str(i) for i in indices]
    assert(sorted(indices) == ['xx', 'xy', 'y'])
    # ...

    # ...
    expr = alpha * dx(dx(F)) + beta * dy(F) + dx(dy(F))

    d = get_max_partial_derivatives(expr, F)
    assert(indices_as_str(d) == 'xxy')

    d = get_max_partial_derivatives(expr)
    assert(indices_as_str(d) == 'xxy')
    # ...

#==============================================================================
@pytest.mark.skip(reason="Must be treated during the Error Messages issue")
def test_derivatives_1():
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
def test_derivatives_2():

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain, kind=None)

    u, v = elements_of(V, names='u, v')

    a = constant('a', dtype=float)

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

#==============================================================================
def test_derivatives_3():

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    V = VectorFunctionSpace('V', domain, kind=None)

    u, v = elements_of(V, names='u, v')

    a = constant('a', dtype=float)

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
@pytest.mark.parametrize('dim', [1, 2, 3])
def test_Grad(dim):

    domain  = Domain('Omega', dim=dim)

    V  = ScalarFunctionSpace('V', domain)
    u1 = element_of(V, name='u1')
    u2 = element_of(V, name='u2')
    u3 = element_of(V, name='u3')

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    # ... constants
    assert grad(1)              == 0  # native int
    assert grad(2.3)            == 0  # native float
    # TODO
#    assert grad(4+5j)           == 0  # native complex
    assert grad(Integer(1))     == 0  # sympy Integer
    assert grad(Float(2.3))     == 0  # sympy Float
    assert grad(Rational(6, 7)) == 0  # sympy Rational

    assert( grad(c1) == 0 )
    assert( grad(1/c1) == 0 )
    # ...

    # ...
    assert grad(u1+u2) == grad(u1) + grad(u2)
    assert grad(c1*u1) == c1*grad(u1)
    assert grad(u1/c1) == grad(u1)/c1
    assert grad(c1*u1 + c2*u2) == c1*grad(u1) + c2*grad(u2)
    assert grad(u1/c1 + u2/c2) == grad(u1)/c1 + grad(u2)/c2

    assert grad(u1*u2) == u1*grad(u2) + u2*grad(u1)
    assert grad(u1/u2) == -u1*grad(u2)/u2**2 + grad(u1)/u2

    assert grad(u1*u2*u3) == u1*u2*grad(u3) + u1*u3*grad(u2) + u2*u3*grad(u1)
    # ...

#==============================================================================
@pytest.mark.parametrize('dim', [2, 3])
def test_Curl(dim):

    domain  = Domain('Omega', dim=dim)

    V  = ScalarFunctionSpace('V', domain)
    W  = VectorFunctionSpace('V', domain)
    u1 = element_of(V, name='u1')
    u2 = element_of(V, name='u2')
    u3 = element_of(V, name='u3')
    F1 = element_of(W, name='F1')
    F2 = element_of(W, name='F2')
    F3 = element_of(W, name='F3')

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    # ... constants
    assert curl(1)              == 0  # native int
    assert curl(2.3)            == 0  # native float
    # TODO
#    assert curl(4+5j)           == 0  # native complex
    assert curl(Integer(1))     == 0  # sympy Integer
    assert curl(Float(2.3))     == 0  # sympy Float
    assert curl(Rational(6, 7)) == 0  # sympy Rational

    assert( curl(c1) == 0 )
    assert( curl(1/c1) == 0 )
    # ...

    # ...
    assert curl(c1*F1) == c1*curl(F1)
    assert curl(F1/c1) == curl(F1)/c1
    assert curl(F1+F2) == curl(F1) + curl(F2)
    assert curl(c1*F1 + c2*F2) == c1*curl(F1) + c2*curl(F2)
    assert curl(F1/c1 + F2/c2) == curl(F1)/c1 + curl(F2)/c2

    assert curl(u1*F1) == u1*curl(F1) + cross(grad(u1), F1)
    assert curl(F1/u1) == curl(F1)/u1 - cross(grad(u1), F1)/u1**2
    # ...

#==============================================================================
@pytest.mark.parametrize('dim', [2, 3])
def test_Div(dim):

    domain  = Domain('Omega', dim=dim)

    V  = ScalarFunctionSpace('V', domain)
    W  = VectorFunctionSpace('V', domain)
    u1 = element_of(V, name='u1')
    u2 = element_of(V, name='u2')
    u3 = element_of(V, name='u3')
    F1 = element_of(W, name='F1')
    F2 = element_of(W, name='F2')
    F3 = element_of(W, name='F3')

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    # ... constants
    assert div(1)              == 0  # native int
    assert div(2.3)            == 0  # native float
    # TODO
#    assert div(4+5j)           == 0  # native complex
    assert div(Integer(1))     == 0  # sympy Integer
    assert div(Float(2.3))     == 0  # sympy Float
    assert div(Rational(6, 7)) == 0  # sympy Rational

    assert( div(c1) == 0 )
    assert( div(1/c1) == 0 )
    # ...

    # ...
    assert div(c1*F1) == c1*div(F1)
    assert div(F1/c1) == div(F1)/c1
    assert div(F1+F2) == div(F1) + div(F2)
    assert div(c1*F1 + c2*F2) == c1*div(F1) + c2*div(F2)
    assert div(F1/c1 + F2/c2) == div(F1)/c1 + div(F2)/c2

    assert div(u1*F1) == u1*div(F1) + dot(grad(u1), F1)
    assert div(F1/u1) == div(F1)/u1 - dot(grad(u1), F1)/u1**2
    # ...

#==============================================================================
def test_Rot():

    dim = 2
    domain  = Domain('Omega', dim=dim)

    V  = ScalarFunctionSpace('V', domain)
    W  = VectorFunctionSpace('V', domain)
    u1 = element_of(V, name='u1')
    u2 = element_of(V, name='u2')
    u3 = element_of(V, name='u3')
    F1 = element_of(W, name='F1')
    F2 = element_of(W, name='F2')
    F3 = element_of(W, name='F3')

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    # ... constants
    assert rot(1)              == 0  # native int
    assert rot(2.3)            == 0  # native float
    # TODO
#    assert rot(4+5j)           == 0  # native complex
    assert rot(Integer(1))     == 0  # sympy Integer
    assert rot(Float(2.3))     == 0  # sympy Float
    assert rot(Rational(6, 7)) == 0  # sympy Rational

    assert( rot(c1) == 0 )
    assert( rot(1/c1) == 0 )
    # ...

    # ...
    assert rot(c1*F1) == c1*rot(F1)
    assert rot(F1/c1) == rot(F1)/c1
    assert rot(F1+F2) == rot(F1) + rot(F2)
    assert rot(c1*F1 + c2*F2) == c1*rot(F1) + c2*rot(F2)
    assert rot(F1/c1 + F2/c2) == rot(F1)/c1 + rot(F2)/c2

    assert rot(u1*F1) == u1*rot(F1) + cross(grad(u1), F1)
    assert rot(F1/u1) == rot(F1)/u1 - cross(grad(u1), F1)/u1**2
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
