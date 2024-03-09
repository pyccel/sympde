# coding: utf-8

import pytest

from sympy import Function
from sympy import Integer, Float, Rational
from sympy import expand

from sympde.core     import constant
from sympde.calculus import grad, dot, inner, outer, cross, rot, curl, div
from sympde.calculus import laplace, hessian, bracket, convect, conv
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

    V  = VectorFunctionSpace('V', domain)
    u1 = element_of(V, name='u1')
    u2 = element_of(V, name='u2')
    u3 = element_of(V, name='u3')

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    vectors = [u1,u2,u3,b1,b2,b3]

    for i1,e1 in enumerate(vectors):
        # Special case: null vector
        assert dot(e1, 0) == 0
        assert dot(0, e1) == 0

#        # Special case: two arguments are the same
#        assert dot(e1, e1).is_real
#        assert dot(e1, e1).is_positive

        for i2,e2 in enumerate(vectors[i1+1:]):
            # Commutativity (symmetry)
            assert dot(e1, e2) == dot(e2, e1)

            # Bilinearity: scalar multiplication
            assert dot(c1 * e1, e2) == c1 * dot(e1, e2)
            assert dot(e1, c2 * e2) == c2 * dot(e1, e2)
            assert dot(c1 * e1, c2 * e2) == c1 * c2 * dot(e1, e2)

            for i3,e3 in enumerate(vectors[i1+1+i2+1:]):
                # Bilinearity: vector addition
                assert dot(e1, e2 + e3) == dot(e1, e2) + dot(e1, e3)
                assert dot(e1 + e2, e3) == dot(e1, e3) + dot(e2, e3)

#==============================================================================
@pytest.mark.parametrize('dim', [1, 2, 3])
def test_Cross(dim):

    domain  = Domain('Omega', dim=dim)

    V  = VectorFunctionSpace('V', domain)
    u1 = element_of(V, name='u1')
    u2 = element_of(V, name='u2')
    u3 = element_of(V, name='u3')

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    vectors = [u1,u2,u3,b1,b2,b3]

    for i1,e1 in enumerate(vectors):
        # Special case: null vector
        assert cross(e1, 0) == 0
        assert cross(0, e1) == 0

        # Special case: two arguments are the same
        assert cross(e1, e1) == 0

        for i2,e2 in enumerate(vectors[i1+1:]):
            # Anti-commutativity (anti-symmetry)
            assert cross(e1, e2) == -cross(e2, e1)

            # Bilinearity: scalar multiplication
            assert cross(c1 * e1, e2) == c1 * cross(e1, e2)
            assert cross(e1, c2 * e2) == c2 * cross(e1, e2)
            assert cross(c1 * e1, c2 * e2) == c1 * c2 * cross(e1, e2)

            for i3,e3 in enumerate(vectors[i1+1+i2+1:]):
                # Bilinearity: vector addition
                assert cross(e1, e2 + e3) == cross(e1, e2) + cross(e1, e3)
                assert cross(e1 + e2, e3) == cross(e1, e3) + cross(e2, e3)

#==============================================================================
@pytest.mark.parametrize('dim', [1, 2, 3])
def test_Inner(dim):

    domain  = Domain('Omega', dim=dim)

    V  = VectorFunctionSpace('V', domain)
    u1 = element_of(V, name='u1')
    u2 = element_of(V, name='u2')
    u3 = element_of(V, name='u3')

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    vectors = [u1,u2,u3,b1,b2,b3]

    for i1,e1 in enumerate(vectors):
        # Special case: null vector
        assert inner(e1, 0) == 0
        assert inner(0, e1) == 0

        for i2,e2 in enumerate(vectors[i1+1:]):
            # Commutativity
            assert inner(e1, e2) == inner(e2, e1)

            # Bilinearity: scalar multiplication
            assert inner(c1 * e1, e2) == c1 * inner(e1, e2)
            assert inner(e1, c2 * e2) == c2 * inner(e1, e2)
            assert inner(c1 * e1, c2 * e2) == c1 * c2 * inner(e1, e2)

            for i3,e3 in enumerate(vectors[i1+1+i2+1:]):
                # Bilinearity: vector addition
                assert inner(e1, e2 + e3) == inner(e1, e2) + inner(e1, e3)
                assert inner(e1 + e2, e3) == inner(e1, e3) + inner(e2, e3)

#==============================================================================
@pytest.mark.parametrize('dim', [1, 2, 3])
def test_Outer(dim):

    domain  = Domain('Omega', dim=dim)

    V  = VectorFunctionSpace('V', domain)
    u1 = element_of(V, name='u1')
    u2 = element_of(V, name='u2')
    u3 = element_of(V, name='u3')

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    vectors = [u1,u2,u3,b1,b2,b3]

    for i1,e1 in enumerate(vectors):
        # Special case: null vector
        assert outer(e1, 0) == 0
        assert outer(0, e1) == 0

        for i2,e2 in enumerate(vectors[i1+1:]):
            # Not commutative
            assert outer(e1, e2) != outer(e2, e1)

            # Bilinearity: scalar multiplication
            assert outer(c1 * e1, e2) == c1 * outer(e1, e2)
            assert outer(e1, c2 * e2) == c2 * outer(e1, e2)
            assert outer(c1 * e1, c2 * e2) == c1 * c2 * outer(e1, e2)

            for i3,e3 in enumerate(vectors[i1+1+i2+1:]):
                # Bilinearity: vector addition
                assert outer(e1, e2 + e3) == outer(e1, e2) + outer(e1, e3)
                assert outer(e1 + e2, e3) == outer(e1, e3) + outer(e2, e3)

#==============================================================================
@pytest.mark.skip(reason="Must be treated during the Error Messages issue")
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
def test_calculus_2d_4():

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
def test_calculus_2d_5():

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
# CLEAN UP SYMPY NAMESPACE
#==============================================================================
def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
