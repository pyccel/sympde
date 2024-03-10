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

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    vectors = [b1,b2,b3]

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

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    vectors = [b1,b2,b3]

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

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    vectors = [b1,b2,b3]

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

    c1 = constant('c1', dtype=float)
    c2 = constant('c2', dtype=float)
    c3 = constant('c3', dtype=float)

    b1 = constant('b1', shape=dim, dtype=float)
    b2 = constant('b2', shape=dim, dtype=float)
    b3 = constant('b3', shape=dim, dtype=float)

    vectors = [b1,b2,b3]

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
# CLEAN UP SYMPY NAMESPACE
#==============================================================================
def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
