# coding: utf-8

import pytest

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import IndexedBase
from sympy import Matrix
from sympy import Function
from sympy import pi, cos, sin
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

from sympde.expr import BilinearForm, LinearForm, Integral
from sympde.expr import atomize
from sympde.expr import evaluate
from sympde.expr import tensorize
from sympde.expr import Mass, Stiffness, Advection, AdvectionT
from sympde.expr import Projection
from sympde.expr import Norm
from sympde.expr import FormCall

from sympde.expr.errors import UnconsistentError
from sympde.expr.errors import UnconsistentLhsError
from sympde.expr.errors import UnconsistentRhsError
from sympde.expr.errors import UnconsistentBCError

from sympde.expr import Equation, EssentialBC


#==============================================================================
def test_essential_bc_1():
    domain = Domain('Omega', dim=2)

    V = FunctionSpace('V', domain)
    W = VectorFunctionSpace('W', domain)

    v = TestFunction(V, name='v')
    w = VectorTestFunction(W, name='w')

    B1 = Boundary(r'\Gamma_1', domain)
    nn = NormalVector('nn')

    # ... scalar case
    bc = EssentialBC(v, 0, B1)

    assert( bc.variable == v )
    assert( bc.order == 0 )
    assert( bc.normal_component == False )
    assert( bc.index_component == None )
    # ...

    # ... scalar case
    bc = EssentialBC(dot(grad(v), nn), 0, B1)

    assert( bc.variable == v )
    assert( bc.order == 1 )
    assert( bc.normal_component == False )
    assert( bc.index_component == None )
    # ...

    # ... vector case
    bc = EssentialBC(w, 0, B1)
    assert( bc.variable == w )
    assert( bc.order == 0 )
    assert( bc.normal_component == False )
    assert( bc.index_component == [0,1] )
    # ...

    # ... vector case
    bc = EssentialBC(dot(w, nn), 0, B1)

    assert( bc.variable == w )
    assert( bc.order == 0 )
    assert( bc.normal_component == True )
    assert( bc.index_component == None )
    # ...

    # ... vector case
    bc = EssentialBC(w[0], 0, B1)
    assert( bc.variable == w )
    assert( bc.order == 0 )
    assert( bc.normal_component == False )
    assert( bc.index_component == [0] )
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
