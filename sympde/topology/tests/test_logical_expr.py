# coding: utf-8

from sympy.core.containers import Tuple
from sympy import Matrix
from sympy.tensor import IndexedBase
from sympy import symbols, simplify, Symbol
from sympy import expand

from sympde.topology import Mapping, DetJacobian
from sympde.topology import Domain

from sympde.core import Constant
from sympde.calculus import grad, dot, inner, cross, rot, curl, div
from sympde.calculus import laplace, hessian, bracket, convect
from sympde.topology import (dx, dy, dz)
from sympde.topology import (dx1, dx2, dx3)
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import Field, VectorField
from sympde.topology import TestFunction
from sympde.topology import VectorTestFunction
from sympde.topology import LogicalExpr
from sympde.topology import SymbolicExpr

# ...
def test_logical_expr_2d_1():
    print('============ test_logical_expr_2d_1 ==============')

    rdim = 2

    M = Mapping('M', rdim)
    domain = Domain('Omega', dim=rdim)

    alpha = Constant('alpha')

    V = FunctionSpace('V', domain)

    u,v = [TestFunction(V, name=i) for i in ['u', 'v']]

    det_M = DetJacobian(M)
    det   = Symbol('det')

    # ...
    expr = 2*u + alpha*v
    expr = LogicalExpr(M, expr)
    print(expr)
    print('')
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(M, expr)
    print(expr.subs(det_M, det))
    print('')
    # ...

    # ...
    expr = dy(u)
    expr = LogicalExpr(M, expr)
    print(expr.subs(det_M, det))
    print('')
    # ...

    # ...
    expr = dx(det_M)
    expr = LogicalExpr(M, expr)
    expr = expr.subs(det_M, det)
    expr = expand(expr)
    print(expr)
    print('')
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(M, expr)
    print(expr.subs(det_M, det))
    print('')
    # ...

# ...
def test_symbolic_expr_2d_1():
    print('============ test_symbolic_expr_2d_1 ==============')

    rdim = 2

    M = Mapping('M', rdim)
    domain = Domain('Omega', dim=rdim)

    alpha = Constant('alpha')

    V = FunctionSpace('V', domain)

    u,v = [TestFunction(V, name=i) for i in ['u', 'v']]

    det_M = DetJacobian(M)
    det_M = SymbolicExpr(det_M)
    print('>>> ', det_M)
    det   = Symbol('det')

    # ...
    expr = u
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    print(expr)
    # ...

    # ...
    expr = dx1(u)
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    print(expr)
    # ...

    # ...
    expr = dx1(dx2(u))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    print(expr)
    # ...

    # ...
    expr = dx1(M[0])
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    print(expr)
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    print(expr)
    # ...

    # ...
    expr = dx(DetJacobian(M))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    print(expand(expr))
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    print(expand(expr))
    # ...

    # ...
    expr = dx(dx(dx(u)))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    print(expand(expr))
    # ...

    # ...
    expr = dx(dx(dx(dx(u))))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    print(expand(expr))
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


#test_logical_expr_2d_1()
test_symbolic_expr_2d_1()
