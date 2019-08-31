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
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import ScalarTestFunction
from sympde.topology import VectorTestFunction
from sympde.topology import LogicalExpr
from sympde.topology import SymbolicExpr
from sympde.topology import IdentityMapping
from sympde.topology import element_of


# ...
def test_logical_expr_1d_1():
    #print('============ test_logical_expr_1d_1 ==============')

    rdim = 1

    M = Mapping('M', rdim)
    domain = Domain('Omega', dim=rdim)

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)

    u,v = [ScalarTestFunction(V, name=i) for i in ['u', 'v']]

    det_M = DetJacobian(M)
    #print('det = ', det_M)
    det   = Symbol('det')

    # ...
    expr = 2*u + alpha*v
    expr = LogicalExpr(M, expr)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(M, expr)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dx(det_M)
    expr = LogicalExpr(M, expr)
    expr = expr.subs(det_M, det)
    expr = expand(expr)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(M, expr)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

# ...
def test_symbolic_expr_1d_1():
    #print('============ test_symbolic_expr_1d_1 ==============')

    rdim = 1

    M = Mapping('M', rdim)
    domain = Domain('Omega', dim=rdim)

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)

    u,v = [ScalarTestFunction(V, name=i) for i in ['u', 'v']]

    det_M = DetJacobian(M)
    det_M = SymbolicExpr(det_M)
    #print('>>> ', det_M)
    det   = Symbol('det')

    # ...
    expr = u
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(u)
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(M[0])
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expr)
    # ...

    # ...
    expr = dx(DetJacobian(M))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(dx(u)))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...



# ...
def test_logical_expr_2d_1():
    #print('============ test_logical_expr_2d_1 ==============')

    rdim = 2

    M = Mapping('M', rdim)
    domain = Domain('Omega', dim=rdim)

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)
    W = VectorFunctionSpace('V', domain)

    u,v = [ScalarTestFunction(V, name=i) for i in ['u', 'v']]
    w = VectorTestFunction(W, name='w')

    det_M = DetJacobian(M)
    #print('det = ', det_M)
    det   = Symbol('det')

    # ...
    expr = 2*u + alpha*v
    expr = LogicalExpr(M, expr)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(M, expr)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dy(u)
    expr = LogicalExpr(M, expr)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dx(det_M)
    expr = LogicalExpr(M, expr)
    expr = expr.subs(det_M, det)
    expr = expand(expr)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(M, expr)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dx(w[0])
    expr = LogicalExpr(M, expr)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

# ...
def test_symbolic_expr_2d_1():
    #print('============ test_symbolic_expr_2d_1 ==============')

    rdim = 2

    M = Mapping('M', rdim)
    domain = Domain('Omega', dim=rdim)

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)

    u,v = [ScalarTestFunction(V, name=i) for i in ['u', 'v']]

    det_M = DetJacobian(M)
    det_M = SymbolicExpr(det_M)
    #print('>>> ', det_M)
    det   = Symbol('det')

    # ...
    expr = u
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(u)
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(dx2(u))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(M[0])
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expr)
    # ...

    # ...
    expr = dx(DetJacobian(M))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(dx(u)))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

# ...
def test_logical_expr_3d_1():
    #print('============ test_logical_expr_3d_1 ==============')

    rdim = 3

    M = Mapping('M', rdim)
    domain = Domain('Omega', dim=rdim)

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)

    u,v = [ScalarTestFunction(V, name=i) for i in ['u', 'v']]

    det_M = DetJacobian(M)
    #print('det = ', det_M)
    det   = Symbol('det')

    # ...
    expr = 2*u + alpha*v
    expr = LogicalExpr(M, expr)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(M, expr)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dy(u)
    expr = LogicalExpr(M, expr)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dx(det_M)
    expr = LogicalExpr(M, expr)
    expr = expr.subs(det_M, det)
    expr = expand(expr)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(M, expr)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

# ...
def test_symbolic_expr_3d_1():
    #print('============ test_symbolic_expr_3d_1 ==============')

    rdim = 3

    M = Mapping('M', rdim)
    domain = Domain('Omega', dim=rdim)

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)

    u,v = [ScalarTestFunction(V, name=i) for i in ['u', 'v']]

    det_M = DetJacobian(M)
    det_M = SymbolicExpr(det_M)
    #print('>>> ', det_M)
    det   = Symbol('det')

    # ...
    expr = u
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(u)
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(dx2(u))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(M[0])
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expr)
    # ...

    # ...
    expr = dx(DetJacobian(M))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(dx(u)))
    expr = LogicalExpr(M, expr)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))

    # ...

#==============================================================================
def test_identity_mapping_2d_1():
    rdim = 2

    x1, x2 = symbols('x1, x2')

    M = IdentityMapping('M', rdim)

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))

    assert(LogicalExpr(M, M[0]) == x1)
    assert(LogicalExpr(M, M[1]) == x2)

    assert(LogicalExpr(M, dx1(M[0])) == 1)
    assert(LogicalExpr(M, dx1(M[1])) == 0)

    assert(LogicalExpr(M, dx2(M[0])) == 0)
    assert(LogicalExpr(M, dx2(M[1])) == 1)

    expected = Matrix([[1, 0], [0, 1]])
    assert(not( M.jacobian == expected))
    assert(LogicalExpr(M, M.jacobian) == expected)


#==============================================================================
def test_identity_mapping_2d_2():
    rdim = 2

    x1, x2 = symbols('x1, x2')

    domain = Domain('Omega', dim=rdim)
    M = IdentityMapping('F', rdim)

    V = ScalarFunctionSpace('V', domain)
    u = element_of(V, name='u')

    # ...
    assert(LogicalExpr(M, dx(u)) == dx1(u))
    assert(LogicalExpr(M, dy(u)) == dx2(u))
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
