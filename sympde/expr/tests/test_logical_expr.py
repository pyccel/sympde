# coding: utf-8

from sympy import Symbol, symbols
from sympy import Matrix
from sympy import expand
from sympy import cos, sin, sqrt, pi

from sympde.core     import constant
from sympde.calculus import grad, dot, inner, rot, div
from sympde.calculus import laplace, bracket, convect
from sympde.calculus import jump, avg, Dn, minus, plus
from sympde.topology import Domain, Mapping, Square
from sympde.topology import dx, dy
from sympde.topology import dx1, dx2, dx3
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import element_of, elements_of
from sympde.topology import LogicalExpr
from sympde.topology import SymbolicExpr
from sympde.topology import IdentityMapping
from sympde.topology import PolarMapping
from sympde.topology import TargetMapping
from sympde.topology import CzarnyMapping
from sympde.topology import CollelaMapping2D
from sympde.topology import TorusMapping
from sympde.topology import TwistedTargetMapping

from sympde.expr     import BilinearForm, integral
from sympde.calculus import grad, div, curl, dot

from sympde.topology.mapping import Jacobian

#==============================================================================
def test_logical_expr_2d_2():
    dim = 2

    A = Square('A')
    B = Square('B')

    M1 = Mapping('M1', dim=dim)
    M2 = Mapping('M2', dim=dim)

    D1 = M1(A)
    D2 = M2(B)

    domains = [D1, D2]
    connectivity = [((0, 0, 1),(1, 0, -1))]
    domain = Domain.join(domains, connectivity, 'domain')


    V1 = ScalarFunctionSpace('V1', domain, kind='h1')
    V2 = VectorFunctionSpace('V2', domain, kind='h1')

    u1,v1 = [element_of(V1, name=i) for i in ['u1', 'v1']]
    u2,v2 = [element_of(V2, name=i) for i in ['u2', 'v2']]

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(domain.interfaces , expr)

    expr = LogicalExpr(int_0(u1*v1), domain)
    assert str(expr.args[0]) == 'Integral(A, u1*v1*sqrt(det(Jacobian(M1).T * Jacobian(M1))))'
    assert str(expr.args[1]) == 'Integral(B, u1*v1*sqrt(det(Jacobian(M2).T * Jacobian(M2))))'

    expr = LogicalExpr(int_1(plus(u1)*plus(v1)), domain)
    assert str(expr) == 'Integral(A|B, PlusInterfaceOperator(u1)*PlusInterfaceOperator(v1)*sqrt(det(Jacobian(M1|M2).T * Jacobian(M1|M2))))'

    expr = LogicalExpr(int_1(minus(u1)*minus(v1)), domain)
    assert str(expr) == 'Integral(A|B, MinusInterfaceOperator(u1)*MinusInterfaceOperator(v1)*sqrt(det(Jacobian(M1|M2).T * Jacobian(M1|M2))))'

    expr = LogicalExpr(int_0(dot(u2,v2)), domain)
    assert str(expr.args[0]) == 'Integral(A, Dot(u2, v2)*sqrt(det(Jacobian(M1).T * Jacobian(M1))))'
    assert str(expr.args[1]) == 'Integral(B, Dot(u2, v2)*sqrt(det(Jacobian(M2).T * Jacobian(M2))))'

    expr = LogicalExpr(int_1(dot(plus(u2),plus(v2))), domain)
    assert str(expr) == 'Integral(A|B, Dot(PlusInterfaceOperator(u2), PlusInterfaceOperator(v2))*sqrt(det(Jacobian(M1|M2).T * Jacobian(M1|M2))))'

    expr = LogicalExpr(int_1(dot(minus(u2), minus(v2))), domain)
    assert str(expr) == 'Integral(A|B, Dot(MinusInterfaceOperator(u2), MinusInterfaceOperator(v2))*sqrt(det(Jacobian(M1|M2).T * Jacobian(M1|M2))))'

    expr = LogicalExpr(int_1(dot(grad(minus(u2)),grad(minus(v2)))), domain)
    assert str(expr) == 'Integral(A|B, Dot((Jacobian(M1)**(-1)).T * Grad(MinusInterfaceOperator(u2)), (Jacobian(M1)**(-1)).T * Grad(MinusInterfaceOperator(v2)))*sqrt(det(Jacobian(M1|M2).T * Jacobian(M1|M2))))'

    expr = LogicalExpr(int_1(dot(grad(plus(u2)),grad(plus(v2)))), domain)
    assert str(expr) == 'Integral(A|B, Dot((Jacobian(M2)**(-1)).T * Grad(PlusInterfaceOperator(u2)), (Jacobian(M2)**(-1)).T * Grad(PlusInterfaceOperator(v2)))*sqrt(det(Jacobian(M1|M2).T * Jacobian(M1|M2))))'

#==============================================================================
def test_logical_expr_2d_3():
    dim = 2

    A = Square('A')
    B = Square('B')

    M1 = Mapping('M1', dim=dim)
    M2 = Mapping('M2', dim=dim)

    D1 = M1(A)
    D2 = M2(B)


    domains = [D1, D2]
    connectivity = [((0, 0, 1),(1, 0, -1))]
    domain = Domain.join(domains, connectivity, 'domain')

    V = VectorFunctionSpace('V', domain, kind='hcurl')

    u,v = [element_of(V, name=i) for i in ['u', 'v']]

    int_0 = lambda expr: integral(domain , expr)

    expr = LogicalExpr(int_0(dot(u,v)), domain)
    assert str(expr.args[0]) == 'Integral(A, Dot((Jacobian(M1)**(-1)).T * u, (Jacobian(M1)**(-1)).T * v)*sqrt(det(Jacobian(M1).T * Jacobian(M1))))'
    assert str(expr.args[1]) == 'Integral(B, Dot((Jacobian(M2)**(-1)).T * u, (Jacobian(M2)**(-1)).T * v)*sqrt(det(Jacobian(M2).T * Jacobian(M2))))'

#==============================================================================
def test_logical_expr_3d_5():

    dim    = 3
    domain = Domain('Omega', dim=dim)
    M      = Mapping('M', dim=dim)

    mapped_domain = M(domain)

    V  = VectorFunctionSpace('V' , domain, kind='hcurl')
    VM = VectorFunctionSpace('VM', mapped_domain, kind='hcurl')

    J   = M.jacobian
    u,v   = elements_of(V,  names='u,v')
    um,vm = elements_of(VM, names='u,v')

    int_md = lambda expr: integral(mapped_domain , expr)
    int_ld = lambda expr: integral(domain , expr)

    am  = BilinearForm((um,vm), int_md(dot(curl(vm),curl(um))))
    a   = LogicalExpr(am, mapped_domain)

    assert a == BilinearForm((u,v), int_ld(sqrt((J.T*J).det())*dot(J/J.det()*curl(u), J/J.det()*curl(v))))

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================
def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
