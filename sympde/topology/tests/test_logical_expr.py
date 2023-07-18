# coding: utf-8

from sympy import Symbol, symbols
from sympy import Matrix
from sympy import expand
from sympy import cos, sin, sqrt, pi

from sympde.core     import Constant
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
def test_logical_expr_1d_1():
    dim = 1

    M = Mapping('M', dim=dim)
    domain = M(Domain('Omega', dim=dim))

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain, kind='h1')

    u,v = [element_of(V, name=i) for i in ['u', 'v']]

    det_M = Jacobian(M).det()
    #print('det = ', det_M)
    det   = Symbol('det')

    # ...
    expr = 2*u + alpha*v
    expr = LogicalExpr(expr, domain)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(expr, domain)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dx(det_M)
    expr = LogicalExpr(expr, domain)
    expr = expr.subs(det_M, det)
    expr = expand(expr)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(expr, domain)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

#==============================================================================
def test_symbolic_expr_1d_1():
    dim = 1

    M = Mapping('M', dim=dim)
    domain = M(Domain('Omega', dim=dim))

    V = ScalarFunctionSpace('V', domain, kind='h1')

    u = element_of(V, name='u')

    det_M = Jacobian(M).det()
    det_M = SymbolicExpr(det_M)
    #print('>>> ', det_M)
    det   = Symbol('det')

    # ...
    expr = u
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(u)
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(M[0])
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expr)
    # ...

    # ...
    expr = dx(Jacobian(M).det())
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(dx(u)))
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

#==============================================================================
def test_logical_expr_2d_1():
    dim = 2

    M = Mapping('M', dim=dim)
    domain = M(Domain('Omega', dim=dim))

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain, kind='h1')
    W = VectorFunctionSpace('V', domain, kind='h1')

    u,v = [element_of(V, name=i) for i in ['u', 'v']]
    w = element_of(W, name='w')

    det_M = Jacobian(M).det()
    #print('det = ', det_M)
    det   = Symbol('det')

    # ...
    expr = 2*u + alpha*v
    expr = LogicalExpr(expr, domain)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(expr, domain)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dy(u)
    expr = LogicalExpr(expr, domain)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dx(det_M)
    expr = LogicalExpr(expr, domain)
    expr = expr.subs(det_M, det)
    expr = expand(expr)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(expr, domain)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dx(w[0])
    expr = LogicalExpr(expr, domain)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

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
def test_symbolic_expr_2d_1():
    dim    = 2
    M      = Mapping('M', dim=dim)
    domain = M(Domain('Omega', dim=dim))

    V = ScalarFunctionSpace('V', domain, kind='h1')

    u = element_of(V, name='u')

    det_M = Jacobian(M).det()
    det_M = SymbolicExpr(det_M)
    #print('>>> ', det_M)
    det   = Symbol('det')

    # ...
    expr = u
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(u)
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(dx2(u))
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(M[0])
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expr)
    # ...

    # ...
    expr = dx(Jacobian(M).det())
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(dx(u)))
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

#==============================================================================
def test_logical_expr_3d_1():
    dim = 3

    M = Mapping('M', dim=dim)
    domain = M(Domain('Omega', dim=dim))

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain, kind='h1')

    u,v = [element_of(V, name=i) for i in ['u', 'v']]

    det_M = Jacobian(M).det()
    #print('det = ', det_M)
    det   = Symbol('det')

    # ...
    expr = 2*u + alpha*v
    expr = LogicalExpr(expr, domain)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(expr, domain)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dy(u)
    expr = LogicalExpr(expr, domain)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

    # ...
    expr = dx(det_M)
    expr = LogicalExpr(expr, domain)
    expr = expr.subs(det_M, det)
    #print(expr)
    #print('')
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(expr, domain)
    #print(expr.subs(det_M, det))
    #print('')
    # ...

def test_logical_expr_3d_2():

    dim   = 3
    domain = Domain('Omega', dim=dim)
    M      = Mapping('M', dim=dim)

    mapped_domain = M(domain)

    V  = ScalarFunctionSpace('V' , domain, kind='h1')
    VM = ScalarFunctionSpace('VM', mapped_domain, kind='h1')

    u,v   = elements_of(V, names='u,v')
    um,vm = elements_of(VM, names='u,v')

    J   = M.jacobian

    a = dot(grad(um),grad(vm))
    e = LogicalExpr(a, mapped_domain)

    assert e == dot(J.inv().T*grad(u), J.inv().T*grad(v))


def test_logical_expr_3d_3():

    dim    = 3
    domain = Domain('Omega', dim=dim)
    M      = Mapping('M', dim=dim)

    mapped_domain = M(domain)

    V  = VectorFunctionSpace('V' , domain, kind='hcurl')
    VM = VectorFunctionSpace('VM', mapped_domain, kind='hcurl')

    u,v   = elements_of(V, names='u,v')
    um,vm = elements_of(VM, names='u,v')

    J   = M.jacobian

    a = dot(curl(um), curl(vm))
    e = LogicalExpr(a, mapped_domain)

    assert e == dot(J/J.det()*curl(u), J/J.det()*curl(v))

def test_logical_expr_3d_4():

    dim    = 3
    domain = Domain('Omega', dim=dim)
    M      = Mapping('M', dim=dim)

    mapped_domain = M(domain)

    V  = VectorFunctionSpace('V' , domain, kind='hdiv')
    VM = VectorFunctionSpace('VM', mapped_domain, kind='hdiv')

    u,v   = elements_of(V, names='u,v')
    um,vm = elements_of(VM, names='u,v')

    J   = M.jacobian

    a = div(um)*div(vm)
    e = LogicalExpr(a, mapped_domain)

    assert e == J.det()**-2*div(u)*div(v)

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
def test_symbolic_expr_3d_1():
    dim    = 3
    M      = Mapping('M', dim=dim)
    domain = M(Domain('Omega', dim=dim))

    V = ScalarFunctionSpace('V', domain, kind='h1')
    u = element_of(V, 'u')

    det_M = Jacobian(M).det()
    det_M = SymbolicExpr(det_M)
    #print('>>> ', det_M)
    det   = Symbol('det')

    # ...
    expr = u
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(u)
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(dx2(u))
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx1(M[0])
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    #print(expr)
    # ...

    # ...
    expr = dx(u)
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expr)
    # ...

    # ...
    expr = dx(Jacobian(M).det())
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(u))
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))
    # ...

    # ...
    expr = dx(dx(dx(u)))
    expr = LogicalExpr(expr, domain)
    expr = SymbolicExpr(expr)
    expr = expr.subs(det_M, det)
    #print(expand(expr))

    # ...

#==============================================================================
def test_identity_mapping_2d_1():
    dim = 2

    x1, x2 = symbols('x1, x2')

    M = IdentityMapping('M', dim=dim)

    domain = M(Domain('Omega', dim=dim))

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))

    assert(LogicalExpr(M[0], domain) == x1)
    assert(LogicalExpr(M[1], domain) == x2)

    assert(LogicalExpr(dx1(M[0]), domain) == 1)
    assert(LogicalExpr(dx1(M[1]), domain) == 0)

    assert(LogicalExpr(dx2(M[0]), domain) == 0)
    assert(LogicalExpr(dx2(M[1]), domain) == 1)

    expected = Matrix([[1, 0], [0, 1]])
    assert( Jacobian(M) == expected )


#==============================================================================
def test_identity_mapping_2d_2():
    dim    = 2
    M      = IdentityMapping('F', dim=dim)
    domain = M(Domain('Omega', dim=dim))


    V = ScalarFunctionSpace('V', domain, kind='h1')
    u = element_of(V, name='u')

    # ...
    assert(LogicalExpr(dx(u), domain) == dx1(u))
    assert(LogicalExpr(dy(u), domain) == dx2(u))
    # ...

#==============================================================================
def test_polar_mapping_2d_1():
    dim = 2

    x1, x2 = symbols('x1, x2')

    constants = ['c1', 'c2', 'rmax', 'rmin']
    c1, c2, rmax, rmin = [Constant(i) for i in constants]

    M = PolarMapping('M', dim=dim)

    domain = M(Domain('Omega', dim=dim))

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))

    assert(LogicalExpr(M[0], domain) == c1 + (rmax*x1 + rmin*(-x1 + 1))*cos(x2))
    assert(LogicalExpr(M[1], domain) == c2 + (rmax*x1 + rmin*(-x1 + 1))*sin(x2))

    assert(LogicalExpr(dx1(M[0]), domain) == (rmax - rmin)*cos(x2))
    assert(LogicalExpr(dx1(M[1]), domain) == (rmax - rmin)*sin(x2))

    expected = -(rmax*x1 + rmin*(-x1 + 1))*sin(x2)
    assert(expand(LogicalExpr(dx2(M[0]), domain)) == expand(expected))
    assert(LogicalExpr(dx2(M[1]), domain) == (rmax*x1 + rmin*(-x1 + 1))*cos(x2))

    expected = Matrix([[(rmax - rmin)*cos(x2), -(rmax*x1 + rmin*(-x1 + 1))*sin(x2)],
                       [(rmax - rmin)*sin(x2), (rmax*x1 + rmin*(-x1 + 1))*cos(x2)]])
    assert(expand(LogicalExpr(Jacobian(M), domain)) == expand(expected))

#==============================================================================
def test_target_mapping_2d_1():
    dim = 2

    x1, x2 = symbols('x1, x2')

    constants = ['c1', 'c2', 'D', 'k']
    c1, c2, D, k = [Constant(i) for i in constants]

    M = TargetMapping('M', dim=dim)

    domain = M(Domain('Omega', dim=dim))

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))

    assert(LogicalExpr(M[0], domain) == -D*x1**2 + c1 + x1*(-k + 1)*cos(x2))
    assert(LogicalExpr(M[1], domain) == c2 + x1*(k + 1)*sin(x2))

    assert(LogicalExpr(dx1(M[0]), domain) == -2*D*x1 + (-k + 1)*cos(x2))
    assert(LogicalExpr(dx1(M[1]), domain) == (k + 1)*sin(x2))

    assert(LogicalExpr(dx2(M[0]), domain) == -x1*(-k + 1)*sin(x2))
    assert(LogicalExpr(dx2(M[1]), domain) == x1*(k + 1)*cos(x2))

    expected = Matrix([[-2*D*x1 + (-k + 1)*cos(x2),
                        -x1*(-k + 1)*sin(x2)],
                       [(k + 1)*sin(x2),
                        x1*(k + 1)*cos(x2)]])

    assert( Jacobian(M) == expected )

#==============================================================================
def test_czarny_mapping_2d_1():
    dim = 2

    x1, x2 = symbols('x1, x2')

    constants = ['c2', 'eps', 'b']
    c2, eps, b = [Constant(i) for i in constants]

    M = CzarnyMapping('M', dim=dim)

    domain = M(Domain('Omega', dim=dim))

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))

    expected =  (-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 1)/eps
    assert(LogicalExpr(M[0], domain) == expected)
    expected =  b*x1*sin(x2)/(sqrt(-eps**2/4 + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)) + c2
    assert(LogicalExpr(M[1], domain) == expected)

    expected =  -cos(x2)/sqrt(eps*(eps + 2*x1*cos(x2)) + 1)
    assert(LogicalExpr(dx1(M[0]), domain) == expected)
    expected =  b*(eps*x1*sin(x2)*cos(x2)/(sqrt(-eps**2/4 + 1)*sqrt(eps*(eps + 2*x1*cos(x2)) + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)**2) + sin(x2)/(sqrt(-eps**2/4 + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)))
    assert((LogicalExpr(dx1(M[1]), domain) - expected).expand() == 0)

    expected =  x1*sin(x2)/sqrt(eps*(eps + 2*x1*cos(x2)) + 1)
    assert(LogicalExpr(dx2(M[0]), domain) == expected)
    expected =  b*x1*(-eps*x1*sin(x2)**2/(sqrt(eps*(eps + 2*x1*cos(x2)) + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)**2) + cos(x2)/(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2))/sqrt(-eps**2/4 + 1)
    assert((LogicalExpr(dx2(M[1]), domain) - expected).expand() == 0)

    expected =  Matrix([[-cos(x2)/sqrt(eps*(eps + 2*x1*cos(x2)) + 1),
                         x1*sin(x2)/sqrt(eps*(eps + 2*x1*cos(x2)) + 1)],
                        [b*(eps*x1*sin(x2)*cos(x2)/(sqrt(-eps**2/4 + 1)*sqrt(eps*(eps + 2*x1*cos(x2)) + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)**2) + sin(x2)/(sqrt(-eps**2/4 + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2))),
                         b*x1*(-eps*x1*sin(x2)**2/(sqrt(eps*(eps + 2*x1*cos(x2)) + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)**2) + cos(x2)/(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2))/sqrt(-eps**2/4 + 1)]])
    assert( all(e.expand().is_zero for e in ( Jacobian(M) - expected)) )

#==============================================================================
def test_collela_mapping_2d_1():
    dim = 2

    x1, x2 = symbols('x1, x2')

    constants = ['eps', 'k1', 'k2']
    eps, k1, k2 = [Constant(i) for i in constants]

    M = CollelaMapping2D('M', dim)

    domain = M(Domain('Omega', dim=dim))

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))

    expected = 2.0*eps*sin(2.0*k1*pi*x1)*sin(2.0*k2*pi*x2) + 2.0*x1 - 1.0
    assert(LogicalExpr(M[0], domain) == expected)
    expected = 2.0*eps*sin(2.0*k1*pi*x1)*sin(2.0*k2*pi*x2) + 2.0*x2 - 1.0
    assert(LogicalExpr(M[1], domain) == expected)

    expected = 4.0*eps*k1*pi*sin(2.0*k2*pi*x2)*cos(2.0*k1*pi*x1) + 2.0
    assert(LogicalExpr(dx1(M[0]), domain) == expected)
    expected = 4.0*eps*k1*pi*sin(2.0*k2*pi*x2)*cos(2.0*k1*pi*x1)
    assert(LogicalExpr(dx1(M[1]), domain) == expected)

    expected = 4.0*eps*k2*pi*sin(2.0*k1*pi*x1)*cos(2.0*k2*pi*x2)
    assert(LogicalExpr(dx2(M[0]), domain) == expected)
    expected = 4.0*eps*k2*pi*sin(2.0*k1*pi*x1)*cos(2.0*k2*pi*x2) + 2.0
    assert(LogicalExpr(dx2(M[1]), domain) == expected)

    expected = Matrix([[4.0*eps*k1*pi*sin(2.0*k2*pi*x2)*cos(2.0*k1*pi*x1) + 2.0,
                        4.0*eps*k2*pi*sin(2.0*k1*pi*x1)*cos(2.0*k2*pi*x2)],
                       [4.0*eps*k1*pi*sin(2.0*k2*pi*x2)*cos(2.0*k1*pi*x1),
                        4.0*eps*k2*pi*sin(2.0*k1*pi*x1)*cos(2.0*k2*pi*x2) + 2.0]])

    assert(( Jacobian(M) == expected ))

#==============================================================================
def test_torus_mapping_3d_1():
    dim = 3

    x1, x2, x3 = symbols('x1, x2, x3')
    R0 = Constant('R0')

    M = TorusMapping('M', dim=dim)

    domain = M(Domain('Omega', dim=dim))

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))
    assert(not( M[2] == x3 ))

    expected = (R0 + x1*cos(x2))*cos(x3)
    assert(LogicalExpr(M[0], domain) == expected)
    expected = (R0 + x1*cos(x2))*sin(x3)
    assert(LogicalExpr(M[1], domain) == expected)
    expected = x1*sin(x2)
    assert(LogicalExpr(M[2], domain) == expected)

    expected = cos(x2)*cos(x3)
    assert(LogicalExpr(dx1(M[0]), domain) == expected)
    expected = sin(x3)*cos(x2)
    assert(LogicalExpr(dx1(M[1]), domain) == expected)
    expected = sin(x2)
    assert(LogicalExpr(dx1(M[2]), domain) == expected)

    expected = -x1*sin(x2)*cos(x3)
    assert(LogicalExpr(dx2(M[0]), domain) == expected)
    expected = -x1*sin(x2)*sin(x3)
    assert(LogicalExpr(dx2(M[1]), domain) == expected)
    expected = x1*cos(x2)
    assert(LogicalExpr(dx2(M[2]), domain) == expected)

    expected = -(R0 + x1*cos(x2))*sin(x3)
    assert(expand(LogicalExpr(dx3(M[0]), domain)) == expand(expected))
    expected = (R0 + x1*cos(x2))*cos(x3)
    assert(LogicalExpr(dx3(M[1]), domain) == expected)
    expected = 0
    assert(LogicalExpr(dx3(M[2]), domain) == expected)

    expected = Matrix([[cos(x2)*cos(x3),-x1*sin(x2)*cos(x3),-(R0+x1*cos(x2))*sin(x3)],
                       [sin(x3)*cos(x2),-x1*sin(x2)*sin(x3),(R0+x1*cos(x2))*cos(x3)],
                       [sin(x2), x1*cos(x2), 0]])

    assert( all(e.expand().is_zero for e in ( Jacobian(M) - expected)) )

#==============================================================================
def test_twisted_target_mapping_3d_1():
    dim = 3

    x1, x2, x3 = symbols('x1, x2, x3')

    constants = ['c1', 'c2', 'c3', 'D', 'k']
    c1, c2, c3, D, k = [Constant(i) for i in constants]

    M = TwistedTargetMapping('M', dim=dim)

    domain = M(Domain('Omega', dim=dim))

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))
    assert(not( M[2] == x3 ))

    expected = -D*x1**2 + c1 + x1*(-k + 1)*cos(x2)
    assert(LogicalExpr(M[0], domain) == expected)
    expected = c2 + x1*(k + 1)*sin(x2)
    assert(LogicalExpr(M[1], domain) == expected)
    expected = c3 + x1**2*x3*sin(2*x2)
    assert(LogicalExpr(M[2], domain) == expected)

    expected = -2*D*x1 + (-k + 1)*cos(x2)
    assert(LogicalExpr(dx1(M[0]), domain) == expected)
    expected = (k + 1)*sin(x2)
    assert(LogicalExpr(dx1(M[1]), domain) == expected)
    expected = 2*x1*x3*sin(2*x2)
    assert(LogicalExpr(dx1(M[2]), domain) == expected)

    expected = -x1*(-k + 1)*sin(x2)
    assert(LogicalExpr(dx2(M[0]), domain) == expected)
    expected = x1*(k + 1)*cos(x2)
    assert(LogicalExpr(dx2(M[1]), domain) == expected)
    expected = 2*x1**2*x3*cos(2*x2)
    assert(LogicalExpr(dx2(M[2]), domain) == expected)

    expected = 0
    assert(expand(LogicalExpr(dx3(M[0]), domain)) == expand(expected))
    expected = 0
    assert(LogicalExpr(dx3(M[1]), domain) == expected)
    expected = x1**2*sin(2*x2)
    assert(LogicalExpr(dx3(M[2]), domain) == expected)

    expected = Matrix([[-2*D*x1 + (-k + 1)*cos(x2), -x1*(-k + 1)*sin(x2), 0],
                       [(k + 1)*sin(x2),             x1*(k + 1)*cos(x2), 0],
                       [2*x1*x3*sin(2*x2),           2*x1**2*x3*cos(2*x2), x1**2*sin(2*x2)]])
    assert( Jacobian(M) == expected )

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================
def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()
