# coding: utf-8

from sympy import Symbol, symbols
from sympy import Matrix
from sympy import expand
from sympy import cos, sin, sqrt, pi

from sympde.core     import Constant
from sympde.topology import Domain, Mapping
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
from sympde.topology import CollelaMapping
from sympde.topology import TorusMapping
from sympde.topology import TwistedTargetMapping

from sympde.calculus import grad, div, curl, dot

from sympde.topology.mapping import Jacobian

#==============================================================================
def test_logical_expr_1d_1():
    rdim = 1

    M = Mapping('M', rdim)
    domain = M(Domain('Omega', dim=rdim))

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)

    u,v = [element_of(V, name=i) for i in ['u', 'v']]

    det_M = Jacobian(M).det()
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

#==============================================================================
def test_symbolic_expr_1d_1():
    rdim = 1

    M = Mapping('M', rdim)
    domain = M(Domain('Omega', dim=rdim))

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)

    u = element_of(V, name='u')

    det_M = Jacobian(M).det()
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
    expr = dx(Jacobian(M).det())
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
def test_logical_expr_2d_1():
    rdim = 2

    M = Mapping('M', rdim)
    domain = M(Domain('Omega', dim=rdim))

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)
    W = VectorFunctionSpace('V', domain)

    u,v = [element_of(V, name=i) for i in ['u', 'v']]
    w = element_of(W, name='w')

    det_M = Jacobian(M).det()
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

#==============================================================================
def test_symbolic_expr_2d_1():
    rdim = 2

    M = Mapping('M', rdim)
    domain = M(Domain('Omega', dim=rdim))

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)

    u = element_of(V, name='u')

    det_M = Jacobian(M).det()
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
    expr = dx(Jacobian(M).det())
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
def test_logical_expr_3d_1():
    rdim = 3

    M = Mapping('M', rdim)
    domain = M(Domain('Omega', dim=rdim))

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)

    u,v = [element_of(V, name=i) for i in ['u', 'v']]

    det_M = Jacobian(M).det()
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

def test_logical_expr_3d_2():

    domain = Domain('Omega', dim=3)
    M      = Mapping('M', 3)

    mapped_domain = M(domain)

    V  = ScalarFunctionSpace('V' , domain)
    VM = ScalarFunctionSpace('VM', mapped_domain)

    u,v   = elements_of(V, names='u,v')
    um,vm = elements_of(VM, names='u,v')

    J   = M.jacobian

    a = dot(grad(um),grad(vm))
    e = LogicalExpr(M, a)

    
    assert e == dot(J.inv().T*grad(u), J.inv().T*grad(v))


def test_logical_expr_3d_3():

    domain = Domain('Omega', dim=3)
    M      = Mapping('M', 3)

    mapped_domain = M(domain)

    V  = VectorFunctionSpace('V' , domain, kind='hcurl')
    VM = VectorFunctionSpace('VM', mapped_domain, kind='hcurl')

    u,v   = elements_of(V, names='u,v')
    um,vm = elements_of(VM, names='u,v')

    J   = M.jacobian

    a = dot(curl(um), curl(vm))
    e = LogicalExpr(M, a)

    assert e == dot(J/J.det()*curl(u), J/J.det()*curl(v))

def test_logical_expr_3d_4():

    domain = Domain('Omega', dim=3)
    M      = Mapping('M', 3)

    mapped_domain = M(domain)

    V  = VectorFunctionSpace('V' , domain, kind='hdiv')
    VM = VectorFunctionSpace('VM', mapped_domain, kind='hdiv')

    u,v   = elements_of(V, names='u,v')
    um,vm = elements_of(VM, names='u,v')

    J   = M.jacobian

    a = div(um)*div(vm)
    e = LogicalExpr(M, a)

    assert e == J.det()**-2*div(u)*div(v)

#==============================================================================
def test_symbolic_expr_3d_1():
    rdim = 3

    M = Mapping('M', rdim)
    domain = M(Domain('Omega', dim=rdim))

    alpha = Constant('alpha')

    V = ScalarFunctionSpace('V', domain)
    u = element_of(V, 'u')

    det_M = Jacobian(M).det()
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
    expr = dx(Jacobian(M).det())
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
    assert(LogicalExpr(M, Jacobian(M)) == expected)


#==============================================================================
def test_identity_mapping_2d_2():
    rdim = 2

    x1, x2 = symbols('x1, x2')

    M      = IdentityMapping('F', rdim)
    domain = M(Domain('Omega', dim=rdim))


    V = ScalarFunctionSpace('V', domain)
    u = element_of(V, name='u')

    # ...
    assert(LogicalExpr(M, dx(u)) == dx1(u))
    assert(LogicalExpr(M, dy(u)) == dx2(u))
    # ...

#==============================================================================
def test_polar_mapping_2d_1():
    rdim = 2

    x1, x2 = symbols('x1, x2')

    constants = ['c1', 'c2', 'rmax', 'rmin']
    c1, c2, rmax, rmin = [Constant(i) for i in constants]

    M = PolarMapping('M', rdim)

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))

    assert(LogicalExpr(M, M[0]) == c1 + (rmax*x1 + rmin*(-x1 + 1))*cos(x2))
    assert(LogicalExpr(M, M[1]) == c2 + (rmax*x1 + rmin*(-x1 + 1))*sin(x2))

    assert(LogicalExpr(M, dx1(M[0])) == (rmax - rmin)*cos(x2))
    assert(LogicalExpr(M, dx1(M[1])) == (rmax - rmin)*sin(x2))

    expected = -(rmax*x1 + rmin*(-x1 + 1))*sin(x2)
    assert(expand(LogicalExpr(M, dx2(M[0]))) == expand(expected))
    assert(LogicalExpr(M, dx2(M[1])) == (rmax*x1 + rmin*(-x1 + 1))*cos(x2))

    expected = Matrix([[(rmax - rmin)*cos(x2), -(rmax*x1 + rmin*(-x1 + 1))*sin(x2)],
                       [(rmax - rmin)*sin(x2), (rmax*x1 + rmin*(-x1 + 1))*cos(x2)]])
    assert(expand(LogicalExpr(M, Jacobian(M))) == expand(expected))

#==============================================================================
def test_target_mapping_2d_1():
    rdim = 2

    x1, x2 = symbols('x1, x2')

    constants = ['c1', 'c2', 'D', 'k']
    c1, c2, D, k = [Constant(i) for i in constants]

    M = TargetMapping('M', rdim)

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))

    assert(LogicalExpr(M, M[0]) == -D*x1**2 + c1 + x1*(-k + 1)*cos(x2))
    assert(LogicalExpr(M, M[1]) == c2 + x1*(k + 1)*sin(x2))

    assert(LogicalExpr(M, dx1(M[0])) == -2*D*x1 + (-k + 1)*cos(x2))
    assert(LogicalExpr(M, dx1(M[1])) == (k + 1)*sin(x2))

    assert(LogicalExpr(M, dx2(M[0])) == -x1*(-k + 1)*sin(x2))
    assert(LogicalExpr(M, dx2(M[1])) == x1*(k + 1)*cos(x2))

    expected = Matrix([[-2*D*x1 + (-k + 1)*cos(x2),
                        -x1*(-k + 1)*sin(x2)],
                       [(k + 1)*sin(x2),
                        x1*(k + 1)*cos(x2)]])

    assert(not( Jacobian(M) == expected))
    assert(expand(LogicalExpr(M, Jacobian(M))) == expand(expected))

#==============================================================================
def test_czarny_mapping_2d_1():
    rdim = 2

    x1, x2 = symbols('x1, x2')

    constants = ['c2', 'eps', 'b']
    c2, eps, b = [Constant(i) for i in constants]

    M = CzarnyMapping('M', rdim)

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))

    expected =  (-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 1)/eps
    assert(LogicalExpr(M, M[0]) == expected)
    expected =  b*x1*sin(x2)/(sqrt(-eps**2/4 + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)) + c2
    assert(LogicalExpr(M, M[1]) == expected)

    expected =  -cos(x2)/sqrt(eps*(eps + 2*x1*cos(x2)) + 1)
    assert(LogicalExpr(M, dx1(M[0])) == expected)
    expected =  b*(eps*x1*sin(x2)*cos(x2)/(sqrt(-eps**2/4 + 1)*sqrt(eps*(eps + 2*x1*cos(x2)) + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)**2) + sin(x2)/(sqrt(-eps**2/4 + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)))
    assert(LogicalExpr(M, dx1(M[1])) == expected)

    expected =  x1*sin(x2)/sqrt(eps*(eps + 2*x1*cos(x2)) + 1)
    assert(LogicalExpr(M, dx2(M[0])) == expected)
    expected =  b*x1*(-eps*x1*sin(x2)**2/(sqrt(eps*(eps + 2*x1*cos(x2)) + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)**2) + cos(x2)/(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2))/sqrt(-eps**2/4 + 1)
    assert(LogicalExpr(M, dx2(M[1])) == expected)

    expected =  Matrix([[-cos(x2)/sqrt(eps*(eps + 2*x1*cos(x2)) + 1),
                         x1*sin(x2)/sqrt(eps*(eps + 2*x1*cos(x2)) + 1)],
                        [b*(eps*x1*sin(x2)*cos(x2)/(sqrt(-eps**2/4 + 1)*sqrt(eps*(eps + 2*x1*cos(x2)) + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)**2) + sin(x2)/(sqrt(-eps**2/4 + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2))),
                         b*x1*(-eps*x1*sin(x2)**2/(sqrt(eps*(eps + 2*x1*cos(x2)) + 1)*(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2)**2) + cos(x2)/(-sqrt(eps*(eps + 2*x1*cos(x2)) + 1) + 2))/sqrt(-eps**2/4 + 1)]])
    assert(not( Jacobian(M) == expected))
    assert(expand(LogicalExpr(M, Jacobian(M))) == expand(expected))

#==============================================================================
def test_collela_mapping_2d_1():
    rdim = 2

    x1, x2 = symbols('x1, x2')

    constants = ['eps', 'k1', 'k2']
    eps, k1, k2 = [Constant(i) for i in constants]

    M = CollelaMapping('M', rdim)

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))

    expected = 2.0*eps*sin(2.0*k1*pi*x1)*sin(2.0*k2*pi*x2) + 2.0*x1 - 1.0
    assert(LogicalExpr(M, M[0]) == expected)
    expected = 2.0*eps*sin(2.0*k1*pi*x1)*sin(2.0*k2*pi*x2) + 2.0*x2 - 1.0
    assert(LogicalExpr(M, M[1]) == expected)

    expected = 4.0*eps*k1*pi*sin(2.0*k2*pi*x2)*cos(2.0*k1*pi*x1) + 2.0
    assert(LogicalExpr(M, dx1(M[0])) == expected)
    expected = 4.0*eps*k1*pi*sin(2.0*k2*pi*x2)*cos(2.0*k1*pi*x1)
    assert(LogicalExpr(M, dx1(M[1])) == expected)

    expected = 4.0*eps*k2*pi*sin(2.0*k1*pi*x1)*cos(2.0*k2*pi*x2)
    assert(LogicalExpr(M, dx2(M[0])) == expected)
    expected = 4.0*eps*k2*pi*sin(2.0*k1*pi*x1)*cos(2.0*k2*pi*x2) + 2.0
    assert(LogicalExpr(M, dx2(M[1])) == expected)

    expected = Matrix([[4.0*eps*k1*pi*sin(2.0*k2*pi*x2)*cos(2.0*k1*pi*x1) + 2.0,
                        4.0*eps*k2*pi*sin(2.0*k1*pi*x1)*cos(2.0*k2*pi*x2)],
                       [4.0*eps*k1*pi*sin(2.0*k2*pi*x2)*cos(2.0*k1*pi*x1),
                        4.0*eps*k2*pi*sin(2.0*k1*pi*x1)*cos(2.0*k2*pi*x2) + 2.0]])
    assert(not( Jacobian(M) == expected))
    assert(expand(LogicalExpr(M, Jacobian(M))) == expand(expected))

#==============================================================================
def test_torus_mapping_3d_1():
    rdim = 3

    x1, x2, x3 = symbols('x1, x2, x3')
    R0 = Constant('R0')

    M = TorusMapping('M', rdim)

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))
    assert(not( M[2] == x3 ))

    expected = (R0 + x1*cos(x2))*cos(x3)
    assert(LogicalExpr(M, M[0]) == expected)
    expected = (R0 + x1*cos(x2))*sin(x3)
    assert(LogicalExpr(M, M[1]) == expected)
    expected = x1*sin(x2)
    assert(LogicalExpr(M, M[2]) == expected)

    expected = cos(x2)*cos(x3)
    assert(LogicalExpr(M, dx1(M[0])) == expected)
    expected = sin(x3)*cos(x2)
    assert(LogicalExpr(M, dx1(M[1])) == expected)
    expected = sin(x2)
    assert(LogicalExpr(M, dx1(M[2])) == expected)

    expected = -x1*sin(x2)*cos(x3)
    assert(LogicalExpr(M, dx2(M[0])) == expected)
    expected = -x1*sin(x2)*sin(x3)
    assert(LogicalExpr(M, dx2(M[1])) == expected)
    expected = x1*cos(x2)
    assert(LogicalExpr(M, dx2(M[2])) == expected)

    expected = -(R0 + x1*cos(x2))*sin(x3)
    assert(expand(LogicalExpr(M, dx3(M[0]))) == expand(expected))
    expected = (R0 + x1*cos(x2))*cos(x3)
    assert(LogicalExpr(M, dx3(M[1])) == expected)
    expected = 0
    assert(LogicalExpr(M, dx3(M[2])) == expected)

    expected = Matrix([[cos(x2)*cos(x3),-x1*sin(x2)*cos(x3),-(R0+x1*cos(x2))*sin(x3)],
                       [sin(x3)*cos(x2),-x1*sin(x2)*sin(x3),(R0+x1*cos(x2))*cos(x3)],
                       [sin(x2), x1*cos(x2), 0]])
    assert(not( Jacobian(M) == expected))
    assert(expand(LogicalExpr(M, Jacobian(M))) == expand(expected))

#==============================================================================
def test_twisted_target_mapping_3d_1():
    rdim = 3

    x1, x2, x3 = symbols('x1, x2, x3')

    constants = ['c1', 'c2', 'c3', 'D', 'k']
    c1, c2, c3, D, k = [Constant(i) for i in constants]

    M = TwistedTargetMapping('M', rdim)

    assert(not( M[0] == x1 ))
    assert(not( M[1] == x2 ))
    assert(not( M[2] == x3 ))

    expected = -D*x1**2 + c1 + x1*(-k + 1)*cos(x2)
    assert(LogicalExpr(M, M[0]) == expected)
    expected = c2 + x1*(k + 1)*sin(x2)
    assert(LogicalExpr(M, M[1]) == expected)
    expected = c3 + x1**2*x3*sin(2*x2)
    assert(LogicalExpr(M, M[2]) == expected)

    expected = -2*D*x1 + (-k + 1)*cos(x2)
    assert(LogicalExpr(M, dx1(M[0])) == expected)
    expected = (k + 1)*sin(x2)
    assert(LogicalExpr(M, dx1(M[1])) == expected)
    expected = 2*x1*x3*sin(2*x2)
    assert(LogicalExpr(M, dx1(M[2])) == expected)

    expected = -x1*(-k + 1)*sin(x2)
    assert(LogicalExpr(M, dx2(M[0])) == expected)
    expected = x1*(k + 1)*cos(x2)
    assert(LogicalExpr(M, dx2(M[1])) == expected)
    expected = 2*x1**2*x3*cos(2*x2)
    assert(LogicalExpr(M, dx2(M[2])) == expected)

    expected = 0
    assert(expand(LogicalExpr(M, dx3(M[0]))) == expand(expected))
    expected = 0
    assert(LogicalExpr(M, dx3(M[1])) == expected)
    expected = x1**2*sin(2*x2)
    assert(LogicalExpr(M, dx3(M[2])) == expected)

    expected = Matrix([[-2*D*x1 + (-k + 1)*cos(x2), -x1*(-k + 1)*sin(x2), 0],
                       [(k + 1)*sin(x2),             x1*(k + 1)*cos(x2), 0],
                       [2*x1*x3*sin(2*x2),           2*x1**2*x3*cos(2*x2), x1**2*sin(2*x2)]])
    assert(not( Jacobian(M) == expected))
    assert(expand(LogicalExpr(M, Jacobian(M))) == expand(expected))


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================
def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
