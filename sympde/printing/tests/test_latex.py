# coding: utf-8

from sympy import sin, cos, pi

from sympde.calculus import grad, dot, inner, rot, div
#from sympde.topology import (dx, dy, dz)
from sympde.topology import Domain, Boundary
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import element_of
from sympde.expr     import BilinearForm, LinearForm, integral
from sympde.exterior import d, wedge, ip, jp, delta, hodge
from sympde.exterior import DifferentialForm
from sympde.printing.latex import latex


#==============================================================================
def test_latex_1d():

    DIM = 1
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    x = V.coordinates

    v = element_of(V, name='v')
    u = element_of(V, name='u')
#    F = element_of(V, name='F')

    int_0 = lambda expr: integral(domain , expr)
    
    assert(latex(grad(v)) == r'\nabla{v}')
    assert(latex(dot(grad(v), grad(u))) == r'\nabla{u} \cdot \nabla{v}')

    a = BilinearForm((v,u), int_0(dot(grad(v), grad(u))))
    print(latex(a))
#    assert(latex(a) == r'\int_{0}^{1} \nabla{v} \cdot \nabla{u} dx')

    b = LinearForm(v, int_0(sin(pi*x)*v))
    print(latex(b))
#    assert(latex(b) == r'\int_{0}^{1} v \sin{\left (\pi x \right )} dx')

#    f = Integral(dx(F)-x, domain)
#    print(latex(f))
##    assert(latex(f) == r'\int_{0}^{1} - x + \partial_{x}F dx')

#==============================================================================
def test_latex_2d_1():

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    x,y = V.coordinates

    v = element_of(V, name='v')
    u = element_of(V, name='u')
#    F = element_of(V, name='F')
    
    int_0 = lambda expr: integral(domain , expr)

    assert(latex(grad(v)) == r'\nabla{v}')
    assert(latex(dot(grad(v), grad(u))) == r'\nabla{u} \cdot \nabla{v}')

    a = BilinearForm((v,u), int_0(dot(grad(v), grad(u))))
    print(latex(a))
#    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1} \nabla{v} \cdot \nabla{u} dxdy')

    b = LinearForm(v, int_0(sin(pi*x)*cos(pi*y)*v))
    print(latex(b))
#    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1} v \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} dxdy')

#    f = Integral(dx(F)-dy(F)-x*y, domain)
#    print(latex(f))
##    assert(latex(f) == r'\int_{0}^{1}\int_{0}^{1} - x y + \partial_{x}F - \partial_{y}F dxdy')

#==============================================================================
def test_latex_2d_2():

    DIM = 2
    domain = Domain('Omega', dim=DIM)

    V = VectorFunctionSpace('V', domain)

    x,y = V.coordinates

    v = element_of(V, name='v')
    u = element_of(V, name='u')
#    F = element_of(V, name='F')

    int_0 = lambda expr: integral(domain , expr)

    assert(latex(v) == r'\mathbf{v}')
    assert(latex(inner(grad(v), grad(u))) == r'\nabla{\mathbf{u}} : \nabla{\mathbf{v}}')

    a = BilinearForm((v,u), int_0(inner(grad(v), grad(u))))
    print(latex(a))
#    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1} \nabla{\mathbf{v}} : \nabla{\mathbf{u}} dxdy')

    b = LinearForm(v, int_0(sin(pi*x)*cos(pi*y)*div(v)))
    print(latex(b))
#    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1} \nabla \cdot \mathbf{v} \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} dxdy')


#==============================================================================
def test_latex_3d_1():

    DIM = 3
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    x,y,z = V.coordinates

    v = element_of(V, name='v')
    u = element_of(V, name='u')
#    F = element_of(V, name='F')

    int_0 = lambda expr: integral(domain , expr)
    
    assert(latex(grad(v)) == r'\nabla{v}')
    assert(latex(dot(grad(v), grad(u))) == r'\nabla{u} \cdot \nabla{v}')

    a = BilinearForm((v,u), int_0(dot(grad(v), grad(u))))
    print(latex(a))
#    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} \nabla{v} \cdot \nabla{u} dxdydz')

    b = LinearForm(v, int_0(sin(pi*x)*cos(pi*y)*cos(2*pi*z)*v))
    print(latex(b))
#    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} v \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} \cos{\left (2 \pi z \right )} dxdydz')

#    f = Integral(dx(F)-dy(F)+dz(F)-x*y*z, domain)
#    print(latex(f))
##    assert(latex(f) == r'\int_{0}^{1}\int_{0}^{1} - x y z + \partial_{x}F - \partial_{y}F + \partial_{z}F dxdy')

#==============================================================================
def test_latex_3d_2():

    DIM = 3
    domain = Domain('Omega', dim=DIM)

    V = VectorFunctionSpace('V', domain)

    x,y,z = V.coordinates

    v = element_of(V, name='v')
    u = element_of(V, name='u')

    int_0 = lambda expr: integral(domain , expr)
    
    assert(latex(v) == r'\mathbf{v}')
    assert(latex(inner(grad(v), grad(u))) == r'\nabla{\mathbf{u}} : \nabla{\mathbf{v}}')

    a = BilinearForm((v,u), int_0(inner(grad(v), grad(u))))
    print(latex(a))
#    assert(latex(a) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} \nabla{\mathbf{v}} : \nabla{\mathbf{u}} dxdydz')

    b = LinearForm(v, int_0(sin(pi*x)*cos(pi*y)*div(v)))
    print(latex(b))
#    assert(latex(b) == r'\int_{0}^{1}\int_{0}^{1}\int_{0}^{1} \nabla \cdot \mathbf{v} \sin{\left (\pi x \right )} \cos{\left (\pi y \right )} dxdydz')



#==============================================================================
def test_latex_2d_3():
    DIM = 2

    domain = Domain('Omega', dim=DIM)

    B1 = Boundary(r'\Gamma_1', domain)
    B2 = Boundary(r'\Gamma_2', domain)
    B3 = Boundary(r'\Gamma_3', domain)

    V = ScalarFunctionSpace('V', domain)

    x = V.coordinates

    v = element_of(V, name='v')
    u = element_of(V, name='u')

    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(B1, expr)
    
    # ...
    expr = dot(grad(v), grad(u))
    a_0 = BilinearForm((v,u), int_0(expr))

    expr = v*u
    a_bnd = BilinearForm((v, u), int_1(expr))

    expr = a_0(v,u) + a_bnd(v,u)
    a = BilinearForm((v,u), expr)
    print(latex(a_0))
    print(latex(a_bnd))
    print(latex(a))
#    print(a)
    print('')
    # ...

#==============================================================================
def test_latex_2d_4():
    DIM = 2

    domain = Domain('Omega', dim=DIM)

    # ... abstract model
    V = VectorFunctionSpace('V', domain)
    W = ScalarFunctionSpace('W', domain)

    v = element_of(V, name='v')
    u = element_of(V, name='u')
    p = element_of(W, name='p')
    q = element_of(W, name='q')

    int_0 = lambda expr: integral(domain , expr)

    a = BilinearForm((v,u), int_0(inner(grad(v), grad(u))))
    b = BilinearForm((v,p), int_0(div(v)*p))
    A = BilinearForm(((v,q),(u,p)), a(v,u) - b(v,p) + b(u,q))
    # ...

    print(latex(A))
#    print(latex(tensorize(A)))
    print('')
    # ...

#==============================================================================
def test_latex_2d_5():
    DIM = 2

    domain = Domain('Omega', dim=DIM)

    # ... abstract model
    W1 = VectorFunctionSpace('W1', domain)

    w1 = element_of(W1, name='w1')
    F  = element_of(W1, 'F')

    int_0 = lambda expr: integral(domain , expr)
    
    # ...
    l1 = LinearForm(w1, int_0(dot(w1, F)))

    print(latex(l1))
    print('')
    # ...

    # ...
    l2 = LinearForm(w1, int_0(rot(w1)*rot(F) + div(w1)*div(F)))

    print(latex(l2))
    print('')
    # ...

#==============================================================================
def test_latex_ec_3d_1():

    n = 3

    # ...
    u_0 = DifferentialForm('u_0', index=0, dim=n)
    v_0 = DifferentialForm('v_0', index=0, dim=n)

    u_1 = DifferentialForm('u_1', index=1, dim=n)
    v_1 = DifferentialForm('v_1', index=1, dim=n)

    u_2 = DifferentialForm('u_2', index=2, dim=n)
    v_2 = DifferentialForm('v_2', index=2, dim=n)

    u_3 = DifferentialForm('u_3', index=3, dim=n)
    v_3 = DifferentialForm('v_3', index=3, dim=n)
    # ...

    # ...
    domain = Domain('Omega', dim=3)
    V = VectorFunctionSpace('V', domain)

    beta = element_of(V, 'beta')
    # ...

    print(latex(u_0))
    print(latex(d(u_0)))
    print(latex(d(delta(u_3))))
    print(latex(d(delta(u_2)) + delta(d(u_2))))
    print(latex(wedge(u_0, u_1)))

    print(latex(ip(beta,u_1)))
    print(latex(hodge(u_1)))
    print(latex(jp(beta,u_1)))

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()

