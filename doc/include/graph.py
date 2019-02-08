# coding: utf-8

import os

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import IndexedBase
from sympy import Matrix
from sympy import Function
from sympy import pi, cos, sin, exp
from sympy import srepr
from sympy.physics.quantum import TensorProduct
from sympy.printing.dot import dotprint

from sympde.core import Constant
from sympde.calculus import grad, dot, inner, cross, rot, curl, div
from sympde.calculus import laplace, hessian, bracket, convect
from sympde.topology import (dx, dy, dz)
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import Field, VectorField
from sympde.topology import ProductSpace
from sympde.topology import TestFunction
from sympde.topology import VectorTestFunction
from sympde.topology import Unknown
from sympde.topology import InteriorDomain, Union
from sympde.topology import Boundary, NormalVector, TangentVector
from sympde.topology import Domain
from sympde.topology import Trace, trace_0, trace_1
from sympde.topology import Mapping
from sympde.topology import Square
from sympde.topology import ElementDomain
from sympde.topology import Area

from sympde.expr.expr import LinearExpr, BilinearExpr
from sympde.expr.expr import LinearForm, BilinearForm
from sympde.expr.expr import DomainIntegral, BoundaryIntegral
from sympde.expr.expr import Functional, Norm
from sympde.expr.expr import linearize
from sympde.expr.evaluation import TerminalExpr


#==============================================================================
def mkdir_p(dir):
    if os.path.isdir(dir):
        return
    os.makedirs(dir)


#==============================================================================
def export(expr, fname):
    """saves the graph as a png/svg file. the extension is eitheir png or svg"""

    name, ext = os.path.splitext(fname)
    ext = ext[1:]

    graph = dotprint(expr)
    f = open('{name}.dot'.format(name=name), 'w')
    f.write(graph)
    f.close()

    cmd = "dot -T{ext} {name}.dot -o{name}.{ext}".format(name=name, ext=ext)
    os.system(cmd)

#==============================================================================
def test_linear_form_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    eps   = Constant('eps', real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    g = Tuple(x**2, y**2)

    # ...
    d_forms = {}

    d_forms['l1'] = LinearForm(v, x*y*v)
    d_forms['l2'] = LinearForm(v, v*trace_1(g, B1))
    d_forms['l3'] = LinearForm(v, v*trace_1(g, B1) + x*y*v)

    # nitsche
    f = Function('f')
    g = Function('g')

    l0   = LinearForm(v, f(x,y)*v)
    l_B1 = LinearForm(v, g(x,y)*trace_0(v, B1)/eps
                         - kappa * g(x,y) * trace_1(grad(v), B1))

    d_forms['l4'] = LinearForm(v, l0(v) + l_B1(v))
    # ...

    # ...
    for name, expr in d_forms.items():
        export(expr, 'linform_2d_{}.png'.format(name))
    # ...

#==============================================================================
def test_bilinear_form_2d_1():

    domain = Domain('Omega', dim=2)
    B1 = Boundary(r'\Gamma_1', domain)

    x,y = domain.coordinates

    kappa = Constant('kappa', is_real=True)
    mu    = Constant('mu'   , is_real=True)
    eps   = Constant('eps', real=True)

    V = FunctionSpace('V', domain)

    u,u1,u2 = [TestFunction(V, name=i) for i in ['u', 'u1', 'u2']]
    v,v1,v2 = [TestFunction(V, name=i) for i in ['v', 'v1', 'v2']]

    # ...
    d_forms = {}

    d_forms['a1'] = BilinearForm((u,v), u*v)
    d_forms['a2'] = BilinearForm((u,v), u*v + dot(grad(u), grad(v)))
    d_forms['a3'] = BilinearForm((u,v), v*trace_1(grad(u), B1))


    # Poisson with Nitsch method
    a0 = BilinearForm((u,v), dot(grad(u),grad(v)))
    a_B1 = BilinearForm((u,v), - kappa * u*trace_1(grad(v), B1)
                               - v*trace_1(grad(u), B1)
                               + trace_0(u, B1) * trace_0(v, B1) / eps)
    a = BilinearForm((u,v), a0(u,v) + a_B1(u,v))

    d_forms['a4'] = a
    # ...

    # ... calls
    d_calls = {}
    for name, a in d_forms.items():
        d_calls[name] = a(u1,v1)
    # ...

    # ... export forms
    for name, expr in d_forms.items():
        export(expr, 'biform_2d_{}.png'.format(name))
    # ...

    # ... export calls
    for name, expr in d_calls.items():
        export(expr, 'biform_2d_call_{}.png'.format(name))
    # ...


#==============================================================================
# stokes
def test_bilinear_form_2d_3():

    domain = Domain('Omega', dim=2)

    x,y = domain.coordinates

    V = VectorFunctionSpace('V', domain)
    W = FunctionSpace('W', domain)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')
    p = TestFunction(W, name='p')
    q = TestFunction(W, name='q')

    a = BilinearForm((u,v), inner(grad(v), grad(u)))
    b = BilinearForm((v,p), div(v)*p)
    A = BilinearForm(((u,p),(v,q)), a(v,u) - b(v,p) + b(u,q))

    export(A, 'stokes_2d.png')

#==============================================================================
def test_linearize_form_2d_3():
    """steady Euler equation."""
    domain = Domain('Omega', dim=2)
    x,y = domain.coordinates

    U = VectorFunctionSpace('U', domain)
    W =       FunctionSpace('W', domain)

    v   = VectorTestFunction(U, name='v')
    phi =       TestFunction(W, name='phi')
    q   =       TestFunction(W, name='q')

    U_0   = VectorField(U, name='U_0')
    Rho_0 = Field(W, name='Rho_0')
    P_0   = Field(W, name='P_0')

    # ...
    expr = div(Rho_0*U_0) * phi
    l1 = LinearForm(phi, expr)

    expr = Rho_0*dot(convect(U_0, grad(U_0)), v) + dot(grad(P_0), v)
    l2 = LinearForm(v, expr)

    expr = dot(U_0, grad(P_0)) * q + P_0 * div(U_0) * q
    l3 = LinearForm(q, expr)
    # ...

    a1 = linearize(l1, [Rho_0, U_0], trials=['d_rho', 'd_u'])
    print(a1)
    print('')

    a2 = linearize(l2, [Rho_0, U_0, P_0], trials=['d_rho', 'd_u', 'd_p'])
    print(a2)
    print('')

    a3 = linearize(l3, [P_0, U_0], trials=['d_p', 'd_u'])
    print(a3)
    print('')

    l = LinearForm((phi, v, q), l1(phi) + l2(v) + l3(q))
    a = linearize(l, [Rho_0, U_0, P_0], trials=['d_rho', 'd_u', 'd_p'])
    print(a)

    export(a, 'steady_euler.png')

############################################
if __name__ == '__main__':

#    test_linear_form_2d_1()
#    test_bilinear_form_2d_1()
#    test_bilinear_form_2d_3()

    test_linearize_form_2d_3()
