# coding: utf-8
from sympy.printing.dot import dotprint
import os

def dotexport(expr, fname):
    txt = str(dotprint(expr))

    name = os.path.splitext(fname)[0]

    f = open('{}.dot'.format(name), 'w')
    f.write(txt)
    f.close()

    cmd = 'dot {name}.dot -Tpng -o {name}.png'.format(name=name)
    os.system(cmd)


# ...
#from sympde import Unknown
#from sympde import grad, div
#
#u = Unknown('u', ldim=2)
#expr = - div(grad(u)) + u
#dotexport(expr, 'graph1.png')
# ...

# ...
#from sympde import Unknown
#from sympde import dx, dy
#
#u = Unknown('u', ldim=2)
#v = Unknown('v', ldim=2)
#
#expr = dx(dy(u*v))
#dotexport(expr, 'graph2.png')
# ...

# ...
#from sympde import Unknown
#from sympde import dx, dy
#
#u = Unknown('u', ldim=2)
#v = Unknown('v', ldim=2)
#
#expr = dy(2*u+3*v)
#print(expr)
# ...

# ...
#from sympde import Unknown, Constant
#from sympde import dx
#
#u = Unknown('u', ldim=1)
#alpha = Constant('alpha')
#
#expr = dx(alpha*u) + dx(dx(2*u))
#print(expr)
# ...

# ...
#from sympde import Constant
#from sympde import dx, dy
#from sympy.abc import x, y
#from sympy import cos, exp
#
#alpha = Constant('alpha')
#
#L = lambda u: -dx(dx(u)) - dy(dy(u)) + alpha * u
#
#expr = L(cos(y)*exp(-x**2))
#print(expr)
# ...

# ...
#from sympde import Constant
#from sympde import dx, dy
#from sympy.abc import x, y
#from sympy import Function
#
#alpha = Constant('alpha')
#f = Function('f')
#
#L = lambda u: -dx(dx(u)) - dy(dy(u)) + alpha * u
#
#expr = L(f(x,y))
#print(expr)
# ...

# ...
#from sympde import grad, dot
#from sympde import FunctionSpace
#from sympde import TestFunction
#from sympde import BilinearForm
#
#V = FunctionSpace('V', ldim=2)
#U = FunctionSpace('U', ldim=2)
#
#v = TestFunction(V, name='v')
#u = TestFunction(U, name='u')
#
#a = BilinearForm((v,u), dot(grad(v), grad(u)) + v*u)
#
#dotexport(a, 'graph_laplace.png')
# ...

# ...
#from sympde import dx
#from sympde import FunctionSpace
#from sympde import TestFunction
#from sympde import BilinearForm
#from sympde import Constant
#
#V = FunctionSpace('V', ldim=1)
#W = FunctionSpace('W', ldim=1)
#
#T = Constant('T', real=True, label='Tension applied to the string')
#rho = Constant('rho', real=True, label='mass density')
#dt = Constant('dt', real=True, label='time step')
#
## trial functions
#u = TestFunction(V, name='u')
#f = TestFunction(W, name='f')
#
## test functions
#v = TestFunction(V, name='v')
#w = TestFunction(W, name='w')
#
#mass = BilinearForm((v,u), v*u)
#adv  = BilinearForm((v,u), dx(v)*u)
#
#expr = rho*mass(v,u) + dt*adv(v, f) + dt*adv(w,u) + mass(w,f)
#a = BilinearForm(((v,w), (u,f)), expr)
#
#print(a)
##dotexport(a, 'graph_wave.png')
# ...

# ...
#from sympde import FunctionSpace
#from sympde import TestFunction
#from sympde import LinearForm
#from sympy import cos
#
#V = FunctionSpace('V', ldim=2)
#
#v = TestFunction(V, name='v')
#
#x,y = V.coordinates
#
#b = LinearForm(v, cos(x-y)*v)
# ...

# ...
#from sympde import grad, div
#from sympde import FunctionSpace
#from sympde import Field
#from sympde import FunctionForm
#from sympy import cos, pi
#
#V = FunctionSpace('V', ldim=1)
#F = Field('F', space=V)
#
#x = V.coordinates
#
#b = FunctionForm(div(grad(F-cos(2*pi*x))))
# ...

# ...
#from sympde import grad, dot
#from sympde import FunctionSpace
#from sympde import TestFunction
#from sympde import BilinearForm
#from sympde import evaluate
#
#V = FunctionSpace('V', ldim=2)
#U = FunctionSpace('U', ldim=2)
#
#v = TestFunction(V, name='v')
#u = TestFunction(U, name='u')
#
#a = BilinearForm((v,u), dot(grad(v), grad(u)) + v*u)
#print(evaluate(a))
# ...

# ...
from sympde import grad, dot
from sympde import FunctionSpace
from sympde import TestFunction
from sympde import BilinearForm
from sympde import atomize

V = FunctionSpace('V', ldim=2)
U = FunctionSpace('U', ldim=2)

v = TestFunction(V, name='v')
u = TestFunction(U, name='u')

a = BilinearForm((v,u), dot(grad(v), grad(u)) + v*u)
print(atomize(a.expr))
# ...
