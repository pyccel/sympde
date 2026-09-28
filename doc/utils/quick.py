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
#from sympde.api import Unknown
#from sympde.api import grad, div
#
#u = Unknown('u', ldim=2)
#expr = - div(grad(u)) + u
#dotexport(expr, 'graph1.png')
# ...

# ...
#from sympde.api import Unknown
#from sympde.api import dx, dy
#
#u = Unknown('u', ldim=2)
#v = Unknown('v', ldim=2)
#
#expr = dx(dy(u*v))
#dotexport(expr, 'graph2.png')
# ...

# ...
#from sympde.api import Unknown
#from sympde.api import dx, dy
#
#u = Unknown('u', ldim=2)
#v = Unknown('v', ldim=2)
#
#expr = dy(2*u+3*v)
#print(expr)
# ...

# ...
#from sympde.api import Unknown, Constant
#from sympde.api import dx
#
#u = Unknown('u', ldim=1)
#alpha = Constant('alpha')
#
#expr = dx(alpha*u) + dx(dx(2*u))
#print(expr)
# ...

# ...
#from sympde.api import Constant
#from sympde.api import dx, dy
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
#from sympde.api import Constant
#from sympde.api import dx, dy
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
#from sympde.api import grad, dot
#from sympde.api import FunctionSpace
#from sympde.api import TestFunction
#from sympde.api import BilinearForm
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
#from sympde.api import dx
#from sympde.api import FunctionSpace
#from sympde.api import TestFunction
#from sympde.api import BilinearForm
#from sympde.api import Constant
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
#from sympde.api import FunctionSpace
#from sympde.api import TestFunction
#from sympde.api import LinearForm
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
#from sympde.api import grad, div
#from sympde.api import FunctionSpace
#from sympde.api import Field
#from sympde.api import FunctionForm
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
#from sympde.api import grad, dot
#from sympde.api import FunctionSpace
#from sympde.api import TestFunction
#from sympde.api import BilinearForm
#from sympde.api import evaluate
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
#from sympde.api import grad, dot
#from sympde.api import FunctionSpace
#from sympde.api import TestFunction
#from sympde.api import BilinearForm
#from sympde.api import atomize
#
#V = FunctionSpace('V', ldim=2)
#U = FunctionSpace('U', ldim=2)
#
#v = TestFunction(V, name='v')
#u = TestFunction(U, name='u')
#
#a = BilinearForm((v,u), dot(grad(v), grad(u)) + v*u)
#print(atomize(a.expr))
# ...

# ...
#from sympde.api import grad, div
#from sympde.api import Unknown
#from sympde.printing import latex
#
#u = Unknown('u', ldim=2)
#
#print(latex(- div(grad(u)) + u))
# ...

# ...
from sympde.api import grad, dot
from sympde.api import FunctionSpace
from sympde.api import TestFunction
from sympde.api import BilinearForm
from sympde.api import atomize
from sympde.printing import latex

V = FunctionSpace('V', ldim=2)
U = FunctionSpace('U', ldim=2)

v = TestFunction(V, name='v')
u = TestFunction(U, name='u')

a = BilinearForm((v,u), dot(grad(v), grad(u)) + v*u)
print(latex(a))
# ...
