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
from sympde import grad, dot
from sympde import FunctionSpace
from sympde import TestFunction
from sympde import BilinearForm
from sympde.core import tensorize
from sympde.printing.latex import latex

V = FunctionSpace('V', ldim=2)
U = FunctionSpace('U', ldim=2)

v = TestFunction(V, name='v')
u = TestFunction(U, name='u')

a = BilinearForm((v,u), dot(grad(v), grad(u)) + v*u)
print('> a            = ', a)
print('> tensorize(a) = ', tensorize(a))

print(latex(tensorize(a)))

#dotexport(a, 'graph_laplace.png')
# ...
