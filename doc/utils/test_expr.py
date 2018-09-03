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
from sympde import Unknown
from sympde import dx, dy

u = Unknown('u', ldim=2)
v = Unknown('v', ldim=2)

expr = dx(dy(u*v))
dotexport(expr, 'graph2.png')
# ...
