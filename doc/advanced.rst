Advanced expression manipulation
********************************

In this section, we shall present different different features for manipulating expressions such us finding the Kronecker representation for bilinear forms.

Tensorization
^^^^^^^^^^^^^

Assume we have the following bilinear form for Laplace

.. code-block:: python

  from sympde import grad, dot
  from sympde import FunctionSpace
  from sympde import TestFunction
  from sympde import BilinearForm
  from sympde.core import tensorize

  V = FunctionSpace('V', ldim=2)
  U = FunctionSpace('U', ldim=2)

  v = TestFunction(V, name='v')
  u = TestFunction(U, name='u')

  a = BilinearForm((v,u), dot(grad(v), grad(u)) + v*u)
  print('> tensorize(a) = ', tensorize(a))

the result is then

.. code-block:: python

  >>> Mass(v1,u1)xMass(v0,u0) + Mass(v1,u1)xStiffness(v0,u0) + Stiffness(v1,u1)xMass(v0,u0) 

which is what we expect, when splitting the integrals with respect to the two directions.
Printing the corresponding latex code gives

.. math::

  {\int_{0}^{1}  u_{1} v_{1} dy}\otimes {\int_{0}^{1}  u_{0} v_{0} dx} + {\int_{0}^{1}  u_{1} v_{1} dy}\otimes {\int_{0}^{1}  u_{0}^\prime v_{0}^\prime dx} + {\int_{0}^{1}  u_{1}^\prime v_{1}^\prime dy}\otimes {\int_{0}^{1}  u_{0} v_{0} dx}
