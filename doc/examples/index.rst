Examples
========

These examples are executable Python scripts. The documentation workflow runs
all of them before building the HTML pages, ensuring that their imports and
demonstrated APIs remain current.

Variational forms
-----------------

This example constructs a scalar Poisson problem and a Nitsche-type form on a
two-patch domain.

.. literalinclude:: variational_forms.py
   :language: python
   :linenos:

Differential calculus
---------------------

This example demonstrates the currently supported symbolic vector-calculus
identities.

.. literalinclude:: differential_calculus.py
   :language: python
   :linenos:

Analytical mappings
-------------------

This example evaluates a built-in analytical mapping and defines a custom
mapping from coordinate expressions.

.. literalinclude:: analytical_mappings.py
   :language: python
   :linenos:

Exterior calculus
-----------------

This example demonstrates differential forms, exterior derivatives, adjoint
operators, products, contractions, and Lie derivatives.

.. literalinclude:: exterior_calculus.py
   :language: python
   :linenos:

Nonlinear linearization
-----------------------

This example constructs a nonlinear Poisson residual and derives its Jacobian
form symbolically.

.. literalinclude:: nonlinear_poisson.py
   :language: python
   :linenos:
