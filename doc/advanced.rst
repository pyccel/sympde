Advanced expression manipulation
********************************

Tensor-product representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``TensorExpr`` lowers an integrated form and expresses it through one-
dimensional symbolic factors. This is useful to inspect the Kronecker
structure consumed by tensor-product discretizations.

.. literalinclude:: examples/tensorization.py
   :language: python

The example is an executable documentation source and is run by the
documentation workflow before Sphinx builds the page.

Logical expressions
^^^^^^^^^^^^^^^^^^^

``LogicalExpr`` pulls expressions on a mapped physical domain back to its
logical domain. The analytical-mapping example demonstrates symbolic and
numerical use of this transformation:

.. literalinclude:: examples/analytical_mappings.py
   :language: python

Nonlinear forms
^^^^^^^^^^^^^^^

``linearize`` differentiates a nonlinear residual with respect to a field and
returns its symbolic Jacobian form:

.. literalinclude:: examples/nonlinear_poisson.py
   :language: python
