Quick start
***********

SymPDE builds symbolic differential operators, variational forms, and PDE
problems on explicitly declared domains and function spaces. User-facing
objects are available from :mod:`sympde.api`.

Domains and fields
^^^^^^^^^^^^^^^^^^

A field belongs to a scalar or vector function space. The following complete
example creates two scalar fields on a square and forms a differential
expression from them:

.. literalinclude:: examples/quickstart.py
   :language: python
   :start-after: # [domains-start]
   :end-before: # [domains-end]

Differential operators remain symbolic, while linearity and product rules are
applied automatically. SymPDE also supports ``curl``, ``laplace``, normal
derivatives, traces, and interface operators.

Variational forms
^^^^^^^^^^^^^^^^^

Integrals identify the domain on which an expression is evaluated. A
``BilinearForm`` declares its trial and test fields, and a ``LinearForm``
declares its test field:

.. literalinclude:: examples/quickstart.py
   :language: python
   :start-after: # [forms-start]
   :end-before: # [forms-end]

The :func:`sympde.api.find` helper combines those forms into a symbolic
equation. Boundary conditions can be supplied with ``EssentialBC`` through
the optional ``bc`` argument.

Expression lowering
^^^^^^^^^^^^^^^^^^^

Generic vector-calculus operators can be lowered to coordinate derivatives
for a specific domain with ``TerminalExpr``:

.. literalinclude:: examples/quickstart.py
   :language: python
   :start-after: # [lowering-start]
   :end-before: # [lowering-end]

Mappings and multipatch domains
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mappings turn logical domains into physical domains. Multiple patches can be
joined by declaring their connected sides and orientation. See the executable
:doc:`examples/index` for analytical mappings, interface forms, exterior
calculus, and nonlinear linearization.

Printing
^^^^^^^^

The public ``latex`` function produces LaTeX for expressions and forms:

.. literalinclude:: examples/quickstart.py
   :language: python
   :start-after: # [printing-start]
   :end-before: # [printing-end]

The documentation workflow executes ``quickstart.py`` before building these
pages, so its imports and demonstrated API are checked continuously.
