SymPDE
======

|CI status|  |binder|  |docs|

**SymPDE** is a symbolic calculus library for partial differential equations and variational forms.
It can be used to have similar capabilities as the fenics_ project, by extending and writing your own *printing* functions.

An example of use can be found in psydac_ or gelato_. 

.. _psydac: https://github.com/pyccel/psydac
.. _gelato: https://github.com/pyccel/gelato
.. _fenics: https://fenicsproject.org/

Install
*******

From PyPi
^^^^^^^^^

Simply run, for a local installation::

  pip3 install --user sympde 

or::

  pip3 install sympde 

for a global installation.

From sources
^^^^^^^^^^^^

* **Standard mode**::

    python3 -m pip install .

* **Development mode**::

    python3 -m pip install --user -e .


For developers
**************

Because many important features of SymPDE are only tested in Psydac, new PRs should also be tested against the test suite of Psydac.
This can be done by opening a PR in Psydac, where the only change consists of installing the corresponding branch of SymPDE.
To achieve this, one just needs to modify the line corresponding to ``sympde`` in the ``pyproject.yaml`` file.

For instance, to test a new SymPDE branch called ``my_feature``, one should write

.. code-block:: python

    # Our packages from PyPi
    'sympde @ https://github.com/pyccel/sympde/archive/refs/heads/my_feature.zip',

Similarly, to test an unreleased version of SymPDE called ``v0.18.4-trunk``, one should write

.. code-block:: python

    # Our packages from PyPi
    'sympde @ https://github.com/pyccel/sympde/archive/refs/tags/v0.18.4-trunk.zip',

Do not forget the comma at the end of the line, as this is an item in a list.
Also, pay attention to the words ``head`` and ``tags`` in the path: the former is used for Git branches, the latter is used for Git tags (which may or may not correspond to GitHub releases).


.. |CI status| image:: https://github.com/pyccel/sympde/actions/workflows/continuous-integration.yml/badge.svg?branch=master&event=push
   :alt: CI status
   :target: https://github.com/pyccel/sympde/actions/workflows/continuous-integration.yml

.. |docs| image:: https://readthedocs.org/projects/sympde/badge/?version=latest
   :alt: Documentation Status
   :target: http://sympde.readthedocs.io/en/latest/?badge=latest

.. |binder| image:: https://mybinder.org/badge_logo.svg
   :alt: Run notebooks in Binder
   :target: https://mybinder.org/v2/gh/pyccel/sympde/master
