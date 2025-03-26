SymPDE
======

|CI status|  |binder|  |docs|

**SymPDE** is a symbolic calculus library for partial differential equations and variational forms.
It can be used to have similar capabilities as the fenics_ project, by extending and writing your own *printing* functions.

An example of use can be found in psydac_ or gelato_. 

.. _psydac: https://github.com/pyccel/psydac
.. _gelato: https://github.com/pyccel/gelato
.. _fenics: https://fenicsproject.org/


Installation
************

Set up a virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We always recommend working in a Python virtual environment.
To create a new one we recommend the venv_ package::

  python3 -m venv <ENV-PATH>

.. _venv: https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment

where ``<ENV-PATH>`` is the location to create the virtual environment.
(A new directory will be created at the required location.)

In order to activate the environment from a new terminal session just run the command ::

  source <ENV-PATH>/bin/activate

Option 1: Install from PyPI
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure that the preferred virtual environment is activated. Then simply run ::

  pip3 install sympde

This will download the correct version of SymPDE from PyPI_ and install it in the virtual environment.

.. _PyPI: https://pypi.org/project/sympde/

Option 2: Install from sources
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, clone the repository with Git to download the source files, and change the current directory::

  git clone https://github.com/pyccel/sympde.git
  cd sympde

To check out a specific branch/tag/commit named ``<TAG>``, just use ``git checkout <TAG>``.

* **Static mode**

  To install the source files in the virtual environment just run::

    python3 -m pip install .

  Further changes to the cloned directory are not reflected in the installed package. This is why we call it a **static** installation.

* **Editable mode**

  In order to make changes to the library, and see these changes when the package is imported, SymPDE should be installed in **editable** mode::

    python3 -m pip install --editable .


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
