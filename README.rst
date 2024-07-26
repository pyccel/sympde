sympde
======

|build-status|  |binder|  |docs|

**sympde** is a Symbolic calculus library for partial differential equations and variational forms. It can be used to have similar capabilities as the fenics_ project, by extending and writing your own *printing* functions.

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

Because many important features of sympde are only tested in psydac, new PRs 
must also be tested against the test suite of psydac. This can be done by opening
a PR in psydac where the only change consists of installing the corresponding 
branch of sympde.

For instance to test a new sympde branch called ``my_feature``, the following
lines can be used in the file ``.github/workflows/continuous-integration.yml`` of psydac:

.. code-block:: YAML

      - name: Download a development version of sympde
        working-directory: /tmp
        run: |
          wget https://github.com/pyccel/sympde/archive/refs/heads/my_feature.zip
          unzip ./my_feature.zip
          python3 -m pip install ./sympde-my_feature        


.. |build-status| image:: https://travis-ci.com/pyccel/sympde.svg?branch=master
    :alt: build status
    :scale: 100%
    :target:  https://travis-ci.com/pyccel/sympde

.. |docs| image:: https://readthedocs.org/projects/sympde/badge/?version=latest
    :alt: Documentation Status
    :scale: 100%
    :target: http://sympde.readthedocs.io/en/latest/?badge=latest

.. |binder| image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/pyccel/sympde/master
