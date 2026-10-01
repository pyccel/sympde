# SymPDE

[![CI status](https://github.com/pyccel/sympde/actions/workflows/testing.yml/badge.svg?branch=master&event=push)](https://github.com/pyccel/sympde/actions/workflows/testing.yml)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pyccel/sympde/master)
[![Documentation Status](https://readthedocs.org/projects/sympde/badge/?version=latest)](http://sympde.readthedocs.io/en/latest/?badge=latest)

**SymPDE** is a symbolic calculus library for partial differential equations and variational forms.
It can be used to provide capabilities similar to the [FEniCS](https://fenicsproject.org/) project by extending and writing your own *printing* functions.

Examples of its use can be found in [Psydac](https://github.com/pyccel/psydac) and [Gelato](https://github.com/pyccel/gelato).

## Installation

### Set up a virtual environment

We always recommend working in a Python virtual environment.
To create a new one, we recommend the [`venv`](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment) package:

```bash
python3 -m venv <ENV-PATH>
```

Here, `<ENV-PATH>` is the location where the virtual environment will be created.
A new directory will be created at that location.

To activate the environment from a new terminal session, run:

```bash
source <ENV-PATH>/bin/activate
```

### Option 1: Install from PyPI

Make sure that the preferred virtual environment is activated, then run:

```bash
pip install sympde
```

This downloads the correct version of SymPDE from [PyPI](https://pypi.org/project/sympde/) and installs it in the virtual environment.

### Option 2: Install from sources

First, clone the repository with Git and change to the repository directory:

```bash
git clone https://github.com/pyccel/sympde.git
cd sympde
```

To check out a specific branch, tag, or commit named `<TAG>`, run `git checkout <TAG>`.

- **Static mode**

  Install the source files in the virtual environment with:

  ```bash
  pip install .
  ```

  Further changes to the cloned directory are not reflected in the installed package. This is why we call it a **static** installation.

- **Editable mode**

  To make changes to the library and see them when the package is imported, install SymPDE in **editable** mode:

  ```bash
  pip install --editable ".[test]"
  ```

### Running the tests

The complete test suite can be run from any directory with:

```bash
pytest -n auto --dist loadgroup --pyargs sympde -ra
```

## For developers

Because many important SymPDE features are only tested in Psydac, new pull requests should also be tested against the Psydac test suite.
This can be done by opening a pull request in Psydac whose only change is to install the corresponding SymPDE branch.
To achieve this, modify the line corresponding to `sympde` in Psydac's `pyproject.toml` file.

For instance, to test a new SymPDE branch called `my_feature`, use:

```python
# Our packages from PyPI
'sympde @ https://github.com/pyccel/sympde/archive/refs/heads/my_feature.zip',
```

Similarly, to test an unreleased version of SymPDE called `v0.18.4-trunk`, use:

```python
# Our packages from PyPI
'sympde @ https://github.com/pyccel/sympde/archive/refs/tags/v0.18.4-trunk.zip',
```

Do not forget the comma at the end of the line, as this is an item in a list.
Also note the words `heads` and `tags` in the paths: the former is used for Git branches, while the latter is used for Git tags, which may or may not correspond to GitHub releases.
