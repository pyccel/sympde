# -*- coding: UTF-8 -*-
#! /usr/bin/python

from pathlib import Path
from setuptools import setup, find_packages

# ...
# Read library version into '__version__' variable
path = Path(__file__).parent / 'sympde' / 'version.py'
exec(path.read_text())
# ...

NAME    = 'sympde'
VERSION = __version__
AUTHOR  = 'Ahmed Ratnani'
EMAIL   = 'ratnaniahmed@gmail.com'
URL     = 'https://github.com/pyccel/sympde'
DESCR   = 'Symbolic calculus for partial differential equations (and variational forms).'
KEYWORDS = ['math']
LICENSE = "LICENSE"

setup_args = dict(
    name                 = NAME,
    version              = VERSION,
    description          = DESCR,
    long_description     = open('README.rst').read(),
    author               = AUTHOR,
    author_email         = EMAIL,
    license              = LICENSE,
    keywords             = KEYWORDS,
    url                  = URL,
#    download_url     = URL+'/tarball/master',
)

# ...
packages = find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"])
# ...

# Dependencies
install_requires = [
    'sympy>=1.2,<1.7',
    'h5py',
    'pytest',
    'pyyaml',
    'yamlloader'
]

def setup_package():

    setup(packages = packages,
          include_package_data = True,
          install_requires = install_requires,
          zip_safe = True,
          **setup_args)


if __name__ == "__main__":
    setup_package()
