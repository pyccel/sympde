# -*- coding: UTF-8 -*-
#! /usr/bin/python

import sys
from setuptools import setup, find_packages

NAME    = 'sympde'
VERSION = '0.9.1'
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

def setup_package():
    if 'setuptools' in sys.modules:
        setup_args['install_requires'] = ['numpy', 'sympy']

    setup(packages = packages, \
          include_package_data = True, \
          zip_safe=True, \
          **setup_args)


if __name__ == "__main__":
    setup_package()
