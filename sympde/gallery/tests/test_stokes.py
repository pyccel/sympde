# coding: utf-8

from sympde.topology import Domain
from sympde.gallery import Stokes

#==============================================================================
def test_stokes_2d():

    domain = Domain(r'\Omega', dim=2)
    model = Stokes(domain=domain)

#    model.preview(outputTexFile='test_stokes_2d.tex')

#==============================================================================
def test_stokes_3d():

    domain = Domain(r'\Omega', dim=3)
    model = Stokes(domain=domain)

#    model.preview(outputTexFile='test_stokes_3d.tex')

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
