# coding: utf-8

# TODO not working yet

from sympde.topology import Domain
from sympde.gallery import Wave

#==============================================================================
def test_wave_1d():

    domain = Domain(r'\Omega', dim=1)
    model = Wave(domain=domain)

    assert(str(model.space) == 'V')

#    model.preview(outputTexFile='test_wave_1d.tex')

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
