# coding: utf-8

from sympde.topology import Domain
from sympde.gallery import Poisson

#==============================================================================
def test_poisson_1d():

    domain = Domain(r'\Omega', dim=1)
    model = Poisson(domain)

    assert(str(model.space) == 'V')

#    model.preview(outputTexFile='test_poisson_1d.tex')

#==============================================================================
def test_poisson_2d():

    domain = Domain(r'\Omega', dim=2)
    model = Poisson(domain)

    assert(str(model.space) == 'V')

#    model.preview(outputTexFile='test_poisson_2d.tex')

#==============================================================================
def test_poisson_3d():

    domain = Domain(r'\Omega', dim=3)
    model = Poisson(domain)

    assert(str(model.space) == 'V')

#    model.preview(outputTexFile='test_poisson_3d.tex')


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
