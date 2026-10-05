# coding: utf-8


from sympy import symbols
from sympy import Tuple
from sympy import Matrix
from sympy import srepr
from sympy import Symbol

from sympde.api import Constant

from sympde.api import d, wedge, ip, jp, delta, hodge
from sympde.api import ZeroForm, OneForm, TwoForm
from sympde.api import ThreeForm, FourForm, FiveForm, SixForm


#==============================================================================
def test_datatype_1():

    assert(ZeroForm < OneForm)
    assert(ZeroForm <= ZeroForm)
    assert(OneForm > ZeroForm)
    assert(ThreeForm <= ThreeForm)


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy.core import cache
    cache.clear_cache()

def teardown_function():
    from sympy.core import cache
    cache.clear_cache()

#test_datatype_1()
