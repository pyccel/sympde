# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy.core.singleton import Singleton
from sympy.core.compatibility import with_metaclass

#==============================================================================
class SpaceType(with_metaclass(Singleton, Basic)):
    """Base class representing function space types"""
    pass

class H1SpaceType(SpaceType):
    pass

class HcurlSpaceType(SpaceType):
    pass

class HdivSpaceType(SpaceType):
    pass

class L2SpaceType(SpaceType):
    pass

H1Space    = H1SpaceType()
HcurlSpace = HcurlSpaceType()
HdivSpace  = HdivSpaceType()
L2Space    = L2SpaceType()

dtype_space_registry = {'h1':    H1Space,
                        'hcurl': HcurlSpace,
                        'hdiv':  HdivSpace,
                        'l2':    L2Space}

