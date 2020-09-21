# coding: utf-8


# TODO - add W(p,2) spaces and Sobolev of higher order => needed for high order
#        derivatives


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
    name = 'h1'

class HcurlSpaceType(SpaceType):
    name = 'hcurl'

class HdivSpaceType(SpaceType):
    name = 'hdiv'

class L2SpaceType(SpaceType):
    name = 'l2'

class UndefinedSpaceType(SpaceType):
    name = 'undefined'

H1Space        = H1SpaceType()
HcurlSpace     = HcurlSpaceType()
HdivSpace      = HdivSpaceType()
L2Space        = L2SpaceType()
UndefinedSpace = UndefinedSpaceType()

dtype_space_registry = {'h1':        H1Space,
                        'hcurl':     HcurlSpace,
                        'hdiv':      HdivSpace,
                        'l2':        L2Space,
                        'undefined': UndefinedSpace}

#==============================================================================
class RegularityType(with_metaclass(Singleton, Basic)):
    """Base class representing the regularity of a space of functions"""
    _index = None

    @property
    def index(self):
        return self._index

    @property
    def name(self):
        if self.index is None:
            return None

        else:
            return 'C^{}'.format(self.index)

    def __str__(self):
        return str(self.name)

    def __lt__(self, other):
        assert(isinstance(other, RegularityType))
        return self.index < other.index

    def __le__(self, other):
        assert(isinstance(other, RegularityType))
        return self.index <= other.index

    def __gt__(self, other):
        assert(isinstance(other, RegularityType))
        return self.index > other.index

    def __ge__(self, other):
        assert(isinstance(other, RegularityType))
        return self.index >= other.index

class H1RegularityType(RegularityType):
    _index = 0
    pass

class HcurlRegularityType(RegularityType):
    _index = -0.3
    pass

class HdivRegularityType(RegularityType):
    _index = -0.6
    pass

class L2RegularityType(RegularityType):
    _index = -1
    pass

H1Regularity        = H1RegularityType()
HcurlRegularity     = HcurlRegularityType()
HdivRegularity      = HdivRegularityType()
L2Regularity        = L2RegularityType()

dtype_regularity_registry = {'h1':    H1Regularity,
                             'hcurl': HcurlRegularity,
                             'hdiv':  HdivRegularity,
                             'l2':    L2Regularity}

