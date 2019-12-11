# coding: utf-8

from sympy.core import Basic
from sympy.core.containers import Tuple
from sympy.core import Expr
from sympy.parsing.sympy_parser import parse_expr

from sympde.topology.basic import Boundary, Union
from sympde.topology.space  import VectorFunctionSpace, ScalarFunctionSpace
from sympde.topology.space import ScalarTestFunction
from sympde.topology.space import VectorTestFunction, IndexedTestTrial
from sympde.topology import Boundary, NormalVector, TangentVector
from sympde.topology import Trace, trace_0, trace_1
from sympde.calculus import grad, dot
from sympde.core.utils import random_string

from .expr import BilinearForm, LinearForm
from .expr import linearize
from .errors import ( UnconsistentLhsError, UnconsistentRhsError,
                      UnconsistentArgumentsError, UnconsistentBCError )


#==============================================================================
class BasicBoundaryCondition(Basic):
    pass

class BasicConstraint(Basic):
    pass
    
#==============================================================================
class EssentialBC(BasicBoundaryCondition):
    _order = None
    _variable = None
    _normal_component = None
    _position = None

    def __new__(cls, lhs, rhs, boundary, position=None, index_component=None):
        # ...
        normal_component = False
        # ...

        # ...
        indexed = list(lhs.atoms(IndexedTestTrial))

        u  = list(lhs.atoms(ScalarTestFunction))

        if not indexed:
            u += list(lhs.atoms(VectorTestFunction))

        else:
            u += indexed
        
        if not( len(u) == 1 ):
            raise ValueError('Expecting one test function')

        u = u[0]
        # ...

        # ... do not allow to use Trace operator
        trace = list(lhs.atoms(Trace))
        if trace:
            raise TypeError('Trace operator is not allowed')
        # ...

        # ...
        order_0_expr = [u]
        order_1_expr = []
        # ...

        # ...
        nn               = list(lhs.atoms(NormalVector))
        normal_component = isinstance(u, VectorTestFunction) and (len(nn) > 0)
        # ...
        # ...
        if nn:
            assert(len(nn) == 1)

            nn = nn[0]
            if isinstance(u.space, VectorFunctionSpace):
                order_0_expr += [dot(u, nn)]

            order_1_expr += [dot(grad(u), nn)]
        # ...

        # ...
        if lhs in order_0_expr:
            order = 0
            if isinstance(u, (IndexedTestTrial)):
                variable = u.base
                index_component = list(u.indices)

            elif isinstance(u, VectorTestFunction) and not normal_component \
                 and isinstance(u.space, VectorFunctionSpace)\
                 and not normal_component:
                variable = u
                index_component = list(range(u.ldim))

            else:
                variable = u
        elif lhs in order_1_expr:
            order = 1
            variable = u

            if isinstance(u, IndexedTestTrial):
                raise NotImplementedError('Indexed case')
        else:
            # TODO change error to unconsistent error
            raise ValueError('Wrong lhs')
        # ...

#        # ... for simple geometries we can compute the indices for the normal
#        # compoenent
#        if normal_component and order == 0:
#            d = boundary.domain.dtype
#            if d:
#                if d['type'] in ['Line', 'Square', 'Cube']:
#                    print('ICI')
#                    index_component = boundary.axis
#                    # TODO shall we use the ext for the sign?
#        # ...

        obj = Basic.__new__(cls, lhs, rhs, boundary)

        obj._order = order
        obj._variable = variable
        obj._normal_component = normal_component
        obj._index_component = index_component
        obj._position = position

        return obj

    @property
    def lhs(self):
        return self._args[0]

    @property
    def rhs(self):
        return self._args[1]

    @property
    def boundary(self):
        return self._args[2]

    @property
    def order(self):
        return self._order

    @property
    def variable(self):
        return self._variable

    @property
    def normal_component(self):
        return self._normal_component

    @property
    def index_component(self):
        return self._index_component

    @property
    def position(self):
        return self._position

    def set_position(self, value):
        self._position = value

#==============================================================================
class Mean(BasicConstraint):

    def __new__(cls, lhs, rhs):
        assert isinstance(lhs,  (ScalarTestFunction, VectorTestFunction))
        rhs = parse_expr(str(rhs))
        if not isinstance(rhs, Expr):
            raise ValueError('only sympy Expr are accepted')

        return Basic.__new__(cls, lhs, rhs)

    @property
    def lhs(self):
        return self._args[0]

    @property
    def rhs(self):
        return self._args[1]

#==============================================================================
# TODO - add check on test/trial functions between lhs/rhs
#      - must be improved: use the bc wrt the variable
#        check that the same boundary is not used in the weak
#        formulation and strong condition
class Equation(Basic):

    def __new__(cls, lhs, rhs, trials, tests, bc=None, constraint=None):
        # ...
        if not isinstance(lhs, BilinearForm):
            raise UnconsistentLhsError('> lhs must be a bilinear')

        if not isinstance(rhs, LinearForm):
            raise UnconsistentRhsError('> rhs must be a linear')
        # ...

        # ...
        _is_test_function = lambda u: isinstance(u, (ScalarTestFunction, VectorTestFunction))
        # ...

        # ...
        if isinstance(tests, (ScalarTestFunction, VectorTestFunction)):
            tests = [tests]

        else:
            assert(isinstance(tests, (list, tuple, Tuple)))
            tests = [*tests]

            assert(all([_is_test_function(i) for i in tests]))
        # ...

        # ...
        if isinstance(trials, (ScalarTestFunction, VectorTestFunction)):
            trials = [trials]

        else:
            assert(isinstance(trials, (list, tuple, Tuple)))
            trials = [*trials]

            assert(all([_is_test_function(i) for i in trials]))
        # ...

#        # ...
#        # find unknowns and tests of the equation
#        # ...
#        tests_lhs, trials_lhs = lhs.variables
#        if isinstance(tests_lhs, (ScalarTestFunction, VectorTestFunction)):
#            tests_lhs = [tests_lhs]
#
#        elif not isinstance(tests_lhs, (list, tuple, Tuple)):
#            msg =  '> Expecting iterable or ScalarTestFunction/VectorTestFunction'
#            raise UnconsistentArgumentsError(msg)
#
#        tests_lhs = Tuple(*tests_lhs)
#
#        if isinstance(trials_lhs, (ScalarTestFunction, VectorTestFunction)):
#            trials_lhs = [trials_lhs]
#
#        elif not isinstance(trials_lhs, (list, tuple, Tuple)):
#            msg =  '> Expecting iterable or ScalarTestFunction/VectorTestFunction'
#            raise UnconsistentArgumentsError(msg)
#
#        trials_lhs = Tuple(*trials_lhs)
#        # ...
#
#        # ... find test functions
#        tests_rhs = rhs.variables
#        if isinstance(tests_rhs, (ScalarTestFunction, VectorTestFunction)):
#            tests_rhs = [tests_rhs]
#
#        elif not isinstance(tests_rhs, (list, tuple, Tuple)):
#            msg =  '> Expecting iterable or ScalarTestFunction/VectorTestFunction'
#            raise UnconsistentArgumentsError(msg)
#
#        tests_rhs = Tuple(*tests_rhs)
#        # ...
#
#        # ...
#        for u_lhs, u_rhs in zip(tests_lhs, tests_rhs):
#            if not( u_lhs is u_rhs ):
#                msg = '> lhs and rhs must have the same test function. '
#                msg += 'given {lhs} & {rhs}'.format(lhs=u_lhs, rhs=u_rhs)
#                raise UnconsistentArgumentsError(msg)
#        # ...

        # ...
        if constraint:
            if isinstance(constraint, BasicConstraint):
                constraint = [constraint]
            elif isinstance(constraint, (list, tuple, Tuple)):
                for i in bc:
                    if not isinstance(i, BasicConstraint):
                        msg = '> Expecting a list of BasicConstraint'
                        raise TypeError(msg)

            else:
                raise TypeError('> Wrong type for bc')

            for i in constraint:
                if not isinstance(i, Mean):
                    raise NotImplementedError('TODO')

            constraint = Tuple(*constraint)

        if bc:
            if isinstance(bc, BasicBoundaryCondition):
                bc = [bc]

            elif isinstance(bc, (list, tuple, Tuple)):
                for i in bc:
                    if not isinstance(i, BasicBoundaryCondition):
                        msg = '> Expecting a list of BasicBoundaryCondition'
                        raise TypeError(msg)

            else:
                raise TypeError('> Wrong type for bc')

            newbc = []
            for i in bc:
                if not isinstance(i, EssentialBC):
                    raise NotImplementedError('')

                if isinstance(i, EssentialBC):
                    if not( i.variable in trials ):
                        msg = 'Essential bc must be on trial functions'
                        raise UnconsistentArgumentsError(msg)

                    else:
                        # TODO treate case of vector test function
                        position = trials.index(i.variable)
                        i.set_position(position)

                if isinstance(i.boundary, Union):
                    if isinstance(i, EssentialBC):
                        newbc += [EssentialBC(i.lhs, i.rhs, j, position=i.position,
                                              index_component=i.index_component)
                                  for j in i.boundary._args]

                else:
                    newbc += [i]

            bc = Tuple(*newbc)
        # ...

        # ... sympify tests/trials
        trials = Tuple(*trials)
        tests  = Tuple(*tests)
        # ...

        return Basic.__new__(cls, lhs, rhs, trials, tests, bc, constraint)

    @property
    def lhs(self):
        return self._args[0]

    @property
    def rhs(self):
        return self._args[1]

    @property
    def trial_functions(self):
        return self._args[2]

    @property
    def test_functions(self):
        return self._args[3]

    @property
    def bc(self):
        return self._args[4]

    @property
    def constraint(self):
        return self._args[5]

#==============================================================================
# TODO must subtitute expr by given args => call then create BasicForm
class NewtonIteration(Equation):

    def __new__(cls, form, fields, bc=None, trials=None):

        assert( isinstance(form, LinearForm) )

        a = linearize(form, fields, trials=trials)

        trials, tests = a.variables

        lhs = a
        rhs = LinearForm(tests, -form.expr)

        return Equation.__new__(cls, lhs, rhs, trials, tests, bc)

#==============================================================================
# user friendly function to create Equation objects
def find(trials, *, forall, lhs, rhs, bc=None, constraint=None):

    tests = forall
    lhs = BilinearForm((trials, tests), lhs)
    rhs =   LinearForm(         tests , rhs)

    return Equation(lhs, rhs, trials, tests, bc=bc, constraint=constraint)
