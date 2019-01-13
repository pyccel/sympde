# coding: utf-8

from sympy.core import Basic
from sympy.core.containers import Tuple
from sympy.core import Expr

from pyccel.ast.core import Nil

from sympde.topology.basic import Boundary, Union
from sympde.topology.space import TestFunction
from sympde.topology.space import VectorTestFunction, IndexedTestTrial
from sympde.topology.space import FunctionSpace
from sympde.topology import Boundary, NormalVector, TangentVector
from sympde.topology import Trace, trace_0, trace_1
from sympde.core import grad, dot
from sympde.core.utils import random_string

from .expr import FormCall, BilinearForm, LinearForm
from .errors import ( UnconsistentLhsError, UnconsistentRhsError,
                      UnconsistentArgumentsError, UnconsistentBCError )

#==============================================================================
class BasicBoundaryCondition(Basic):
    pass

#==============================================================================
class DirichletBC(BasicBoundaryCondition):

    def __new__(cls, boundary, value=None):
        return Basic.__new__(cls, boundary, value)

    @property
    def boundary(self):
        return self._args[0]

    @property
    def value(self):
        return self._args[1]

#==============================================================================
class EssentialBC(BasicBoundaryCondition):
    _order = None
    _variable = None
    _normal_component = None
    _position = None

    def __new__(cls, lhs, rhs, boundary, position=None):
        # ...
        if not( rhs == 0 ):
            raise NotImplementedError('Only homogeneous case is available')
        # ...

        # ...
        normal_component = False
        index_component = None
        # ...

        # ...
        indexed = list(lhs.atoms(IndexedTestTrial))

        u  = list(lhs.atoms(TestFunction))
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
            if isinstance(u, VectorTestFunction):
                order_0_expr += [dot(u, nn)]

            order_1_expr += [dot(grad(u), nn)]
        # ...

        # ...
        if lhs in order_0_expr:
            order = 0
            if isinstance(u, IndexedTestTrial):
                variable = u.base
                index_component = list(u.indices)

            elif isinstance(u, VectorTestFunction) and not normal_component:
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
            print(order_1_expr)
            print(lhs)
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
# TODO add check on test/trial functions between lhs/rhs
class Equation(Basic):
    def __new__(cls, lhs, rhs, bc=None):
        # ...
        if not isinstance(lhs, FormCall):
            raise UnconsistentLhsError('> lhs must be a call')

        if not isinstance(rhs, FormCall):
            raise UnconsistentRhsError('> rhs must be a call')
        # ...

        # ...
        if not isinstance(lhs.expr, BilinearForm):
            raise UnconsistentLhsError('> lhs must be a bilinear')

        if not isinstance(rhs.expr, LinearForm):
            raise UnconsistentRhsError('> rhs must be a linear')
        # ...

        # ...
        # find unknowns and tests of the equation
        # ...
        tests_lhs, trials_lhs = lhs.arguments
        if isinstance(tests_lhs, (TestFunction, VectorTestFunction)):
            tests_lhs = [tests_lhs]

        elif not isinstance(tests_lhs, (list, tuple, Tuple)):
            msg =  '> Expecting iterable or TestFunction/VectorTestFunction'
            raise UnconsistentArgumentsError(msg)

        tests_lhs = Tuple(*tests_lhs)

        if isinstance(trials_lhs, (TestFunction, VectorTestFunction)):
            trials_lhs = [trials_lhs]

        elif not isinstance(trials_lhs, (list, tuple, Tuple)):
            msg =  '> Expecting iterable or TestFunction/VectorTestFunction'
            raise UnconsistentArgumentsError(msg)

        trials_lhs = Tuple(*trials_lhs)
        # ...

        # ... find test functions
        tests_rhs = rhs.arguments
        if isinstance(tests_rhs, (TestFunction, VectorTestFunction)):
            tests_rhs = [tests_rhs]

        elif not isinstance(tests_rhs, (list, tuple, Tuple)):
            msg =  '> Expecting iterable or TestFunction/VectorTestFunction'
            raise UnconsistentArgumentsError(msg)

        tests_rhs = Tuple(*tests_rhs)
        # ...

        # ...
        for u_lhs, u_rhs in zip(tests_lhs, tests_rhs):
            if not( u_lhs is u_rhs ):
                msg = '> lhs and rhs must have the same test function. '
                msg += 'given {lhs} & {rhs}'.format(lhs=u_lhs, rhs=u_rhs)
                raise UnconsistentArgumentsError(msg)
        # ...

        # ...
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
                if not isinstance(i, (DirichletBC, EssentialBC)):
                    raise NotImplementedError('')

                if isinstance(i, EssentialBC):
                    if not( i.variable in trials_lhs ):
                        msg = 'Essential bc must be on trial functions'
                        raise UnconsistentArgumentsError(msg)

                    else:
                        # TODO treate case of vector test function
                        position = trials_lhs.index(i.variable)
                        i.set_position(position)

                if isinstance(i.boundary, Union):
                    if isinstance(i, DirichletBC):
                        newbc += [DirichletBC(j) for j in i.boundary._args]

                    if isinstance(i, EssentialBC):
                        newbc += [EssentialBC(i.lhs, i.rhs, j, position=i.position)
                                  for j in i.boundary._args]

                else:
                    newbc += [i]

            bc = Tuple(*newbc)

#            # TODO must be improved: use the bc wrt the variable
#            # ... check that the same boundary is not used in the weak
#            #     formulation and strong condition
#            lhs_bnd = lhs.atoms(Boundary)
#            rhs_bnd = rhs.atoms(Boundary)
#            bc_bnd  = bc.atoms(Boundary)
#
#            if lhs_bnd & bc_bnd:
#                msg = '> {} used for lhs and Dirichlet'.format(lhs_bnd & bc_bnd)
#                raise UnconsistentBCError(msg)
#
#            if rhs_bnd & bc_bnd:
#                msg = '> {} used for rhs and Dirichlet'.format(rhs_bnd & bc_bnd)
#                raise UnconsistentBCError(msg)
#            # ...
        # ...

        return Basic.__new__(cls, lhs, rhs, trials_lhs, bc)

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
    def bc(self):
        return self._args[3]

    @property
    def test_functions(self):
        return self.rhs.arguments

    @property
    def is_undefined(self):
        return isinstance(self.lhs, Nil) or isinstance(self.rhs, Nil)

#==============================================================================
class LambdaEquation(Equation):

    @property
    def variables(self):
        return self.trial_functions

#==============================================================================
class Projection(LambdaEquation):
    def __new__(cls, expr, space, kind='l2', bc=None, name=None):
        # ...
        tests = expr.atoms((TestFunction, VectorTestFunction))
        if tests or not isinstance(expr, Expr):
            msg = '> Expecting an Expression without test functions'
            raise UnconsistentArgumentsError(msg)
        # ...

        # ...
        if not isinstance(space, FunctionSpace):
            raise UnconsistentArgumentsError('> Expecting a FunctionSpace')
        # ...

        # ...
        if not(kind in ['l2']):
            raise ValueError('> Only L2 projector is available')
        # ...

        # ... defining the lhs and rhs
        V = space
        if kind == 'l2':
            tag = random_string( 3 )
            v_name = 'v_{}'.format(tag)
            u_name = 'u_{}'.format(tag)
            lhs_name = 'lhs_{}'.format(tag)
            rhs_name = 'rhs_{}'.format(tag)
            if V.shape == 1:
                v = TestFunction(V, name=v_name)
                u = TestFunction(V, name=u_name)

                expr_lhs = v*u
                expr_rhs = expr*v

            else:
                v = VectorTestFunction(V, name=v_name)
                u = VectorTestFunction(V, name=u_name)

                expr_lhs = Dot(v,u)
                expr_rhs = Dot(expr, v)

            lhs = BilinearForm((v,u), expr_lhs, name=lhs_name)
            rhs = LinearForm(v, expr_rhs, name=rhs_name)
        # ...

        obj = Equation.__new__(cls, lhs(v,u), rhs(v), bc=bc)
        obj._name = name

    @property
    def name(self):
        return self._name

#==============================================================================
class Interpolation(LambdaEquation):
    def __new__(cls, expr, space, kind='nodal', name=None):
        raise NotImplementedError('TODO')

        # ...
        tests = expr.atoms((TestFunction, VectorTestFunction))
        if tests or not isinstance(expr, Expr):
            msg = '> Expecting an Expression without test functions'
            raise UnconsistentArgumentsError(msg)
        # ...

        # ...
        if not isinstance(space, FunctionSpace):
            raise UnconsistentArgumentsError('> Expecting a FunctionSpace')
        # ...

        # ...
        if not(kind in ['nodal']):
            raise ValueError('> Only nodal interpolation is available')
        # ...

        # ... defining the lhs and rhs
        V = space
        # ...

    @property
    def name(self):
        return self._name

