# coding: utf-8

from sympy.core import Basic
from sympy.core.containers import Tuple
from sympy.core import Expr

from pyccel.ast.core import Nil

from sympde.topology.basic import Boundary
from sympde.topology.space import TestFunction
from sympde.topology.space import VectorTestFunction
from sympde.topology.space import FunctionSpace

from .utils import random_string
from .expr import FormCall, BilinearForm, LinearForm
from .errors import ( UnconsistentLhsError, UnconsistentRhsError,
                      UnconsistentArgumentsError, UnconsistentBCError )

#==============================================================================
class BasicBoundaryCondition(Basic):

    def __new__(cls, boundary, value=None):
        return Basic.__new__(cls, boundary, value)

    @property
    def boundary(self):
        return self._args[0]

    @property
    def value(self):
        return self._args[1]

#==============================================================================
class DirichletBC(BasicBoundaryCondition):
    pass


#==============================================================================
# TODO add check on test/trial functions between lhs/rhs
class Equation(Basic):
    def __new__(cls, lhs, rhs, bc=None):
        # ...
        if not isinstance(lhs, FormCall):
            raise UnconsistentLhsError('> lhs must be a call')
        # ...

        # ...
        if not isinstance(rhs, FormCall):
            raise UnconsistentRhsError('> rhs must be a call')
        # ...

        # find unknowns of the equation
        # ...
        if isinstance(lhs.expr, BilinearForm):
            tests_lhs, trials_lhs = lhs.arguments

            # ...
            if isinstance(tests_lhs, (TestFunction, VectorTestFunction)):
                tests_lhs = [tests_lhs]

            elif not isinstance(tests_lhs, (list, tuple, Tuple)):
                msg =  '> Expecting iterable or TestFunction/VectorTestFunction'
                raise UnconsistentArgumentsError(msg)

            tests_lhs = Tuple(*tests_lhs)
            # ...

            # ...
            if isinstance(trials_lhs, (TestFunction, VectorTestFunction)):
                trials_lhs = [trials_lhs]

            elif not isinstance(trials_lhs, (list, tuple, Tuple)):
                msg =  '> Expecting iterable or TestFunction/VectorTestFunction'
                raise UnconsistentArgumentsError(msg)

            trials_lhs = Tuple(*trials_lhs)
            # ...

        else:
            raise UnconsistentLhsError('> lhs must be a bilinear')

        if isinstance(rhs.expr, LinearForm):
            tests_rhs = rhs.arguments

            # ...
            if isinstance(tests_rhs, (TestFunction, VectorTestFunction)):
                tests_rhs = [tests_rhs]

            elif not isinstance(tests_rhs, (list, tuple, Tuple)):
                msg =  '> Expecting iterable or TestFunction/VectorTestFunction'
                raise UnconsistentArgumentsError(msg)

            tests_rhs = Tuple(*tests_rhs)
            # ...

        else:
            raise UnconsistentRhsError('> rhs must be a linear')
        # ...

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

            bc = Tuple(*bc)

            # ... check that the same boundary is not used in the weak
            #     formulation and strong condition
            lhs_bnd = lhs.atoms(Boundary)
            rhs_bnd = rhs.atoms(Boundary)
            bc_bnd  = bc.atoms(Boundary)

            if lhs_bnd & bc_bnd:
                msg = '> {} used for lhs and Dirichlet'.format(lhs_bnd & bc_bnd)
                raise UnconsistentBCError(msg)

            if rhs_bnd & bc_bnd:
                msg = '> {} used for rhs and Dirichlet'.format(rhs_bnd & bc_bnd)
                raise UnconsistentBCError(msg)
            # ...
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
    def __new__(cls, expr, space, kind='l2', mapping=None, bc=None, name=None):
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

            lhs = BilinearForm((v,u), expr_lhs, mapping=mapping, name=lhs_name)
            rhs = LinearForm(v, expr_rhs, mapping=mapping, name=rhs_name)
        # ...

        obj = Equation.__new__(cls, lhs(v,u), rhs(v), bc=bc)
        obj._name = name

    @property
    def name(self):
        return self._name

#==============================================================================
class Interpolation(LambdaEquation):
    def __new__(cls, expr, space, kind='nodal', mapping=None, name=None):
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

