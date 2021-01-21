# coding: utf-8

from sympy.core import Expr
from sympy.core.containers import Tuple

from sympde.core.basic     import Constant
from sympde.topology.space import ScalarFunctionSpace,VectorFunctionSpace
from sympde.topology.space import ScalarTestFunction
from sympde.topology.space import VectorTestFunction
from sympde.topology.space import IndexedTestTrial
from sympde.topology.space import ScalarField, VectorField

#==============================================================================
def _sanitize_arguments(arguments, is_bilinear=False, is_linear=False):

    # ...
    if is_bilinear or is_linear:

        if is_bilinear:
            test_functions = arguments[1]

        elif is_linear:
            test_functions = arguments

        if isinstance(test_functions, (ScalarTestFunction, VectorTestFunction)):
            test_functions = [test_functions]

        elif isinstance(test_functions, (tuple, list, Tuple)):
            are_valid = [isinstance(i, (ScalarTestFunction, VectorTestFunction)) for i in test_functions]

            if not all(are_valid):
                raise TypeError('> Wrong arguments for test functions')

        else:
            msg = 'Wrong type for test function(s). given {}'.format(type(test_functions))
            raise TypeError(msg)

        test_functions = Tuple(*test_functions)
    # ...

    if is_bilinear:

        trial_functions = arguments[0]
        if isinstance(trial_functions, (ScalarTestFunction, VectorTestFunction)):
            trial_functions = [trial_functions]

        elif isinstance(trial_functions, (tuple, list, Tuple)):
            are_valid = [isinstance(i, (ScalarTestFunction, VectorTestFunction)) for i in trial_functions]
            if not all(are_valid):
                raise TypeError('> Wrong arguments for trial functions')
        else:
            msg = 'Wrong type for trial function(s). given {}'.format(type(trial_functions))
            raise TypeError(msg)

        trial_functions = Tuple(*trial_functions)
    # ...

    args = Tuple(trial_functions, test_functions) if is_bilinear else Tuple(*test_functions)
    return args

#==============================================================================
class BasicExpr(Expr):
    is_Function   = True
    is_linear     = False
    is_bilinear   = False
    is_functional = False

    @property
    def fields(self):
        atoms  = self.expr.atoms(ScalarField, VectorField, ScalarTestFunction, VectorTestFunction)
        if self.is_bilinear or self.is_linear:
            args = self.variables
            if self.is_bilinear:
                args = args[0]+args[1]
            fields = tuple(atoms.difference(args))
        else:
            fields = tuple(atoms)
        return fields

    # TODO use .atoms
    @property
    def constants(self):
        ls = self.expr.atoms(Constant)
        # no redanduncy
        return tuple(ls)

#==============================================================================
# TODO check unicity of domain in __new__
class BasicForm(Expr):
    is_Function = True
    is_linear   = False
    is_bilinear = False
    is_functional = False
    is_norm       = False
    _domain     = None
    _ldim        = None
    _body       = None
    _kwargs     = None

    @property
    def fields(self):
        atoms  = self.expr.atoms(ScalarField, VectorField, ScalarTestFunction, VectorTestFunction)
        if self.is_bilinear or self.is_linear:
            args = self.variables
            if self.is_bilinear:
                args = args[0]+args[1]
            fields = tuple(i for i in atoms if i not in args)
        else:
            fields = tuple(atoms)
        return fields

    # TODO use .atoms
    @property
    def constants(self):
        ls = self.expr.atoms(Constant)
        # no redanduncy
        return tuple(ls)

    @property
    def domain(self):
        return self._domain

    @property
    def ldim(self):
        return self._ldim

    def get_free_variables(self):
        if self._kwargs is None:
            fields = self.fields
            constants = self.constants

            _kwargs = {}
            for i in fields + constants:
                _kwargs[i.name] = i

            self._kwargs = _kwargs

        return self._kwargs

    def _update_free_variables(self, **kwargs):

        expr = self.expr

        if not kwargs:
            return expr

        # ... use free variables if given and available
        _kwargs = self.get_free_variables()
        _kwargs_names = list(_kwargs.keys())
        for name, v in kwargs.items():
            if not(name in _kwargs_names):
                raise ValueError('{} is not a free variable'.format(name))

            var = _kwargs[name]
            expr = expr.xreplace({var: v})
        # ...

        return expr
