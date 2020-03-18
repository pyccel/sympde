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
            test_functions = arguments[0]

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

    # ...
    if is_bilinear:

        trial_functions = arguments[1]
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

    if is_bilinear:
        args = [test_functions, trial_functions]
        args = Tuple(*args)

    else:
        args = Tuple(*test_functions)

    return args

#==============================================================================
class BasicExpr(Expr):
    is_Function   = True
    is_linear     = False
    is_bilinear   = False
    is_functional = False
    is_annotated  = False

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

    def annotate(self):
    
        if self.is_annotated:
            return self

        if self.is_bilinear or self.is_linear:

            fields          = self.fields
            scalar_fields   = [f for f in fields if isinstance(f.space, ScalarFunctionSpace)]
            vector_fields   = [f for f in fields if isinstance(f.space, VectorFunctionSpace)]

            new_scalar_fields   = [f.space.field(f.name) for f in scalar_fields]
            new_vector_fields   = [f.space.field(f.name) for f in vector_fields]
            
            indexed             = self.expr.atoms(IndexedTestTrial)
            indexed             = [f for f in indexed if f.base in vector_fields]
            new_indexed         = [VectorField(f.base.space, f.base.name)[f.indices[0]] for f in indexed]
            
            expr = self.subs(zip(indexed, new_indexed))
            expr = expr.subs(zip(vector_fields, new_vector_fields))
            expr = expr.subs(zip(scalar_fields, new_scalar_fields))
            
        elif self.is_functional:
            domain = self.domain
            expr   = self.expr
            fields = self.fields
            scalar_fields   = [f for f in fields if isinstance(f.space, ScalarFunctionSpace)]
            vector_fields   = [f for f in fields if isinstance(f.space, VectorFunctionSpace)]

            new_scalar_fields   = [f.space.field(f.name) for f in scalar_fields]
            new_vector_fields   = [f.space.field(f.name) for f in vector_fields]

            indexed             = self.expr.atoms(IndexedTestTrial)
            indexed             = [f for f in indexed if f.base in vector_fields]
            new_indexed         = [VectorField(f.base.space, f.base.name)[f.indices[0]] for f in indexed]
            
            expr = expr.subs(zip(indexed, new_indexed))
            expr = expr.subs(zip(vector_fields, new_vector_fields))
            expr = expr.subs(zip(scalar_fields, new_scalar_fields))
            expr = self.func(expr, self.domain, eval=False)
            if self.is_norm:
                expr._exponent = self._exponent

        expr.is_annotated = True
        return expr

#==============================================================================
# TODO check unicity of domain in __new__
class BasicForm(Expr):
    is_Function = True
    is_linear   = False
    is_bilinear = False
    is_functional = False
    is_annotated  = False
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
        
    def annotate(self):
    
        if self.is_annotated:
            return self

        if self.is_bilinear or self.is_linear:

            fields          = self.fields
            scalar_fields   = [f for f in fields if isinstance(f.space, ScalarFunctionSpace)]
            vector_fields   = [f for f in fields if isinstance(f.space, VectorFunctionSpace)]

            new_scalar_fields   = [f.space.field(f.name) for f in scalar_fields]
            new_vector_fields   = [f.space.field(f.name) for f in vector_fields]
            
            indexed             = self.expr.atoms(IndexedTestTrial)
            indexed             = [f for f in indexed if f.base in vector_fields]
            new_indexed         = [VectorField(f.base.space, f.base.name)[f.indices[0]] for f in indexed]
            
            expr = self.subs(zip(indexed, new_indexed))
            expr = expr.subs(zip(vector_fields, new_vector_fields))
            expr = expr.subs(zip(scalar_fields, new_scalar_fields))
            
        elif self.is_functional:
            domain = self.domain
            expr   = self.expr
            fields = self.fields
            scalar_fields   = [f for f in fields if isinstance(f.space, ScalarFunctionSpace)]
            vector_fields   = [f for f in fields if isinstance(f.space, VectorFunctionSpace)]

            new_scalar_fields   = [f.space.field(f.name) for f in scalar_fields]
            new_vector_fields   = [f.space.field(f.name) for f in vector_fields]

            indexed             = self.expr.atoms(IndexedTestTrial)
            indexed             = [f for f in indexed if f.base in vector_fields]
            new_indexed         = [VectorField(f.base.space, f.base.name)[f.indices[0]] for f in indexed]
            
            expr = expr.subs(zip(indexed, new_indexed))
            expr = expr.subs(zip(vector_fields, new_vector_fields))
            expr = expr.subs(zip(scalar_fields, new_scalar_fields))
            expr = self.func(expr, self.domain, eval=False)
            if self.is_norm:
                expr._exponent = self._exponent

        expr.is_annotated = True
        return expr
