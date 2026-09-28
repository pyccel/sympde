from pathlib import Path


def test_non_api_package_initializers_are_empty():
    package_dir = Path(__file__).parents[2]
    initializers = package_dir.rglob('__init__.py')

    for initializer in initializers:
        if initializer.parent != package_dir / 'api':
            assert initializer.read_text() == ''


def test_api_has_only_explicit_exports():
    import sympde.api as api

    assert len(api.__all__) == len(set(api.__all__))
    assert all(hasattr(api, name) for name in api.__all__)
    assert {name for name in vars(api) if not name.startswith('_')} == {
        name for name in api.__all__ if not name.startswith('_')
    }

    namespace = {}
    exec('from sympde.api import *', namespace)
    assert set(namespace).difference({'__builtins__'}) == set(api.__all__)


def test_api_objects_come_from_defining_modules():
    from sympde.api import BilinearForm, Constant, Cube, Mapping, elements_of
    from sympde.api import ScalarFunctionSpace, VectorFunctionSpace
    from sympde.core.basic import Constant as CoreConstant
    from sympde.expr.expr import BilinearForm as ExprBilinearForm
    from sympde.topology.domain import Cube as TopologyCube
    from sympde.topology.mapping import Mapping as TopologyMapping
    from sympde.topology.space import elements_of as topology_elements_of
    from sympde.topology.space import ScalarFunctionSpace as TopologyScalarSpace
    from sympde.topology.space import VectorFunctionSpace as TopologyVectorSpace

    assert BilinearForm is ExprBilinearForm
    assert Constant is CoreConstant
    assert Cube is TopologyCube
    assert Mapping is TopologyMapping
    assert elements_of is topology_elements_of
    assert ScalarFunctionSpace is TopologyScalarSpace
    assert VectorFunctionSpace is TopologyVectorSpace
