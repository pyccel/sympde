import warnings
from typing import Protocol, runtime_checkable

from .mapping import DefinedMapping, PointEvaluableMapping, _DefinedMappingEvaluator

__all__ = ('BasicCallableMapping', 'CallableMapping')

_LEGACY_CALLABLE_MAPPING_MSG = (
    'sympde.topology.callable_mapping is deprecated and will be removed in a '
    'future release. Use DefinedMapping.get_callable_mapping() and '
    'is_point_evaluable_mapping from sympde.topology.mapping.'
)


@runtime_checkable
class BasicCallableMapping(PointEvaluableMapping, Protocol):
    """Deprecated compatibility alias for point-evaluable mappings."""
    pass


class CallableMapping:
    """Deprecated compatibility wrapper around a DefinedMapping evaluator."""

    def __init__(self, mapping, **kwargs):
        warnings.warn(_LEGACY_CALLABLE_MAPPING_MSG, DeprecationWarning, stacklevel=2)

        if not isinstance(mapping, DefinedMapping):
            raise TypeError(
                f'mapping must be a DefinedMapping, got {type(mapping)} instead'
            )

        self._mapping = mapping
        if mapping.expressions is None:
            if kwargs:
                raise ValueError(
                    'Cannot override symbolic constants for a mapping without '
                    'analytical expressions.'
                )
            self._impl = mapping.get_callable_mapping()
        else:
            self._impl = _DefinedMappingEvaluator(mapping, **kwargs)

    def __call__(self, *eta):
        return self._impl(*eta)

    def jacobian(self, *eta):
        return self._impl.jacobian(*eta)

    def jacobian_inv(self, *eta):
        return self._impl.jacobian_inv(*eta)

    def metric(self, *eta):
        return self._impl.metric(*eta)

    def metric_det(self, *eta):
        return self._impl.metric_det(*eta)

    @property
    def ldim(self):
        return self._mapping.ldim

    @property
    def pdim(self):
        return self._mapping.pdim

    @property
    def params(self):
        return getattr(self._impl, 'params', {})

    @property
    def symbolic_mapping(self):
        return getattr(self._impl, 'symbolic_mapping', self._mapping)