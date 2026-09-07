# The point-evaluable analytic wrapper ``CallableMapping`` was removed in
# WP06d-2: an ``AnalyticMapping`` lambdifies its own expressions (WP06c) and is
# its own callable mapping. This module now only re-exports
# ``BasicCallableMapping`` (the plain-Python point-evaluation ABC) for the
# psydac call sites that still import it from here; those migrate to
# ``DefinedMapping`` in WP06d-3.

from .mapping import BasicCallableMapping

__all__ = ('BasicCallableMapping',)


def __getattr__(name):
    if name == 'CallableMapping':
        raise AttributeError(
            "CallableMapping was removed in sympde WP06d-2. Analytic mappings "
            "are point-evaluable directly: subclass AnalyticMapping (an "
            "AnalyticMapping instance is its own callable mapping), or attach a "
            "callable with Mapping.set_callable_mapping().")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
