from importlib import import_module

from .utils import lambdify_sympde


_PLOTTING_EXPORTS = {
    'TopologyVertexIncidence',
    'TopologyVertex',
    'collect_topology_vertices',
    'format_topology',
    'print_topology',
    'plot_domain',
}


__all__ = (
    'TopologyVertexIncidence',
    'TopologyVertex',
    'collect_topology_vertices',
    'format_topology',
    'print_topology',
    'lambdify_sympde',
    'plot_domain',
)


def __getattr__(name):
    if name not in _PLOTTING_EXPORTS:
        raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

    value = getattr(import_module('.plotting', __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | _PLOTTING_EXPORTS)
