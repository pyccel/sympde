from importlib import import_module
import warnings

from .basic              import *
from .derivatives        import *
from .datatype           import *
from .domain             import *
from .mapping            import *
from .measure            import *
from .space              import *
from .analytical_mapping import *

_LEGACY_CALLABLE_NAMES = {'CallableMapping', 'BasicCallableMapping'}


def __getattr__(name):
	if name in _LEGACY_CALLABLE_NAMES:
		warnings.warn(
			f'Importing {name} from sympde.topology is deprecated and will be removed '
			'in a future release. Import from sympde.topology.mapping instead.',
			DeprecationWarning,
			stacklevel=2,
		)
		legacy_module = import_module('.callable_mapping', __name__)
		return getattr(legacy_module, name)

	raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
