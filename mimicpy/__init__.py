from .shell import shell
from ._version import __version__
from ._authors import __authors__
from ._global import _Global as _global

_global.host = shell.Local('.', None)
