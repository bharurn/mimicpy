import logging
from .shell import shell
from ._version import __version__
from ._authors import __authors__
from ._global import _Global as gbl

gbl.host = shell.Local('.', None)
logging.basicConfig(filename='mimicpy.log',
                    format='%(asctime)s %(levelname)s %(message)s',
                    filemode='w',
                    level=logging.INFO)
