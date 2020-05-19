from .system.ligand import NonStdLigand, StdLigand
from .system.protein import Protein
from .core import prepare, simulate
from .shell import shell
from ._global import _Global as gbl
from .utils.logger import Logger, StdOut
from .utils.errors import MiMiCPyError

def setHost(dirc='.', *args, path=None):
    closeHost()
    
    if ':' not in dirc:
        gbl.host = shell.Local(dirc, path, *args)
    else:
        gbl.host = shell.Remote(dirc, path, *args)

def getHost(): return gbl.host

def setEnv(**kwargs):
    for k,v in kwargs.items():
        if hasattr(gbl, k) and (k != 'host' or k != 'logger'):
            setattr(gbl, k, v)
        else:
            raise MiMiCPyError(f"{k} is not an enviornment executable/path!")

def setLogger(level, redirect=StdOut):
    if level == 0:
        gbl.logger.info = None
        gbl.logger.debug = None
        gbl.logger.debug2 = None
    elif level == 1:
        gbl.logger.info = redirect
        gbl.logger.debug = None
        gbl.logger.debug2 = None
    elif level == 2:
        gbl.logger.info = redirect
        gbl.logger.debug = redirect
        gbl.logger.debug2 = None
    elif level == 3:
        gbl.logger.info = redirect
        gbl.logger.debug = redirect
        gbl.logger.debug2 = redirect
        
def closeHost(): gbl.host.__del__()

gbl.host = shell.Local('.', None)
gbl.logger = Logger(info=StdOut(), debug=StdOut(), debug2=None)