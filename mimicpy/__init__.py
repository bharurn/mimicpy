# This is part of MiMiCPy

"""

This initializer script of mimicpy

"""

from .system.ligand import NonStdLigand, StdLigand
from .system.protein import Protein
from .core import prepare, simulate
from .shell import shell
from ._global import _Global as gbl
from .utils.logger import Logger, StdOut
from .utils.errors import MiMiCPyError

def setHost(dirc='.', *args, path=None):
    """Wrapper function to set-up host"""

    closeHost()
    
    if ':' not in dirc:
        gbl.host = shell.Local(dirc, path, *args)
    else:
        gbl.host = shell.Remote(dirc, path, *args)

def getHost():
    """Convenience function to return host"""
    return gbl.host

def setEnv(**kwargs):
    """Set all executable paths/names in _Global"""
    for k,v in kwargs.items():
        if hasattr(gbl, k) and (k != 'host' or k != 'logger'):
            setattr(gbl, k, v)
        else:
            raise MiMiCPyError(f"{k} is not an enviornment executable/path!")

def setLogger(level, redirect=StdOut):
    """Set the logger level and stream"""
    if level == 0: # no output
        gbl.logger.info = None
        gbl.logger.debug = None
        gbl.logger.debug2 = None
    elif level == 1: # only info
        gbl.logger.info = redirect
        gbl.logger.debug = None
        gbl.logger.debug2 = None
    elif level == 2: # some debug commands
        gbl.logger.info = redirect
        gbl.logger.debug = redirect
        gbl.logger.debug2 = None
    elif level == 3: # verbose debug commands
        gbl.logger.info = redirect
        gbl.logger.debug = redirect
        gbl.logger.debug2 = redirect
        
def closeHost():
    """Convenience function to close host"""
    gbl.host.__del__()
    
gbl.host = shell.Local('.', None)
gbl.logger = Logger(info=StdOut(), debug=StdOut(), debug2=None)