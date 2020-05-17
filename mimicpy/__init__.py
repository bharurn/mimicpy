from .system.ligand import NonStdLigand, StdLigand
from .system.protein import Protein
from .core import prepare, simulate
from .shell import shell
import scripts
from ._global import _Global as _global
import sys

def setHost(dirc='.', *args, path=None):
    closeHost()
    
    if ':' not in dirc:
        _global.host = shell.Local(dirc, path, *args)
    else:
        _global.host = shell.Remote(dirc, path, *args)

def getHost(): return _global.host

def setEnv(**kwargs):
    for k,v in kwargs.items(): exec(f'_global.{k}="{v}"')

def setLogger(level, redirect=sys.stdout):
    if level == 0:
        _global.logger.info = None
        _global.logger.debug = None
        _global.logger.debug2 = None
    elif level == 1:
        _global.logger.info = redirect
        _global.logger.debug = None
        _global.logger.debug2 = None
    elif level == 2:
        _global.logger.info = redirect
        _global.logger.debug = redirect
        _global.logger.debug2 = None
    elif level == 3:
        _global.logger.info = redirect
        _global.logger.debug = redirect
        _global.logger.debug2 = redirect
        
def closeHost(): _global.host.__del__()
