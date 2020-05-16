from .system.ligand import NonStdLigand, StdLigand
from .system.protein import Protein
from .core import prepare, simulate
from .shell import shell
import scripts
from ._global import _Global as _global

def setHost(dirc='.', *args, path=None):
    closeHost()
    
    if ':' not in dirc:
        _global.host = shell.Local(dirc, path, *args)
    else:
        _global.host = shell.Remote(dirc, path, *args)

def getHost(): return _global.host

def setEnv(**kwargs):
    for k,v in kwargs.items(): exec(f'_global.{k}="{v}"')
    
def closeHost(): _global.host.__del__()
