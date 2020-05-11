from .system.ligand import NonStdLigand, StdLigand
from .system.protein import Protein
from .core import prepare, simulate
from . import analysis
from .analysis import dashboard
from .utils import shell
from .utils import scripts
import mimicpy._global as _global

def setHost(dirc, **kwargs):
    closeHost()
    
    if ':' not in dirc:
        _global.host = shell.Local(dirc)
    else:
        if 'ssh_config' not in kwargs:
            from os.path import expanduser
            ssh_config = expanduser("~")+'/.ssh/config'
        else:
            ssh_config = kwargs['ssh_config']
            del kwargs['ssh_config']
            
        _global.host = shell.Remote(dirc, ssh_config, **kwargs)

def getHost(): return _global.host

def setEnv(**kwargs):
    for k,v in kwargs.items(): exec(f'_global.{k}="{v}"')
    
def closeHost(): _global.host.close()