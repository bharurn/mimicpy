from . import shell
from .utils.logger import Logger
import sys

class _Global:
    host = shell.Local('.', None)
    logger = Logger(info=sys.stdout, debug=sys.stdout, debug2=None)
    gmx = 'gmx'
    cpmd = 'cpmd.x'
    cpmd_pp = None
    gauss = 'g09'
    obabel = 'obabel'
    antechamber = 'antechamber'
    parmchk = 'parmchk2'
    tleap = 'tleap'
    acpype = 'acpype'