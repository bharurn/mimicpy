from . import shell

class _Global:
    host = shell.Local('.', None)
    gmx = 'gmx'
    cpmd = 'cpmd.x'
    cpmd_pp = None
    gauss = 'g09'
    obabel = 'obabel'
    antechamber = 'antechamber'
    parmchk = 'parmchk2'
    tleap = 'tleap'
    acpype = 'acpype'