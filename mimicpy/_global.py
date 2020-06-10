# This is part of MiMiCPy

"""

This module contains the _Global which holds all global variables

"""

class _Global:
    """
    Class to hold all module scope varibale
    All variables are static members, so changes in their values are reflected everywhere
    """
    host = None
    logger = None
    gmx = None
    # mimic_path
    # from this get gmx and cpmd bath
    # source GMXRC.bash from this
    # get gmx force field path
    # get cpmd path
    cpmd = 'cpmd.x'
    cpmd_pp = None