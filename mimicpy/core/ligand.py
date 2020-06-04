# This is part of MiMiCPy

"""

This module contains NonStdLigand and StdLigand classes to load ligands and protnate/get topology

"""

from ..parsers import pdb as hpdb
from .._global import _Global as _global
from ..utils.errors import EnvNotSetError

class NonStdLigand: 
    
    def __init__(self, elems, atm_types):
        self.elems = elems
        self.atm_types = atm_types