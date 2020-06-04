# This is part of MiMiCPy

"""

This module contains NonStdLigand and StdLigand classes to load ligands and protnate/get topology

"""

from ..parsers import pdb as hpdb
from .._global import _Global as _global
from ..utils.errors import EnvNotSetError

class NonStdLigand: 
    
    def __init__(self, pdb, itp, posre, chains, name, elems, atm_types):
        self.pdb = pdb
        self.itp = itp
        self.posre = posre
        self.chains = chains
        self.name = name
        self.elems = elems # for MPT
        self.atm_types = atm_types # for MPT