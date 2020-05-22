# This is part of MiMiCPy

"""

This module contains the Section, Atom and Input classes that
allows for pythonic creation/manipulation of CPMD scripts

"""

from collections import OrderedDict 
from .base import Script

class Section(Script):
   
    def __str__(self):
        
        val = ''
        
        for d in self.params():
            if getattr(self, d) == None:
                continue
            
            d_ = d.replace('_','-').replace('--', ' ').upper()
            v = str(getattr(self, d)).strip()
            if v == '':
                val += d_+'\n'
            else:
                val += f"{d_}\n{v}\n"
        
        return val
    
class Atom:
    def __init__(self, coords=[], lmax='s', pp='MT_BLYP', labels=''):
        self.coords = coords
        self.pp = pp
        self.labels = labels
        self.lmax = lmax
    
    def __str__(self):
        if not self.pp.startswith('_'): self.pp = '_' + self.pp
        if not self.labels.startswith(' ') and self.labels != '': self.labels = ' ' + self.labels
        val = f'{self.pp}{self.labels}\n'
        val += f'   LMAX={self.lmax.upper()}\n'
        val += f'    {len(self.coords)}\n'
        for d in self.coords:
            val += f'  {d[0]}   {d[1]}   {d[2]}\n'
        val += '\n'
        return val
    
class Input(Script):
    def __init__(self, *args):
        super().__init__()
        self.atoms = OrderedDict()
        self.info = 'MiMiC Run'
        for val in args:
            setattr(self, val, Section())
    
    def checkSection(self, section):
        return self.hasparam(section)
    
    def __str__(self):
        val = ''
        
        for d in self.params():
            if getattr(self, d) == None:
                continue
            elif d.upper() == 'INFO':
                info = f'\n&INFO\n{getattr(self, d)}\n&END\n'
            elif d.upper() == 'ATOMS':
                atoms = '\n&ATOMS\n'
                for k, v in getattr(self, d).items():
                    atoms += f'*{k.replace("*", "")}{str(v)}'
                atoms += '&END\n'
            else:
                v = str(getattr(self, d))
                val += f"\n&{d.upper()}\n{v}&END\n"
        
        return info+val+atoms