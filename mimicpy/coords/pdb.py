"""Module for gro files"""

import pandas as pd
from .base import Coords

class Pdb(Coords):
    """reads pdb files"""

    def __read_line(line):
        vals = {}
        vals['record'] = line[:6].strip()
        
        if vals['record'] != 'HETATM' and vals['record'] != 'ATOM':
            vals['content'] = line[6:]
            return vals
        
        vals['id'] = line[6:11].strip()
        vals['name'] = line[12:16].strip()
        vals['altLoc'] = line[16]
        vals['resname'] = line[17:20].strip()
        vals['chainID'] = line[21]
        vals['resSeq'] = line[22:26].strip()
        vals['iCode'] = line[26]
        vals['x'] = line[30:38].strip()
        vals['y'] = line[38:46].strip()
        vals['z'] = line[46:54].strip()
        vals['occupancy'] = line[54:60].strip()
        vals['tempFactor'] = line[60:66].strip()
        vals['element'] = line[76:78].strip()
        vals['charge'] = line[78:80].strip()
        
        return vals
        
    def _read(self):
        # ATOM/HETATM line is always 78 bytes/chars
        self.file .buffer = 78*self.buffer
        
        pdb_lst = []
        
        for chunk in self.file:
            
            for line in chunk.splitlines():
                try:
                    vals = self.__read_line(line)
                except: # if only part of line was read
                    self.file.f.seek(-len(line), 1) # push back the file pointer to start of line
                
                if vals['record'] == 'ATOM' or vals['record'] == 'HETATM':
                    pdb_lst.append(vals)
                    
        self._coords = pd.DataFrame(pdb_lst)
        
        dims = [0, 0, 0]
        for i, r in enumerate(['x', 'y', 'z']):
            self._coords[r] /= 10 # convert ang to nm
            dims[i] = abs(max(self._coords[r]) - min(self._coords[r])) # find box size
        self._box = dims
        
    def _write(self):
        pass