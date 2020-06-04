# This is part of MiMiCPy

"""

This module contains the helper functions to efficiently
convert a pdb line to dict, or convert the whole file to a dataframe

"""

from .._global import _Global as gbl
import pandas as pd

keys = ['record', 'serial', 'name', 'altLoc', 'resName', 'chainID', 'resSeq', 'iCode', 'x', 'y', 'z',\
        'occupancy', 'tempFactor', 'element', 'charge']

def readLine(line):
    vals = {}
    vals['record'] = line[:6].strip()
    
    if vals['record'] != 'HETATM' and vals['record'] != 'ATOM':
        vals['content'] = line[6:]
        return vals
    
    vals['serial'] = line[6:11].strip()
    vals['name'] = line[12:16].strip()
    vals['altLoc'] = line[16]
    vals['resName'] = line[17:20].strip()
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

def parseFile(self, pdb, lines=1000):
    f = gbl.host.open(pdb, 'rb')
    
    pdb_lst = []
    
    while True:
        # ATOM/HETATM line is always 78 bytes/chars
        bytes_data = f.read(78*lines)
        if bytes_data == b'': break
    
        for line in bytes_data.decode().splitlines():
            try:
                vals = readLine(line)
            except: # if only part of line was read
                f.seek(-len(line), 1) # push back the file pointer to start of line
            
            if vals['record'] == 'ATOM' or vals['record'] == 'HETATM':
                pdb_lst.append(vals)
                
    return pd.DataFrame(pdb_lst)