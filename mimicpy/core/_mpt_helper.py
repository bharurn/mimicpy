# This is part of MiMiCPy

"""

This module contains the helper functions and classes
to handle MPT files

"""

from .._global import _Global as _global
import re
import pandas as pd
from ._constants import bohr_rad
import pickle
from ..utils.errors import ParserError, asserter

def _do_combine(line, x, y, z):
    """Code for adding x,y,z to list from gro file line"""
    if line.strip() == '': return False
    
    coords = line[20:].split()
    x.append(coords[0])
    y.append(coords[1])
    z.append(coords[2])
    return True
   
def _combine_remote(no, gro):
    """Code for reading large gro file on remote server"""
    # head always gives better performance than python, esp on remote server
    txt = _global.host.run(f'head -n {no} {gro}') # requires UNIX shell
    lines = (x.group(0) for x in re.finditer(r"^(.*?)\n", txt, re.MULTILINE)) # find \n faster than splitlines
    
    x = []
    y = []
    z = []
        
    for i, line in enumerate(lines):
        if i<=1: continue
        if not _do_combine(line, x, y, z):
            break
    
    return (x,y,z)

def _combine_local(no, gro):
    """Code for reading large gro file on local server"""
    x = []
    y = []
    z = []
    
    with open(gro) as file: # iterate without loading in memeory
        for i in range(no):
            line = next(file)
            if i<=1: continue
            if not _do_combine(line, x, y, z):
                break
        
    return (x,y,z)

def combine(mpt, gro):
    """Combine MPT with gro"""
    # we assume that waters are present after all protein/ligands
    # read gro file until waters, to decrease time of exec
    no = len(mpt)+2 # +2 for the first two lines in gro file
    
    # return list of xyz, open gro depending on local or remote
    if _global.host.isLocal():
        x,y,z = _combine_local(no, gro)
    else:
        x,y,z = _combine_remote(no, gro)
        
    mpt['x'] = pd.Series(x, index=mpt.index)
    mpt['y'] = pd.Series(y, index=mpt.index)
    mpt['z'] = pd.Series(z, index=mpt.index)
    
    last = _global.host.run(f'tail -n 1 {gro}') # fast access of last line, requires UNIX shell
    
    return mpt, [float(v)/bohr_rad for v in last.split()]

def read(mpt, gro=None):
    """Convenience function to unpickle MPT and combine with coords"""
    
    _global.logger.write('debug', f"Reading topology from {mpt}..")
    df_f = _global.host.open(mpt, 'rb') # open as bytes
    df = pickle.load(df_f) # unpickle
    
    # combine with gro file, if it was passed
    if gro:
        _global.logger.write('info', f"Combining with latest coordinates data from {gro}..")
        return combine(df, gro)
    else:
        return df

class MPTWriter:
    """Convenience class to right DF and pickle as MPT"""
    def __init__(self):
        """Class constructor"""
        self.i = 0 # atom number
        self.df_lst = [] # list to be converted to dataframe
    
    def write_row(self, vals):
        """Write a single row of dataframe list"""
        
        self.i += 1
        
        vals.update({'number':self.i}) # write atoms number as number
        vals.update({'resNo': int(vals['resSeq'])}) # rename resSeq to resNo, make int
        #vals.update({'charge': 0}) --> collect charge info from pdb, usually its wrong so no point
        
        self.df_lst.append(vals)
        
    def write(self, name, dirc):
        """Convert list to dataframe and pickle"""
        df = pd.DataFrame(self.df_lst)
        lst = ['serial','record','altLoc','resSeq','iCode','occupancy','x','y','z','tempFactor']
        # error checking, before dropping check if those columns exists
        # if not mean PDB was not parsed correctly
        asserter(set(lst) <= set(df.columns), ParserError, "PDB")
        
        df = df.drop(lst, axis=1) # drop unneeded rows from PDB list
        df = df.set_index(['number'])
        
        MPTWriter.dump(df, name, dirc) # pickle file
        
    @staticmethod
    def dump(df, name, dirc):
        """Convenience function to pickle dataframe"""
        _global.logger.write('debug', f"Saving MiMiCPy topology as {name}..")
        f = _global.host.open(f'{dirc}/{name}', 'wb')
        pickle.dump(df, f)