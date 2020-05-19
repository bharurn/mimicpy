from .._global import _Global as _global
import re
import pandas as pd
from ._constants import bohr_rad
import pickle

def _do_combine(line, x, y, z):
    if line.strip() == '': return False
    
    coords = line[20:].split()
    x.append(coords[0])
    y.append(coords[1])
    z.append(coords[2])
    return True
   
def _combine_remote(no, gro):
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
    no = len(mpt)+2 # +2 for the first two lines in gro file
    
    if _global.host.name == 'localhost':
        x,y,z = _combine_local(no, gro)
    else:
        x,y,z = _combine_remote(no, gro)
        
    mpt['x'] = pd.Series(x, index=mpt.index)
    mpt['y'] = pd.Series(y, index=mpt.index)
    mpt['z'] = pd.Series(z, index=mpt.index)
    
    last = _global.host.run(f'tail -n 1 {gro}') # fast access of last line, requires UNIX shell
    
    return mpt, [float(v)/bohr_rad for v in last.split()]

def read(mpt, gro=None):
    _global.logger.write('debug', f"Reading topology from {mpt}..")
    df_f = _global.host.open(mpt, 'rb')
    df = pickle.load(df_f)
    if gro:
        _global.logger.write('info', f"Combining with latest coordinates data from {gro}..")
        return combine(df, gro)
    else:
        return df

class MPTWriter:
    def __init__(self):
        self.i = 0
        self.df_lst = []
    
    def write_row(self, vals):
        self.i += 1
        vals.update({'number':self.i})
        vals.update({'resNo': int(vals['resSeq'])})
        vals.update({'charge': 0})
        
        self.df_lst.append(vals)
        
    def write(self, name, dirc):
        df = pd.DataFrame(self.df_lst)
        df = df.drop(['serial','record','altLoc','resSeq','iCode','occupancy','x','y','z','tempFactor'], axis=1)
        df = df.set_index(['number'])
        
        MPTWriter.dump(df, name, dirc)
        
    @staticmethod
    def dump(df, name, dirc):
        _global.logger.write('debug', f"Saving MiMiCPy topology as {name}..")
        f = _global.host.open(f'{dirc}/{name}', 'wb')
        pickle.dump(df, f)