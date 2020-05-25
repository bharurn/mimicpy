import pickle
from .._global import _Global as gbl
from . import _mpt_writer
import pandas as pd

def write(inp, out, nonstd_atm_types={}, buff=1000):
    tail = gbl.host.run(f'tail -n 30 {inp}')
    mols = _mpt_writer.molecules(tail)
        
    file = gbl.host.open(inp, 'rb')
    atm_types_to_symb = _mpt_writer.atomtypes(file, buff)
    # extend atm_types_to_symb with nonstdligands
    # should be dict of atom type --> symb
    atm_types_to_symb.extend(nonstd_atm_types)
    ap = _mpt_writer.AtomsParser(file, mols, atm_types_to_symb, buff)
    
    pkl = gbl.host.open(out, 'wb')
    pickle.dump(ap.mol_df, pkl)
    
class Reader:
    
    columns = ['type	', 'resNo','resName','name','charge','element',	'mass',	'mol']
    
    def __init__(self, file):
        pkl = gbl.host.open(file, 'rb') # open as bytes
        self.mpt = pickle.load(pkl) # unpickle
        pkl.close()
    
    def selectAtom(self, idx, mol=None, relative=False):
        orig_id = idx
        if not mol:
            mol = ''
            for k,v in self.mpt.items():
                if idx <= v[2]: break
                else: mol = k 
        
        df = self.mpt[mol][3]
        if not relative:
            atms_before = self.mpt[mol][2]
            idx = idx-atms_before
        
        natms = self.mpt[mol][1]
        
        if natms > 1:
            idx -= 1
            col = idx % natms
            idx = col+1
        #if idx > natms:
        #    idx = idx%natms
        #    if idx == 0:
        #        idx = natms
        
        srs = df.loc[idx]
        return srs.append(pd.Series({'mol':mol, 'id':orig_id}))
    
    def selectAtoms(self, ids):
        s = [self.selectAtom(i) for i in ids]
        return pd.concat(s, axis=1).T