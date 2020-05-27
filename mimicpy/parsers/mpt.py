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
    atm_types_to_symb.update(nonstd_atm_types)
    ap = _mpt_writer.AtomsParser(file, mols, atm_types_to_symb, buff)
    
    df = ap.mol_df
    
    # replace repeating dataframes with the string name of prev mol
    keys = list(df.keys())
    vals = []
    for i in range(len(keys)):
        key_i = keys[i]
        for j in range(i+1, len(keys)):
            key_j = keys[j]
        
            try:
                if all(df[key_i][3] == df[key_j][3]): 
                    vals.append((key_i, key_j))
            except ValueError:
                continue
    
    for k,v in vals: df[v][3] = k
    
    pkl = gbl.host.open(out, 'wb')
    pickle.dump(df, pkl)
    
class Reader:
    
    def __init__(self, file):
        pkl = gbl.host.open(file, 'rb') # open as bytes
        self.mpt = pickle.load(pkl) # unpickle
        pkl.close()
    
    def _get_df(self, mol):
        third_val = self.mpt[mol][3]
        
        if isinstance(third_val, str):
            return self.mpt[third_val][3]
        else:
            return third_val
    
    def selectByID(self, idx, mol=None, relative=False):
        orig_id = idx
        if not mol:
            mol = ''
            for k,v in self.mpt.items():
                if idx <= v[2]: break
                else: mol = k 
        
        df = self._get_df(mol)
        if not relative:
            atms_before = self.mpt[mol][2]
            idx = idx-atms_before
        
        natms = self.mpt[mol][1]
        
        if natms > 1:
            idx -= 1
            col = idx % natms
            idx = col+1
            
        srs = df.loc[idx]
        
        return srs.append(pd.Series({'mol':mol, 'id':orig_id}))
    
    def selectByIDs(self, ids):
        s = [self.selectAtom(i) for i in ids]
        return pd.concat(s, axis=1).T
    
    def getFull(self):
        # not complete!!
        df = None
        
        for mol in self.mpt:    
            _df = self._get_df(mol)
            if df is None:
                df = _df
            else:
                df = df.append(_df)
        
        return df
    
    # TO DO: add func to get each column seperately
    
    # TO DO: need to make a selection language for current mpt format