import pickle
from .._global import _Global as gbl
from . import _mpt_writer
import pandas as pd

def write(inp, out, nonstd_atm_types={}, buff=1000, guess_elems=True):
    tail = gbl.host.run(f'tail -n 30 {inp}')
    mols = _mpt_writer.molecules(tail)
        
    file = gbl.host.open(inp, 'rb')
    atm_types_to_symb = _mpt_writer.atomtypes(file, buff)
    # extend atm_types_to_symb with nonstdligands
    # should be dict of atom type --> symb
    atm_types_to_symb.update(nonstd_atm_types)
    ap = _mpt_writer.AtomsParser(file, mols, atm_types_to_symb, buff, guess_elems)
    
    df = ap.mol_df
    
    # replace repeating dataframes with the string name of prev mol 
    # I strongly belive that Python handles this on its own
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
    
    previous_atoms = 0
    for i in range(len(mols)):
        mols[i][1] = [mols[i][1], df[mols[i][0]][1], previous_atoms, df[mols[i][0]][3]]
        previous_atoms += mols[i][1][0] * mols[i][1][1]

    pkl = gbl.host.open(out, 'wb')
#    pickle.dump(df, pkl)
    pickle.dump(mols, pkl)
    
class Reader:
    
    def __init__(self, file):
        pkl = gbl.host.open(file, 'rb') # open as bytes
        self.mpt = pickle.load(pkl) # unpickle
        pkl.close()
    
    def selectByID(self, idx, mol=None, relative=False):
        orig_id = idx
        if not mol:
            mol = ''
            resn = 0
            curr_res = 0
            for k,v in self.mpt.items():
                if idx <= v[2]:
                    break
                else:
                    mol = k
                    # to get correct resid
                    curr_res = self._get_df(k)['resid'].iloc[-1]
                    resn += curr_res # no. of res till now
        
        resn -= curr_res # resn includes no. of res for current mol, remove it
        
        df = self._get_df(mol)
        if not relative:
            atms_before = self.mpt[mol][2]
            idx = idx-atms_before
        
        ## Accounting for multiple molecules
        natms = self.mpt[mol][1]
        
        idx -= 1
        
        if self.mpt[mol][0] > 1: # to get correct resid
            n_res_before = idx//natms
            resn += n_res_before
        
        # atom id for multiple molecules case
        col = idx % natms
        idx = col+1
            
        srs = df.loc[idx]
        resn += srs['resid']
        # drop resid and assing it again, to avoid pandas warning
        srs = srs.drop(labels=['resid'])
        
        return srs.append(pd.Series({'mol':mol, 'id':orig_id, 'resid': resn}))
    
    def selectByIDs(self, ids):
        s = [self.selectAtom(i) for i in ids]
        return pd.concat(s, axis=1).T
    
    def r(self, a, no):
        """Function to keep track of res counter in getDF()"""
        if self._res_i%no == 0:
            self._res_before += 1
        self._res_i += 1
        return a+self._res_before
    
    
    def getDF(self):
        df = None
        resn = 0
        for mol in self.mpt:
#            _df = self._get_df(mol)
            _df = mol[1][3]
#            no = self.mpt[mol][0]
            no = mol[1][0]
            _df = _df.loc[_df.index.repeat(no)].reset_index(drop=True)
            _df['resid'] += resn
            
            if no > 1:
                self._res_i = 0 # res counter
                self._res_before = -1 # residues before current one
                
  #              _df['resid'] = _df['resid'].apply(self.r, args=[self.mpt[mol][1]])
                _df['resid'] = _df['resid'].apply(self.r, args=[mol[1][1]])
            
            resn = _df['resid'].iloc[-1]
            if df is None:
                df = _df
            else:
                df = df.append(_df)
        
    	# atom id is automatically generated when multipling df
    	# but resid in not, TO DO: resid handling
        df['id'] = df.index+1

        return df.set_index(['id'])
    
    def getProperty(self, prop):
        df = None
        
        for mol in self.mpt:    
            _df = self._get_df(mol)[prop]
            no = self.mpt[mol][0]
            if df is None:
                df = pd.concat([_df]*no)
            else:
                df = df.append(pd.concat([_df]*no))
        
        return df.to_list()
