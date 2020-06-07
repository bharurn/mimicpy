import pickle
from .._global import _Global as gbl
from . import top
import pandas as pd
import numpy as np

def write(topol, mpt, nonstd_atm_types={}, buff=1000, guess_elems=True):
    mols, df = top.read(topol, nonstd_atm_types, buff, guess_elems)

    # replace repeating dataframes with the string name of prev mol 
    
    # equality b/w dataframes checking seems to be incorrect
    # seems to be only checking the size of df
    keys = list(df.keys())
    vals = []
    for i in range(len(keys)):
        key_i = keys[i]
        for j in range(i+1, len(keys)):
            key_j = keys[j]
        
            try:
                if all(df[key_i][1] == df[key_j][1]): 
                    vals.append((key_i, key_j))
            except ValueError:
                continue
    
    for k,v in vals: df[v][1] = k

    pkl = gbl.host.open(mpt, 'wb')
    pickle.dump(mols, pkl)
    pickle.dump(df, pkl)
    
class Reader:
    
    def __init__(self, file):
        pkl = gbl.host.open(file, 'rb') # open as bytes
        self.mol_names = pickle.load(pkl) # unpickle
        self.atom_info = pickle.load(pkl)
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
        print("a :" + str(a))
        print("no:" + str(no))
        if self._res_i%no == 0:
            self._res_before += 1
        self._res_i += 1
        return a+self._res_before
                
    def buildSystemTopology(self):
        molecule_topology = pd.DataFrame()
        for mol, n_mols in self.mol_names:
            _df = self.atom_info[mol][1]
            # repeat the molecule topology n_mol times and preserve the atom order
            _df = pd.DataFrame(np.tile(_df.values, (n_mols, 1)), columns = _df.columns)
            # reset index to consecutive numbering
            _df = _df.reset_index(drop=True)
            
            molecule_topology = molecule_topology.append(_df,  ignore_index=True)

    	# atom id is automatically generated when multipling df
    	# but resid in not, TO DO: resid handling
        molecule_topology['id'] = molecule_topology.index+1

        return molecule_topology.set_index(['id'])
    
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
