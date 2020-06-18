from .._global import _Global as gbl
from .mpt_xdr import pack_strlist, pack_topol_dict, unpack_strlist, unpack_topol_dict
from . import top
import pandas as pd
import numpy as np
import xdrlib

class MPT:
    def __init__(self, mols, topol_dict):
        self.mol_list = mols
        self.topol_dict = topol_dict
    
    @classmethod
    def fromTop(cls, topol, nonstd_atm_types={}, buff=1000, guess_elems=True):
        mols, topol_dict = top.read(topol, nonstd_atm_types, buff, guess_elems)
        return cls(mols, topol_dict)
    
    def write(self, fname):
        """Function to write mpt file
        based on XDR format. Format given below:
        ##Header
         mol names from mol_list
         mol nos. from mol_list
        ##TopolDict
         repeating dict keys
         repeating dict values
         mol name of first entry in dict_df
         col1 of df in mol name
         col2 of df in mol name
         col3 ....
         col4 ....
         ... continue for all columns of df
         mol name of second entry in dict_df
         .... continue for all entries in dict_df
        ##End
        """
        
        packer = xdrlib.Packer()
        
        mol_names, no = list(zip(*self.mol_list)) # unzip list of tuples to get mol_name and nums
        pack_strlist(packer, mol_names) #pack mol names as string list
        packer.pack_list(no, packer.pack_int) #pack num of mols as list of ints
        pack_topol_dict(packer, self.topol_dict) #pack topol dict object
        gbl.host.write(packer.get_buffer(), fname, asbytes=True)
    
    @classmethod
    def fromFile(cls, file):
        unpacker = xdrlib.Unpacker(gbl.host.read(file, asbytes=True)) # open as bytes
        
        # unpack mol list
        mol_names = unpack_strlist(unpacker) # unpack mol names
        nos = unpacker.unpack_list(unpacker.unpack_int) # unpack num of mols
        mol_list = list(zip(mol_names, nos)) # zip together
        topol_dict = unpack_topol_dict(unpacker) # unpack topol_dict
        return cls(mol_list, topol_dict)
    
    ##rewrite all this!!
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
