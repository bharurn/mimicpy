from .._global import _Global as gbl
from .mpt_xdr import pack_strlist, pack_topol_dict, unpack_strlist, unpack_topol_dict
from ..utils.errors import SelectionError, MiMiCPyError
from . import top
import pandas as pd
import numpy as np
import xdrlib

columns = top.top_reader.ITPParser.columns.copy() # copy it, otherwise editing will change both
columns.remove('number') # remove 'number'
columns.append('mol')

class MPT:    
    def __init__(self, mols, topol_dict, mode='r'):
        self.__mol_list = mols
        self.__topol_dict = topol_dict
        self.__all_data = None
        self.__columns = columns
        
        if mode == 'r':
            # preload data into memory
            self.__generateData()
        elif mode != 'w':
            raise MiMiCPyError(f"{mode} not a mode. Only r/w can be used")
    
    @classmethod
    def fromTop(cls, topol, nonstd_atm_types={}, buff=1000, guess_elems=True, mode='r'):
        mols, topol_dict = top.read(topol, nonstd_atm_types, buff, guess_elems)
        return cls(mols, topol_dict, mode)
    
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
        
        mol_names, no = list(zip(*self.__mol_list)) # unzip list of tuples to get mol_name and nums
        pack_strlist(packer, mol_names) #pack mol names as string list
        packer.pack_list(no, packer.pack_int) #pack num of mols as list of ints
        pack_topol_dict(packer, self.__topol_dict) #pack topol dict object
        gbl.host.write(packer.get_buffer(), fname, asbytes=True)
    
    @classmethod
    def fromFile(cls, file, mode='r'):
        unpacker = xdrlib.Unpacker(gbl.host.read(file, asbytes=True)) # open as bytes
        
        # unpack mol list
        mol_names = unpack_strlist(unpacker) # unpack mol names
        nos = unpacker.unpack_list(unpacker.unpack_int) # unpack num of mols
        mol_list = list(zip(mol_names, nos)) # zip together
        topol_dict = unpack_topol_dict(unpacker) # unpack topol_dict
        return cls(mol_list, topol_dict, mode)
    
    def __getitem__(self, key):
        """Select an atom by passing the atom ID to key
            it can be a single int, list or a slice
            The indexing starts from 1
            If a string is passed as key then that property is returned
        """
        ## If string --> return property
        if isinstance(key, str): return self.__getProperty(key)
        
        ##Else it is atom ID that is requested
        
        elif isinstance(key, int): key = [key] 
        elif isinstance(key, slice): key = list(range(key.stop)[key])
        
        return self.__selectbyID(key)
    
    def __generateData(self):
        """Generate list of all data from topol_dict"""
        # generate everything as python lists, much faster than np or df
        # generating resid is slow, otherwise everthing is fast enough
        self.__all_data = [self.__getProperty(i) for i in self.__columns]
        
    def __selectbyID(self, ids):
        if self.__all_data == None: self.__generateData() # generate data if not already done
        
        # select necessary rows; keep everything as python list, much faster than df or np
        # i-1 to convert gromacs/mpt id to list indexing
        data_list = [[row[i-1] for row in self.__all_data] for i in ids]
        
        df = pd.DataFrame(data_list, columns=self.__columns)
        df['id'] = ids
        return df.set_index(['id'])
    
    def close(self):
        self.__all_data = None
        self.__topol_dict = None
    
    def __getProperty(self, prop):
        if self.__topol_dict == None: raise MiMiCPyError("MPT file is closed")
        if self.__all_data != None: return self.__all_data[self.__columns.index(prop)]
        
        if prop == 'resid': return self.__getResID()
        
        prop_list = []
        
        if prop == 'mol':
            for mol, n_mols in self.__mol_list:
                prop_list += [mol]*len(self.__topol_dict[mol])*n_mols
        
        else:
            for mol, n_mols in self.__mol_list:
                prop_list += self.__topol_dict[mol][prop].to_list()*n_mols
        
        return prop_list 
    
    def __getResID(self):
        if self.__topol_dict == None: raise MiMiCPyError("MPT file is closed")
        resn_so_far = 0
        resn_list = []
        for mol, n_mols in self.__mol_list:
            for n in range(n_mols): # this part makes its slow
                lst = self.__topol_dict[mol]['resid'].to_numpy()+resn_so_far
                resn_list += lst.tolist()
                resn_so_far = lst[-1]
        
        return resn_list
    
    def __translate(self, selection):
        """Translates selection langauge to numpy boolean
           selection eg., resname is SER and id < 25 and mol not Protein_chain_B
           will be translated to np_vals['resname'] == 'SER' and np_vals['id'] < 25 and np_vals['mol'] != 'Protein_chain_B'
         """
         
        selection = selection.replace('(', ' ( ').replace(')', ' ) ') # put space between brackets
         
        ev = '' # converted string
        i = 0 # counter to keep track of word position
        n_brack = 0 # brackets counter
        keys = []
        for s in selection.split():
            if i == 0: # if starting of set
                if s in self.__columns or s == 'id':
                    ev += f"(np_vals['{s}']"
                    keys.append(s)
                elif s == '(':
                    ev += s
                    i = -1 # don't increment
                    n_brack += 1
                else:
                    raise SelectionError(f"{s} is not a valid selection keyword")
                    
            elif i == 1:
                if s == 'is':
                    ev += '=='
                elif s == 'not':
                    ev += '!='
                elif s == '>' or s == '>=' or s == '<' or s == '<=':
                    ev += s
                else:
                    raise SelectionError(f"{s} is not a valid logical operator")
            elif i == 2: # parse everything else, meant for the third word
                if s.isnumeric():
                    ev += f"{s})"
                else:
                    ev += f"'{s}')"
            elif i == 3:
                # if and/or encountered, reset i to -1 (will become 0 due to i+= 1 at end)
                # so we can start parsing again
                # if ) encpuntered, set i to 2 so we can parse and/or again
                if s == 'or':
                    ev += f' | '
                    i = -1
                elif s == 'and':
                    ev += f' & '
                    i = -1
                elif s == ')':
                    ev += ')'
                    i = 2
                    n_brack -= 1
                else:
                    raise SelectionError(f"{s} is not a valid boolean operator")
    
            i += 1
        
        # check brackets
        if n_brack > 0:
            raise SelectionError("Missing closing bracket is selection")
        elif n_brack < 0:
            raise SelectionError("Missing open bracket is selection")
    
        # return the translated string and the keywords used
        return ev, keys
    
    def select(self, selection):
        if selection is None or selection.strip() == '':
            raise SelectionError("The selection cannot be empty")
        
        if self.__all_data == None: self.__generateData() # generate data if not already done
        
        np_str, vals = self.__translate(selection) # get transaltes str and keywords
        
        # the keywords 'vals' has to added to np_vals dict for np_str to be executed
        np_vals = {}
        for i in vals:
            if i == 'id':
                natms = len(self.__all_data[0])
                arr = np.array(list(range(natms)))+1
            else:
                arr = np.array(self.__getProperty(i))
            np_vals[i] = arr
        
        ids = (np.where(eval(np_str))[0]+1).tolist()
        
        if ids == []:
            raise SelectionError("The selection did not return any atoms")
        
        return self.__selectbyID(ids)
    