from .._global import _Global as gbl
from .mpt_xdr import pack_strlist, pack_topol_dict, unpack_strlist, unpack_topol_dict
from ..utils.errors import SelectionError, MiMiCPyError


import xdrlib

import pandas as pd
import numpy as np
from .top import Top


def _get_columns():

    from .itp import Itp

    columns = Itp.columns.copy()
    columns.remove('number')
    columns.append('mol')

    return columns


class Mpt:
    
    columns = _get_columns()

    def __init__(self, molecules, topol_dict, mode='r'):
        self.molecules = molecules
        self.topol_dict = topol_dict
        self.__all_data = None  # What is all data ?
        self.natms = None  # Does this mean number of atoms ?
        if mode == 'r':
            self.__generate_data()
        elif mode == 'w':
            pass
        else:
            raise MiMiCPyError(f"{mode} not a mode. Only r or w can be used")


    @classmethod
    def from_top(cls, top_file, mode='r', buffer=1000, nonstandard_atom_types=None):
        top = Top(top_file, mode=mode, buffer=buffer, nonstandard_atom_types=nonstandard_atom_types)
        molecules = top.molecules
        topol_dict = top.topol_dict
        return cls(molecules, topol_dict)


    @classmethod
    def from_mpt(cls, file):
        unpacker = xdrlib.Unpacker(gbl.host.read(file, asbytes=True)) # open as bytes
        molecule_names = unpack_strlist(unpacker)
        number_of_molecules = unpacker.unpack_list(unpacker.unpack_int)
        molecules = list(zip(molecule_names, number_of_molecules))
        topol_dict = unpack_topol_dict(unpacker)
        return cls(molecules, topol_dict)


    def write(self, file_name):
        """ writes mpt file based on XDR format. Format given below:
        ##Header
         molecule names from self.molecules
         number of each molecule from self.molecules
        ##TopolDict
         repeating dictionary keys
         repeating dictionary values
         molecule name of first entry in dictionary of dataframes
         col1 of dataframe in molecule name
         col2 of dataframe in molecule name
         ... continue for all columns of df
         molecule name of second entry in dictionary of dataframes
         .... continue for all entries in dict_df
        ##End
        """

        packer = xdrlib.Packer()
        molecule_names, number_of_molecules = list(zip(*self.molecules))

        pack_strlist(packer, molecule_names)
        packer.pack_list(number_of_molecules, packer.pack_int)
        pack_topol_dict(packer, self.topol_dict)

        gbl.host.write(packer.get_buffer(), file_name, asbytes=True)




    def __getitem__(self, key):
        """Select an atom by passing the atom ID to key
            it can be a single int, list or a slice
            The indexing starts from 1
            If a string is passed as key then that property is returned
        """
        ## If string --> return property
        if isinstance(key, str):
            return self.__getProperty(key)

        ##Else it is atom ID that is requested

        elif isinstance(key, int): key = [key]
        elif isinstance(key, slice): key = list(range(key.stop)[key])

        return self.__selectbyID(key)

    def __generate_data(self):
        """Generate list of all data from topol_dict"""
        # generate everything as python lists, much faster than np or df
        # generating resid is slow, otherwise everthing is fast enough
        self.__all_data = [self.__getProperty(i) for i in self.columns]
        self.natms = len(self.__all_data[0])

    def __selectbyID(self, ids):
        if self.__all_data == None: self.__generate_data() # generate data if not already done

        # select necessary rows; keep everything as python list, much faster than df or np
        # i-1 to convert gromacs/mpt id to list indexing
        data_list = [[row[i-1] for row in self.__all_data] for i in ids]

        df = pd.DataFrame(data_list, columns=self.columns)
        df['id'] = ids
        return df.set_index(['id'])

    def close(self):
        self.__all_data = None
        self.topol_dict = None

    def __getProperty(self, prop):
        if self.topol_dict == None: raise MiMiCPyError("MPT file is closed")
        if self.__all_data != None: return self.__all_data[self.columns.index(prop)]

        if prop == 'resid': return self.__getResID()

        prop_list = []

        if prop == 'mol':
            for mol, n_mols in self.molecules:
                prop_list += [mol]*len(self.topol_dict[mol])*n_mols

        else:
            for mol, n_mols in self.molecules:
                prop_list += self.topol_dict[mol][prop].to_list()*n_mols

        return prop_list

    def __getResID(self):
        if self.topol_dict == None: raise MiMiCPyError("MPT file is closed")
        resn_so_far = 0
        resn_list = []
        for mol, n_mols in self.molecules:
            for n in range(n_mols): # this part makes its slow
                lst = self.topol_dict[mol]['resid'].to_numpy()+resn_so_far
                resn_list += lst.tolist()
                resn_so_far = lst[-1]

        return resn_list

    def __translate(self, selection):
        """Translates selection langauge to numpy boolean
           selection eg., resname is SER and id < 25 and mol not Protein_chain_B
           will be translated to np_vals['resname'] == 'SER' and np_vals['id'] < 25 and np_vals['mol'] != 'Protein_chain_B'
         """

        add_space = lambda string, a: string.replace(a, f' {a} ') # put space between and after character

        selection = add_space(add_space(selection,'('),')') # put space before/after brackets

        ev = '' # converted string
        i = 0 # counter to keep track of word position
        n_brack = 0 # brackets counter
        keys = []
        for s in selection.split():
            if i == 0: # if starting of set
                if s in self.columns or s == 'id':
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

        if self.__all_data == None: self.__generate_data() # generate data if not already done

        np_str, vals = self.__translate(selection) # get transaltes str and keywords

        # the keywords 'vals' has to added to np_vals dict for np_str to be executed
        np_vals = {}
        for i in vals:
            if i == 'id':
                arr = np.array(list(range(self.natms)))+1
            else:
                arr = np.array(self.__getProperty(i))
            np_vals[i] = arr

        ids = (np.where(eval(np_str))[0]+1).tolist()

        if ids == []:
            raise SelectionError("The selection did not return any atoms")

        return self.__selectbyID(ids)
    