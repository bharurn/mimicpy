""" Module for MiMiCPy-specific topology """

import xdrlib
import pandas as pd
import numpy as np
from .top import Top
from .itp import Itp
from .topol_dict import TopolDict
from .._global import _Global as gbl
from ..utils.errors import SelectionError, MiMiCPyError


ENCODER = 'utf-8'


def _get_itp_columns():
    columns = Itp.columns.copy()
    return columns


def _get_mpt_columns():
    columns = _get_itp_columns().copy()
    columns.remove('number')
    columns.append('mol')
    return columns


class Mpt:
    """ provides static methods for topology-specific xdr packing/unpacking,
        class methods to create new Mpt objects from top or mpt files,
        public methods for atom selection, writing and closing mpt files
    """

    columns = _get_mpt_columns()

    def __init__(self, molecules, topol_dict, mode='r'):
        self.molecules = molecules
        self.topol_dict = topol_dict
        self._expanded_data = None
        self._number_of_atoms = None

        if mode == 'r':
            self.__generate_data()
        elif mode == 'w':
            pass
        else:
            raise MiMiCPyError(f"{mode} not a mode. Only r or w can be used.")


    @staticmethod
    def __pack_strlist(packer, lst):
        s = ",".join(lst).encode(ENCODER)
        packer.pack_string(s)


    @staticmethod
    def __pack_df(packer, df):
        packer.pack_list(df.index.to_list(), packer.pack_int)  # Pack index
        for i, col in enumerate(df.columns):
            lst = df[col].to_list()
            if i in [0, 2, 3, 5]:  # Pack string columns, .i.e, atom name, type, resname, etc.
                Mpt.__pack_strlist(packer, lst)
            elif i == 1:  # Pack resid column as list of ints
                packer.pack_list(lst, packer.pack_int)
            elif i in [4, 6]:  # Pack charge and mass as list of floats
                packer.pack_list(lst, packer.pack_float)


    @staticmethod
    def __pack_topol_dict(packer, topol_dict):
        #  Pack repeating dict
        Mpt.__pack_strlist(packer, topol_dict.repeating.keys())
        Mpt.__pack_strlist(packer, topol_dict.repeating.values())

        # Pack dict of dataframes
        for k, v in topol_dict.dict_df.items():
            packer.pack_string(k.encode(ENCODER))
            Mpt.__pack_df(packer, v)


    @staticmethod
    def __unpack_strlist(unpacker):
        return unpacker.unpack_string().decode(ENCODER).split(',')


    @staticmethod
    def __unpack_df(unpacker):
        atom_numbers = unpacker.unpack_list(unpacker.unpack_int)
        atom_types = Mpt.__unpack_strlist(unpacker)
        residue_ids = unpacker.unpack_list(unpacker.unpack_int)
        residue_names = Mpt.__unpack_strlist(unpacker)
        atom_names = Mpt.__unpack_strlist(unpacker)
        charges = unpacker.unpack_list(unpacker.unpack_float)
        elements = Mpt.__unpack_strlist(unpacker)
        masses = unpacker.unpack_list(unpacker.unpack_float)
        df = pd.DataFrame([atom_numbers, atom_types,
                           residue_ids, residue_names,
                           atom_names, charges,
                           elements, masses]).T
        df.columns = _get_itp_columns()
        return df.set_index(df.columns[0])


    @staticmethod
    def __unpack_topol_dict(unpacker):
        # Unpack repeating dict
        repeating_keys = Mpt.__unpack_strlist(unpacker)
        repeating_vals = Mpt.__unpack_strlist(unpacker)
        repeating = dict(zip(repeating_keys, repeating_vals))
        if repeating == {'': ''}:
            repeating = {}

        # Unpack dataframe dict until EOF
        dict_df = {}
        while True:
            try:
                mol = unpacker.unpack_string().decode(ENCODER)
            except EOFError:
                break
            df = Mpt.__unpack_df(unpacker)
            dict_df[mol] = df

        return TopolDict(dict_df, repeating)


    @classmethod
    def from_top(cls, top_file, mode='r', buffer=1000, nonstandard_atom_types=None):
        """ Create a Mpt object from a top file """
        top = Top(top_file, mode=mode, buffer=buffer, nonstandard_atom_types=nonstandard_atom_types)
        molecules = top.molecules
        topol_dict = top.topol_dict
        return cls(molecules, topol_dict)


    @classmethod
    def from_mpt(cls, file):
        """ Create a Mpt object from an excisting mpt file """
        unpacker = xdrlib.Unpacker(gbl.host.read(file, asbytes=True))
        molecule_names = Mpt.__unpack_strlist(unpacker)
        number_of_molecules = unpacker.unpack_list(unpacker.unpack_int)
        molecules = list(zip(molecule_names, number_of_molecules))
        topol_dict = Mpt.__unpack_topol_dict(unpacker)
        return cls(molecules, topol_dict)


    def __getitem__(self, key):
        """ Select an atom by passing the atom ID to key.
            Atom ID can be a single int, list, or a slice. The indexing starts from 1.
            If a string is passed as key, then that property is returned.
        """

        if isinstance(key, str):
            return self.__get_property(key)

        if isinstance(key, int):
            key = [key]
        elif isinstance(key, slice):
            key = list(range(key.stop)[key])

        return self.__select_by_id(key)


    def __generate_data(self):
        # Generate everything as python lists, much faster than np or df
        # Generating resid is slow, otherwise everthing is fast enough
        self._expanded_data = [self.__get_property(i) for i in self.columns]
        self._number_of_atoms = len(self._expanded_data[0])


    def __select_by_id(self, ids):
        if self._expanded_data is None:
            self.__generate_data()

        # Select necessary rows; keep everything as python list, much faster than df or np
        # i-1 to convert gromacs/mpt id to list indexing
        data_list = [[row[i-1] for row in self._expanded_data] for i in ids]

        df = pd.DataFrame(data_list, columns=self.columns)
        df['id'] = ids
        return df.set_index(['id'])


    def __get_property(self, prop):
        if self._expanded_data is not None:
            return self._expanded_data[self.columns.index(prop)]

        if prop == 'resid':
            return self.__get_residue_id()

        prop_list = []
        if prop == 'mol':
            for mol, n_mols in self.molecules:
                prop_list += [mol]*len(self.topol_dict[mol]) * n_mols
        else:
            for mol, n_mols in self.molecules:
                prop_list += self.topol_dict[mol][prop].to_list() * n_mols

        return prop_list


    def __get_residue_id(self):
        resn_so_far = 0
        resn_list = []
        for mol, n_mols in self.molecules:
            for _ in range(n_mols):  # This part makes its slow
                lst = self.topol_dict[mol]['resid'].to_numpy() + resn_so_far
                resn_list += lst.tolist()
                resn_so_far = lst[-1]

        return resn_list


    def __translate(self, selection):  # Maybe use state pattern
        """ Tanslate selection langauge to numpy boolean.
            Selection 'resname is SER and id < 25 and mol not Protein_chain_B'
            will be translated to
            np_vals['resname'] == 'SER' and np_vals['id'] < 25 and np_vals['mol'] != 'Protein_chain_B'
         """

        add_space = lambda string, a: string.replace(a, f' {a} ')
        selection = add_space(add_space(selection, '('), ')')

        np_selection_expression = ''
        word_position = 0
        open_brackets = 0
        selectors = []
        for s in selection.split():
            if word_position == 0:
                if s in self.columns or s == 'id':
                    np_selection_expression += f"(np_vals['{s}']"
                    selectors.append(s)
                elif s == '(':
                    np_selection_expression += s
                    word_position = -1
                    open_brackets += 1
                else:
                    raise SelectionError(f"{s} is not a valid selection keyword")

            elif word_position == 1:
                if s == 'is':
                    np_selection_expression += '=='
                elif s == 'not':
                    np_selection_expression += '!='
                elif s == '>' or s == '>=' or s == '<' or s == '<=':
                    np_selection_expression += s
                else:
                    raise SelectionError(f"{s} is not a valid logical operator")

            elif word_position == 2:
                if s.isnumeric():
                    np_selection_expression += f"{s})"
                else:
                    np_selection_expression += f"'{s}')"

            elif word_position == 3:
                # If and' or 'or' encountered, reset word_position to -1
                # if ) encpuntered, reset word_position to 2 to parse 'and' or 'or' again
                if s == 'or':
                    np_selection_expression += f' | '
                    word_position = -1
                elif s == 'and':
                    np_selection_expression += f' & '
                    word_position = -1
                elif s == ')':
                    np_selection_expression += ')'
                    word_position = 2
                    open_brackets -= 1
                else:
                    raise SelectionError(f"{s} is not a valid boolean operator")

            word_position += 1

        if open_brackets > 0:
            raise SelectionError("Missing closing bracket is selection")
        if open_brackets < 0:
            raise SelectionError("Missing open bracket is selection")

        return np_selection_expression, selectors


    def select(self, selection):  # Maybe move to a new module
        """ Select atoms based on selection language expression """
        if selection is None or selection.strip() == '':
            raise SelectionError("The selection cannot be empty")

        if self._expanded_data is None:
            self.__generate_data()

        np_str, vals = self.__translate(selection)
        np_vals = {}

        for i in vals:
            if i == 'id':
                arr = np.array(list(range(self._number_of_atoms)))+1
            else:
                arr = np.array(self.__get_property(i))
            np_vals[i] = arr

        ids = (np.where(eval(np_str))[0]+1).tolist()

        if ids == []:
            raise SelectionError("The selection did not return any atoms")

        return self.__select_by_id(ids)


    def write(self, file_name):
        """ Write mpt file based on XDR format. Format given below:
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

        Mpt.__pack_strlist(packer, molecule_names)
        packer.pack_list(number_of_molecules, packer.pack_int)
        Mpt.__pack_topol_dict(packer, self.topol_dict)

        gbl.host.write(packer.get_buffer(), file_name, asbytes=True)
